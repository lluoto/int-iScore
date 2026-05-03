"""
Non-reference AlphaFold3 analysis mode (without DockQ).
"""

import os
import csv
import json
import string
from pathlib import Path
from tqdm import tqdm

from ..core.metrics import (
    count_clashes,
    convert_and_soap,
    dp2_and_cpscore,
    procheck_analysis,
)
from ..core.parser import parse_nonref_args, parse_chain_specification
from ..utils import calculate_sc


def extract_summary_and_estimate(json_file):
    """Extract ranking scores and check for issues."""
    with open(json_file, "r") as f:
        data = json.load(f)
    
    ranking_score = data["ranking_score"]
    iptm = data["iptm"]
    ptm = data["ptm"]
    has_clash = data.get("has_clash", 0)
    
    notification = {"iptm": [], "ptm": []}
    alphabet_uppercase = list(string.ascii_uppercase)
    
    for index, chain_iptm in enumerate(data["chain_iptm"]):
        if chain_iptm < 0.6:
            notification["iptm"].append((alphabet_uppercase[index], chain_iptm))
    
    for index, chain_ptm in enumerate(data["chain_ptm"]):
        if chain_ptm < 0.5:
            notification["ptm"].append((alphabet_uppercase[index], chain_ptm))
    
    mark_summary = 0
    if ranking_score < 0 or iptm < 0.6 or ptm < 0.5:
        mark_summary = 1
    if has_clash != 0:
        mark_summary = 1
    
    return ranking_score, iptm, ptm, notification, mark_summary


def calculate_average_plddt(json_file):
    """Calculate average pLDDT from JSON."""
    with open(json_file, "r") as f:
        data = json.load(f)
    
    plddt_scores = data["atom_plddts"]
    return sum(plddt_scores) / len(plddt_scores)


def process_nonref_file(file_path, chain_info, args, rank=0):
    """Process a single non-reference AlphaFold3 model."""
    from Bio import PDB
    
    pdbfile_parser = PDB.PDBParser()
    
    if file_path.endswith(".cif"):
        soap, model_chain, ref_chains, new_name, dope = convert_and_soap(file_path)
        pdb_path = new_name
    else:
        pdb_path = file_path
        soap, model_chain, ref_chains, new_name, dope = convert_and_soap(file_path)
        new_name = file_path
    
    working_dir = os.path.dirname(os.path.abspath(pdb_path))
    cp_score, interface_residues, pairs = dp2_and_cpscore(pdb_path, chain_info, working_dir)
    
    structure = pdbfile_parser.get_structure(rank, pdb_path)
    number_contact, bsa = count_clashes(structure, chain_info)
    rama = procheck_analysis(structure)
    
    query_chains, partner_chains = chain_info
    sc_results = calculate_sc(pdb_path, query_chains, partner_chains)
    sc_mean = sc_results.get('sc_mean', 0.0)
    
    return [
        os.path.basename(file_path),
        dope,
        soap,
        sc_mean,
        cp_score,
        bsa,
        interface_residues,
        pairs,
        rama,
    ]


def run_nonref_analysis(cli_args=None):
    """Main entry point for non-reference AF3 mode."""
    args = cli_args if cli_args is not None else parse_nonref_args()
    
    SCRIPT_DIR = Path(__file__).parent.parent.parent
    os.chdir(SCRIPT_DIR)
    
    query_chains, partner_chains = parse_chain_specification(args.chains)
    chain_info = (query_chains, partner_chains)
    
    os.makedirs(args.temp_dir, exist_ok=True)
    
    output_file = os.path.join(args.output, f"{args.name}_results.csv")
    
    with open(output_file, "w", newline="") as out_file:
        w = csv.writer(out_file)
        w.writerow([
            "filename", "dope", "soap", "sc", "cpscore",
            "bsa", "interface_residues", "pairs", "ramachandran"
        ])
    
    extension = ".cif" if args.format == "cif" else ".pdb"
    file_list = [
        os.path.join(args.input, f)
        for f in os.listdir(args.input)
        if f.endswith(extension)
    ]
    
    for file_path in tqdm(file_list, desc="Processing AF3 models"):
        try:
            result = process_nonref_file(file_path, chain_info, args)
            with open(output_file, "a", newline="") as out_file:
                w = csv.writer(out_file)
                w.writerow(result)
            
            if file_path.endswith(".cif"):
                pdb_path = file_path.replace(".cif", ".pdb")
                if os.path.exists(pdb_path):
                    import shutil
                    shutil.copy(pdb_path, args.temp_dir)
        
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    run_nonref_analysis()
