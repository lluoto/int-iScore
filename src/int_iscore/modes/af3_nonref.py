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
from ..utils.external_sc_ec import calculate_ec_external


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
    
    # Initialize default values
    cp_score = 0
    interface_residues = []
    pairs = 0
    sc_mean = 0
    query_chains = list(chain_info[0])
    partner_chains = list(chain_info[1])
    
    try:
        if file_path.endswith(".cif"):
            soap, model_chain, ref_chains, new_name, dope = convert_and_soap(file_path, reference=None)
            pdb_path = new_name
        else:
            pdb_path = file_path
            soap, model_chain, ref_chains, new_name, dope = convert_and_soap(file_path, reference=None)
            new_name = file_path
    except Exception as e:
        print(f"Error in convert_and_soap: {e}")
        soap, model_chain, ref_chains, new_name, dope = 0, [], [], "", 0
        pdb_path = file_path
    
    working_dir = os.path.dirname(os.path.abspath(pdb_path))
    if not working_dir.endswith('/'):
        working_dir += '/'
    pdb_basename = os.path.basename(pdb_path)
    old_cwd = os.getcwd()
    os.chdir(working_dir)
    try:
        cp_score, interface_residues, pairs = dp2_and_cpscore(pdb_basename, (query_chains.copy(), partner_chains.copy()))
        # Set to defaults if empty
        if not interface_residues:
            interface_residues = []
            cp_score = 0
    except Exception as e:
        print(f"Error in dp2_and_cpscore: {e}")
        cp_score = 0
        interface_residues = []
        pairs = 0
    finally:
        os.chdir(old_cwd)
    
    try:
        structure = pdbfile_parser.get_structure(rank, pdb_path)
    except Exception as e:
        print(f"Error parsing PDB: {e}")
        structure = None
    
    try:
        if structure is not None:
            number_contact, bsa = count_clashes(structure, (query_chains.copy(), partner_chains.copy()))
        else:
            number_contact, bsa = 0, 0
    except Exception as e:
        print(f"Error in count_clashes: {e}")
        number_contact, bsa = 0, 0
    
    try:
        if structure is not None:
            rama = procheck_analysis(structure)
        else:
            rama = {}
    except Exception as e:
        print(f"Error in procheck: {e}")
        rama = {}
    
    try:
        sc_results = calculate_sc(pdb_path, query_chains, partner_chains)
        sc_mean = sc_results.get('sc_mean', 0.0)
    except Exception as e:
        print(f"Error in calculate_sc: {e}")
        sc_mean = 0

    try:
        ec_external = calculate_ec_external(
            pdb_path,
            interface_residues,
            chimera_path=getattr(args, "chimera_path", None),
            vmd_path=getattr(args, "vmd_path", None),
            pdb2pqr_path=getattr(args, "pdb2pqr_path", None),
            apbs_path=getattr(args, "apbs_path", None),
        )
        ec = ec_external if ec_external is not None else 0
    except Exception as e:
        print(f"Error in calculate_ec: {e}")
        ec = 0
    
    return [
        os.path.basename(file_path),
        dope,
        soap,
        sc_mean,
        ec,
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
    
    os.makedirs(args.output, exist_ok=True)
    output_file = os.path.join(args.output, f"{args.name}_results.csv")
    
    with open(output_file, "w", newline="") as out_file:
        w = csv.writer(out_file)
        w.writerow([
            "filename", "dope", "soap", "sc", "ec", "cpscore",
            "bsa", "interface_residues", "pairs", "ramachandran"
        ])
    
    candidate_extensions = [".cif", ".pdb"] if args.format == "auto" else [f".{args.format}"]
    file_list = [
        os.path.join(args.input, f)
        for f in sorted(os.listdir(args.input))
        if any(f.endswith(ext) for ext in candidate_extensions)
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
