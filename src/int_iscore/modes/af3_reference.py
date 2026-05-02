"""
AlphaFold3 analysis mode with reference structure for DockQ calculation.
"""

import os
import json
import time
import csv
import string
from pathlib import Path
from tqdm import tqdm

from ..core.metrics import (
    count_clashes,
    convert_and_soap,
    calculate_dockq,
    dp2_and_cpscore,
    procheck_analysis,
)
from ..core.parser import parse_af3_args, parse_chain_specification
from ..utils import calculate_sc
from ..utils.external_sc_ec import calculate_ec_external


def time_comsuming_printer(func):
    """Decorator to print function execution time."""
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__}, cost {end_time - start_time:.2f}s")
        return result
    return wrapper


@time_comsuming_printer
def extract_summary_and_estimate(json_file):
    """Extract ranking scores and metrics from AlphaFold3 JSON."""
    with open(json_file, "r") as f:
        data = json.load(f)
    
    ranking_score = data["ranking_score"]
    iptm = data["iptm"]
    ptm = data["ptm"]
    
    notification = {"iptm": [], "ptm": []}
    alphabet_uppercase = list(string.ascii_uppercase)
    
    for index, chain_iptm in enumerate(data["chain_iptm"]):
        if chain_iptm < 0.6:
            notification["iptm"].append((alphabet_uppercase[index], chain_iptm))
    
    for index, chain_ptm in enumerate(data["chain_ptm"]):
        if chain_ptm < 0.5:
            notification["ptm"].append((alphabet_uppercase[index], chain_ptm))
    
    return ranking_score, iptm, ptm, notification


@time_comsuming_printer
def calculate_average_plddt(json_file, server_output=False):
    """Calculate average pLDDT from JSON file."""
    with open(json_file, "r") as f:
        data = json.load(f)
    
    key = "full_data" if server_output else "confidences"
    json_path = json_file.replace("summary_confidences", key)
    
    with open(json_path, "r") as f:
        plddt_data = json.load(f)
    
    return sum(plddt_data["atom_plddts"]) / len(plddt_data["atom_plddts"])


def process_single_file(file_path, chain_info, args, rank=0):
    """Process a single AlphaFold3 result file."""
    from Bio import PDB
    import numpy as np
    
    pdbfile_parser = PDB.PDBParser()
    query_chains = list(chain_info[0])
    partner_chains = list(chain_info[1])
    
    # Initialize default values
    cp_score = 0
    interface_residues = []
    pairs = 0
    interface_plddt = 0
    sc_mean = 0
    
    try:
        ranking_score, iptm, ptm, notification = extract_summary_and_estimate(file_path)
    except Exception as e:
        print(f"Error extracting summary: {e}")
        ranking_score, iptm, ptm, notification = 0, 0, 0, {}
    
    try:
        average_plddt = calculate_average_plddt(
            file_path.replace("summary_confidences", "confidences"),
            args.server_output
        )
    except Exception as e:
        print(f"Error calculating plddt: {e}")
        average_plddt = 0
    
    try:
        cif_file = file_path.replace("summary_confidences", "model")[:-5] + ".cif"
        reference_name = Path(args.input).name.upper()
        
        if args.reference:
            reference = args.reference
        else:
            reference = os.path.join(args.reference_dir, f"{reference_name}_ignorechain.pdb")
            if not os.path.exists(reference) and not args.no_clean_pdb:
                os.system(f"python clean_pdb.py {reference_name} ignorechain")
                generated_reference = f"{reference_name}_ignorechain.pdb"
                if os.path.exists(generated_reference):
                    reference = generated_reference
    except Exception as e:
        print(f"Error setting up reference: {e}")
        reference = None
    
    try:
        soap, model_chain, refe_chain, new_name, dope = convert_and_soap(cif_file, reference)
    except Exception as e:
        print(f"Error in convert_and_soap: {e}")
        soap, model_chain, refe_chain, new_name, dope = 0, [], [], "", 0
    
    try:
        dockq = calculate_dockq(new_name, reference, model_chain, refe_chain)
    except Exception as e:
        print(f"Error calculating DockQ: {e}")
        dockq = 0
    
    try:
        cp_score, interface_residues, pairs = dp2_and_cpscore(new_name, (query_chains.copy(), partner_chains.copy()))
        # Set to defaults if empty
        if not interface_residues:
            interface_residues = []
            cp_score = 0
    except Exception as e:
        print(f"Error in dp2_and_cpscore: {e}")
        cp_score = 0
        interface_residues = []
        pairs = 0
    
    try:
        structure = pdbfile_parser.get_structure(rank, new_name)
    except Exception as e:
        print(f"Error parsing PDB: {e}")
        structure = None
    
    try:
        sc_results = calculate_sc(new_name, query_chains, partner_chains)
        sc_mean = sc_results.get('sc_mean', 0.0)
    except Exception as e:
        print(f"Error calculating SC: {e}")
        sc_mean = 0
    
    # Calculate interface_plddt only if interface_residues exists
    residues_plddt = []
    if structure is not None and interface_residues:
        try:
            model = structure[0]
            for residue in interface_residues:
                residue_plddt = []
                parts = residue.rsplit("_", 1)
                if len(parts) == 2:
                    res_num, chain_id = parts
                    try:
                        par_chain = model[chain_id]
                        par_residue = par_chain[int(res_num)]
                        for atom in par_residue:
                            residue_plddt.append(atom.bfactor)
                    except (KeyError, ValueError):
                        continue
                if residue_plddt:
                    residues_plddt.append(np.mean(residue_plddt))
            
            if residues_plddt:
                interface_plddt = np.mean(residues_plddt)
            else:
                interface_plddt = 0
        except Exception as e:
            print(f"Error calculating interface plddt: {e}")
            interface_plddt = 0
    
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
        ec_external = calculate_ec_external(
            new_name,
            interface_residues,
            chimera_path=getattr(args, "chimera_path", None),
            vmd_path=getattr(args, "vmd_path", None),
            pdb2pqr_path=getattr(args, "pdb2pqr_path", None),
            apbs_path=getattr(args, "apbs_path", None),
        )
        ec = ec_external if ec_external is not None else 0
    except Exception as e:
        print(f"Error calculating EC: {e}")
        ec = 0
    
    return [
        os.path.basename(new_name) if new_name else os.path.basename(file_path),
        average_plddt,
        ranking_score,
        iptm,
        ptm,
        dope,
        soap,
        dockq,
        sc_mean,
        cp_score,
        bsa,
        interface_plddt,
        notification,
        interface_residues,
        pairs,
        rama,
        ec,
    ]


def run_af3_analysis(cli_args=None):
    """Main entry point for AlphaFold3 reference mode."""
    args = cli_args if cli_args is not None else parse_af3_args()
    
    SCRIPT_DIR = Path(__file__).parent.parent.parent
    os.chdir(SCRIPT_DIR)
    
    query_chains, partner_chains = parse_chain_specification(args.chains)
    chain_info = (query_chains, partner_chains)
    
    os.makedirs(args.output, exist_ok=True)
    output_file = os.path.join(args.output, f"{args.name}_results.csv")
    
    with open(output_file, "w", newline="") as out_file:
        w = csv.writer(out_file)
        w.writerow([
            "filename", "average_plddt", "ranking_score", "iptm", "ptm",
            "dope", "soap", "dockq", "sc", "cpscore", "bsa",
            "interface_plddt", "notification", "interface_residues",
            "pairs", "ramachandran", "ec"
        ])
    
    file_list = []
    for dirpath, dirnames, filenames in os.walk(args.input):
        for filename in filenames:
            if filename.endswith(".json") and "summa" in filename:
                file_list.append(os.path.join(dirpath, filename))
    
    if args.mpi:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        
        file_list = comm.bcast(file_list if rank == 0 else None, root=0)
        
        files_per_process = len(file_list) // size
        start_index = rank * files_per_process
        end_index = start_index + files_per_process if rank != size - 1 else len(file_list)
        
        local_files = file_list[start_index:end_index]
    else:
        rank = 0
        local_files = file_list
    
    for file_path in tqdm(local_files, desc="Processing AF3 models"):
        try:
            result = process_single_file(file_path, chain_info, args, rank)
            with open(output_file, "a", newline="") as out_file:
                w = csv.writer(out_file)
                w.writerow(result)
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    if args.mpi:
        MPI.Finalize()
    
    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    run_af3_analysis()
