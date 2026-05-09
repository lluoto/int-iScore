"""
Molecular dynamics snapshots analysis mode with clustering.
"""

import os
import csv
import time
import numpy as np
from pathlib import Path
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ..core.metrics import (
    count_clashes,
    convert_and_soap,
    dp2_and_cpscore,
    procheck_analysis,
    compute_frustration_score,
)
from ..core.parser import parse_md_args, parse_chain_specification
from ..utils import calculate_sc
from ..utils.external_sc_ec import calculate_ec_external


def calculate_rmsd(structure1, structure2):
    """Calculate RMSD between two structures using CA atoms."""
    from Bio import PDB
    
    atoms1 = [atom for atom in structure1.get_atoms() if atom.get_name() == "CA"]
    atoms2 = [atom for atom in structure2.get_atoms() if atom.get_name() == "CA"]
    
    if len(atoms1) != len(atoms2):
        return float("inf")
    
    sup = PDB.Superimposer()
    sup.set_atoms(atoms1, atoms2)
    return sup.rms


def select_representative_spectral(affinity_matrix, member_indices):
    """Select representative structure from spectral cluster."""
    sub_affinity = affinity_matrix[np.ix_(member_indices, member_indices)]
    total_affinity = np.sum(sub_affinity, axis=1)
    rep_idx_within_cluster = np.argmax(total_affinity)
    return member_indices[rep_idx_within_cluster]


def spectral_clustering_rmsd(rmsd_matrix, valid_files, n_clusters=20):
    """Perform spectral clustering on RMSD matrix and select representatives."""
    from sklearn.cluster import SpectralClustering
    from sklearn.manifold import MDS
    
    sigma = np.median(rmsd_matrix[rmsd_matrix > 0])
    similarity_matrix = np.exp(-rmsd_matrix**2 / (2 * sigma**2))
    np.fill_diagonal(similarity_matrix, 1)
    
    spectral = SpectralClustering(
        n_clusters=n_clusters,
        affinity="precomputed",
        random_state=42,
        n_init=10
    )
    
    labels = spectral.fit_predict(similarity_matrix)
    representatives = []
    str_clu_map = {}
    
    for cluster_id in np.unique(labels):
        member_indices = np.where(labels == cluster_id)[0]
        rep_idx = select_representative_spectral(spectral.affinity_matrix_, member_indices)
        representatives.append(rep_idx)
        str_clu_map[cluster_id] = {
            "size": len(member_indices),
            "indices": member_indices.tolist(),
        }
    
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
    coords_2d = mds.fit_transform(rmsd_matrix)
    
    np.save("spe_clumap.npy", str_clu_map)
    
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(coords_2d[:, 0], coords_2d[:, 1],
                         c=labels, cmap="tab20", s=100, alpha=0.8)
    plt.colorbar(scatter, label="Cluster")
    plt.title(f"Spectral Clustering (n_clusters={len(set(labels))})")
    plt.xlabel("MDS Component 1")
    plt.ylabel("MDS Component 2")
    plt.grid(True, alpha=0.3)
    plt.savefig("spectral_clustering_mds.png")
    plt.close()
    
    return [valid_files[i] for i in representatives]


def RMSD_calculate(pdb_files, cache_dir="cache", n_clusters=20):
    """Calculate RMSD matrix and cluster structures."""
    from Bio import PDB
    
    os.makedirs(cache_dir, exist_ok=True)
    cache_file = os.path.join(cache_dir, "rmsd_matrix.npy")
    
    if os.path.exists(cache_file):
        rmsd_matrix = np.load(cache_file)
    else:
        pdbfile_parser = PDB.PDBParser()
        models = []
        valid_files = []
        
        for f in tqdm(pdb_files, desc="Parsing PDB files"):
            structure = pdbfile_parser.get_structure("temp", f)
            first_model = list(structure.get_models())[0]
            models.append(first_model)
            valid_files.append(f)
        
        n = len(models)
        if n == 0:
            return []
        
        rmsd_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                rmsd = calculate_rmsd(models[i], models[j])
                rmsd_matrix[i, j] = rmsd
                rmsd_matrix[j, i] = rmsd
        
        np.save(cache_file, rmsd_matrix)
    
    return spectral_clustering_rmsd(rmsd_matrix, pdb_files, n_clusters)


def process_md_snapshot(file_path, chain_info, args, rank=0):
    """Process a single MD snapshot."""
    from Bio import PDB
    
    pdbfile_parser = PDB.PDBParser()
    
    # Initialize default values
    cp_score = 0
    interface_residues = []
    pairs = 0
    sc_mean = 0
    interface_plddt = 0
    query_chains = list(chain_info[0])
    partner_chains = list(chain_info[1])
    
    try:
        soap, model_chain, ref_chains, new_name, dope = convert_and_soap(file_path, reference=None)
    except Exception as e:
        print(f"Error in convert_and_soap {os.path.basename(file_path)}: {e}")
        return None
    
    try:
        cp_score, interface_residues, pairs = dp2_and_cpscore(new_name, (query_chains.copy(), partner_chains.copy()))
        # Set to defaults if empty
        if not interface_residues:
            interface_residues = []
            cp_score = 0
    except Exception as e:
        print(f"Error in dp2_and_cpscore {os.path.basename(file_path)}: {e}")
        cp_score = 0
        interface_residues = []
        pairs = 0
    
    try:
        structure = pdbfile_parser.get_structure(rank, file_path)
    except Exception as e:
        print(f"Error in PDB parsing {os.path.basename(file_path)}: {e}")
        return None
    
    try:
        number_contact, bsa = count_clashes(structure, (query_chains.copy(), partner_chains.copy()))
    except Exception as e:
        print(f"Error in count_clashes {os.path.basename(file_path)}: {e}")
        number_contact, bsa = 0, 0
    
    try:
        rama = procheck_analysis(structure)
    except Exception as e:
        print(f"Error in procheck {os.path.basename(file_path)}: {e}")
        rama = {}
    
    try:
        sc_results = calculate_sc(new_name, query_chains, partner_chains)
        sc_mean = sc_results.get('sc_mean', 0.0)
    except Exception as e:
        print(f"Error in calculate_sc {os.path.basename(file_path)}: {e}")
        sc_mean = 0

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
        print(f"Error in calculate_ec {os.path.basename(file_path)}: {e}")
        ec = 0
    
    return [
        os.path.basename(file_path),
        dope,
        soap,
        compute_frustration_score(file_path, partner_chains, mode="configurational", electrostatics_k=0),  # frustration_score
        cp_score,
        bsa,
        interface_residues,
        rama,
        pairs,
        sc_mean,
        ec,
    ]


def run_md_analysis(cli_args=None):
    """Main entry point for MD snapshots mode."""
    args = cli_args if cli_args is not None else parse_md_args()
    
    SCRIPT_DIR = Path(__file__).parent.parent.parent
    os.chdir(SCRIPT_DIR)
    
    query_chains, partner_chains = parse_chain_specification(args.chains)
    chain_info = (query_chains, partner_chains)
    
    os.makedirs(args.output, exist_ok=True)
    os.makedirs(args.output, exist_ok=True)
    output_file = os.path.join(args.output, f"{args.name}_results.csv")
    
    with open(output_file, "w", newline="") as out_file:
        w = csv.writer(out_file)
        w.writerow([
            "filename", "dope", "soap", "frustration_score", "cpscore",
            "bsa", "interface_residues", "ramachandran", "pairs", "sc", "ec"
        ])
    
    pdb_files = [
        os.path.join(args.input, f)
        for f in os.listdir(args.input)
        if f.endswith(".pdb")
    ]
    
    if args.cluster_method == "spectral":
        representative_pdbs = RMSD_calculate(
            pdb_files,
            cache_dir=args.cache_dir,
            n_clusters=args.n_clusters
        )
    else:
        representative_pdbs = pdb_files
    
    if args.mpi:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        
        representative_pdbs = comm.bcast(
            representative_pdbs if rank == 0 else None, root=0
        )
    else:
        rank = 0
        size = 1
    
    files_per_process = len(representative_pdbs) // size if args.mpi else len(representative_pdbs)
    start_index = rank * files_per_process if args.mpi else 0
    end_index = (start_index + files_per_process if rank != size - 1 else len(representative_pdbs)
                 if args.mpi else len(representative_pdbs))
    
    local_files = representative_pdbs[start_index:end_index]
    
    for file_path in tqdm(local_files, desc="Processing MD snapshots"):
        try:
            result = process_md_snapshot(file_path, chain_info, args, rank)
            if result is not None:
                with open(output_file, "a", newline="") as out_file:
                    w = csv.writer(out_file)
                    w.writerow(result)
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    if args.mpi:
        MPI.Finalize()
    
    print(f"Results saved to {output_file}")
    print(f"Processed {len(representative_pdbs)} representative structures")


if __name__ == "__main__":
    run_md_analysis()
