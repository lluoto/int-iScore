import os
import json
import string
import csv
from Bio import PDB
import intercaat_functions as icaat
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from scipy.stats import zscore
from sklearn.cluster import DBSCAN,AffinityPropagation
from mpi4py import MPI
import subprocess
from modeller import *
from modeller import soap_pp
import math
import itertools
import numpy as np
import csv
import ast
import subprocess
import time
from tqdm import tqdm
from pathlib import Path
temp_dir='int_iscore_temp'
os.makedirs(temp_dir,exist_ok=True)
onset_time=time.time()
SCRIPT_DIR = Path(__file__).parent.absolute()
os.chdir(SCRIPT_DIR)
rankling_list=[]
alphabet_uppercase = list(string.ascii_uppercase)
reference_list=[]
pdbfile_parser = PDB.PDBParser()

    # Write the configuration to a temporary file
# Atomic radii for clash detection 
atom_radii = {
#    "H": 1.20,  # Who cares about hydrogen??
    "C": 1.70, 
    "N": 1.55, 
    "O": 1.52,
    "S": 1.80,
    "F": 1.47, 
    "P": 1.80, 
    "CL": 1.75, 
    "MG": 1.73,
    'ILE  CD ':1.70,
    'TRP  OT1': 1.52,
    'TRP  OT2': 1.52,
    'CD':1.70,
    'OT1':1.52,
    'OT2':1.52,
    'CD ':1.70,
}
def select_representative_spectral(affinity_matrix, member_indices):
    """
    Select representative in spectral embedding space.
    This is the most natural choice for spectral clustering.
    """
    # Get submatrix for this cluster
    sub_affinity = affinity_matrix[np.ix_(member_indices, member_indices)]
    
    # Find structure with maximum total affinity to cluster members
    total_affinity = np.sum(sub_affinity, axis=1)
    rep_idx_within_cluster = np.argmax(total_affinity)
    rep_idx = member_indices[rep_idx_within_cluster]
    
    return rep_idx
def count_clashes(structure, chain_info, clash_cutoff=0.4):
    """Calculate atomic clashes and buried surface area using freesasa."""
    query_chains = [chain_info[-1]]
    partner_chains = list(chain_info[:-1])

    q_atoms = []
    for q_id in query_chains:
        if q_id in structure[0]:
            q_atoms.extend(extract_atoms(structure[0][q_id]))

    p_atoms = []
    for p_id in partner_chains:
        if p_id in structure[0]:
            p_atoms.extend(extract_atoms(structure[0][p_id]))

    sasa_q = calculate_sasa(q_atoms) if q_atoms else 0
    sasa_p = calculate_sasa(p_atoms) if p_atoms else 0
    complex_atoms = q_atoms + p_atoms
    sasa_complex = calculate_sasa(complex_atoms) if complex_atoms else 0
    bsa = (sasa_q + sasa_p - sasa_complex) / 2

    clashes = []

    clash_cutoffs = {i + "_" + j: (clash_cutoff * (atom_radii[i] + atom_radii[j])) for i in atom_radii for j in atom_radii}

    selected_ids = {id(atom) for atom in complex_atoms}
    atoms = [x for x in structure.get_atoms() if id(x) in selected_ids and x.element in atom_radii]
    coords = np.array([a.coord for a in atoms], dtype="d")

    # Build a KDTree (speedy!!!)
    kdt = PDB.kdtrees.KDTree(coords)

    

    # Iterate through all atoms
    for atom_1 in atoms:
        # Find atoms that could be clashing
        kdt_search = kdt.search(np.array(atom_1.coord, dtype="d"), max(clash_cutoffs.values()))

        # Get index and distance of potential clashes
        potential_clash = [(a.index, a.radius) for a in kdt_search]

        for ix, atom_distance in potential_clash:
            atom_2 = atoms[ix]

            # Exclude clashes from atoms in the same residue
            if atom_1.parent.id == atom_2.parent.id:
                continue

            # Exclude clashes from peptide bonds
            elif (atom_2.name == "C" and atom_1.name == "N") or (atom_2.name == "N" and atom_1.name == "C"):
                continue

            # Exclude clashes from disulphide bridges
            elif (atom_2.name == "SG" and atom_1.name == "SG") and atom_distance > 1.88:
                continue

            if atom_distance < clash_cutoffs[atom_2.element + "_" + atom_1.element]:
                clashes.append((atom_1, atom_2))
    # bsa=(sasa_chain_A+sum(part_b_list)-result.totalArea())/2
    return len(clashes) // 2,bsa
def contact_detect(file_name,chain_info):
    """Detect inter-chain contacts using Voronoi tessellation"""
    if os.path.exists(f'{file_name}_interface_residues.npy') and os.path.exists(f'{file_name}_pair_info.npy'):
        print(f"Loading precomputed contacts for {file_name}")
        interface_residues = np.load(f'{file_name}_interface_residues.npy', allow_pickle=True)
        pairs = np.load(f'{file_name}_pair_info.npy', allow_pickle=True)
        return pairs, interface_residues
    qc=[chain_info[-1]]
    ic=chain_info[:-1]
    dirs='/'.join(file_name.split('/')[:-1])
    pdb = icaat.parse(file_name, (qc+ic+[]), dirs)
    coordinates = []
    match = []
    for line in tqdm(pdb,desc=f'parsing {file_name} in contact_detect'):
        coordinates.append([line[8], line[9], line[10]])

    # Creates 3D voronoi diagram and returns indices of neighboring atoms
    contacts = icaat.run_voro(coordinates)
    # Creates a list (pdbAtomClass) that contains classes for all atoms analayzed 
    pdbAtomClass = icaat.appendAtomClasses(pdb)

    # Creates list of all atomic interactions
    count1 = 0
    count2 = 0
    for buddies in contacts:
        while count2 < len(buddies):
            # XYZ = atom, X, Y, Z 
            XYZ1 = [pdb[count1][2][0:2], float(pdb[count1][8]), float(pdb[count1][9]), float(pdb[count1][10])]
            XYZ2 = [pdb[buddies[count2]][2][0:2], \
            float(pdb[buddies[count2]][8]), float(pdb[buddies[count2]][9]), float(pdb[buddies[count2]][10])]
            # Returns the distance between two atoms and min distance reqeuired for solvent molecule
            Ad, Vd = icaat.inter(XYZ1,XYZ2,1.4)
            # Finds class of atom1 and atom2
            class1 = pdbAtomClass[count1]
            class2 = pdbAtomClass[buddies[count2]]
            # atom = residue, residue #, chain, atom
            atom1 = '{0:<3} {1:>5} {2} {3:<4}'.format(pdb[count1][4], pdb[count1][6], pdb[count1][5], pdb[count1][2])
            atom2 = '{0:<3} {1:>5} {2} {3:<4}'.format(pdb[buddies[count2]][4], pdb[buddies[count2]][6],\
                    pdb[buddies[count2]][5], pdb[buddies[count2]][2])
            Line = '{0} | {1} | {2:<4} |    {3}   {4}'.format(atom1, atom2, str(round(Ad,2)), str(class1), str(class2))
            # Only appends if classes are compatible or the user inputs that class compatibility does not matter
            # if class compatibility is unknown, the atomic interaction will be shown 
            if icaat.compatible(class1,class2) == True:
                # line 1: appends if distance between atoms is within bound and accounts for occupancy < 1
                # line 2: appends only for specific chain
                if Ad < Vd and any((atom1 + atom2) in sub for sub in match) == False and Ad < 5 \
                and pdb[count1][5] == qc[0]:
                    # appends all residues from specified chain
                    if len(qc) == 1:
                        # appends only against queried neighbor chains
                        if pdb[buddies[count2]][5] in ic:
                            match.append(Line)
                    #appends only specified residues from specified chain
                    elif pdb[count1][6] in qc:
                        # appends only against queried neighbor chains
                        if pdb[buddies[count2]][5] in ic:
                            match.append(Line)
            count2 += 1
        count2 = 0
        count1 += 1

    # Filters the list generated above (match) based on arg2 and arg4
    # Also displays interactions matrix based on arg5
    newMatch = icaat.filterMatch(match, pdb, qc, 4, 'no')
    pairs=[]
    interface_residues=[]
    for pair in newMatch:
        curr_pre=pair.split(' ')
        curr=[element for element in curr_pre if element!='' and element!='|']
        pairs.append((curr[0],curr[4]))
        if f'{curr[1]}_{curr[2]}' not in interface_residues:
            interface_residues.append(f'{curr[1]}_{curr[2]}')
        elif f'{curr[5]}_{curr[6]}' not in interface_residues:
            interface_residues.append(f'{curr[5]}_{curr[6]}')
    np.save(f'{file_name}_interface_residues.npy',interface_residues)
    np.save(f'{file_name}_pair_info.npy',pairs)
    return pairs,interface_residues

def calculate_rmsd(structure1, structure2):
    """Calculate RMSD between two structures using Biopython"""
    atoms1 = [atom for atom in structure1.get_atoms() if atom.get_name() == 'CA']
    atoms2 = [atom for atom in structure2.get_atoms() if atom.get_name() == 'CA']
    
    if len(atoms1) != len(atoms2):
        return float('inf')  # Structures have different number of CA atoms
    
    sup = PDB.Superimposer()
    sup.set_atoms(atoms1, atoms2)
    return sup.rms
def spectral_clustering_rmsd(rmsd_matrix,vaild, n_clusters=20):
    """Spectral clustering often works better than hierarchical for RMSD."""
    
    from sklearn.cluster import SpectralClustering
    from sklearn.manifold import MDS
    
    # Convert RMSD to similarity (higher RMSD = lower similarity)
    # Use Gaussian kernel: similarity = exp(-RMSD² / (2*sigma²))
    
    # Estimate sigma from RMSD distribution
    sigma = np.median(rmsd_matrix[rmsd_matrix > 0])
    similarity_matrix = np.exp(-rmsd_matrix**2 / (2 * sigma**2))
    np.fill_diagonal(similarity_matrix, 1)  # Self-similarity = 1
    
        
    # Try Spectral Clustering
    spectral = SpectralClustering(
        n_clusters=n_clusters,
        affinity='precomputed',
        random_state=42,
        n_init=10
    )
    representatives = []
    str_clu_map={}
    labels = spectral.fit_predict(similarity_matrix)
    unique_clusters = np.unique(labels)
    for cluster_id in unique_clusters:
        
        member_indices = np.where(labels == cluster_id)[0]
        rep_idx = select_representative_spectral(
                    spectral.affinity_matrix_, 
                    member_indices)
        representatives.append(rep_idx)
        str_clu_map[cluster_id] = {
            'size': len(member_indices),
            'indices': member_indices.tolist()}
    # Visualize with MDS
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    coords_2d = mds.fit_transform(rmsd_matrix)
    np.save('spe_clumap.npy',str_clu_map)
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(coords_2d[:, 0], coords_2d[:, 1], 
                         c=labels, cmap='tab20', s=100, alpha=0.8)
    plt.colorbar(scatter, label='Cluster')
    plt.title(f'Spectral Clustering (n_clusters={len(set(labels))})')
    plt.xlabel('MDS Component 1')
    plt.ylabel('MDS Component 2')
    plt.grid(True, alpha=0.3)
    # plt.show()
    plt.savefig('spectral_clustering_mds.png')
    Valid_files=[vaild[i] for i in representatives]
    return Valid_files

# Try spectral clustering

def RMSD_calculate(pdb_files, eps=0.5, min_samples=5):
    structures = []
    models = []  # Store models instead of structures
    valid_files = pdb_files
    # Read structures and get first model
    
    if os.path.exists('rmsd_matrix.npy'):
        rmsd_matrix = np.load('rmsd_matrix.npy')

    else:
        for f in tqdm(pdb_files,desc=f'parsering pdb files by rank {rank}'):
            structure = pdbfile_parser.get_structure('temp', f)
            structures.append(structure)
            
            # Get the first model
            first_model = list(structure.get_models())[0]
            models.append(first_model)
            
            valid_files.append(f)
        n = len(models)

        if n == 0:
            return []
        
        # Compute RMSD matrix
        rmsd_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                rmsd = calculate_rmsd(models[i], models[j])  # Pass models, not structures
                rmsd_matrix[i, j] = rmsd
                rmsd_matrix[j, i] = rmsd
        np.save('rmsd_matrix.npy',rmsd_matrix)
    # DBSCAN clustering
    # condensed=squareform(rmsd_matrix)
    # Z = linkage(condensed, method='average')
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    # t_cutoff = np.percentile(condensed, 95)  # Start with 95th percentile
    
    # clusters = fcluster(Z, t=t_cutoff, criterion='distance')
    # unique_clusters = len(set(clusters))
    # # 1. Plot dendrogram
    # dendrogram(Z, ax=ax1, truncate_mode='lastp', p=min(15, len(rmsd_matrix)))
    # ax1.set_title(f'hierarchical Linkage Dendrogram')
    # ax1.set_xlabel('Structure Index')
    # ax1.set_ylabel('Distance')
    
    # # Add cutoff line for clusters
    # cutoff = max(Z[-unique_clusters, 2], Z[-unique_clusters+1, 2]) if unique_clusters < len(Z) else Z[-1, 2]
    # ax1.axhline(y=cutoff, color='r', linestyle='--', alpha=0.5)
    # ax1.text(0.5, cutoff, f' {unique_clusters} clusters cutoff', 
    #          transform=ax1.get_yaxis_transform(), color='r')
    
    # # 2. Plot cluster size distribution
    # cluster_sizes = np.bincount(clusters)  # Count structures per cluster
    # cluster_sizes = cluster_sizes[1:]  # Remove index 0 (no cluster)
    
    # ax2.bar(range(1, unique_clusters + 1), cluster_sizes)
    # ax2.set_title(f'Cluster Size Distribution ({unique_clusters} clusters)')
    # ax2.set_xlabel('Cluster ID')
    # ax2.set_ylabel('Number of Structures')
    # ax2.set_xticks(range(1, unique_clusters + 1))
    
    # plt.tight_layout()
    # plt.show()
    # Try to get exactly 20 clusters
    # You might need to adjust the t parameter
    
    # rep_indices=[]
    # for cid in sorted(set(clusters)):
    #     indice=np.where(clusters==cid)[0]
    #     sub_rmsd=rmsd_matrix[np.ix_(indice,indice)]
    #     avg_rmsd=sub_rmsd.mean(axis=1)
    #     rep_idx=indice[np.argmin(avg_rmsd)]
    #     rep_indices.append(rep_idx)
    

    # The medoids are the representative structures
    # clustering = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed').fit(rmsd_matrix)
    # labels = clustering.labels_
    
    # # Form clusters
    # clusters = {}
    # noise = []
    # for i, label in enumerate(labels):
    #     if label == -1:
    #         noise.append(i)
    #     else:
    #         clusters.setdefault(label, []).append(i)
    # print()
    # all_clusters = list(clusters.values()) + [[i] for i in noise]
    # all_clusters_sorted = sorted(all_clusters, key=len, reverse=True)
    # print(all_clusters_sorted,'check dbscan data')
    # Select up to 20 representatives
    # representatives = []
    # for cluster in all_clusters_sorted[:20]:
        # Select the first structure in each cluster
        # representatives.append(cluster[0])
      # or set a fixed number
    # return quick_rmsd_clustering_top20(rmsd_matrix,n_top=20,labels=valid_files)
    spectral_reps = spectral_clustering_rmsd(rmsd_matrix,valid_files, n_clusters=20)

    return spectral_reps
    # return [valid_files[i] for i in representatives]
    # repre_hi=[valid_files[i] for i in rep_indices]
    # return repre_hi

def kabsch_rmsd_centered_fast(P, Q):
    """Compute RMSD between centered structures using Kabsch algorithm (optimized)."""
    H = P.T @ Q
    try:
        U, S, Vt = np.linalg.svd(H)
    except np.linalg.LinAlgError:
        return 0.0
    det = np.linalg.det(U) * np.linalg.det(Vt)
    d = -1 if det < 0 else 1
    trace = S.sum() if d == 1 else S[0] + S[1] - S[2]
    rmsd_sq = (np.sum(P**2) + np.sum(Q**2) - 2 * trace) / len(P)
    return np.sqrt(max(rmsd_sq, 0))



rankling_list=[]
alphabet_uppercase = list(string.ascii_uppercase)
reference_list=[]
one_to_three_letter = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
    'C': 'CYS', 'E': 'GLU', 'Q': 'GLN', 'G': 'GLY',
    'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
    'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
    'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
}

def procheck_analysis(structure):
    procheck=[ # 4 regions defined by PROCHECK
        "AFFFFFFFFFFAAAGGDDDDDGGGGGDDDDDGGGGA", # F - Favored
        "AFFFFFFFFFFFAAAGGDDDDDDDDDDDDDDGGGAA", # A - Additional allowed
        "AFFFFFFFFFFFAAAGGGGDDDDDDDDDDDDGGAAA", # G - Generously allowed
        "AAFFFFFFFFFFFAAAAGGDDDDDDDDDDDDGGGAA", # D - Disallowed
        "AAFFFFFFFFFFFFAAGGGDDDDDDDDDDDDGGGGA",
        "AAFFFFFFFFFFFAAAGGGDDDDDDDDDDDDDGGGA",
        "AAAFFFFFFFFFAAAAGGGDDDGGGGGGDDDDDGGA",
        "AAAAFFFFFFFAAAAAAGGGGGGGGGGGDDDDDGGA",
        "AAAAAFAAFAAAAAAAAGGGGGGGAAGGDDDDDGGA",
        "AAAAAAAAAAAAAAGGGGGGGAAAAGGGDDDDDGGG",
        "AAAAAAAAAAAAGGGGGGGGGAAAAGGGDDDDDGGG",
        "GAAAAAAAAAAAAGGDDDDGGGAAAGGGGDDDDGGG",
        "GGAAAAAAAAAAAGGDDDDGGAAAAAAGGDDDDGGG",
        "GGAAAAAAAAAAGGGDDDDGGAAFAAGGGDDDDGGG",
        "GAAAAAAAAAAAGGGDDDDGGGAFFAGGGDDDDGGG",
        "GAAAAAAFAAAAAGGGDDDGGGAAAAAGGDDDDGGG",
        "GAAAAAFFFFAAAGGGGDDDGGGAAAGGGDDDDGGG",
        "GAAAAAFFFFFAAAGGGDDDGGGGAAAGGDDDDGGG",
        "GAAAAFFFFFFFAAAGGGDDGGGAGAAGGDDDDGGG",
        "GAAAAAFFFFFFFAAGGGGDGGGGGGGGGDDDDGGG",
        "GGAAAAFFFFFFFFAAGGGDGGGGGGGGGDDDDGGG",
        "GGAAAAAFFFFFFFAAAGGGDDDDDDDDDDDDDGGG",
        "GGGAAAAAFFFFFFFAAGGGDDDDDDDDDDDDDGGG",
        "GAAAAAAAAFFFFFFAAAGGGDDDDDDDDDDDDGGG",
        "AAGAAAAAAAAFFFFAAAGGGDDDDDDDDDDDDGGG",
        "GGGAAAAAAAAAAAAAAAAGGDDDDDDDDDDDDGGG",
        "GGGGGAAAAAAAAAAAAAGGGDDDDDDDDDDDDGGG",
        "DGGGGAAAAAAAAAAGGGGGGDDDDDDDDDDDDDDD",
        "DDDGGGGGGGGAGGGGGGGGDDDDDDDDDDDDDDDD",
        "DDDGGGGAAGGGGGGGGDDDDDDDDDDDDDDDDDDD",
        "GGGGGAAGGGGGGGDDDDDDDGGGGGDDDDDDDDDD",
        "GGGGGGAAAAGGGGDDDDDDDGGGGGDDDDDDDDDD",
        "GAAAGAAAAAGGGGGDDDDDDGGAGGGDDDDDDDDD",
        "GAAAAAAAAAAGGGGDDDDDDGGAGGGDDDDDDGGG",
        "GAAAAAAAAAAAAGGDDDDDDGGAAGGDDDDDDGGG",
        "AAAAAAAAAAAAAGGDDDDDDGGGGGGDDDDDDGGA"]
    residues_count = {"F": 0,       # Favored
                                "A": 0,       # Additional allowed
                                "G": 0,       # Generously allowed
                                "D": 0,       # Disallowed
                                "gly": 0,     # Glycines
                                "pro": 0,     # Prolines
                                "end_res": 0, # end-residues
                                "total": 0,   # total number of residues
                                "T": 0,
                                }
    residues_tags = []
    color_type={'F':'gray','A':'lightskyblue','G':'lightgreen','D':'lightcoral'}

    plot_data_x=[]
    plot_data_y=[]
    color_data=[]
    for model in structure:
        for chain in model:
            polypeptides = PDB.CaPPBuilder().build_peptides(chain)
            for poly_idx,poly in enumerate(polypeptides):
                phi_psi = poly.get_phi_psi_list()
                for res_idx, residue in enumerate(poly):
                    residue_tags = {}
                    phi, psi = phi_psi[res_idx]
                    if not (phi and psi):
                        residue_tags["region"] = "end_res"
                        residues_tags.append(residue_tags)
                        continue
                    het, resseq, icode = residue.id
                    phi_str = "%.2f" % (phi/math.pi*180)
                    psi_str = "%.2f" % (psi/math.pi*180)
                    plot_data_x.append(float(phi_str))
                    plot_data_y.append(float(psi_str))
                    # key = residue.resname + str(resseq)
                    residue_tags.update({"resname": residue.resname,
                                            "position": str(resseq),
                                            "phi": phi, "psi": psi,
                                            "phi_str": phi_str, "psi_str": psi_str})
                    # Glycines.
                    if residue.resname == "GLY":
                        residue_tags["region"] = "gly"
                        color_data.append('gray')
                    # Prolines.
                    elif residue.resname == "PRO":
                        residue_tags["region"] = "pro"
                        color_data.append('gray')
                    # Other residues.
                    else:
                        region = procheck[int(18-psi/math.pi*18)][int(phi/math.pi*18+18)]
                        residue_tags["region"] = region
                        color_data.append(color_type[region])
                    residues_tags.append(residue_tags)

    for res_idx, res_tags in enumerate(residues_tags):
        residues_count["total"] += 1
        if res_tags["region"] == "end_res":
            residues_count["end_res"] +=1
            continue
        # Glycines.
        if res_tags["resname"] == "GLY":
            residues_count["gly"] += 1
        # Prolines.
        elif res_tags["resname"] == "PRO":
            residues_count["pro"] += 1
        # Other residues.
        else:
            if res_tags["region"] in residues_count:
                residues_count[res_tags["region"]] += 1
    residues_count["T"] = sum([residues_count[k] for k in ("F", "A", "G", "D")])
    residues_count["F"] = (residues_count["F"],f'{residues_count["F"]/residues_count["T"]*100:.2f}%')
    residues_count["A"] = (residues_count["A"],f'{residues_count["A"]/residues_count["T"]*100:.2f}%')
    residues_count["G"] = (residues_count["G"],f'{residues_count["G"]/residues_count["T"]*100:.2f}%')
    residues_count["D"] = (residues_count["D"],f'{residues_count["D"]/residues_count["T"]*100:.2f}%')
    if float(residues_count["F"][-1][:-1])<95:
        residues_count["warning"] = 'True'
    else:
        residues_count["warning"] = 'False'
    return residues_count

# Accessing elements using two indices
# Example: Accessing the value at index (row_index1, row_index2)
def _read_contpref(contpref_file):
    aa=['ILE','VAL','LEU','PHE','CYS','MET','ALA','GLY','THR','SER','TRP','TYR','PRO','HIS','GLU','GLN','ASP','ASN','LYS','ARG']
    mat=np.loadtxt(contpref_file)
    return(aa,mat)
def dp2_and_cpscore(file_name,chain_info):
    all_pairs,interface_residues=contact_detect(file_name,chain_info)
    aa,mat=_read_contpref('svmlight/contpref.mat')
    pairs={}
    unique_residue_pairs = [tuple(sorted(pair)) for pair in all_pairs]
    # Print the unique pairs
    for pair in unique_residue_pairs:
        if not pairs.get(f'{pair[0]}:{pair[1]}'):
            pairs[f'{pair[0]}:{pair[1]}']=1
        else:
            pairs[f'{pair[0]}:{pair[1]}']+=1
    sort_index=1
    svm_input=f'{file_name}.svm'
    with open(svm_input,'w') as svm_file:
        svm_file.write('0.0 ')
        for pair in unique_residue_pairs:
            probability=mat[aa.index(pair[0])][aa.index(pair[1])]
            count_pair=pairs[f'{pair[0]}:{pair[1]}']
            cpcount=float(count_pair)*probability/len(all_pairs)
            current_pair=f'{sort_index}:{cpcount} '
            svm_file.write(current_pair)
            sort_index+=1
        svm_file.write('\n')
    svm_file.close()
    preds=[]
    svm_outputs=[]
    cmds=[]
    for i,svm_model in enumerate(os.listdir('svmlight/SVMmodels.CPscore')):
        svm_output=f'{svm_input}.{i}'
        cp_model=os.path.join('svmlight/SVMmodels.CPscore',svm_model)
        cmd=f'svmlight/svm_classify {svm_input} {cp_model} {svm_output} > /dev/null &'
        cmds.append(cmd)
        svm_outputs.append(svm_output)
    cmds.append('wait')
    os.system("".join(cmds))
    for svm_output in svm_outputs:
        with open(svm_output) as f:
            pred=f.read().rstrip().split()[0]
            preds.append(float(pred))
    try:
        os.remove(f'{file_name}*.svm*')
    except IsADirectoryError:
        pass  # Skip directories silently
    except FileNotFoundError:
        pass  # File already deleted (race condition)
    
    CPscore=np.mean(preds)
    return CPscore,interface_residues,len(unique_residue_pairs)

def convert_and_soap(model):
    env = Environ()
    path=model.split('/')    

    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    # Read a model previously generated by Modeller's AutoModel class
    mdl = Model(env, file=model)
    chains=mdl.chains
    dope_scores=[]
    for chain_id in range(len(chains)):
        sel = Selection(mdl.chains[chain_id])  # 创建当前链的选择对象
        temp_filename = f'{model}_chain_{chain_id}.pdb'  # 临时文件名
        sel.write(file=temp_filename)  # 将当前链写入单独的 PDB 文件
        singlemodel=Model(env, file=temp_filename)  # 读取单独的链作为新的模型
        dope = singlemodel.assess_normalized_dopehr()  # 评估归一化 DOPE
        dope_scores.append(dope)
        try:
            os.remove(temp_filename)  # 删除临时文件
        except IsADirectoryError:
            pass  # Skip directories silently
        except FileNotFoundError:
            pass 

    dope=np.sum(dope_scores)


# Write the structure to a PDB file
    sp = soap_pp.PairScorer()

    atmsel = Selection(mdl.chains[0])

    # Assess with the above Scorer
    try:
        score = atmsel.assess(sp)

    except ModellerError:
        print("The SOAP-Protein-OD library file is not included with MODELLER.")
        print("Please get it from https://salilab.org/SOAP/.")

    chains=mdl.chains[0:]
    model_chains=[single.name for single in chains]
    return score,model_chains,dope
   
def calculate_average_plddt(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    plddt_scores = data['atom_plddts']
    total_plddt = sum(plddt_scores)
    average_plddt = total_plddt / len(plddt_scores)
    return average_plddt
def process_json_files(output_file):

    # Each process processes its assigned files
    for file_path in tqdm(file_list[start_index:end_index]):
        split_file=file_path.split('/')
        soap,model_chain,dope=convert_and_soap(file_path)
        cp_score,interface_residues,pairs=dp2_and_cpscore(file_path,model_chain)

        structure = pdbfile_parser.get_structure('0',file_path)
        number_contact,bsa=count_clashes(structure,model_chain)
        frustration_score=0

        
        rama=procheck_analysis(structure)
        ec=0
        sc=0
        with open(output_file, 'a',newline='') as out_file:
            w=csv.writer(out_file)

            output_line=[split_file[-1],dope,soap,frustration_score,cp_score,bsa,interface_residues,rama,pairs,sc,ec,f'rank {0} do this']
            w.writerow(output_line)
        print(f'Finished processing {file_path}')


output_file = 'md_3070_0082_hi_1_3.csv'
with open(output_file, 'w',newline='') as out_file:
    w=csv.writer(out_file)
    output_line=['filename','dope','soap','frustration_score','cpscore','bsa','interface_residues','rama','pairs','ec']
    w.writerow(output_line)

# Gather all files in the directory

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Gather all files in the directory
input_dir='/mnt/sdb2/lluoto/split_frames'
if rank==0:
    if os.path.exists('affinity_str.npy'):
        # pdb_files = [os.path.join(input_dir, file) for file in os.listdir(input_dir) 
                        # if '3070_0038' not in file and 'test' not in file and file.endswith('.pdb')]
        # pdb_files=np.load('affinity_str.npy',allow_pickle=True).tolist()
        # representative_pdbs = RMSD_calculate(pdb_files, eps=0.5, min_samples=5)
        # print(representative_pdbs,'check representative pdbs',len(representative_pdbs))
        # for i in representative_pdbs:
        #     os.system(f'cp {i} /mnt/sdb2/lluoto/3070_0082_repre')
        
        # representative_pdbs=np.load('affinity_str.npy',allow_pickle=True).tolist()
        # representative_pdbs=[i.replace('/mnt/sdb2/lluoto/', '/mnt/sdb2/') for i in representative_pdbs ]
        # representative_pdb2=np.load('3070_0082_representative_hi.npy',allow_pickle=True).tolist()
        # print(representative_pdb1,representative_pdb2)
        # for pdb in representative_pdb1:
        representative_pdbs=[os.path.join(input_dir, file) for file in  '''3070_0082-repeat_3_0975.pdb
3070_0082-repeat_3_0898.pdb
3070_0082-repeat_3_0829.pdb
3070_0082-repeat_3_0774.pdb
3070_0082-repeat_3_0549.pdb
3070_0082-repeat_3_0418.pdb
3070_0082-repeat_3_0310.pdb
3070_0082-repeat_3_0233.pdb
3070_0082-repeat_3_0178.pdb
3070_0082-repeat_3_0127.pdb
3070_0082-repeat_1_0991.pdb
3070_0082-repeat_1_0800.pdb
3070_0082-repeat_1_0731.pdb
3070_0082-repeat_1_0618.pdb
3070_0082-repeat_1_0527.pdb
3070_0082-repeat_1_0458.pdb
3070_0082-repeat_1_0337.pdb
3070_0082-repeat_1_0277.pdb
3070_0082-repeat_1_0232.pdb
3070_0082-repeat_1_0107.pdb'''.split('\n')]
        #     representative_pdb2.append(pdb)
        # representative_pdbs=list(set(representative_pdb2))
        # print("Loaded representative structures from file:", representative_pdbs)
        
        # pdb_files = [os.path.join(input_dir, file) for file in os.listdir(input_dir) 
        #                 if '3070_0038' not in file and 'test' not in file and file.endswith('.pdb')]
        np.save('affinity_str.npy',representative_pdbs)
    else:
        pdb_files = [os.path.join(input_dir, file) for file in os.listdir(input_dir) 
                        if '3070_0038' not in file and 'test' not in file and file.endswith('.pdb')]
        representative_pdbs = RMSD_calculate(pdb_files, eps=0.5, min_samples=5)
        print("Representative structures:", representative_pdbs)
        # np.save('3070_0082_representative_hi.npy',representative_pdbs)
else:
    representative_pdbs = None
# Broadcast the file list from rank 0 to all processes
file_list = comm.bcast(representative_pdbs, root=0)

# Now all processes have the same file_list
# Broadcast the file list to all processes


# Distribute the work among processes
files_per_process = len(file_list) // size

start_index = rank * files_per_process
end_index = start_index + files_per_process if rank != size - 1 else len(file_list)
process_json_files(output_file)


comm.Barrier()

# Finalize MPI
MPI.Finalize()

print((representative_pdbs))


# rankling_list.sort(key=lambda x: x[1],reverse=True)  # Sort by the first value of each dictionary
# print(rankling_list)
