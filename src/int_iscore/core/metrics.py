"""
Core metrics and calculation functions for protein interface analysis.
"""

import os
import math
import string
import tempfile
import numpy as np
import itertools
from pathlib import Path
from Bio import PDB
import freesasa
from modeller import *
from modeller import soap_pp
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces

# Import intercaat functions
import intercaat_functions as icaat

# Atomic radii for clash detection 
atom_radii = {
    "C": 1.70, 
    "N": 1.55, 
    "O": 1.52,
    "S": 1.80,
    "F": 1.47, 
    "P": 1.80, 
    "CL": 1.75, 
    "MG": 1.73,
    'ILE  CD ': 1.70,
    'TRP  OT1': 1.52,
    'TRP  OT2': 1.52,
    'CD': 1.70,
    'OT1': 1.52,
    'OT2': 1.52,
    'CD ': 1.70,
}

def procheck_analysis(structure):
    """Perform Procheck-like Ramachandran analysis."""
    procheck = [
        "AFFFFFFFFFFAAAGGDDDDDGGGGGDDDDDGGGGA",
        "AFFFFFFFFFFFAAAGGDDDDDDDDDDDDDDGGGAA",
        "AFFFFFFFFFFFAAAGGGGDDDDDDDDDDDDGGAAA",
        "AAFFFFFFFFFFFAAAAGGDDDDDDDDDDDDGGGAA",
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
        "AAAAAAAAAAAAAGGDDDDDDGGGGGGDDDDDDGGA"
    ]
    residues_count = {
        "F": 0, "A": 0, "G": 0, "D": 0,
        "gly": 0, "pro": 0, "end_res": 0,
        "total": 0, "T": 0,
    }
    residues_tags = []
    color_type = {'F': 'gray', 'A': 'lightskyblue', 'G': 'lightgreen', 'D': 'lightcoral'}

    for model in structure:
        for chain in model:
            polypeptides = PDB.CaPPBuilder().build_peptides(chain)
            for poly in polypeptides:
                phi_psi = poly.get_phi_psi_list()
                for res_idx, residue in enumerate(poly):
                    residue_tags = {}
                    phi, psi = phi_psi[res_idx]
                    if not (phi and psi):
                        residue_tags["region"] = "end_res"
                        residues_tags.append(residue_tags)
                        continue
                    
                    phi_str = "%.2f" % (phi / math.pi * 180)
                    psi_str = "%.2f" % (psi / math.pi * 180)
                    residue_tags.update({
                        "resname": residue.resname,
                        "position": str(residue.id[1]),
                        "phi": phi, "psi": psi,
                        "phi_str": phi_str, "psi_str": psi_str
                    })
                    
                    if residue.resname == "GLY":
                        residue_tags["region"] = "gly"
                    elif residue.resname == "PRO":
                        residue_tags["region"] = "pro"
                    else:
                        region = procheck[int(18 - psi / math.pi * 18)][int(phi / math.pi * 18 + 18)]
                        residue_tags["region"] = region
                    residues_tags.append(residue_tags)

    for res_tags in residues_tags:
        residues_count["total"] += 1
        if res_tags["region"] == "end_res":
            residues_count["end_res"] += 1
            continue
        if res_tags.get("resname") == "GLY":
            residues_count["gly"] += 1
        elif res_tags.get("resname") == "PRO":
            residues_count["pro"] += 1
        else:
            region = res_tags["region"]
            if region in residues_count:
                residues_count[region] += 1
                
    residues_count["T"] = sum([residues_count[k] for k in ("F", "A", "G", "D")])
    if residues_count["T"] > 0:
        for k in ("F", "A", "G", "D"):
            count = residues_count[k]
            residues_count[k] = (count, f'{count / residues_count["T"] * 100:.2f}%')
    
    return residues_count

def count_clashes(structure, chain_info, clash_cutoff=0.4):
    """Calculate atomic clashes and buried surface area."""
    query_chains = chain_info[0]
    partner_chains = chain_info[1]
    
    # Handle both list and string inputs for chains
    if isinstance(query_chains, str): query_chains = [query_chains]
    if isinstance(partner_chains, str): partner_chains = [partner_chains]
    
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
    
    clash_cutoffs = {i + "_" + j: (clash_cutoff * (atom_radii[i] + atom_radii[j])) 
                     for i in atom_radii for j in atom_radii}

    selected_ids = {id(atom) for atom in complex_atoms}
    atoms = [x for x in structure.get_atoms() if id(x) in selected_ids and x.element in atom_radii]
    coords = np.array([a.coord for a in atoms], dtype="d")
    kdt = PDB.kdtrees.KDTree(coords)

    clashes = []
    for atom_1 in atoms:
        kdt_search = kdt.search(np.array(atom_1.coord, dtype="d"), max(clash_cutoffs.values()))
        potential_clash = [(a.index, a.radius) for a in kdt_search]

        for ix, atom_distance in potential_clash:
            atom_2 = atoms[ix]
            if atom_1.parent.id == atom_2.parent.id:
                continue
            elif (atom_2.name == "C" and atom_1.name == "N") or (atom_2.name == "N" and atom_1.name == "C"):
                continue
            elif (atom_2.name == "SG" and atom_1.name == "SG") and atom_distance > 1.88:
                continue

            if atom_distance < clash_cutoffs[atom_2.element + "_" + atom_1.element]:
                clashes.append((atom_1, atom_2))

    return len(clashes) // 2, bsa

def calculate_sasa(atoms):
    """Calculate SASA for a list of atoms."""
    if not atoms:
        return 0
    atom_coords = []
    curr_radii = []
    for atom in atoms:
        atom_coords += [atom.coord[0], atom.coord[1], atom.coord[2]]
        curr_radii.append(atom_radii.get(atom.element, 1.70))
    result = freesasa.calcCoord(atom_coords, curr_radii)
    return result.totalArea()

def extract_atoms(chain):
    """Extract all atoms from a chain."""
    atoms = []
    for residue in chain:
        for atom in residue:
            atoms.append(atom)
    return atoms

def _read_contpref(contpref_file):
    """Read contact preference matrix."""
    aa = ['ILE','VAL','LEU','PHE','CYS','MET','ALA','GLY','THR','SER','TRP','TYR','PRO','HIS','GLU','GLN','ASP','ASN','LYS','ARG']
    if not os.path.exists(contpref_file):
        # Default to ones if file missing, though in practice it should be there
        return aa, np.ones((20, 20))
    mat = np.loadtxt(contpref_file)
    return aa, mat

def contact_detect(file_name, chain_info):
    """Detect inter-chain contacts using Voronoi tessellation."""
    qc = chain_info[0]
    ic = chain_info[1]
    
    if isinstance(qc, str): qc = [qc]
    if isinstance(ic, str): ic = [ic]
    
    pdb = icaat.parse(file_name, (qc + ic), './')
    coordinates = []
    for line in pdb:
        coordinates.append([line[8], line[9], line[10]])

    contacts = icaat.run_voro(coordinates)
    pdbAtomClass = icaat.appendAtomClasses(pdb)

    match = []
    count1 = 0
    for buddies in contacts:
        for buddy_idx in buddies:
            XYZ1 = [pdb[count1][2][0:2], float(pdb[count1][8]), float(pdb[count1][9]), float(pdb[count1][10])]
            XYZ2 = [pdb[buddy_idx][2][0:2], float(pdb[buddy_idx][8]), float(pdb[buddy_idx][9]), float(pdb[buddy_idx][10])]
            Ad, Vd = icaat.inter(XYZ1, XYZ2, 1.4)
            class1 = pdbAtomClass[count1]
            class2 = pdbAtomClass[buddy_idx]
            
            atom1 = '{0:<3} {1:>5} {2} {3:<4}'.format(pdb[count1][4], pdb[count1][6], pdb[count1][5], pdb[count1][2])
            atom2 = '{0:<3} {1:>5} {2} {3:<4}'.format(pdb[buddy_idx][4], pdb[buddy_idx][6], pdb[buddy_idx][5], pdb[buddy_idx][2])
            Line = '{0} | {1} | {2:<4} |    {3}   {4}'.format(atom1, atom2, str(round(Ad, 2)), str(class1), str(class2))
            
            if icaat.compatible(class1, class2):
                if Ad < Vd and any((atom1 + atom2) in sub for sub in match) == False and Ad < 5 \
                and pdb[count1][5] in qc:
                    if pdb[buddy_idx][5] in ic:
                        match.append(Line)
        count1 += 1

    newMatch = icaat.filterMatch(match, pdb, qc, 4, 'no')
    pairs = []
    interface_residues = []
    for pair in newMatch:
        curr_pre = pair.split(' ')
        curr = [element for element in curr_pre if element != '' and element != '|']
        pairs.append((curr[0], curr[4]))
        res1 = f'{curr[1]}_{curr[2]}'
        res2 = f'{curr[5]}_{curr[6]}'
        if res1 not in interface_residues: interface_residues.append(res1)
        if res2 not in interface_residues: interface_residues.append(res2)
        
    return pairs, interface_residues

def dp2_and_cpscore(file_name, chain_info):
    """Calculate CPscore using SVMlight."""
    all_pairs, interface_residues = contact_detect(file_name, chain_info)
    
    # Try to find contpref.mat in expected locations
    contpref_path = 'svmlight/contpref.mat'
    if not os.path.exists(contpref_path):
        # Try relative to script
        base_dir = Path(__file__).parent.parent.parent.parent
        contpref_path = str(base_dir / 'svmlight' / 'contpref.mat')
        
    aa, mat = _read_contpref(contpref_path)
    
    pairs = {}
    unique_residue_pairs = [tuple(sorted(pair)) for pair in all_pairs]
    for pair in unique_residue_pairs:
        key = f'{pair[0]}:{pair[1]}'
        pairs[key] = pairs.get(key, 0) + 1
        
    svm_input = f'{file_name}.svm'
    with open(svm_input, 'w') as svm_file:
        svm_file.write('0.0 ')
        sort_index = 1
        for pair in unique_residue_pairs:
            try:
                probability = mat[aa.index(pair[0])][aa.index(pair[1])]
            except ValueError:
                probability = 1.0
            count_pair = pairs[f'{pair[0]}:{pair[1]}']
            cpcount = float(count_pair) * probability / len(all_pairs) if all_pairs else 0
            svm_file.write(f'{sort_index}:{cpcount} ')
            sort_index += 1
        svm_file.write('\n')
        
    preds = []
    models_dir = 'svmlight/SVMmodels.CPscore'
    if not os.path.exists(models_dir):
        base_dir = Path(__file__).parent.parent.parent.parent
        models_dir = str(base_dir / 'svmlight' / 'SVMmodels.CPscore')
        
    svm_classify_bin = 'svmlight/svm_classify'
    if not os.path.exists(svm_classify_bin):
        base_dir = Path(__file__).parent.parent.parent.parent
        svm_classify_bin = str(base_dir / 'svmlight' / 'svm_classify')

    svm_outputs = []
    cmds = []
    if os.path.exists(models_dir):
        for i, svm_model in enumerate(os.listdir(models_dir)):
            svm_output = f'{svm_input}.{i}'
            cp_model = os.path.join(models_dir, svm_model)
            cmd = f'{svm_classify_bin} {svm_input} {cp_model} {svm_output} > /dev/null &'
            cmds.append(cmd)
            svm_outputs.append(svm_output)
        
        if cmds:
            cmds.append('wait')
            os.system("".join(cmds))
            
        for svm_output in svm_outputs:
            if os.path.exists(svm_output):
                with open(svm_output) as f:
                    content = f.read().rstrip().split()
                    if content:
                        preds.append(float(content[0]))
                os.remove(svm_output)
    
    if os.path.exists(svm_input):
        os.remove(svm_input)
        
    CPscore = np.mean(preds) if preds else 0.0
    return CPscore, interface_residues, len(unique_residue_pairs)

def calculate_dockq(model_name, reference, model_chains, reference_chains):
    """Calculate DockQ score."""
    model = load_PDB(model_name)
    native = load_PDB(reference)
    
    chain_pairs = []
    for m_chain in model_chains:
        temp = []
        for r_chain in reference_chains:
            temp.append({r_chain: m_chain})
        chain_pairs.append(temp)
        
    combined_list = []
    for combination in itertools.product(*chain_pairs):
        combined_dict = {}
        for element in combination:
            combined_dict.update(element)
        combined_list.append(combined_dict)
        
    combined_list = [c for c in combined_list if len(c) == 2]
    
    dockq_values = []
    for chain_map in combined_list:
        results = run_on_all_native_interfaces(model, native, chain_map=chain_map)
        if results and results[0]:
            for key in results[0]:
                val = results[0][key].get('DockQ', 0)
                dockq_values.append(val)
                
    return max(dockq_values) if dockq_values else 0.0

def convert_and_soap(model, reference):
    """Convert structure and calculate SOAP/DOPE scores."""
    env = Environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    
    # Create output directory if it doesn't exist
    if not os.path.exists('all_pdb'):
        os.makedirs('all_pdb')
        
    model_path = Path(model)
    new_filename = f"{model_path.parent.name}_{model_path.stem}.pdb"
    new_path = os.path.join('all_pdb', new_filename)
    
    mdl = Model(env, file=model)
    dope_scores = []
    for chain in mdl.chains:
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as handle:
            temp_chain_path = handle.name
        try:
            Selection(chain).write(file=temp_chain_path)
            single_chain_model = Model(env, file=temp_chain_path)
            dope_scores.append(single_chain_model.assess_normalized_dopehr())
        finally:
            if os.path.exists(temp_chain_path):
                os.remove(temp_chain_path)
    dope = float(np.sum(dope_scores)) if dope_scores else 0.0
    mdl.write(file=new_path)
    
    sp = soap_pp.PairScorer()
    atmsel = Selection(mdl.chains[0])
    try:
        score = atmsel.assess(sp)
    except:
        score = 0.0
        
    model_chains = [c.name for c in mdl.chains]

    ref_chains = []
    if reference is not None and os.path.exists(reference):
        ref_mdl = Model(env, file=reference)
        ref_chains = [c.name for c in ref_mdl.chains]

    return score, model_chains, ref_chains, new_path, dope
