import os
import json
import string
import csv
from mpi4py import MPI
import tqdm
import subprocess
import math
import time
from modeller import *
from modeller import soap_pp
from DockQ.DockQ import load_PDB,run_on_all_native_interfaces
import itertools
import numpy as np
import freesasa
import argparse
from Bio import PDB
import intercaat_functions as icaat
from pathlib import Path

# --- Global Parameters ---
onset_time=time.time()
SCRIPT_DIR = Path(__file__).parent.absolute()
os.chdir(SCRIPT_DIR)
rankling_list=[]
alphabet_uppercase = list(string.ascii_uppercase)
reference_list=[]
pdbfile_parser = PDB.PDBParser()
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
    "MG": 1.73

}
# --- Function Definitions ---

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


# [Previous imports remain the same...]
def time_comsuming_printer(func):
    def wrapper(*args,**kwargs):
        start_time=time.time()
        result=func(*args,**kwargs)
        end_time=time.time()
        print(f'{func.__name__},cost {end_time-start_time} s, already spent')
        return result
    return wrapper
def parse_arguments():
    
    parser = argparse.ArgumentParser(
        description="Protein complex evaluation data extraction and computation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '-i', '--input', 
        required=True,
        help="Input path of AI generated models",
        type=str
    )
    parser.add_argument(
        '-o', '--output',
        default=".",
        help="Output perliminary csv results and pdb file path for further calculation",
        type=str
    )
    parser.add_argument(
        '-c', '--composition',
        default=["A","B"],
        help="The output model chain format, default as ['A','B']",
        type=str
    )
    parser.add_argument(
        '-n', '--name',
        default='result',
        help="The output csv file prefix",
        type=str
    )
    return parser.parse_args()
parser=parse_arguments()
@time_comsuming_printer
def contact_detect(file_name,chain_info):
    """Detect inter-chain contacts using Voronoi tessellation"""
    
    qc=chain_info[0]
    ic=chain_info[1]
    pdb = icaat.parse(file_name, (qc+ic+[]), './')
    coordinates = []
    match = []
    for line in pdb:
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
    return pairs,interface_residues
@time_comsuming_printer
def count_clashes(structure,chain_info, clash_cutoff=0.4):
    """Calculate atomic clashes and buried surface area using AXVG-compatible metrics"""
    # Set what we count as a clash for each pair of atoms
    result, sasa_classes = freesasa.calcBioPDB(structure)
    # Select specific chains (e.g., chain A and chain B)
    chain_A = structure[0][chain_info[0]]  # Chain A
    chain_B = structure[0][chain_info[1]]
    # Extract atoms for chains
    chain_A_atoms = extract_atoms(chain_A)
    chain_B_atoms = extract_atoms(chain_B)
    # Calculate SASA for individual chains
    sasa_chain_A = calculate_sasa(chain_A_atoms)
    sasa_chain_B = calculate_sasa(chain_B_atoms)
    bsa=(sasa_chain_A+sasa_chain_B-result.totalArea())/2
    clash_cutoffs = {i + "_" + j: (clash_cutoff * (atom_radii[i] + atom_radii[j])) for i in atom_radii for j in atom_radii}

    # Extract atoms for which we have a radii
    atoms = [x for x in structure.get_atoms() if x.element in atom_radii]
    coords = np.array([a.coord for a in atoms], dtype="d")

    # Build a KDTree (speedy!!!)
    kdt = PDB.kdtrees.KDTree(coords)

    # Initialize a list to hold clashes
    clashes = []

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

    return len(clashes) // 2,bsa
@time_comsuming_printer
def calculate_sasa(atoms):
    atom_coords=[]
    curr_radii=[]
    for atom in atoms:
        atom_coords += [atom.coord[0], atom.coord[1], atom.coord[2]]
        curr_radii.append(atom_radii[atom.element])
    result = freesasa.calcCoord(atom_coords,curr_radii)
    return result.totalArea()
def extract_atoms(chain):
    atoms = []
    for residue in chain:
        for atom in residue:
            atoms.append(atom)
    return atoms
def _read_contpref(contpref_file):
    aa=['ILE','VAL','LEU','PHE','CYS','MET','ALA','GLY','THR','SER','TRP','TYR','PRO','HIS','GLU','GLN','ASP','ASN','LYS','ARG']
    mat=np.loadtxt(contpref_file)
    return(aa,mat)
def calculate_dockq(model_name,reference,model_chains,reference_chains):
    model=load_PDB(model_name)
    native=load_PDB(reference)
    dockq=[]
    chain_pairs=[]
    for model_chain in model_chains:
        temp=[]
        for reference_chain in reference_chains:
            temp.append({reference_chain:model_chain})
        chain_pairs.append(temp)
    combined_list = []
    for combination in itertools.product(*chain_pairs):
        combined_dict = {}
        for element in combination:
            combined_dict.update(element)  # Combine key-value pairs into one dictionary
        combined_list.append(combined_dict)
    combined_list=[file for file in combined_list if len(file)==2]
    dockq=[]
    dockqvalue=0
    for chain_map in combined_list:
        dock_key=run_on_all_native_interfaces(model, native, chain_map=chain_map)[0].keys()
        for key in dock_key:
            dockq_pre=run_on_all_native_interfaces(model, native, chain_map=chain_map)[0].get(key,'0')
            dockq+=[dockq_pre['DockQ']]
        dockqvalue= max(dockq) if len(dockq)!=0 else 0
    return dockqvalue
@time_comsuming_printer
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

    os.system(f'rm {file_name}*.svm*')
    CPscore=np.mean(preds)
    return CPscore,interface_residues,len(unique_residue_pairs)
@time_comsuming_printer
def convert_and_soap(model,reference):
    """Convert structure and calculate SOAP/DOPE scores using MODELLER, with gather reference model chain information for calculating DockQ score"""

    env = Environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    path=model.split('/')
    print(path,'check dir rebuild component')
    combine_path=[path[1],'_',path[2],path[-1][:-4],'.pdb']
    new_file=''.join(combine_path)
    pre_path=['all_pdb',new_file]
    new_path='/'.join(pre_path)
    # Read a model previously generated by Modeller's AutoModel class
    mdl = Model(env, file=model)
    dope=mdl.assess_normalized_dopehr()
    # Write the structure to a PDB file
        
    mdl.write(file=new_path)
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
    mdl = Model(env,file=reference)
    chains=mdl.chains[0:]
    reference_chains=[single.name for single in chains]
    return score,model_chains,reference_chains,new_path,dope
@time_comsuming_printer
def extract_summary_and_estimate(json_file,mark_summary):
    with open(json_file, 'r') as f:
        data = json.load(f)
    ranking_score = data['ranking_score']
    iptm = data['iptm']
    ptm = data['ptm']
    notification={'iptm':[],'ptm':[]}
    rankling_list.append([ranking_score,iptm,ptm])
    for index,chain_iptm in enumerate(data['chain_iptm']):
        if chain_iptm<0.6:
            notification['iptm']+=[(alphabet_uppercase[index],chain_iptm)]
    for index,chain_ptm in enumerate(data['chain_ptm']):
        if chain_ptm<0.5:
            notification['ptm']+=[(alphabet_uppercase[index],chain_ptm)]
    mark_summary=0
    return ranking_score,iptm,ptm,notification,mark_summary
@time_comsuming_printer
def calculate_average_plddt(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    plddt_scores = data['atom_plddts']
    average_plddt=np.mean(plddt_scores)
    return average_plddt
@time_comsuming_printer
def process_json_files(output_file,server_output):
    """Main processing function adapted for AI generated models"""
    for file_path in tqdm(file_list[start_index:end_index]):
        mark_summary=0
        print(f'{file_path} running')
        ranking_score,iptm,ptm,notification,mark_summary=extract_summary_and_estimate(file_path,mark_summary)
        average_plddt = calculate_average_plddt(file_path.replace('summary_confidences','confidences'))
        if server_output:
            average_plddt = calculate_average_plddt(file_path.replace('summary_confidences','full_data'))
        split_path=file_path.split('/')
        cif_file=file_path.replace('summary_confidences','model')[:-5]+'.cif'
        reference_name=split_path[1].upper()
        reference=f'{reference_name}_ignorechain.pdb'
        if os.path.exists(reference)==False:
            os.system(f'./clean_pdb.py {reference_name} ignorechain')
        soap,model_chain,refe_chain,new_name,dope=convert_and_soap(cif_file,reference)
        dockq=calculate_dockq(new_name,reference,model_chain,refe_chain)
        cp_score,interface_residues,pairs=dp2_and_cpscore(new_name,chain_info)
        structure = pdbfile_parser.get_structure(rank,new_name)
        model=structure[0]
        residues_plddt=[]
        for residue in interface_residues:
            residue_plddt=[]
            chain_and_id=residue.split('_')
            par_chain=model[chain_and_id[-1]]
            par_residue=par_chain[int(chain_and_id[0])]
            for atom in par_residue:
                residue_plddt.append(atom.bfactor)
            residues_plddt.append(np.mean(residue_plddt))
        interface_plddt=np.mean(residues_plddt)
        number_contact,bsa=count_clashes(structure,chain_info)
        frustration_score=0
        rama=procheck_analysis(structure)

        with open(output_file, 'a',newline='') as out_file:
            w=csv.writer(out_file)
            final_name=os.path.split(new_name)[:-4]
            output_line=[final_name,average_plddt,ranking_score,iptm,ptm,dope,soap,dockq,frustration_score,cp_score,bsa,interface_plddt,notification,interface_residues,pairs,rama,f'rank {rank} do this']
            w.writerow(output_line)
chain_info=parser.chain
input_file_path=parser.input
prefix=parser.name
output_file = f'{prefix}_1_3.csv'
output_dir=parser.output
reference_dir='reference'
final_output_file = f'{prefix}_3_3.csv'
with open(final_output_file, 'w',newline='') as out_file:
    w=csv.writer(out_file)
    output_line=['filename','average_plddt','ranking_score','iptm','ptm','dope','soap','dockq','frustrate','cpscore','interface_plddt','bsa','sc','ec','notification','interface_residue','ramachandran']
    w.writerow(output_line)
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Gather all files in the directory
if rank == 0:
    file_list = []
    for dirpath, dirnames, filenames in os.walk(input_file_path):

        for filename in filenames:
            if filename.endswith('.json') and 'summa' in filename and 'seed' in dirpath:
                file_list.append(os.path.join(dirpath, filename))
else:
    file_list = None

# Broadcast the file list to all processes
file_list = comm.bcast(file_list, root=0)

# Distribute the work among processes
files_per_process = len(file_list) // size

start_index = rank * files_per_process
end_index = start_index + files_per_process if rank != size - 1 else len(file_list)
process_json_files(output_file)

comm.Barrier()

# Finalize MPI
MPI.Finalize()
