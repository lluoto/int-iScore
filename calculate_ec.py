

import numpy as np
import ast
import argparse
import csv
import os
import subprocess

# --- Function Definitions ---
def get_args():
    parser = argparse.ArgumentParser(description='calculate the ec on data')
    parser.add_argument('--model','-m',   type=str, default='', help='apbs outputed dx file')
    parser.add_argument('--sc','-q',   type=str, default='', help='acquire sc score')
    parser.add_argument('--output','-o',   type=str, default='', help='accepted sc score')
    parser.add_argument('--name','-n',   type=str, default='', help='file prefix of csv file')
    parser.add_argument('--chimera_path','-c',   type=str, default='', help='Chimera or ChimeraX for patch mapping')

    return parser.parse_args()

def read_pqr(pqr_file_path):
    """processing pqr data for setting steric distribution"""
    # Open and read the file
    with open(pqr_file_path, "r") as file:
        lines = file.readlines()

    # Parse the file line by line
    atoms = []
    for line in lines:
        if line.startswith('TER'):
            atoms.append('boundary')
            continue
        elif line.startswith("ATOM") :  # Identify atom records
            record = {
                "atom_name": line[12:16].strip(),
                "res_id": int(line[22:26].strip()),
                "xyz": [float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip())],
                "radius": float(line[62:70].strip())
            }
            if record["atom_name"].startswith('H'):
                pass
            else:
                atoms.append(record)

    # Print parsed atoms
    return atoms

def map_atoms(pqr_file,interface_residue,chimera_path):
    '''map spatial electrostatic potential patches to it's atom/residue'''
    result=subprocess.run([f'{chimera_path}/bin/chimera','--nogui','--script',"potential_patch.py -i {} -d {}.dx".format(pqr_file,pqr_file)],stdout=subprocess.PIPE)
    run_output=result.stdout.decode('utf-8').split('\n')
    vertex_value=[]
    for line in run_output:
        if line.startswith('['):
            vertex_value.append(ast.literal_eval(line))
    atom_coords=[]
    atom_radis={}
    pqr_atoms=read_pqr(pqr_file)
    index_boundary=pqr_atoms.index('boundary')

    for residue in interface_residue:
        current_atom=residue.split('_')
        if current_atom[-1]=='A':
            for position in pqr_atoms[:index_boundary]:
                if position['res_id']==int(current_atom[0]) and 'H' not in position['atom_name']:# # Build a KD-tree for atom coordinates
                    atom_coords.append(position['xyz'])
                    atom_radis[str(position['xyz'])]=[residue,position['radius']]
        elif current_atom[-1]=='B':
            for position in pqr_atoms[index_boundary+1:-1]:
                if position['res_id']==int(current_atom[0]) and 'H' not in position['atom_name']:# # Build a KD-tree for atom coordinates
                    atom_coords.append(position['xyz'])
                    atom_radis[str(position['xyz'])]=[residue,position['radius']]
    potential=[]
    vertex=[]
    for line in vertex_value:
        vertex.append(line[:-1])
        potential.append(line[-1])

    vertex=np.array(vertex)
    atom_coords=np.array(atom_coords)
    if vertex.size == 0:
        vertex = vertex.reshape(0, 3)
    else:
        # Move dimension of size 3 to last position if exists
        if vertex.shape[-1] != 3 and 3 in vertex.shape:
            pos = list(vertex.shape).index(3)
            vertex = np.moveaxis(vertex, pos, -1)
        
        # Ensure last dimension has exactly 3 coordinates
        if vertex.shape[-1] < 3:
            raise ValueError("Vertex array has < 3 coordinates")
        elif vertex.shape[-1] > 3:
            vertex = vertex[..., :3]  # Truncate extra coordinates
        
        # Flatten to 2D array
        if vertex.ndim > 2:
            vertex = vertex.reshape(-1, 3)
        elif vertex.ndim == 1:
            vertex = vertex.reshape(1, 3)
    
    # Compute distances (now safe with standardized shapes)
    diff = atom_coords[:, np.newaxis, :] - vertex[np.newaxis, :, :]
    distances = np.linalg.norm(diff, axis=2)
    residue_potentials_dict= {key: [] for key in interface_residue}
    # Map each vertex to its closest atom
    for index_i in range(len(atom_coords)):
        for index_j in range(len(vertex)):
            coord='['+', '.join(map(str, atom_coords[index_i]))+']'
            if atom_radis[coord][-1]+1.3<distances[index_i,index_j]<atom_radis[coord][-1]+1.5:
                residue_potentials_dict[atom_radis[coord][0]].append(potential[index_j])
    return residue_potentials_dict

def calculate_residue_scores(residue_potentials_dict):
    """
    Calculate the score for each residue based on potentials.
    
    Parameters:
    - residue_potentials_dict: Dictionary where keys are residue indices and values are lists 
      of potentials for each residue.

    Returns:
    - scores: List of calculated scores for each residue.
    """
    all_numerator=0
    all_denominator=0
    # Extract potentials and keys (residue indices) from dictionary

    residue_potentials = list(residue_potentials_dict.values())
    
    # Iterate through each residue in the dictionary
    for  check,resi_own in enumerate(residue_potentials):

        if len(resi_own)==0 or len(resi_own)==1:
            pass
        else:
        # Convert phi_own to a numpy array
            for index,phi_current in enumerate(resi_own):
                phi_own = np.array(phi_current)
            
                other_potentials = [resi_own[i] for i in range(len(resi_own)) if i != index]
                

                phi_rest = np.array(other_potentials)

                phi_mean = np.mean(resi_own)
                phi_prime_mean = np.mean(phi_rest)

                numerator = (phi_own - phi_mean) * (np.sum(phi_rest) - phi_prime_mean)
                # print(check,index,numerator)
                denominator = np.sqrt(np.sum((phi_own - phi_mean) ** 2) * np.sum((phi_rest - phi_prime_mean) ** 2))
                all_numerator+=numerator
                all_denominator+=denominator
    Em = -all_numerator / all_denominator if all_denominator != 0 else 0
    return Em
# --- Global Parameters ---
args=get_args()
pqr_file=args.model
all_lines=[]
file_dir=args.output
chimera_path=args.chimera_path
prefix=args.name
with open(os.path.join(file_dir,f'{prefix}_1_3.csv'),'r') as file:
    re=csv.reader(file)
    for line in re:
        all_lines.append(line)
#write into file right now

output_sheet= os.path.join(file_dir,f'{prefix}_2_3.csv')
chimera_path
for index,line in enumerate(all_lines):
    if line[0]==os.path.split(pqr_file)[:-4]:
        print(ast.literal_eval)
        residue_potentials_dict=map_atoms(pqr_file,ast.literal_eval(line[-3]),chimera_path)
        score=calculate_residue_scores(residue_potentials_dict)
        sc_score=args.sc
        new_line=line[:-4]+[sc_score,score]+line[-4:]
        with open(output_sheet,'a',newline='') as file:
            write=csv.writer(file)
            write.writerow(new_line)
