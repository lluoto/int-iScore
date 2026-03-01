#
# Write the surface vertices of a molecular surface together with the chain
# id and residue number of the closest residue.
#
# To use this script, create the molecular surface, select it, and then
# open the script (File / Open).
#

# Get selected molecular surface models
from chimera import selection
from chimera import runCommand as rc
import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--input','-i', type=str, help='Input file name')
args = parser.parse_args()
filename=args.input
rc("open " + filename)
rc('select :.A')
rc('zone sel 5')
another_chain=selection.currentAtoms()
residues=[]
plddts=[]
pairs=[]
positions=[]
for atom in another_chain:
    plddt=atom.bfactor
    plddts.append(plddt)

    # break
    re_ch=str(atom.residue).split(' ')
    positions.append(('{}@{}'.format(re_ch[-1],atom.name),[atom.coord().x,atom.coord().y,atom.coord().z]))
    # print(re_ch)
    if re_ch[-1] not in residues:
        residues.append(re_ch[-1])
        # print(re_ch[-1])
        rc('zone {} 5'.format(re_ch[-1]))
        selected=selection.currentAtoms()
        for candidate in selected:
            if str(candidate.residue).split('.')[-1] == 'A':
                
                plddts.append(atom.bfactor)
                pair_resi=str(candidate.residue).split(' ')
                if pair_resi[-1] not in residues:
                    residues.append(pair_resi[-1])
                pair=(atom.residue.type,candidate.residue.type)
                if pair not in pairs:
                    pairs.append(pair)
                current_coord='{}@{}'.format(pair_resi[-1],candidate.name)
                if current_coord not in positions:
                    positions.append((current_coord,[candidate.coord().x,candidate.coord().y,candidate.coord().z]))

# help(atom)
try:
    #rc('surface')
    #rc('measure bur :.A :.B')
    #print('residues:',residues)
    #print(pairs)
    print(sum(plddts)/len(plddts))
    # print(positions)
except ZeroDivisionError: 
    print('no interact')
#rc('close all ')
