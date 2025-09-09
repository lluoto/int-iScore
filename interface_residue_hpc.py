import os
import string
import csv
import ast
import logging
from multiprocessing import Pool 
import subprocess
import numpy as np
alphabet_uppercase = list(string.ascii_uppercase)
output_sheet= 'new_bench_2_3.csv'
sheet_path='new_bench_1_3.csv'
all_lines=[]
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

# Display the DataFrame

one_to_three_letter = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
    'C': 'CYS', 'E': 'GLU', 'Q': 'GLN', 'G': 'GLY',
    'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
    'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
    'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
}
# Accessing elements using two indices
# Example: Accessing the value at index (row_index1, row_index2)
with open(sheet_path, 'r') as out_file:
    w=csv.reader(out_file)
    for line in w:
        all_lines.append(line)
def _read_contpref(contpref_file):
    aa=['ILE','VAL','LEU','PHE','CYS','MET','ALA','GLY','THR','SER','TRP','TYR','PRO','HIS','GLU','GLN','ASP','ASN','LYS','ARG']
    mat=np.loadtxt(contpref_file)
#    print(mat.shape)
    return(aa,mat)
def dp2_and_cpscore(file_name):

    pairs={}
    interface_residues=''
    result=subprocess.run(['chimera','--nogui','--nostatus','--script',f"surfvert.py -i {file_name}"],stdout=subprocess.PIPE)
    # aj=result.stdout.decode('utf-8').split('\n')
    if result.stdout.decode('utf-8')=='no interact':
        return False,False,False,False
    else:
        run_output=result.stdout.decode('utf-8').split('\n')
        print(run_output)
        all_pairs=ast.literal_eval(run_output[-3])
        # print(count)
    interface_residues=run_output[-4]
    interface_average_plddt=run_output[-2]
    bsa_index=run_output.index('Buried solvent accessible surface area')+1
    final_bsa=run_output[bsa_index].split('=')[-1].strip(' ')
    aa,mat=_read_contpref('svmlight/contpref.mat')
    pairs={}
    unique_residue_pairs = [tuple(sorted(pair)) for pair in all_pairs]
    # Convert back to a list if needed

    # Print the unique pairs
    for pair in unique_residue_pairs:

        if not pairs.get(f'{pair[0]}:{pair[1]}'):
            pairs[f'{pair[0]}:{pair[1]}']=1
        else:
            pairs[f'{pair[0]}:{pair[1]}']+=1
    sort_index=1
    svm_input=f'{file_name}.svm'
    # svm_input='test.txt'
    with open(svm_input,'a') as svm_file:
        svm_file.write('0.0 ')
        # svm_file.write("\n")
        for pair in unique_residue_pairs:
            # try:
            probability=mat[aa.index(pair[0])][aa.index(pair[1])]
            count_pair=pairs[f'{pair[0]}:{pair[1]}']
            cpcount=float(count_pair)*probability/len(all_pairs)
            current_pair=f'{sort_index}:{cpcount} '
            svm_file.write(current_pair)
            # svm_file.write('\n')
            sort_index+=1
        svm_file.write('\n')

       
        svm_file.close()
        preds=[]
        svm_outputs=[]
        cmds=[]
        for i,svm_model in enumerate(os.listdir('svmlight/SVMmodels.CPscore')):
            #print(svm_model)
            svm_output=f'{svm_input}.{i}'
            cp_model=os.path.join('svmlight/SVMmodels.CPscore',svm_model)
            cmd=f'svmlight/svm_classify {svm_input} {cp_model} {svm_output} > /dev/null &'
            #print(cmd)
            cmds.append(cmd)
            svm_outputs.append(svm_output)

        cmds.append('wait')
        # print(cmds)
        os.system("".join(cmds))
        for svm_output in svm_outputs:
            with open(svm_output) as f:
                pred=f.read().rstrip().split()[0]
                #pred=subprocess.check_output(f'cat {svm_output}', shell=True,stderr=subprocess.STDOUT).decode('UTF-8').rstrip().split()[0]
                preds.append(float(pred))

            
        #os.system(f'cp {tmpdir}/* /home/x_bjowa/proj/local/ProQDock/bar10/')

        os.system(f'rm {file_name}*.svm*')
        CPscore=np.mean(preds)

    
    return CPscore,interface_residues,interface_average_plddt,final_bsa,len(unique_residue_pairs)

    

pdb_path='output/'
# pdb_path='PBS/1zli_dock'
file_list=[]
for filename in os.listdir(pdb_path):
    if filename.endswith('.pdb'):
        file_list.append(os.path.join(pdb_path, filename))
def main(file_list):
    for filename in file_list:

        # filename='chan_g/benchmark/benchmark_relax/fold_2024_09_11_19_48_2ido_20_model_1_0001.pdb'
        log.info(f'Processing file: {filename}')
        pdb_name=filename.split('/')[-1][:-4]
        cp_score,interface_residues,average_plddt,final_bsa,pair_number=dp2_and_cpscore(filename)
        frustration_score=0
        name=filename.split[-1]
        os.system(f'cp {filename} frustratometer2-master/ &wait')
        os.system(f'perl frustratometer2-master/RunFrustratometer.pl {name} configurational 0')
        with open(f'{name}.done/{name}_configurational','r') as file:
            interface_frustraction=[]
            f=file.readlines()
            for i in f:
                if 'A B' in i:
                    current_line=i.split(' ')
                    interface_frustraction.append(float(current_line[-3]))
                    frustration_score=np.sum(interface_frustraction)
        
        os.system(f'rm -rf frustratometer2-master/{name}*')
        for index,line in enumerate(all_lines[1:]):
            if line[0]==pdb_name:
            # if line[0]==pdb_name[:-5]:

                all_lines[index]=line[:-1]+[cp_score,average_plddt,final_bsa,frustration_score]+[line[-1]]+[interface_residues,pair_number]+[f'rank {rank} do this']
                with open(output_sheet,'a',newline='') as file:
                    write=csv.writer(file)
                    write.writerow(all_lines[index])
file_list=[file_list]
pool = Pool(8)
pool.map(main, file_list) 
pool.close()
pool.join()