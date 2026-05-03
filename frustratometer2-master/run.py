import numpy as np
import os
import csv
from mpi4py import MPI
import subprocess
all_lines=[]
with open('/mnt/sdb2/lluoto/int-iScore/md_3070_0082_hi_2_3.csv','r') as file:
    re=csv.reader(file)
    for line in re:
        all_lines.append(line)
input_dir='/mnt/sdb2/lluoto/3070_0082_repre'

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
work_list=[line[0] for line in all_lines]
# Gather all files in the directory
if rank == 0:
    file_list = []
    for input_file in os.listdir(input_dir):
        if input_file.endswith('.pdb'):
            file_list.append(os.path.join(input_dir, input_file))
else:
    file_list = None

# Broadcast the file list to all processes
file_list = comm.bcast(file_list, root=0)

# Distribute the work among processes
files_per_process = len(file_list) // size
start_index = rank * files_per_process
end_index = start_index + files_per_process if rank != size - 1 else len(file_list)
output_path='.'
for input_file in file_list[start_index:end_index]:
        
        os.system(f'cp {input_file} {output_path} &wait')
        input_file=input_file.split('/')[-1]
        name=input_file[:-4]
        os.system(f'perl RunFrustratometer.pl {input_file} configurational 0')
        with open(f'{input_file}.done/{input_file}_configurational','r') as file:
            interface_frustraction=[]
            f=file.readlines()
            for i in f:
                if 'A E' in i or 'B E' in i or 'C E' in i or 'D E' in i:
                    current_line=i.split(' ')
                    interface_frustraction.append(float(current_line[-3]))
            frustration_score=np.sum(interface_frustraction)
            for line in all_lines:
                print(line,input_file,'check')
                if input_file[:-4] in line[0]:
                    new_line=line
                    print(frustration_score,'score')
                    new_line+=[frustration_score]
                    with open('/mnt/sdb2/lluoto/int-iScore/md_3070_0082_hi_3_3.csv','a') as new_file:
                        re=csv.writer(new_file)
                        re.writerow(new_line)
        
        #os.system(f'rm -rf {name}*')


