import os
import intercaat_functions as icaat

qc=['E']
ic_list=['A','B','C','D']
# ic='B'
pdb = icaat.parse('/media/cuixi/data0/lluoto/split_frames/3070_0082-repeat_3_0969.pdb', (qc+ic_list+[]), '')
coordinates = []
match = []
for line in pdb:
    coordinates.append([line[8], line[9], line[10]])
print('get done')
contacts = icaat.run_voro(coordinates)
print('done')