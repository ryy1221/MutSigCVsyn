import shutil
import os

source = '/gpfs/scratch/yur97/compressed_coverage_tracks/'
dest1 = '/gpfs/scratch/yur97/coverage_tracks_greylist/'

list_white = []
with open('/gpfs/scratch/yur97/MutSigCV/PCAWG_sample_info.txt', 'r') as f:
    next(f)
    for lines in f:
        lines_split = lines.split('\t')
        tumor_id = lines_split[0]
        list_white.append(tumor_id)

files = os.listdir(source)

for f in files:
    f_name = f.split('.')[0]
    try:
        if f_name not in list_white:
            shutil.move(source+f, dest1)
    except shutil.Error:
        pass



