## This notebook simply get patient 0 positions
import os
from utils import *
from tqdm import tqdm
from multiprocessing import Pool
import pickle

dir_tracks = '/gpfs/scratch/yur97/compressed_coverage_tracks/'
dir_out = '/gpfs/scratch/yur97/patient_zero_position'
ltrack =  os.listdir(dir_tracks)

ltrack = [os.path.join(dir_tracks,i) for i in ltrack]

def get_zero_position(coverage_file):
    id_aliquot = coverage_file.split('.')[0].split('/')[-1]
    dict_patient = {}
    dict_patient[id_aliquot] = {}
    with open(coverage_file, 'r+b') as wig_f:
        mwig = mmap.mmap(wig_f.fileno(), 0, prot=mmap.PROT_READ)
        itmwig = iter(mwig.readline, b"")
        next(itmwig)
        for lines in itmwig:
            lines = lines.decode("utf-8") 
            if lines.startswith('fixed'):
                line_list = re.findall(r'\d+', lines)
                chr_n = line_list[0]
                start_pos = int(line_list[1])
                position = start_pos - 1
            else:
                position += 1
            if chr_n not in dict_patient[id_aliquot]:
                dict_patient[id_aliquot][chr_n] = []
            if lines == '0\n':
                dict_patient[id_aliquot][chr_n].append(position)
        mwig.close()
    
    pickle.dump(dict_patient, open(os.path.join(dir_out,id_aliquot+'.pkl'),'wb'))

p = Pool(20)
p.map(get_zero_position, ltrack)
p.close()
p.join()