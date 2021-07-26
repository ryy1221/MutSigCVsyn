### This notebook run coverage calculation
# Author: Yiyun

import pickle, itertools, time
import pandas as pd
from utils import *
import numpy as np
from multiprocessing import Pool
import os
import argparse

parser = argparse.ArgumentParser(description='calculate coverage for input patients')
parser.add_argument('patient_list_path', metavar='LIST_P', type=str,
                    help='The path for patient list')
arg = parser.parse_args()
lp_path = vars(arg)['patient_list_path']
lp = pickle.load(open(lp_path,'rb'))
print(lp)

### Read fasta, annotation file and name file
dict_record = pickle.load(open('/gpfs/group/epo2/default/yur97/gitlab/pcawg-to-mutsigcv/anno_ref/proc_refs/dict_record_new.pkl','rb'))
dict_transcript_info = pickle.load(open('/gpfs/group/epo2/default/yur97/gitlab/pcawg-to-mutsigcv/anno_ref/proc_refs/dict_transcript_info_062121.pkl','rb'))
dict_name_forcov = pickle.load(open('/gpfs/group/epo2/default/yur97/gitlab/pcawg-to-mutsigcv/anno_ref/proc_refs/dict_name_forcov_062121.pkl', 'rb'))

# define calculation function
def coverage_calculation(params):  # params:[transcript, patient]
    global dict_name_forcov, dict_transcript_info, pat

#     print(f'START {params[0]}')
    # get transcript sequence, cdns sequence and positions
    strand = dict_transcript_info[params[0]]['strand']
    chromosome = dict_transcript_info[params[0]]['chr'].strip('chr')
    list_cds = get_mrna_position(params[0], dict_transcript_info, strand)
    seq_transcript = get_transcript_sequence(params[0], dict_transcript_info, dict_record, strand)
    seq_cds = get_cdna_sequence(list_cds, strand, seq_transcript)
    transcript_start_pos = dict_transcript_info[params[0]]['transcript'][0]
    dict_position_cov, flag_lp = calculate_coverage(seq_transcript, seq_cds, list_cds, transcript_start_pos)
    
    # Initialize empty dataframe
    name_gene = dict_name_forcov[params[0]];zone = ['nonsilent','flank','silent'];categ = [1,2,3,4,5,6,7]
    idx = pd.MultiIndex.from_tuples(list(itertools.product([name_gene],zone,categ)))
    df_out = pd.DataFrame(index=idx, columns = [pat])

    for keys in dict_position_cov:
        coverage_patient = calculate_patient(dict_position_cov[keys], patientp[params[1]], chromosome)
        for i in range(1,8):
            df_out.loc[(name_gene,keys,i),:] = coverage_patient[i-1]
    
    return df_out

### Patient dictionary 
dir_pos = '/gpfs/scratch/yur97/patient_zero_position'
dir_out = '/gpfs/scratch/yur97/cov_out'

for p in lp:
    pat = p.split('.')[0]
    patientp = pickle.load(open(os.path.join(dir_pos,p),'rb'))
    list_patient = list(patientp.keys())
    paramlist = list(itertools.product(dict_name_forcov.keys(), list_patient))
    
    # Start processing
    start1 = time.time()
    p = Pool(20)
    res = p.map(coverage_calculation, paramlist)
    p.close()
    p.join()
    end1 = time.time()
    
    # Finish processing, save file
    cov = pd.concat(res)
    cov.to_csv(os.path.join(dir_out,pat+'.csv'), sep = '\t')
    
    print(f"finish: {pat}")
    print(f'time used: {end1 - start1}')