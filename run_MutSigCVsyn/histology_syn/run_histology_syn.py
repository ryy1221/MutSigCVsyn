### This script run MutSigCVsyn -- histology
import matlab.engine
import os
import time
import argparse
import os
import pandas as pd
eng = matlab.engine.start_matlab()

# add histology option parser and get histology type
feature_type = 'histology'
parser = argparse.ArgumentParser(description='Run MutSigCVsyn for input histology')
parser.add_argument('histology_type', metavar='HIS', type=str,
                    help='The histology that will run MutSigCV on')
arg = parser.parse_args()
feat = vars(arg)['histology_type']
print(feat)

# Relative path to the mutsigsyn folder
covariates = '../anno_ref/MutSigCVsyn_inputs/gene.covariates.new.txt'
mut_dic = '../anno_ref/MutSigCVsyn_inputs/mutation_type_dictionary_file_nbg.txt'
chain_f = '../anno_ref/MutSigCVsyn_inputs/chr_files_hg19'
dir_cov = '../cov_out/cov_cohorts_072221'
dir_maf = '../maf_out/maf_cohorts_060121'

def runMutSigCV_syn(feat):
    global feature_type,covariates, mut_dic, chain_f,dir_cov,dir_maf,out_dir
    
    maf = os.path.join(dir_maf, feature_type,feat+'.csv')
    cov = os.path.join(dir_cov, feature_type,feat+'.csv')
    outp = os.path.join('../mutsig_out/syn/cohort_090821/', feature_type, feat, feat)

    eng.cd(r'../../mutsigsyn/', nargout = 0)
    eng.MutSigCVsyn(maf,cov, covariates, outp, mut_dic, chain_f, nargout = 0)
    end = time.time()

runMutSigCV_syn(feat)


