### This script run MutSigCVsyn -- histology
import matlab.engine
import os
import time
import argparse
import os
import pandas as pd
eng = matlab.engine.start_matlab()

# add histology option parser and get histology type
feature_type = 'organ'
parser = argparse.ArgumentParser(description='Run MutSigCVsyn for input feature')
parser.add_argument('feature', metavar='HIS', type=str,
                    help='The histology that will run MutSigCV on')
arg = parser.parse_args()
feat = vars(arg)['feature']
print(feat)

# Relative path to the mutsigsyn folder
covariates = '../anno_ref/gene.covariates.converted.txt'
mut_dic = '../anno_ref/mutation_type_dictionary_file.txt'
chain_f = '../anno_ref/chr_files_hg19'
dir_cov = '../cov_out/cov_cohorts_060921_filtered'
dir_maf = '../maf_out/maf_cohorts_060121'

def runMutSigCV_syn(feat):
    global feature_type,covariates, mut_dic, chain_f,dir_cov,dir_maf,out_dir
    
    maf = os.path.join(dir_maf, feature_type,feat+'.csv')
    cov = os.path.join(dir_cov, feature_type,feat+'.csv')
    outp = os.path.join('../mutsig_out/syn/cohort_060921/', feature_type, feat, feat)

    eng.cd(r'../../mutsigsyn/', nargout = 0)
    eng.MutSigCVs_nbg(maf,cov, covariates, outp, mut_dic, chain_f, nargout = 0)
    end = time.time()

runMutSigCV_syn(feat)


