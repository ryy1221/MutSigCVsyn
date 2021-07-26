### This script run MutSigCVsyn -- histology
import matlab.engine
import os
import shutil
import time

eng = matlab.engine.start_matlab()

covariates = '../anno_ref/gene.covariates.converted.txt'
mut_dic = '../anno_ref/mutation_type_dictionary_file.txt'
chain_f = '../anno_ref/chr_files_hg19'
dir_cov = '../cov_out/cov_cohorts_060121'
dir_maf = '../maf_out/maf_cohorts_060121'
dir_cohorts = '../anno_ref/cohorts/'


feature_type = 'organ'
lfeat = []
with open(os.path.join(dir_cohorts,feature_type+'.csv')) as f:
    next(f)
    for lines in f:
        lfeat.append(lines.split(',')[-1].strip('\n'))
lfeat = list(set(lfeat))
print lfeat

def runMutSigCV_syn(feat):
    global feature_type,covariates, mut_dic, chain_f,dir_cov,dir_maf,out_dir
    
    maf = os.path.join(dir_maf, feature_type,feat+'.csv')
    cov = os.path.join(dir_cov, feature_type,feat+'.csv')
    outp = os.path.join('../mutsig_out/nsyn/cohort/', feature_type, feat, feat)

    eng.cd(r'../mutsig/', nargout = 0)
    eng.MutSigCV_nbg(maf,cov, covariates, outp, mut_dic, chain_f, nargout = 0)
    end = time.time()
    
for feat in lfeat:
    runMutSigCV_syn(feat)

feature_type = 'system'
lfeat = []
with open(os.path.join(dir_cohorts,feature_type+'.csv')) as f:
    next(f)
    for lines in f:
        lfeat.append(lines.split(',')[-1].strip('\n'))
lfeat = list(set(lfeat))
print lfeat

def runMutSigCV_syn(feat):
    global feature_type,covariates, mut_dic, chain_f,dir_cov,dir_maf,out_dir
    
    maf = os.path.join(dir_maf, feature_type,feat+'.csv')
    cov = os.path.join(dir_cov, feature_type,feat+'.csv')
    outp = os.path.join('../mutsig_out/nsyn/cohort/', feature_type, feat, feat)

    eng.cd(r'../mutsig/', nargout = 0)
    eng.MutSigCV_nbg(maf,cov, covariates, outp, mut_dic, chain_f, nargout = 0)
    end = time.time()
    
for feat in lfeat:
    runMutSigCV_syn(feat)

