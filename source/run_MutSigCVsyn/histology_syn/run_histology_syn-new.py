### This script run MutSigCVsyn -- histology
import matlab.engine
import os,time,argparse
import pandas as pd
eng = matlab.engine.start_matlab()

# add histology option parser and get histology type
feature_type = 'histology'
parser = argparse.ArgumentParser(description='Run MutSigCVsyn for input histology')
parser.add_argument('histology_type', metavar='HIS', type=str,
                    help='The histology that will run MutSigCV on')
parser.add_argument('--c', help='Use old covariate file or new re-annotated covariate file, enter old/new ')
parser.add_argument('--hyperm', help='exclude or not exclude hypermutator, enter yes/no')

arg = parser.parse_args()
feat = vars(arg)['histology_type'];cova = vars(arg)['c'];hyper = vars(arg)['hyperm']

print(feat);print(cova);print(hyper)

# Relative path to the mutsigsyn folder
# 1. covariate file
if cova == 'old':
    covariates = '../data/MutSigCVsyn_inputs/gene.covariates.txt'
elif cova == 'new':
    covariates = '../data/MutSigCVsyn_inputs/gene.covariates.new.txt'

# 2. mutation dictionary and chain files -- nochange
mut_dic = '../data/MutSigCVsyn_inputs/mutation_type_dictionary_file_nbg.txt'
chain_f = '../data/MutSigCVsyn_inputs/chr_files_hg19'

# 3. mutation/coverage files:
if hyper in ['no','n']:
    dir_cov = '../data/cov/histology_new/'
    dir_maf = '../data/maf/histology_new/'
elif hyper in ['yes','y']:
    dir_cov = '../data/cov/histology_new_nohypermutator/'
    dir_maf = '../data/maf/histology_new_nohypermutator/'

# 4. output directory name
out_dir = 'cohort_new'
if cova == 'new': out_dir = out_dir+'_newcova'
if hyper in ['yes','y']: out_dir = out_dir + '_nohypermutator'
print(out_dir)


# define function to run msigsyn
def runMutSigCV_syn(feat):
    global feature_type,covariates, mut_dic, chain_f,dir_cov,dir_maf,out_dir
    
    maf = os.path.join(dir_maf,feat+'.csv')
    cov = os.path.join(dir_cov,feat+'.csv')
    outp = os.path.join('../out/syn/', out_dir ,feat, feat)

    eng.cd(r'../MutSigCVsyn/', nargout = 0)
    eng.MutSigCVsyn(maf,cov, covariates, outp, mut_dic, chain_f, nargout = 0)
    end = time.time()

runMutSigCV_syn(feat)


