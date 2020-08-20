from os import listdir
import os
import shutil
from os.path import join
import time
import matlab.engine
eng = matlab.engine.start_matlab()

covariates = '/Users/stella/Downloads/MutSig/MutSigCV_example_data.1.0.1/gene.covariates.txt'
mut_dic = '/Users/stella/Downloads/MutSig/mutation_type_dictionary_file.txt'
chain_f = '/Users/stella/Downloads/MutSig/chr_files_hg19'

### ORGAN
org_cov_maf = '/Users/stella/Downloads/MutSigrelated_local_input/org_cov_maf'
list_org = listdir(org_cov_maf)
for folders in list_org:
    start = time.time()
    organ = folders
    print "START ORGAN: "+ organ
    out_p = join(org_cov_maf, folders) + '/' + organ + '.out.txt'
    list_f = listdir(join(org_cov_maf, folders))
    for fs in list_f:
        if fs.endswith('maf'):  # maf file
            maf_input = join(org_cov_maf, folders, fs)
        elif fs.startswith('merged'):
            cov_input = join(org_cov_maf, folders, fs)
    if not maf_input and not cov_input:
        print "NO COVERAGE OR NO MAF FILE: " + organ
        continue
    # run matlab code
    eng.cd(r'/Users/stella/Downloads/MutSig/MutSigCV_1.41/', nargout = 0)
    eng.MutSigCV(maf_input,cov_input,covariates, out_p, mut_dic, chain_f, nargout = 0)
    end = time.time()
    print "TIME USED FOR "+ organ+ " is "+ str(end-start)


### HISTOLOGY
# his_cov_maf = '/Users/stella/Downloads/MutSigrelated_local_input/his_cov_maf'
# list_his = listdir(his_cov_maf)
# for folders in list_his:
#     start = time.time()
#     histology = folders
#     print "START HISTOLOGY: "+ histology
#     out_p = join(his_cov_maf, folders) + '/'+ histology + '.out.txt'
#     list_f = listdir(join(his_cov_maf, folders))
#     for fs in list_f:
#         if fs.endswith('maf'):  # maf file
#             maf_input = join(his_cov_maf, folders, fs)
#         elif fs.startswith('merged'):
#             cov_input = join(his_cov_maf, folders, fs)
#     if not maf_input and not cov_input:
#         print "NO COVERAGE OR NO MAF FILE: " + histology
#         continue
#     # run matlab code
#     eng.cd(r'/Users/stella/Downloads/MutSig/MutSigCV_1.41/', nargout = 0)
#     eng.MutSigCV(maf_input,cov_input,covariates, out_p, mut_dic, chain_f, nargout = 0)
#     end = time.time()
#     print "TIME USED FOR "+ histology+ " is "+ str(end-start)


#######The following code run put coverage file and mutation file in a same folder
# his_cov = '/Users/stella/Downloads/MutSigrelated_local_input/merged_his_cov/'
# his_cov_list = listdir(his_cov)
# tum_cov = '/Users/stella/Downloads/MutSigrelated_local_input/merged_tumor_cov/'
# tum_cov_list = listdir(tum_cov)
# his_maf = '/Users/stella/Downloads/MutSigrelated_local_input/his_maf'
# his_maf_list = listdir(his_maf)
# tum_maf = '/Users/stella/Downloads/MutSigrelated_local_input/tum_maf'
# tum_maf_list = listdir(tum_maf)
#
# # create histology foler
# for cf in his_cov_list:
#     histology = cf.split('_')[1].split('.')[0]
#     # create directory for histology
#     os.mkdir('/Users/stella/Downloads/MutSigrelated_local_input/his_cov_maf/'+histology)
#     # move file
#     cdest = '/Users/stella/Downloads/MutSigrelated_local_input/his_cov_maf/'+histology+'/'+cf
#     csource = join(his_cov, cf)
#     shutil.move(csource, cdest)
#     for mf in his_maf_list:
#         if histology in mf:
#             msource = join(his_maf, mf)
#             mdest = '/Users/stella/Downloads/MutSigrelated_local_input/his_cov_maf/'+histology+'/'+mf
#             shutil.move(msource,mdest)
#
#
# # create organ folder
# for co in tum_cov_list:
#     organ =co.split('.')[0].split('_')[-1]
#     # create directory for organ
#     os.mkdir('/Users/stella/Downloads/MutSigrelated_local_input/org_cov_maf/' + organ)
#     # move file
#     codest = '/Users/stella/Downloads/MutSigrelated_local_input/org_cov_maf/' + organ+'/'+co
#     cosourse = join(tum_cov, co)
#     shutil.move(cosourse, codest)
#     for mo in tum_maf_list:
#         if organ in mo:
#             mosource = join(tum_maf, mo)
#             modest = '/Users/stella/Downloads/MutSigrelated_local_input/org_cov_maf/' + organ+'/'+mo
#             shutil.move(mosource, modest)
