# This script is for running MutSig
import matlab.engine
import os
import time
import matlab.engine
eng = matlab.engine.start_matlab()

covariates = '/storage/home/yur97/yur97/gitlab/pcawg-to-mutsigcv/anno_ref/gene.covariates.converted.txt'
mut_dic = '/storage/home/yur97/yur97/gitlab/pcawg-to-mutsigcv/anno_ref/mutation_type_dictionary_file.txt'
chain_f = '/storage/home/yur97/yur97/gitlab/pcawg-to-mutsigcv/anno_ref/chr_files_hg19'
org_dir = '/storage/home/yur97/yur97/gitlab/pcawg-to-mutsigcv/org_cov_maf'
his_dir = '/storage/home/yur97/yur97/gitlab/pcawg-to-mutsigcv/his_cov_maf'
out_dir = '/storage/home/yur97/yur97/gitlab/pcawg-to-mutsigcv/mutsig_out/'


def runMutSig(org, syn):
    global covariates, mut_dic, chain_f, org_dir, his_dir, out_dir

    # if organ mode than use organ directory, if not then histology directory
    if org:
        tdir = org_dir
    else:
        tdir = his_dir
    listf = os.listdir(tdir)

    # get input maf and coverage files
    for folders in listf:
        start = time.time()
        organ = folders

        if os.path.isdir(os.path.join(tdir, folders)):
            # get input file in folder, set maf and coverage input to false
            list_input = os.listdir(os.path.join(tdir, folders))
            maf_input = False
            cov_input = False
            # find inputs
            for fs in list_input:
                if fs.endswith('maf'):  # maf file
                    maf_input = os.path.join(tdir, folders, fs)
                elif fs == 'merged_' + organ + '.txt':
                    cov_input = os.path.join(tdir, folders, fs)
                elif fs.startswith('.'):
                    pass

            # check if both input exists
            if not (maf_input and cov_input):
                print "NO COVERAGE OR NO MAF FILE: " + organ
                continue

            # if syn
            if syn:
                out_p = os.path.join(out_dir, 'syn', 'out08262020', organ, organ)
                try:
                    list_out = os.listdir(os.path.join(out_dir, 'syn', 'out08262020', organ))
                    if organ + '.sig_genes.txt' in list_out:
                        print "NOT RUNNING DONE FOR: " + organ
                        continue
                except OSError:
                    pass
                # run matlab code
                eng.cd(r'../mutsigsyn/', nargout=0)
                eng.MutSigCVs_nbg(maf_input, cov_input, covariates, out_p, mut_dic, chain_f, nargout=0)
                end = time.time()
                print "TIME USED FOR " + organ + " is " + str(end - start)

            else:
                out_p = os.path.join(out_dir, 'nsyn', 'out08262020', organ, organ)
                try:
                    list_out = os.listdir(os.path.join(out_dir, 'nsyn', 'out08262020', organ))
                    if organ + '.sig_genes.txt' in list_out:
                        print "NOT RUNNING DONE FOR: " + organ
                        continue
                except OSError:
                    pass
                # run matlab code
                eng.cd(r'../mutsig/', nargout=0)
                eng.MutSigCV_nbg(maf_input, cov_input, covariates, out_p, mut_dic, chain_f, nargout=0)
                end = time.time()
                print "TIME USED FOR " + organ + " is " + str(end - start)

runMutSig(org=False, syn = True)