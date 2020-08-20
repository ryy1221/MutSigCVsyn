# This code is for process PCAWG ICGC consensus somatic pass only maf files. The output files are inputs for MutSigCV
# Author: Yiyun Rao
# Date last edit:05222020

from os import listdir
from os.path import join

def get_tumor_list(list_f):
    list_tumor = []
    for f in list_f:
        organ_type = f.split('-')[0].split('_')[-1]
        if organ_type not in list_tumor:
            list_tumor.append(organ_type)

    return list_tumor


maf_dir = '/Users/stella/PycharmProjects/MutSig/merged_tumor_type_maf/categ_out/'
merged_maf_list = listdir(maf_dir)
tumor_list = get_tumor_list(merged_maf_list)


for his in tumor_list:
    list_done = []
    f_name = '/Users/stella/PycharmProjects/MutSig/merged_organ_maf/merged_'+his+'.maf'
    tmp = '/Users/stella/PycharmProjects/MutSig/tmp.txt'
    with open(tmp, 'wb') as f:
        for his_maf in merged_maf_list:
            if his in his_maf and his_maf not in list_done:
                with open(join(maf_dir,his_maf),'rb') as in_f:
                        f.write(in_f.read())
                        list_done.append(his_maf)
    tmp_f =open(tmp,'r')
    lines = tmp_f.readlines()
    tmp_f.close()

    with open(f_name, 'wb') as f:
        n=0
        for line in lines:
            n+=1
            if n==1:
                f.write(line)
            if n >1 and not line.startswith('Hugo'):
                f.write(line)

