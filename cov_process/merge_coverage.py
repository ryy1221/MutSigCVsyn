#This script is for merging the coverage file of patient to generate a consensus coverage file

import os
from os import listdir
from os.path import isfile, join
import re
import time
import shutil
import os.path
import csv

# get all patient information
# this function parse the patient information file
def parse_info_file(info_file):
    donor_dict = {}
    with open(info_file, 'r') as pcawg_info:
        next(pcawg_info)
        for lines in pcawg_info:
            line_split = lines.split('\t')
            tumor_aliquot_id, normal_aliquot_id, donor_id, sample_id, specimen_id, project_code, gender, age \
                = line_split[0:8]
            if donor_id not in donor_dict:
                donor_dict[donor_id] = []
            donor_dict[donor_id] = tumor_aliquot_id

    return donor_dict

dict_donor = parse_info_file('/gpfs/scratch/yur97/coverage_run/PCAWG_sample_info.txt')

coverage_dir = '/gpfs/scratch/yur97/merge_coverage_test/coverage_files/'
list_file = listdir(coverage_dir)
print list_file

list_patient_id = []
for files in list_file:
    patient_id = files.split('.')[0]
    list_patient_id.append(patient_id)
print list_patient_id


# this function create a dictionary to store patient coverage values
# dict_coverage[gene][zone][categ][patient1] = x
def create_coverage_dict(patient_file, dict_coverage):
    with open(patient_file, 'r') as read_patient_file:
        next(read_patient_file)
        patient_csv = csv.reader(read_patient_file, delimiter='\t')
        for row in patient_csv:
            gene, zone, patient, categ1, categ2, categ3, categ4, categ5, categ6, categ7 = row[0:10]
            donor = dict_donor[patient]
            if gene not in dict_coverage:
                dict_coverage[gene] = {}
            if zone not in dict_coverage[gene]:
                dict_coverage[gene][zone] = {}
                dict_coverage[gene][zone][1] = {}
                dict_coverage[gene][zone][2] = {}
                dict_coverage[gene][zone][3] = {}
                dict_coverage[gene][zone][4] = {}
                dict_coverage[gene][zone][5] = {}
                dict_coverage[gene][zone][6] = {}
                dict_coverage[gene][zone][7] = {}

            dict_coverage[gene][zone][1][donor] = categ1
            dict_coverage[gene][zone][2][donor] = categ2
            dict_coverage[gene][zone][3][donor] = categ3
            dict_coverage[gene][zone][4][donor] = categ4
            dict_coverage[gene][zone][5][donor] = categ5
            dict_coverage[gene][zone][6][donor] = categ6
            dict_coverage[gene][zone][7][donor] = categ7

    return dict_coverage


coverage_dir = '/gpfs/scratch/yur97/merge_coverage_test/coverage_files/'
list_file = listdir(coverage_dir)
coverage_dict = {}
for f in list_file:
    absolute_f = join(coverage_dir, f)
    coverage_dict = create_coverage_dict(absolute_f, coverage_dict)


with open('/gpfs/scratch/yur97/merge_coverage_test/coverage_merged.txt','wb') as merge_f:
    merge_f.write('gene'+'\t'+'zone'+'\t'+'categ')
    for p in list_patient_id:
        merge_f.write('\t'+p)
        if p == list_patient_id[-1]:
            merge_f.write('\n')
    for genes in coverage_dict:
        for zones in coverage_dict[genes]:
            for i in coverage_dict[genes][zones]:
                merge_f.write(genes+'\t' + zones +'\t' + str(i))
                for patients in coverage_dict[genes][zones][i]:
                    merge_f.write('\t'+str(coverage_dict[genes][zones][i][patients]))
                    if patients == coverage_dict[genes][zones][i].keys()[-1]:
                        merge_f.write('\n')