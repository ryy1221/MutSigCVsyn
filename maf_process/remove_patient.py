from os import listdir
from os.path import join
import time

his_cov = '/Users/stella/Downloads/merged_his_cov/'
his_cov_list = listdir(his_cov)
# tum_cov = '/Users/stella/Downloads/merged_tumor_cov/'
# tum_cov_list = listdir(tum_cov)
his_maf = '/Users/stella/PycharmProjects/MutSig/merged_tumor_type_maf/categ_out/'
his_maf_list = listdir(his_maf)
# tum_maf = '/Users/stella/PycharmProjects/MutSig/merged_organ_maf'
# tum_maf_list = listdir(tum_maf)


def get_cov_patient_list(cov_f):
    list_p = []
    with open(cov_f,'r') as f:
        list_p = next(f).split('\t')[3:]

    return list_p


for his_file in his_cov_list:
    his_path = join(his_cov,his_file)
    list_patient = get_cov_patient_list(his_path)
    out = '/Users/stella/Downloads/MutSigrelated_local_input/his_maf/'+his_file+'_filtered.maf'
    for his_maff in his_maf_list:
        histology = his_maff.split('_')[6]
        if histology in his_file:
            with open(join(his_maf, his_maff),'r') as m:
                with open(out,'wb') as out_f:
                    out_f.write(next(m))
                    for lines in m:
                        patient = lines.split('\t')[12]
                        # patient = patient.replace('-','')
                        if patient in list_patient:
                            out_f.write(lines)

# for tum_file in tum_cov_list:
#     start1 = time.time()
#     his_path = join(tum_cov,tum_file)
#     list_patient = get_cov_patient_list(his_path)
#     print list_patient
#     out = '/Users/stella/Downloads/MutSigrelated_local_input/tum_maf/'+tum_file+'_filtered.maf'
#     for tum_maff in tum_maf_list:
#         organ = tum_maff.split('_')[1].split('.')[0]
#         if organ in tum_file:
#             with open(join(tum_maf, tum_maff),'r') as m:
#                 with open(out,'wb') as out_f:
#                     out_f.write(next(m))
#                     for lines in m:
#                         patient = lines.split('\t')[12]
#                         # patient = patient.replace('-','')
#                         if patient in list_patient:
#                             out_f.write(lines)
#     end1 = time.time()
#     print 'finish '+tum_file + '. Time used: '+ str(end1-start1)


