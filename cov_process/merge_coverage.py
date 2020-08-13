#This script is for merging the coverage file of patient to generate a consensus coverage file

from os import listdir
from os.path import isfile, join
from optparse import OptionParser
import csv

def parse_arguments():
    parser = OptionParser("Usage: %prog -i <inputpath>")
    parser.add_option("-i", "--info_file",
                      dest="inputPath",
                      help="The full absolute path of the PCAWG ICGC sample info file")

    parser.add_option("-d", "--cov_dir",
                     dest="inputPath1",
                     help="The full absolute path of the coverage directory")

    parser.add_option("-o", "--output_path",
                      dest="inputPath1",
                      help="The full absolute path of the output folder")

    (options, arg) = parser.parse_args()
    if not options.inputPath:
        parser.error("Requires at least one input file, the path of this script. Run with --help for help.")
        quit(1)
    return options

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
            histology_abbreviation = line_split[12]
            if donor_id not in donor_dict:
                donor_dict[donor_id] = []
            donor_dict[donor_id] = tumor_aliquot_id

    return donor_dict

# get patient list in organ
def get_patient_list(coverage_organ_dir):
    list_file = listdir(coverage_organ_dir)
    list_patient_id = []
    for files in list_file:
        patient_id = files.split('.')[0]
        list_patient_id.append(patient_id)

    return list_patient_id


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


if __name__ == '__main__':
    arguments = parse_arguments()
    info_f = arguments.inputPath
    cov_dir = arguments.inputPath1
    out_p = arguments.inputPath2

    dict_donor = parse_info_file(info_f)
    list_dirs = listdir(cov_dir)
    for dirs in list_dirs:
        organ_folder_path = join(cov_dir, dirs)
        organ_type = dirs
        # get patient list for this organ
        patient_list = get_patient_list(organ_folder_path)


        # get file list in this organ type
        organ_f_list = listdir(organ_folder_path)
        coverage_dict = {}
        for f in organ_f_list:
            p_f = join(organ_folder_path, f)
            coverage_dict = create_coverage_dict(p_f, coverage_dict)


        # creat output file
        with open(out_p+'merged'+organ_type+'.txt', 'wb') as merge_f:
            merge_f.write('gene'+'\t'+'zone'+'\t'+'categ')
            for p in patient_list:
                merge_f.write('\t'+p)
                if p == patient_list[-1]:
                    merge_f.write('\n')
            for genes in coverage_dict:
                for zones in coverage_dict[genes]:
                    for i in coverage_dict[genes][zones]:
                        merge_f.write(genes+'\t' + zones +'\t' + str(i))
                        for patients in coverage_dict[genes][zones][i]:
                            merge_f.write('\t'+str(coverage_dict[genes][zones][i][patients]))
                            if patients == coverage_dict[genes][zones][i].keys()[-1]:
                                merge_f.write('\n')


