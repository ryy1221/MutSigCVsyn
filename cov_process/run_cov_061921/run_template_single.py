# This script is the pipeline of running coverage calculation
# Author Yiyun
import re
from utils import *
import numpy as np
import pickle
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import itertools
from multiprocessing import Pool
from os.path import join
from os import listdir
import time
from optparse import OptionParser


# Use the optional parse function for input files.
def parse_arguments():
    parser = OptionParser("Usage: %prog -i <inputpath>")
    parser.add_option("-c", "--coverage_file",
                      dest="inputPath1",
                      help="The full absolute path of the coveragefile")

    (options, arg) = parser.parse_args()
    if not options.inputPath1:
        parser.error("Requires coverage file. Run with --help for help.")
        quit(1)
    return options


# define calculation function
def coverage_calculation(params):  # params:[transcript, patient]
    global dict_name, dict_transcript_info, file_out
    try:
        name_gene = dict_name[params[0]]
        # print "***START CALCULATE: " + name_gene + "***"
        # get transcript sequence, cdns sequence and positions
        strand = dict_transcript_info[params[0]]['strand']
        chromosome = dict_transcript_info[params[0]]['chr'].strip('chr')
        list_cds = get_mrna_position(params[0], dict_transcript_info, strand)
        seq_transcript = get_transcript_sequence(params[0], dict_transcript_info, dict_record, strand)
        seq_cds = get_cdna_sequence(list_cds, strand, seq_transcript)
        transcript_start_pos = dict_transcript_info[params[0]]['transcript'][0]
        dict_position_cov, flag_lp = calculate_coverage(seq_transcript, seq_cds, list_cds, transcript_start_pos)

        for keys in dict_position_cov:
            coverage_patient = calculate_patient(dict_position_cov[keys], patientp[params[1]], chromosome)

            with open(file_out, 'a') as f_cov:
                f_cov.write(name_gene + '\t' + keys + '\t' + params[1] + '\t' + str(coverage_patient[0]) + '\t' + str(
                    coverage_patient[1]) + '\t' + \
                            str(coverage_patient[2]) + '\t' + str(coverage_patient[3]) + '\t' + str(
                    coverage_patient[4]) + '\t' + str(coverage_patient[5]) \
                            + '\t' + str(coverage_patient[6]) + '\n')
                # writer = csv.writer(f_cov, delimiter='\t')
                # writer.writerows(zip(*list_print))
    except KeyError:
        pass


if __name__ == '__main__':
    arguments = parse_arguments()
    cov = arguments.inputPath1

    # output folder
    dir_out = '/gpfs/scratch/yur97/yur97/gitlab/pcawg-to-mutsigcv/cov_out/add_coverage/'
    list_out = listdir(dir_out)
    pat = cov.split('/')[-1].split('.')[0]
    pat = pat+'.csv'
    if pat in list_out:
        print "ALREADY EXIST:" +pat
        exit()

    # Load fasta file
    dict_record = SeqIO.to_dict(SeqIO.parse('/storage/home/yur97/scratch/yur97/gitlab/pcawg-to-mutsigcv/anno_ref'
                                            '/gencode_v19/GRCh37.p13.genome.fa', 'fasta'))
    # Load annotation data
    dict_name = pickle.load(open('/storage/home/yur97/scratch/yur97/gitlab/pcawg-to-mutsigcv/proc_09152020/dict_name'
                                 '.pkl', 'rb'))
    dict_transcript_info = pickle.load(open('/storage/home/yur97/scratch/yur97/gitlab/pcawg-to-mutsigcv/proc_09152020'
                                            '/dict_transcript_info.pkl', 'rb'))
    # Load annotation data
    dict_name = pickle.load(open('/storage/home/yur97/scratch/yur97/gitlab/pcawg-to-mutsigcv/proc_09152020/dict_name'
                                 '.pkl', 'rb'))
    dict_transcript_info = pickle.load(open('/storage/home/yur97/scratch/yur97/gitlab/pcawg-to-mutsigcv/proc_09152020'
                                            '/dict_transcript_info.pkl', 'rb'))
    list_addt = pickle.load(open('/storage/home/yur97/scratch/yur97/gitlab/pcawg-to-mutsigcv/proc_09152020/list_addt'
                                 '.pkl', 'rb'))

    # Read patient file
    patientp = get_zero_position(cov)
    list_patient = patientp.keys()
    paramlist = list(itertools.product(list_addt, list_patient))

    # initialize output files
    file_out = list_patient[0] + '.csv'
    file_out = join(dir_out, file_out)
    f_cov = open(file_out, "wb")
    f_cov.write(
        'gene_name' + '\t' + 'zone' + '\t' + 'patient' + '\t' + 'Categ1' + '\t' + 'Categ2' + '\t' + 'Categ3' + '\t' + \
        'Categ4' + '\t' + 'Categ5' + '\t' + 'Categ6' + '\t' + 'Categ7' + '\n')
    f_cov.close()

    # calculate patient coverage
    start1 = time.time()
    p = Pool(4)
    p.map(coverage_calculation, paramlist)
    p.close()
    p.join()
    end1 = time.time()
    print "finish: " + str(cov)
    print 'time used: ' + str(end1 - start1)
