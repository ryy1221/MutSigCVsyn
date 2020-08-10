from optparse import OptionParser
import multiprocessing
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import time
import re
from functools import partial
from contextlib import contextmanager
from os import listdir
from os.path import isfile, join
import csv
from bisect import bisect_left
import itertools

# nucleotides and condon table
nuc = ["A", "T", "C", "G"]
table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}
gc_di_nuc = "CG"  # definition of CG dinucleotide, from 5' to 3'
stop_aa = ["TAA", "TAG", "TGA"]
start_aa = "ATG"


# Error definition
class FracError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return 'FracError, {0} '.format(self.message)
        else:
            return 'FracError has been raised'


class CategError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return 'CategError, {0} '.format(self.message)
        else:
            return 'CategError has been raised'


class AnnotationError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return 'AnnotationError, {0} '.format(self.message)
        else:
            return 'AnnotationError has been raised'


# Use the optional parse function for input files.
def parse_arguments():
    parser = OptionParser("Usage: %prog -i <inputpath>")
    parser.add_option("-a", "--annotation_file",
                      dest="inputPath1",
                      help="The full absolute path of the GENCODEv19 annotation file")
    parser.add_option("-f", "--fasta_file",
                      dest="inputPath2",
                      help="The full absolute path of the GENCODEv19 genome fasta file")
    parser.add_option("-i", "--patient_info_file",
                      dest="inputPath3",
                      help="The full absolute path of the PCAWG reference file")
    parser.add_option("-o", "--output", dest="outputPath", default=False,
                      help="The full absolute path of the output csv file")
    parser.add_option("-c", "--coverage_file", dest="inputPath4", default=False,
                      help="If use -c, will enter full coverage mode")

    (options, arg) = parser.parse_args()
    if not options.inputPath1 or not options.inputPath2 or not options.inputPath3:
        parser.error("Requires annotation file, fasta file and patient information. Run with --help for help.")
        quit(1)
    return options

#
# def grouper(iterable, n, fillvalue=None):
#     "Collect data into fixed-length chunks or blocks"
#     # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
#     args = [iter(iterable)] * n
#     return izip_longest(fillvalue=fillvalue, *args)


# Determine if the mutation is a transition mutation.
def is_transition(ori_allele, alt_allele):
    p = False
    transition = [("A", "G"), ("G", "A"), ("T", "C"), ("C", "T")]
    transversion = [("C", "A"), ("C", "G"), ("G", "T"), ("T", "G"), ("G", "C"), ("A", "C"), ("A", "T"), ("T", "A")]
    if (ori_allele, alt_allele) in transition:
        p = True
    elif (ori_allele, alt_allele) in transversion:
        p = False
    else:
        print ori_allele, alt_allele
        raise CategError("NEITHER TRANSITION NOR TRANSVERSION MUTATION")

    return p


# if base in gc dinucleotide
def is_cpg(base, base_before, base_after):
    flag = False
    if (base + base_after == gc_di_nuc) or (base_before + base == gc_di_nuc):
        flag = True

    return flag


# if base is A:T pair
def is_at(base):
    flag = False
    if base == "A" or base == "T":
        flag = True

    return flag


# if a mutation cause silent codon
def is_silent(before_context, after_context):
    before_aa = table[before_context]
    after_aa = table[after_context]
    flag = False
    if before_aa == after_aa:
        flag = True

    return flag


# just add count to whatever dictionary passed
def add_count(base, base_before, base_after, mut_bases, dict_cov):
    dict_cov[6] += float(1) / 3
    if is_at(base):
        if is_transition(base, mut_bases):
            dict_cov[4] += float(1) / 3
        else:
            dict_cov[5] += float(1) / 3
    elif is_cpg(base, base_before, base_after):  # This is a CG dinucleotide
        if is_transition(base, mut_bases):
            dict_cov[0] += float(1) / 3
        else:
            dict_cov[1] += float(1) / 3
    else:
        if is_transition(base, mut_bases):
            dict_cov[2] += float(1) / 3
        else:
            dict_cov[3] += float(1) / 3

    return dict_cov


# Transform all the elements in the input list into numeric.
def list_tonumeric(test_list):
    test_list = list(map(int, test_list))
    return test_list


# Group the positions in a list into position start and end pairs in a gene.
def list_to_pairs(list_old):
    new_list = []
    for i in range(0, len(list_old), 2):
        new_list.append(list_old[i:i + 2])
    return new_list


# This code is for process sequences on negative strand For the sequences in list, this code will get the reverse
# complement for each sequence and then reverse the sequence order.
def negative_process(list_sequence):
    new_list = []
    for sequences in list_sequence:
        sequences = sequences.reverse_complement()
        new_list.append(sequences)
    new_list = new_list[::-1]

    return new_list


# Find out if the base fall into non-coding region or coding region.
def get_region(region_list, base_pos, region_f):
    r_flag = "nc"
    j_flag = 0
    for regions in region_list:
        if regions[1] - 1 >= base_pos >= regions[0]:
            r_flag = region_f
            if r_flag == "coding":
                if base_pos == regions[1] - 1 or base_pos == regions[1] - 2 or base_pos == regions[0] or base_pos == \
                        regions[0] + 1:
                    j_flag = 1
            break

    return r_flag, j_flag


def get_annotation(annotation_file_path):
    name_dict = {}  # {transcript:[gene_name]}
    transcript_dict = {}  # {'gene_name':['ENST00000']}
    transcript_info_dict = {}  # {transcriptID:'exon'[]; 'CDS'[];'gene'[];'UTR

    with open(annotation_file_path, 'r') as annotation_f:
        for lines in annotation_f:
            if lines.startswith("chr"):
                line_split = lines.split('\t')
                chr_n, source, feature, start_pos, end_pos = line_split[0:5]
                strand, frame, attribute = line_split[6:9]
                position_list = [int(start_pos), int(end_pos)]

                # feature = line_split[2] # 'exon' / 'CDS'/ 'UTR'/ 'start_codon'/ 'stop_codon'/ 'transcript'/ 'gene'

                # if this is a protein coding gene and is the primary transcript:
                if 'protein_coding' in attribute and "KNOWN" in attribute and 'appris_principal' in attribute:
                    # parse the attributes
                    attribute_split = attribute.split(';')
                    gene_id_col, transcript_id_col, gene_type_col, gene_status_col, gene_name_col, transcript_type_col, \
                    transcript_status_col, transcript_name_col = attribute_split[0:8]
                    # not include the version number
                    gene_id = re.findall(r'(ENSG\d+|ENSGR\d+)', gene_id_col)[0]
                    transcript_id = re.findall(r'(ENST\d+|ENSTR\d+)', transcript_id_col)[0]
                    gene_status = re.findall(r'\"(.*?)\"', gene_status_col)[0]  # KNOWN or not
                    gene_name = re.findall(r'\"(.*?)\"', gene_name_col)[0]
                    transcript_status = re.findall(r'\"(.*?)\"', transcript_status_col)[0]  # KNOWN or not

                    if gene_name not in transcript_dict:
                        transcript_dict[gene_name] = []
                    if transcript_id not in transcript_dict[gene_name]:
                        transcript_dict[gene_name].append(transcript_id)

                    if transcript_id not in transcript_info_dict:
                        transcript_info_dict[transcript_id] = {}
                        transcript_info_dict[transcript_id]['strand'] = strand
                        transcript_info_dict[transcript_id]['chr'] = chr_n
                        transcript_info_dict[transcript_id]['exon'] = []
                        transcript_info_dict[transcript_id]['CDS'] = []
                        transcript_info_dict[transcript_id]['UTR'] = []
                        transcript_info_dict[transcript_id]['transcript'] = []

                        # if the transcript and genes are known, parse according to feature type
                    if feature == 'exon':
                        # exon_n = re.findall(r'\d+', attribute_split[8])
                        transcript_info_dict[transcript_id]['exon'].append(position_list)
                    elif feature == 'CDS':
                        transcript_info_dict[transcript_id]['CDS'].append(position_list)
                    elif feature == 'UTR':
                        transcript_info_dict[transcript_id]['UTR'].append(position_list)
                    elif feature == 'transcript':
                        transcript_info_dict[transcript_id]['transcript'] = position_list

                    if transcript_id not in name_dict:
                        name_dict[transcript_id] = gene_name

    # delete the transcript record that have 2 or more principle transcripts, only keep the longest transcript
    for names in transcript_dict:
        if len(transcript_dict[names]) > 1:
            store_max_t = 0
            for t in transcript_dict[names]:
                len_t = abs(transcript_info_dict[t]['transcript'][1] - transcript_info_dict[t]['transcript'][0])
                if len_t >= store_max_t:
                    store_max_t = len_t
                else:
                    del transcript_info_dict[t]

    return name_dict, transcript_info_dict


# get mrna positions in gene, then delete the utr positions
# the returned result is the cds positions
def get_mrna_position(transcript, info_dict, strand_gene):
    list_all_pos = info_dict[transcript]['UTR'] + info_dict[transcript]['CDS']
    list_all_pos.sort()

    # get the index of utr pair
    utr_index_list = []
    for pairs in info_dict[transcript]['UTR']:
        utr_idx = list_all_pos.index(pairs)
        utr_index_list.append(utr_idx)

    # if negative strand, the start position is the biggest position
    list_gene_position = []
    if strand_gene == '-':
        transcript_start = list_all_pos[-1][1]
        for pairs in list_all_pos:
            for positions in pairs:
                gene_position = -(positions - transcript_start)
                list_gene_position.append(gene_position)
    else:
        transcript_start = list_all_pos[0][0]
        for pairs in list_all_pos:
            for positions in pairs:
                gene_position = positions - transcript_start
                list_gene_position.append(gene_position)

    # put the mrna positions into pairs
    list_gene_position.sort()
    mrna_list = list_to_pairs(list_gene_position)
    # negative strand need to turn the utr index list around
    if strand_gene == '-':
        exon_n = len(mrna_list) - 1
        utr_index_list = [abs(i - exon_n) for i in utr_index_list]

    # remove utr positions according to the utr index
    for utr_index in sorted(utr_index_list,
                            reverse=True):  # reverse the list so that index in list won't change after deletion
        del mrna_list[utr_index]

    return mrna_list


# get the transcript sequence
def get_transcript_sequence(transcript, info_dict, fasta_dict, strand_gene):
    chromosome = info_dict[transcript]['chr']
    transcript_start_pos = info_dict[transcript]['transcript'][0]
    transcript_end_pos = info_dict[transcript]['transcript'][1]
    if strand_gene == '-':
        transcript_seq = fasta_dict[chromosome][transcript_start_pos - 1:transcript_end_pos].reverse_complement()
    else:
        transcript_seq = fasta_dict[chromosome][transcript_start_pos - 1:transcript_end_pos]

    return transcript_seq


# get the cdna sequence from the transcript sequence, inputs are cdna positions within the gene
def get_cdna_sequence(cds_list, strand_gene, sequence_transcript):
    list_cds_seq = []
    for cds_pairs in cds_list:
        cds_start = cds_pairs[0]
        cds_end = cds_pairs[1]
        if strand_gene == '-':
            list_cds_seq.append(sequence_transcript[cds_start:cds_end + 1])
        else:
            list_cds_seq.append(sequence_transcript[cds_start:cds_end + 1])

    concatenated_cds = Seq("", generic_dna)
    for s in list_cds_seq:
        concatenated_cds += s

    return concatenated_cds.seq


# this code loop through all the bases and get the coverage.
def calculate_coverage(transcript_seq, cds_seq, mrna_list, start_pos):
    # create dictionary to store position values, use genome positions as key
    cov_dict = {'nonsilent': {}, 'flank': {}, 'silent': {}}

    # Initialize the base value
    before_base = "N"
    after_base = "N"
    base_context = None
    context_position = None
    lp_flag = 0
    region_flag = 'nc'  # 'nc','coding'
    base_transcirpt_pos = 0  # position in transcript
    base_genome_pos = start_pos - 1  # position in gene, use genome position as the key

    n_context = 0

    for bases in transcript_seq:
        base_genome_pos += 1
        base_transcirpt_pos += 1
        for idx in cov_dict:
            if base_genome_pos not in cov_dict[idx]:
                cov_dict[idx][base_genome_pos] = [0, 0, 0, 0, 0, 0, 0]

        # skip the base if it's N, and give a low quality flag to the gene
        if bases == "N":
            lp_flag = 1
            continue

        # Determine before base and after base
        if base_transcirpt_pos == 1:
            after_base = transcript_seq[base_transcirpt_pos]
        elif base_transcirpt_pos == len(transcript_seq):
            before_base = transcript_seq[base_transcirpt_pos - 2]
        else:
            before_base = transcript_seq[base_transcirpt_pos - 2]
            after_base = transcript_seq[base_transcirpt_pos]

        # Determine regions and if the base is at junction position
        for regions in mrna_list:
            if regions[1] + 1 >= base_transcirpt_pos >= regions[0] + 1:
                region_flag = 'coding'
        if region_flag == "nc":
            other_nuc = filter(lambda i: i != bases, nuc)
            for mut_bases in other_nuc:
                cov_dict['flank'][base_genome_pos] = add_count(bases, before_base, after_base, mut_bases,
                                                               cov_dict['flank'][base_genome_pos])

        elif region_flag == "coding":
            n_context += 1
            if n_context <= 3:
                junction_flag = 1
            if n_context % 3 != 0:
                context_position = n_context % 3 - 1
                codon_n = (n_context // 3) * 3
                base_context = cds_seq[codon_n:codon_n + 3]
            else:
                context_position = 2
                codon_n = n_context
                base_context = cds_seq[codon_n - 3:codon_n]

            other_nuc = filter(lambda i: i != bases, nuc)
            for mut_bases in other_nuc:
                context_before = base_context  # This is the original tri-nucleotide context
                context_after = base_context[:context_position] + mut_bases + base_context[context_position + 1:]
                try:
                    if is_silent(context_before, context_after):
                        cov_dict['silent'][base_genome_pos] = add_count(bases, before_base, after_base, mut_bases,
                                                                        cov_dict['silent'][base_genome_pos])
                    else:
                        cov_dict['nonsilent'][base_genome_pos] = add_count(bases, before_base, after_base, mut_bases,
                                                                           cov_dict['nonsilent'][base_genome_pos])
                except KeyError:
                    lp_flag = 1
                    continue

        region_flag = "nc"

    return cov_dict, lp_flag


# this function parse the patient information file
def parse_info_file(info_file):
    donor_dict = {}
    pcawg_info = open(info_file, 'r')
    next(pcawg_info)
    for lines in pcawg_info:
        line_split = lines.split('\t')
        tumor_aliquot_id, normal_aliquot_id, donor_id, sample_id, specimen_id, project_code, gender, age \
            = line_split[0:8]
        if tumor_aliquot_id not in donor_dict:
            donor_dict[tumor_aliquot_id] = []
        donor_dict[tumor_aliquot_id] = donor_id

    return donor_dict


# this function get the file list in the coverage directory
def get_file_list(coverage_dir, donor_dict):
    patient_list = []
    file_list = [f for f in listdir(coverage_dir) if isfile(join(coverage_dir, f))]
    for file_names in file_list:
        aliquot_id = file_names.split('.')[0]
        try:
            id_donor = donor_dict[aliquot_id]
            patient_list.append(id_donor)
        except KeyError:
            pass

    return file_list, patient_list


# this function integrates bisect to find the leftmost value exactly equal to x. Could be only used for sorted list
def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return True
    else:
        return False


# this function is for appending the 0 positions in patient wig files then return a list
def get_zero_position(coverage_file):
    id_aliquot = coverage_file.split('.')[0].split('/')[-1]
    try:
        id_donor = dict_donor[id_aliquot]
        if id_donor not in dict_patient:
            dict_patient[id_donor] = {}
        with open(coverage_file, 'r') as wig_f:
            next(wig_f)
            for lines in wig_f:
                if lines.startswith('fixed'):
                    line_list = re.findall(r'\d+', lines)
                    chr_n = line_list[0]
                    start_pos = int(line_list[1])
                    position = start_pos - 1
                else:
                    position += 1
                if chr_n not in dict_patient[id_donor]:
                    dict_patient[id_donor][chr_n] = []
                if lines == '0\n':
                    dict_patient[id_donor][chr_n].append(position)
    except KeyError:
        pass

    return dict_patient


# this function is used when parallelzing the combination of transcripts and patient
def calculate_patient(position_dict_func, patient_dict, chr):
    patient_categ_list = [0, 0, 0, 0, 0, 0, 0]
    for positions, categs in position_dict_func.iteritems():
        try:
            if not index(patient_dict[chr], positions):
                for i in range(0, 7):
                    patient_categ_list[i] += categs[i]
        except KeyError:
            for i in range(0, 7):
                patient_categ_list[i] += categs[i]
    round_patient_list = [round(num) for num in patient_categ_list]

    return round_patient_list


# @contextmanager
# def poolcontext(*args, **kwargs):
#     pool = multiprocessing.Pool(*args, **kwargs)
#     yield pool
#     pool.terminate()


# this function include the above functions to calculate coverage and write into a csv file
# every combination of transcript and patient will be passed as arguments. The process is parallelized
def coverage_calculation(params):  # params:[transcript, patient]
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
        coverage_patient = calculate_patient(dict_position_cov[keys], dict_patient[params[1]], chromosome)

        with open(file_out, 'a') as f_cov:
            f_cov.write(name_gene + '\t' + keys + '\t' + params[1] + '\t' + str(coverage_patient[0]) + '\t' + str(
                coverage_patient[1]) + '\t' + \
                        str(coverage_patient[2]) + '\t' + str(coverage_patient[3]) + '\t' + str(
                coverage_patient[4]) + '\t' + str(coverage_patient[5]) \
                        + '\t' + str(coverage_patient[6]) + '\n')
            # writer = csv.writer(f_cov, delimiter='\t')
            # writer.writerows(zip(*list_print))


if __name__ == '__main__':
    arguments = parse_arguments()
    file_anno = arguments.inputPath1
    file_fasta = arguments.inputPath2
    file_info = arguments.inputPath3
    file_coverage = arguments.inputPath4
    file_out = arguments.outputPath

    # first get the donor id according to the information files
    dict_donor = parse_info_file(file_info)

    # open the file and write the first few columns of the header line, this also prevents the file already exisits
    f_cov = open(file_out, "wb")
    f_cov.write('gene_name' + '\t' + 'zone' + '\t' + 'patient' + '\t' + 'Categ1'+ '\t' + 'Categ2'+'\t' + 'Categ3'+'\t' + \
                'Categ4'+'\t' + 'Categ5'+'\t' + 'Categ6'+'\t' + 'Categ7'+'\n')
    f_cov.close()

    print "START READING PATIENT COVERAGE"
    start1 = time.time()
    if file_coverage:
        dict_patient = {}
        dict_patient = get_zero_position(file_coverage)

        list_patient = dict_patient.keys()
    end1 = time.time()
    print "FINISH READING PATIENT COVERAGE: " + '\t' + str(end1 - start1)

    # initialize coverage dictionary
    zone_list = ["s", "ns", "nc"]
    categ_list = list(range(1, 8))
    # print "***START CALCULATING COVERAGE: " + str(name_transcripts)
    start1 = time.time()

    # parse annotaiton file
    # dict_name: {gene_name: transcript ID}
    # dict_transcript_info{transcript ID: UTR, CDS, exon, transcript, chr, strand}
    print "***START PARSING ANNOTATION***"
    dict_name, dict_transcript_info = get_annotation(file_anno)
    print "***FINISH PARSING ANNOTATION***"

    # store genome fasta into dictionary
    dict_record = SeqIO.to_dict(SeqIO.parse(file_fasta, 'fasta'))  # the keys are ['chrN']

    start1 = time.time()
    list_transcript = dict_transcript_info.keys()
    paramlist = list(itertools.product(list_transcript, list_patient))
    p = Pool(2)
    p.map(coverage_calculation, paramlist)
    p.close()
    p.join()
    end1 = time.time()
    print "TIME FOR ALL TRNASCRIPTS: " + str(end1 - start1)
