# This script is for assigning the mutation category to mutations in MAF file
# Author: Yiyun
# Date: Jul-21-2020

from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import re
import time
import os
from os import listdir

# Global variable nucleotides and condon table
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
stop_codon = ["TAA", "TAG", "TGA"]
null_list = ["Nonstop_Mutation", "Nonsense_Mutation", "Splice_Site", "DEL", "INS"]


# Use the optional parse function for input files.
def parse_arguments():
    parser = OptionParser("Usage: %prog -i <inputpath>")
    parser.add_option("-m", "--mutation_dir",
                      dest="inputPath1",
                      help="The full absolute path of the mutation file")
    parser.add_option("-f", "--fasta_file",
                      dest="inputPath2",
                      help="The full absolute path of the GENCODEv19 genome fasta file")

    (options, arg) = parser.parse_args()
    if not options.inputPath1 or not options.inputPath2:
        parser.error("Requires mutation file and fasta file. Run with --help for help.")
        quit(1)
    return options


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


# this function determines transition mutation
def is_transition(ori_allele, alt_allele):
    p = False
    transition = [("A", "G"), ("G", "A"), ("T", "C"), ("C", "T")]
    transversion = [("C", "A"), ("C", "G"), ("G", "T"), ("T", "G"), ("G", "C"), ("A", "C"), ("A", "T"), ("T", "A")]
    if (ori_allele, alt_allele) in transition:
        p = True
    elif (ori_allele, alt_allele) in transversion:
        p = False
    else:
        raise CategError("NEITHER TRANSITION NOR TRANSVERSION MUTATION")

    return p


# this function calculate the contribution of the coverage of a certain base in CpG or not to different categs
def validate_categ(original_allele, altered_allele,
                   tri_nucleotide):  # This before and after base is not same with tri-nucleotide context. It's not in codon context.
    base_before = tri_nucleotide[0]
    base_after = tri_nucleotide[2]

    # Initialize the categ fractions
    categ = 0
    # If a C:G basepairs mutation
    # If mutation in CpG dinucleotide
    if original_allele == "A" or original_allele == "T":
        if is_transition(original_allele, altered_allele):
            categ = 5
        else:
            categ = 6
    elif (original_allele + base_after == gc_di_nuc) or (
            base_before + original_allele == gc_di_nuc):  # This is a CG dinucleotide
        if is_transition(original_allele, altered_allele):
            categ = 1
        else:
            categ = 2
    else:
        if is_transition(original_allele, altered_allele):
            categ = 3
        else:
            categ = 4

    return categ


# this function create a DNA instances
def into_dna(sequence):
    dna_sequence = Seq(sequence, generic_dna)

    return dna_sequence


# this function get all the start position in the mutation file and return a dictionary
def get_mutation_position(maf_file):
    mut_dict = {}
    maf_f = open(maf_file, "r")
    next(maf_f)
    for lines in maf_f:
        line_split = lines.split("\t")
        chr_n = line_split[1]
        pos_start = int(line_split[2])

        # Creat dictionary:
        if chr_n not in mut_dict:
            mut_dict[chr_n] = {}
        # first only validate point mutation
        if pos_start not in mut_dict[chr_n]:
            mut_dict[chr_n][
                pos_start] = []  # This empty list will store the reference base, the base before and after the
            # reference base
    maf_f.close()

    return mut_dict


# this function use the mutation dictionary as input and get the reference base as well as base before and after
def get_base(mutation_dict):
    # open the fasta file
    for seq_record in dict_record:
        for keys in mutation_dict:
            if seq_record == 'chr' + keys:
                dict_chr_sorted = sorted(mutation_dict[keys])  # sort the mutation dictionary for the mapped chromosome
                for chromosome_position in dict_chr_sorted:  # for the positions, find the base and before and after base
                    base = dict_record[seq_record][chromosome_position - 1]
                    b_base = dict_record[seq_record][chromosome_position - 2]
                    a_base = dict_record[seq_record][chromosome_position]
                    mutation_dict[keys][chromosome_position] = [base, b_base, a_base]

    return mutation_dict


# this function identify and assgin categ
def assign_categ(maf_file, base_dict):
    # open the ouput file
    output_maf = open(maf_file + '.output.txt', 'wb')
    maf = open(maf_file, "r")

    header = '{}\t{}\n'.format(next(maf).rstrip(), 'categ')
    output_maf.write(header)

    for lines in maf:
        line_split = lines.split("\t")
        gene_name = line_split[0]
        chr_n, pos_start, pos_end, strand, variant, classification = line_split[1:7]
        genome_change = line_split[14]
        pos_start = int(pos_start)
        if variant in null_list or classification in null_list:
            predicted_categ = 7

        else:  # for SNP, DNP, TNP, the mutations are defined according to the first base
            # print genome_change
            alleles = re.findall(r'[A-Z]+>[A-Z]+', genome_change)[0].split('>')
            ref_allele = alleles[0]
            alter_allele = alleles[1]
            # print genome_change + '\t' + ref_allele + '\t' + alter_allele

            # Get reference base
            [genome_base, base_b, base_a] = base_dict[chr_n][pos_start]
            base_b = into_dna(base_b);
            genome_base = into_dna(genome_base);
            base_a = into_dna(base_a)
            tri_nuc = base_b + genome_base + base_a
            ori_allele = ref_allele[0]
            mut_allele = alter_allele[0]

            # For those mutations in exon, check if the reference allele is the same with the mutation allele
            # if cdna_change and re.findall(r'[A-Z]', cdna_change):
            #     ori_allele = re.findall(r'[A-Z]', cdna_change)[0]

            if not str(genome_base) == ori_allele:
                tri_nuc = tri_nuc.reverse_complement()
                ori_allele = into_dna(ref_allele).complement()[0]
                mut_allele = into_dna(alter_allele).complement()[0]
            predicted_categ = validate_categ(ori_allele, mut_allele, tri_nucleotide=tri_nuc)

        # write the cate to the file
        new_line = '{}\t{}\n'.format(lines.rstrip(), str(predicted_categ))
        output_maf.write(new_line)

        # if int(predicted_categ) != int(categ):
        #     print 'n'
    maf.close()
    output_maf.close()


if __name__ == '__main__':
    arguments = parse_arguments()
    dir_mut = arguments.inputPath1
    file_fasta = arguments.inputPath2
    # second, store the genome file into a dictionary
    dict_record = SeqIO.to_dict(SeqIO.parse(file_fasta, 'fasta'))  # the keys are ['chrN']
    print "***FINISH READING FASTA FILE"

    list_f = listdir(dir_mut)
    for f in list_f:
        file_mut = os.path.join(dir_mut,f)
        # first, get the mutation position
        start1 = time.time()
        dict_mut = get_mutation_position(file_mut)
        end1 = time.time()
        print "***FINISH GETTING POSITION: "+ str(end1-start1)
        # third, according to the positions find before and after bases
        start2 = time.time()
        dict_mut = get_base(dict_mut)
        end2 = time.time()
        print "***FINISH GETTING GENOME BASE: "+ str(end2-start2)
        # finally, assign categ and write to new mutation file
        start3 = time.time()
        assign_categ(file_mut, dict_mut)
        end3 = time.time()
        print "***FINISH ASSIGNING CATEG: "+ str(end3-start3)
