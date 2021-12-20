from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import re,time,os,pickle
from tqdm import tqdm
from multiprocessing import Pool
import multiprocessing as mp
from functools import partial


### reference table
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
maf_cols = pickle.load(open('../data/proc_refs/maf_cols.pkl', 'rb')) #maf columns 

### This function split the splitted raw maf files into individual ones
def split_patient(fsplit, dir_maf = '../maf_raw/maf_split', dir_out = '../data/maf/intermediate/individual/'):
    df_maf = pd.read_csv(os.path.join(dir_maf,fsplit), sep = '\t', header = None)
    lpatient = df_maf.iloc[:,12].unique()
    for patient in lpatient:
        outpath = os.path.join(dir_out,patient+'.to_merge.csv')
        df_patient = df_maf[df_maf.iloc[:,12] == patient]
        if os.path.exists(outpath):
            df_patient.to_csv(outpath, mode ='a', header = None)
        else:
            df_patient.to_csv(outpath)
            
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
#         print("NEITHER TRANSITION NOR TRANSVERSION MUTATION")
        pass

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

def allele_categ_assign(idx,df):
    genome_change = df.loc[idx,'Genome_Change']
    ref_contxt = df.loc[idx,'ref_context']
    if '>' in genome_change:
        alleles = re.findall(r'[A-Z]+>[A-Z]+', genome_change)[0].split('>')
        ref_allele = alleles[0]
        alter_allele = alleles[1]
        mid_pos = int((len(ref_contxt)+1)/2)
        bbase = ref_contxt[mid_pos-1]
        abase = ref_contxt[mid_pos+1]
        tri_nuc = bbase+ref_allele+abase
        categ = validate_categ(ref_allele, alter_allele, tri_nucleotide=tri_nuc)
    else:
        categ = 7
    df.loc[idx,'categ'] = int(categ)


### This function assing categ to patients
def categ_assign(patient, dir_ind = '../data/maf/intermediate/individual/split',dir_categ_out = \
                 '../data/maf/intermediate/individual/categ'):
    
    # check if file exists, not run if exists
    outf = os.path.join(dir_categ_out, patient+'.to_merge.categ.csv')
    outfzip = os.path.join(dir_categ_out, patient+'.to_merge.categ.csv.gz')
    # if os.path.exists(outf) or os.path.exists(outfzip):
    #     return
    # print(f'start assigning categ for {outf}...')
    
    df_maf = pd.read_csv(os.path.join(dir_ind,patient+'.to_merge.csv'), index_col = 0)
    df_maf = df_maf.reset_index(drop = True)
    df_maf.columns = maf_cols
    
    # Multiprocess the mutations
    for idx in df_maf.index:
        genome_change = df_maf.loc[idx,'Genome_Change']
        ref_contxt = df_maf.loc[idx,'ref_context']
        if '>' in genome_change:
            alleles = re.findall(r'[A-Z]+>[A-Z]+', genome_change)[0].split('>')
            ref_allele = alleles[0]
            alter_allele = alleles[1]
            mid_pos = int((len(ref_contxt)+1)/2)
            bbase = ref_contxt[mid_pos-1]
            abase = ref_contxt[mid_pos+1]
            tri_nuc = bbase+ref_allele+abase
            categ = validate_categ(ref_allele, alter_allele, tri_nucleotide=tri_nuc)
        else:
            categ = 7
        df_maf.loc[idx,'categ'] = int(categ)

    print(f'saving {patient} maf file...')
    # df_maf.to_csv(os.path.join(dir_categ_out, patient+'.to_merge.categ.csv'))
    df_maf.to_csv(os.path.join(dir_categ_out, patient+'.to_merge.categ.csv.gz')\
         ,chunksize=100000,compression='gzip',encoding='utf-8')
    print(f'finish saving {patient} maf file...')

