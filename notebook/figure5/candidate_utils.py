import pandas as pd
import os, pickle
from multiprocessing import Pool
import numpy as np
from collections import Counter
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns


#-------------Read Data------------------

# Read the significant gene dataframe after FDR calculation
dir_out = '../figure4/'
feature_type = 'histology';syn_nsyn = 'syn';run = 'cohort_072221';threshold = 1;
df_syn = pd.read_csv(os.path.join(dir_out,feature_type+'.syn.df_all_forheatmap.'+run+'.'+str(threshold)+'.csv'),index_col = 0)
df_syn = df_syn.set_index('gene')

# Maf file
dir_maf = '../../maf_out/maf_cohorts_060121'

# The sample information file
dir_refs = '../../anno_ref/proc_refs'
df_info = pd.read_csv(os.path.join(dir_refs, 'PCAWG_sample_info.txt'), sep = '\t')

# The copynumber file
df_cnv = pd.read_csv('../../anno_ref/ICGC/consensus_cnv/consensus_CN.by_gene.170214.txt', sep = '\t')
df_cnv = df_cnv.set_index(['Gene Symbol', 'Locus ID', 'Cytoband'])

# Read expression data
exp_dir = '../../anno_ref/ICGC/pcawg_rnaseq/'
gene_tophatuq = 'tophat_star_fpkm_uq.v2_aliquot_gl.tsv'
gene_tophat = 'tophat_star_fpkm.v2.aliquot_gl.tsv'
# Read aliquot id information
df_exp_info = pd.read_csv(os.path.join(exp_dir,'rnaseq.metadata.tsv'), sep = '\t')
# Read expression information
# df_exp_uq = pd.read_csv(join(exp_dir,gene_tophatuq),sep = '\t', index_col = 0)
df_exp = pd.read_csv(os.path.join(exp_dir,gene_tophat),sep = '\t', index_col = 0)
# Read the significant gene-id dataframe
fsig_name = 'sig_gene_name_id.csv'
df_nsig = pd.read_csv(os.path.join(dir_out, fsig_name))
df_nsig.columns = ['name', 'id']

# Columns needed in the information dataframe
num_cols = ['SV.events','Coding.SNVs', 'Non.coding.SNVs','CNA.events..do_not_use.', 'Retrotransposon.insertions',\
            'Mitochondrial.mutations','all.Indels']
str_cols = ['dcc_specimen_type','histology_tier3','histology_tier4',
        'tumour_histological_type',\
       'tumour_stage', 'tumour_histological_comment', 'specimen_donor_treatment_type']

#-------------Syn Patient sub histology and mutation info------------------
def get_num_str(gene):
    global df_syn, dir_maf
    histology = df_syn.loc[gene,'feature'] #histology type
    df_maf = pd.read_csv(os.path.join(dir_maf,feature_type, histology+'.csv'), sep = '\t')
    df_maf = df_maf.set_index(['Hugo_Symbol','Variant_Classification','Donor_ID' ])
    df_gene_mut = df_maf.loc[pd.IndexSlice[gene, 'Silent',:],:]### All patients' specific histology in this histology type
    ldonor = df_gene_mut.index.get_level_values('Donor_ID').unique().tolist()
    
    df_histology = df_info[df_info['histology_abbreviation'] == histology]
    df_histology_syn = df_histology[df_histology['icgc_donor_id'].isin(ldonor)]

    ### The tumor patient mutation information for synonymous patients and all patients
    df_num_summ = pd.concat([pd.concat({'syn': df_histology_syn[num_cols].describe()}, names=['mut_status'])\
    ,pd.concat({'all': df_histology[num_cols].describe()}, names=['mut_status'])])

    ### The tumor histology subtype information
    df_str_info = pd.DataFrame(columns = str_cols, index = ['syn', 'all'])
    for columns in str_cols:
        if columns not in ['icgc_donor_id','histology_abbreviation','histology_tier2']:
            df_str_info.loc['all', columns] = str(dict(Counter(df_histology[columns].tolist())))
            df_str_info.loc['syn', columns] = str(dict(Counter(df_histology_syn[columns].tolist())))
    return df_str_info, df_num_summ

#-------------Syn Patient Mutations------------------
def get_mut(gene):
    global df_syn
    histology = df_syn.loc[gene,'feature'] #histology type
    df_maf = pd.read_csv(os.path.join(dir_maf,feature_type, histology+'.csv'), sep = '\t')
    df_maf = df_maf.set_index(['Hugo_Symbol','Variant_Classification','Donor_ID' ])
    ldonor_all = df_maf.index.get_level_values('Donor_ID').unique().tolist()
    print(f'Total {len(ldonor_all)} patients')
    df_gene_mut = df_maf.loc[pd.IndexSlice[gene, 'Silent',:],:]### All patients' specific histology in this histology type
    ldonor = df_gene_mut.index.get_level_values('Donor_ID').unique().tolist()
    print(f'{len(ldonor)} patients have synonymous mutations, they are: {ldonor}')
    for muts in df_gene_mut['Genome_Change'].tolist():
        chromosome = muts.split(':')[0].split('.')[-1].split('chr')[-1]
        change = muts.split(':')[1]
        mutation = chromosome+':g.'+change
        print(mutation)

    return df_maf, df_gene_mut, ldonor

#-------------Syn Patient copy number------------------
def get_cnv(ldonor, gene):
    for donors in ldonor:
        idx = df_info[df_info['icgc_donor_id'] == donors].index
        aliquot_id = df_info.loc[idx, 'tumour_specimen_aliquot_id']
        cnv = df_cnv.loc[pd.IndexSlice[gene,:,:],aliquot_id]
        print(cnv)
    return cnv

#-------------Syn Patient Expression------------------
def get_gene_exp(gene_name, df_expression):
    global df_nsig
    idx = df_nsig[df_nsig['name'] == gene_name].index
    gene_id = df_nsig.loc[idx, 'id'].values[0]
    df = df_expression.loc[df_expression.index.str.contains(rf'{gene_id}'),:]
    return df

def get_syn_mut(gene_name):
    global df_syn
    histology = df_syn.loc[gene_name,'feature']
    df_maf = pd.read_csv(os.path.join(dir_maf,feature_type, histology+'.csv'), sep = '\t')
    df_maf = df_maf.set_index(['Hugo_Symbol','Variant_Classification','Donor_ID' ])
    df_silent = df_maf.loc[pd.IndexSlice[gene_name, 'Silent',:],:]
    patient = df_silent.index.get_level_values('Donor_ID').unique().tolist()
    
    return df_silent, patient

def get_patient_id(gene_name, patients):
    global df_exp_info, df_syn
    histology = df_syn.loc[gene_name,'feature']
    df = df_exp_info[df_exp_info['histology_abbreviation'] == histology]
    
    # Get tumor, syn and normal patient aliquot id
    normal_id = df[df['tumor.normal'] == 'normal']['aliquot_id']
    df_tumor = df[df['tumor.normal'] == 'tumor']
    tumor_syn_id = df_tumor[df_tumor['icgc_donor_id'].isin(patients)]['aliquot_id']
    tumor_other_id = df_tumor[~df_tumor['icgc_donor_id'].isin(patients)]['aliquot_id']
    
    return normal_id, tumor_syn_id, tumor_other_id

def get_patient_exp(ids, df_gene_exp, tissue_type = None):
    df_for_test = df_gene_exp[ids].transpose()
    df = df_gene_exp[ids].transpose().reset_index()
    df['tumor.normal'] = tissue_type
    df.columns = ['id','exp','tumor.normal']
    
    return df, df_for_test

def get_expression(gene):
    histology = df_syn.loc[gene,'feature']

    df_exp_gene = get_gene_exp(gene, df_exp)
    df_synmut, synp = get_syn_mut(gene)
    id_normal, id_syn, id_other =get_patient_id(gene, synp)

    df_normal, normal_test = get_patient_exp(id_normal, df_exp_gene, 'normal')
    df_tsyn, syn_test = get_patient_exp(id_syn, df_exp_gene,'tumor_syn')
    df_tother, other_test = get_patient_exp(id_other, df_exp_gene, 'tumor_other')
    df_all = pd.concat([df_normal,df_tsyn,df_tother])

    nnorm = len(df_normal['id'].unique())
    print(f'Number of normal patient: {nnorm}')
    nsyn = len(df_tsyn['id'].unique())
    print(f'Number of synonymous patients: {nsyn}')
    nother = len(df_tother['id'].unique())
    print(f'Number of other tumor patients:{nother}')

    fig,ax = plt.subplots(figsize=(5,3))
    ax = sns.boxplot(x = 'tumor.normal', y = 'exp', data = df_all)

    mannwhitneyu(syn_test, other_test)
    text = f'test statistic:{round(mannwhitneyu(syn_test, other_test)[0],4)},p-value:{round(mannwhitneyu(syn_test, other_test)[1],4)}'
    print(text)
    ax.set_title(histology+'_'+gene+'_'+text)
    # plt.savefig('./res/'+organ_type+'_'+gene_name+'.png')
    plt.show()

###-------------------CERES------------------------
def get_lineage(gene, histology):
    ### Read file
    dir_depmap = './depmap'
    df_depmap = pd.read_csv(os.path.join(dir_depmap, gene+'_21Q2.csv'))
    print(df_depmap[df_depmap['Lineage'] == histology]['Lineage Subtype'].unique())
    
def get_CERES(gene, lineage):
    dir_depmap = './depmap'
    df_depmap = pd.read_csv(os.path.join(dir_depmap, gene+'_21Q2.csv'))
    df_depmap.columns = ['ID', 'CERES', 'Name', 'Primary Disease', 'Lineage', 'Lineage Subtype', 'Expression', 'Mutation']

    df_gene = df_depmap[df_depmap['Lineage Subtype'] == lineage]
    print(f'{len(df_gene)} in histology subtype')
    df_other = df_depmap[df_depmap['Lineage Subtype'] != lineage]
    print(f'{len(df_other)} not in histology subtype')
    df_all = pd.DataFrame({'inlineage':df_gene['CERES'],
                          'other':df_other['CERES']})
    df_all = df_all.melt()

    fig,ax = plt.subplots(figsize=(5,3))
    ax = sns.boxplot(x = 'variable', y = 'value', data = df_all)

    mannwhitneyu(df_gene['CERES'], df_other['CERES'])
    text = f'test statistic:{round(mannwhitneyu(df_gene["CERES"], df_other["CERES"])[0],4)},\
    p-value:{round(mannwhitneyu(df_gene["CERES"], df_other["CERES"])[1],4)}'
    print(text)
    ax.set_title('CERES:'+lineage+'_'+gene+'_'+text)
    
    # Expression
    df_all_exp = pd.DataFrame({'inlineage':df_gene['Expression'],
                      'other':df_other['Expression']})
    df_all_exp = df_all_exp.melt()
    fig,ax = plt.subplots(figsize=(5,3))
    ax = sns.boxplot(x = 'variable', y = 'value', data = df_all_exp)

    mannwhitneyu(df_gene['Expression'], df_other['Expression'])
    text = f'test statistic:{round(mannwhitneyu(df_gene["Expression"], df_other["Expression"])[0],4)},\
    p-value:{round(mannwhitneyu(df_gene["Expression"], df_other["Expression"])[1],4)}'
    print(text)
    ax.set_title('DepMap Expression:'+lineage+'_'+gene+'_'+text)