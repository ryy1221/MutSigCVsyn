import pandas as pd
import os,pickle
from tqdm import tqdm
from utils import *
from multiprocessing import Pool
from datetime import date

today = date.today()

dir_cohorts = '../anno_ref/cohorts'
dir_maf = '../maf_out/maf_cohorts_060121'
feature_type = 'histology'
df = pd.read_csv(os.path.join(dir_cohorts,feature_type+'.csv'))
lfeat = df[feature_type].unique()

### Load annotation data
dict_transcript_info = pickle.load(open('../proc_09152020/dict_transcript_info_062121.pkl','rb'))
dict_record = pickle.load(open('../proc_09152020/dict_record_new.pkl','rb'))

### Parameter card
chromosome = '18'
feat = 'Skin-Melanoma'
start_range = 27; end_range=77

interval = 1000000
### Create breakpoints for chr1
l_bk = []
breakpoints = list(range(0, len(dict_record['chr'+chromosome]), interval))
for i in range(0, len(breakpoints)-1):
    group = [breakpoints[i], breakpoints[i+1]]
    l_bk.append(group)


### Define function to calculate
def get_chr_mut(i):
    global chromosome, feat
    
    df = pd.DataFrame(columns = ['Donor','chr-pos','n_syn','n_nsyn','n_coding','n_all'])
    df_maf = pd.read_csv(os.path.join(dir_maf,feature_type, feat+'.csv'), sep = '\t')
    ldonor = df_maf['Donor_ID'].unique().tolist()

    for donors in tqdm(ldonor):
#             # Initialize a donor dataframe
        df_donor = df_maf[df_maf['Donor_ID'] == donors]
        df_donor = df_donor[df_donor['Chromosome'] == chromosome]
        df_donor = df_donor[df_donor['Variant_Type'] == 'SNP']
        df_donor = df_donor[(df_donor['Start_position'] < l_bk[i][1]) & \
                            (df_donor['Start_position'] > l_bk[i][0])]
        n_all = len(df_donor)
        n_syn = len(df_donor[df_donor['Variant_Classification'] == 'Silent'])
        n_nsyn = len(df_donor[df_donor['Variant_Classification'].isin(['Missense_Mutation', 'Nonsense_Mutation','Splice_Site', 'Nonstop_Mutation','Start_Codon_SNP'])])
        l_donor = [donors, '1-'+str(i), n_syn, n_nsyn, n_syn+n_nsyn, n_all]
        df_temp = pd.DataFrame([l_donor], columns = ['Donor','chr-pos','n_syn','n_nsyn','n_coding','n_all'])
        df = pd.concat([df,df_temp], ignore_index = True)
        
    return df

### Run get mutation, we don't get the full range
p = Pool(15)
res = p.map(get_chr_mut, range(start_range, end_range))
p.close()
p.join()

### Concat the results
df_all = pd.concat(res, ignore_index = True)

### Save results
fname = feat+'_chr'+chromosome+'_'+str(interval)+'_'+str(start_range)+'-'+str(end_range)+'_'+today.strftime("%b-%d-%Y")+'.csv'
df_all.to_csv(os.path.join('./figure1',fname))