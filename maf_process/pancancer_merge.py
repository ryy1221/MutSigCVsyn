### This script merge the individual coverage files into the cohorts
### Author: Yiyun
import pandas as pd
import os
from tqdm import tqdm
import pickle

### File path and read directory
dir_cohorts = '../anno_ref/cohorts/'
dir_maf_ind = '../maf_out/maf_ind/'
dir_out = '../maf_out/maf_cohorts_060121'

feature_type = 'pancancer'
df = pd.read_csv(os.path.join(dir_cohorts,feature_type+'.csv'))
lfeat = df[feature_type].unique()
print(lfeat)
# merge patients
def merge_maf(feature):
    global df, feature_type
    ldf = []
    # Initialize dataframe and output path
    path_out = os.path.join(dir_out,feature_type)
    
    if os.path.exists(os.path.join(path_out,feature+'.csv')):
        print(f'exists{feature}')
        return
    
    # Get patient list for histology
    df_filtered = df[df[feature_type] == feature]
    lp = df_filtered['tumor_aliquot_id'].unique()

    #read patient file
    for p in tqdm(lp):
        df_maf_ind = pd.read_csv(os.path.join(dir_maf_ind,p+'.to_merge.categ.csv'), index_col = 0)
        ldf.append(df_maf_ind)
#     pickle.dump(ldf, open(os.path.join(path_out,feature+'.pkl'), 'wb'))  
#     print(f'finish dump {feature}')
    df_maf_merged = pd.concat(ldf, axis = 0)
    df_maf_merged.to_csv(os.path.join(path_out,feature+'.csv'),sep='\t', index = False)
    print(f'finish{feature}')

merge_maf('pancancer')