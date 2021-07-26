import pandas as pd
import os
from tqdm import tqdm
import pickle

### File path and read directory
dir_cohorts = '../anno_ref/cohorts/'
dir_cov_ind = '../cov_out/cov_ind/'
dir_out = '../cov_out/cov_cohorts_060121'

l_carcinoma = pickle.load(open('../cov_out/cov_cohorts_060121/origin/carcinoma.pkl','rb'))
df_cov_merged = pd.concat(l_carcinoma, axis = 1)
df_cov_merged.to_csv('../cov_out/cov_cohorts_060121/origin/carcinoma.csv',sep='\t')