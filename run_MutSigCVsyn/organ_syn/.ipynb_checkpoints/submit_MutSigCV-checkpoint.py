#!/usr/bin/env python
import subprocess
import pandas as pd
import os
# Get feat
dir_cohorts = '../../anno_ref/cohorts/'
feature_type = 'organ'

df = pd.read_csv(os.path.join(dir_cohorts,feature_type+'.csv'))
lfeat = df[feature_type].unique()

for feat in lfeat:
    qsub_command =f'qsub -v FEAT_NAME={feat} -A wff3_a_g_hc_default run_organ_syn.pbs'
    print(qsub_command) # Uncomment this line when testing to view the qsub command

# Comment the following 3 lines when testing to prevent jobs from being submitted
#     exit_status = subprocess.call(qsub_command, shell=True)
#     if exit_status is 1:  # Check to make sure the job submitted
#         print(f'Job {qsub_command} failed to submit')
        
print("Done submitting jobs!")