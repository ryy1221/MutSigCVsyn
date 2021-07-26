#!/usr/bin/env python
import subprocess
import pandas as pd
import os
# Get feat
dir_cohorts = '../../anno_ref/cohorts/'
feature_type = 'histology'

df = pd.read_csv(os.path.join(dir_cohorts,feature_type+'.csv'))
lfeat = df[feature_type].unique()

for feat in lfeat:
    qsub_command =f'qsub -v HISTOLOGY_NAME={feat} -A epo2_f_g_sc_default run_histology_syn.pbs'
    print(qsub_command) # Uncomment this line when testing to view the qsub command

# Comment the following 3 lines when testing to prevent jobs from being submitted
#     exit_status = subprocess.call(qsub_command, shell=True)
#     if exit_status is 1:  # Check to make sure the job submitted
#         print "Job {0} failed to submit".format(qsub_command)
        
print("Done submitting jobs!")