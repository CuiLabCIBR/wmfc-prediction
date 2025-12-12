#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os

sublist_file = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/PNC/batchcode/sublist_new.txt" 
base_path = "/ibmgpfs/cuizaixu_lab/xuhaoshu/WM_prediction/datasets/PNC/fmriprep"  
output_base_path = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/PNC/results/custom_confounds"
# run_ids = ["rest_run-1", "rest_run-2", "rest_run-3", "rest_run-4"] #EFNY
run_ids = ["rest"] #PNC
# run_ids = ["rest_run-01","rest_run-02"] #CCNP

if not os.path.exists(sublist_file):
    print(f"error: could not find file {sublist_file}")
    exit(1)


with open(sublist_file, 'r') as f:
    subjects = [line.strip() for line in f if line.strip()]

sub_num = len(subjects)
print(f"read {sub_num} subjects.")


if not os.path.exists(output_base_path):
    os.makedirs(output_base_path)


columns_to_extract = [
    'csf', 'csf_derivative1', 'csf_derivative1_power2', 'csf_power2',
    'global_signal', 'global_signal_derivative1', 'global_signal_derivative1_power2', 'global_signal_power2']


for subid in subjects:
    sub_output_path = os.path.join(output_base_path)
    if not os.path.exists(sub_output_path):
        os.makedirs(sub_output_path)

    for run in run_ids:
        # tsv_file = os.path.join(base_path, subid, "func", f"{subid}_task-{run}_desc-confounds_timeseries.tsv") #EFNY
        # tsv_file = os.path.join(base_path, subid, "ses-01/func", f"{subid}_ses-01_task-{run}_desc-confounds_timeseries.tsv")   #CCNP
        tsv_file = os.path.join(base_path, subid, "ses-PNC1/func", f"{subid}_ses-PNC1_task-{run}_acq-singleband_desc-confounds_timeseries.tsv")   #pnc
        # tsv_file = os.path.join(base_path, subid, "ses-PNC1/func", f"{subid}_ses-PNC1_task-{run}_acq-100_desc-confounds_timeseries.tsv")   #pnc

        if not os.path.exists(tsv_file):
            print(f"warning: could not find file {tsv_file}，skipping...")
            continue
        
        # read TSV
        df = pd.read_csv(tsv_file, sep='\t')
        
        # extract 
        extracted_df = df[columns_to_extract].copy()
        extracted_df.fillna(0, inplace=True)
        

        # output_file = os.path.join(sub_output_path, f"{subid}_task-{run}_desc-confounds_timeseries.tsv")   #EFNY
        # output_file = os.path.join(sub_output_path, f"{subid}_ses-01_task-{run}_desc-confounds_timeseries.tsv")  #CCNP
        output_file = os.path.join(sub_output_path, f"{subid}_ses-PNC1_task-{run}_acq-singleband_desc-confounds_timeseries.tsv")  #pnc
        # output_file = os.path.join(sub_output_path, f"{subid}_ses-PNC1_task-{run}_acq-100_desc-confounds_timeseries.tsv")  #pnc
        extracted_df.to_csv(output_file, sep='\t', index=False)
        print(f"save to {output_file}")
