#coding: utf-8
import scipy.io as sio
import numpy as np
import pandas as pd
import os
import sys
sys.path.append('/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/code/4th_prediction/s02_prediction/FC_nosmooth');
import PLSr1_CZ_Random_RegressCovariates

targetStr = 'age'
outFolder = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/code/4th_prediction/s02_prediction/FC_nosmooth/age_withHaufe_Schaefer100_522/brainproject522_'+ targetStr

# Import data
# 1. atlas loading
GG_datapath = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/FC_merge_nosmooth/total_FCvector_GG_n522.txt';
GW_datapath = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/FC_merge_nosmooth/total_FCvector_GW_60p_n522.txt';
WW_datapath = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/FC_merge_nosmooth/total_FCvector_WW_60p_n522.txt';

GG_data_files_all = np.genfromtxt(GG_datapath, filling_values=0.0)
GG_data_files_all = np.float32(GG_data_files_all) 
  
GW_data_files_all = np.genfromtxt(GW_datapath, filling_values=0.0)
GW_data_files_all = np.float32(GW_data_files_all)
  
WW_data_files_all = np.genfromtxt(WW_datapath, filling_values=0.0)
WW_data_files_all = np.float32(WW_data_files_all)

SubjectsData = []
SubjectsData.append(GG_data_files_all)
SubjectsData.append(GW_data_files_all)
SubjectsData.append(WW_data_files_all)

# 2. subject label: prediction score
#labelpath =  '/ibmgpfs/cuizaixu_lab/xuhaoshu/WM_prediction/datasets/brainproject/code_WMpost/step03_getFCmatrix/getDemo/CCNP_only/CCNP_used_demo_n209.csv';
labelpath = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/code/4th_prediction/s01_preparedata/subid_meanFD_age_sex_522.csv';
label_files_all = pd.read_csv(labelpath)
dimention = targetStr 
label = label_files_all[dimention]
y_label = np.array(label)
OverallPsyFactor = y_label

# 3. covariates  
#covariatespath = '/ibmgpfs/cuizaixu_lab/xuhaoshu/WM_prediction/datasets/brainproject/code_WMpost/step03_getFCmatrix/getDemo/CCNP_only/CCNP_used_demo_n209.csv';
covariatespath = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/code/4th_prediction/s01_preparedata/subid_meanFD_age_sex_522.csv';
Covariates = pd.read_csv(covariatespath, header=0)
Covariates = Covariates.values
Covariates = Covariates[:, [2,3]].astype(float) # sex, motion
# subID,age,sex,meanFD
# Range of parameters
ComponentNumber_Range = np.arange(10) + 1
FoldQuantity = 5
Parallel_Quantity = 1
CVtimes = 101

# # Predict
# ResultantFolder = outFolder + '/RegressCovariates_RandomCV'
# PLSr1_CZ_Random_RegressCovariates.PLSr1_KFold_RandomCV_MultiTimes(SubjectsData, OverallPsyFactor, Covariates, FoldQuantity, ComponentNumber_Range, CVtimes, ResultantFolder, Parallel_Quantity, 0)

# Permutation
ResultantFolder = outFolder + '/RegressCovariates_RandomCV_Permutation';
PLSr1_CZ_Random_RegressCovariates.PLSr1_KFold_RandomCV_MultiTimes(SubjectsData, OverallPsyFactor, Covariates, FoldQuantity, ComponentNumber_Range, 1000, ResultantFolder, Parallel_Quantity, 1)


