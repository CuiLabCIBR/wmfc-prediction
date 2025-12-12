#coding: utf-8
import scipy.io as sio
import numpy as np
import pandas as pd
import os
import sys
sys.path.append('/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/code/5th_prediction/PLS_prediction/networkFC/cognition_nosmooth_motion2runFD');
import PLSr1_CZ_Random_RegressCovariates

targetStr_list = ['nihtbx_fluidcomp', 'nihtbx_cryst', 'nihtbx_totalcomp']
           
FCpath =  '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/results/FC_nosmooth_merge/'
demopath = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/code/5th_prediction/prepare_data/cognition/'
 
for targetStr in targetStr_list:
    outFolder = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/code/5th_prediction/PLS_prediction/networkFC/cognition_nosmooth_motion2runFD/'+ targetStr
    # Import data
    # 1. atlas loading
    GG_datapath = FCpath + 'total_FCvector_GG_n4388_cognition.txt';
    GW_datapath = FCpath + 'total_FCvector_GW_60p_n4388_cognition.txt';
    WW_datapath = FCpath + 'total_FCvector_WW_60p_n4388_cognition.txt';

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
    labelpath =  demopath + 'abcd_mergedTable_with_motion2run_deNaN_from_n5959_QCpassed_n4388_label.csv';
    label_files_all = pd.read_csv(labelpath)
    dimention = targetStr 
    label = label_files_all[dimention]
    y_label = np.array(label)
    OverallPsyFactor = y_label

    # 3. covariates  
    covariatespath = demopath + 'abcd_mergedTable_with_motion2run_deNaN_from_n5959_QCpassed_n4388_covariates.csv';
    Covariates = pd.read_csv(covariatespath, header=0)
    Covariates = Covariates.values
    Covariates = Covariates[:, 1:].astype(float) # age, sex, motion, site
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


