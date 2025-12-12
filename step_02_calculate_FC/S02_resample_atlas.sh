#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH -p q_fat_c


module load afni
base="/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/xcpd/sub-THU20231119140XJR/func/sub-THU20231119140XJR_task-rest_run-1_space-MNI152NLin6Asym_res-2_desc-denoised_bold.nii.gz"
GM_atlas="/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/code/atlas/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"
WM_atlas="/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/code/atlas/rICBM_DTI_81_WMPM_60p_FMRIB58.nii.gz"
save_dir="/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/code/3rd_cal_FC"

# resample GM atlas
3dresample -master $base -input $GM_atlas -prefix ${save_dir}/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm_resample.nii.gz
# resample WM atlas
3dresample -master $base -input $WM_atlas -prefix ${save_dir}/rICBM_DTI_81_WMPM_60p_FMRIB58_resample.nii.gz