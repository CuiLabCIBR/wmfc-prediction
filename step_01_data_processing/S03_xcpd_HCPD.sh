#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH -p q_fat


module load singularity
subj=$1
PublicDataPath_rsfMRI=/ibmgpfs/cuizaixu_lab/Public_Data/HCPD/processed_data/rsfMRI/fmriprep_rest_no_MSM
PublicDataPath_sMRI=/ibmgpfs/cuizaixu_lab/Public_Data/HCPD/rawdata/sMRI_processed_20210514

Path=/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/HCPD/data

fslic=/ibmgpfs/cuizaixu_lab/xulongzhou/tool/freesurfer
templateflow=/ibmgpfs/cuizaixu_lab/xulongzhou/tool/templateflow

### 1. construct bids-format file ###
bids=${Path}/bids

## 1.1 copy func data from public path but without surface data(91k)
mkdir -p ${bids}/sub-${subj}/func
rsync -av --exclude='*91k*' ${PublicDataPath_rsfMRI}/sub-${subj}/func/ /${bids}/sub-${subj}/func
# rename file with 2009cAsym to 6Asym, actually is 6Asym!
for f in ${bids}/sub-${subj}/func/*2009cAsym*; do
    newname="${f//2009cAsym/6Asym}"
    mv "$f" "$newname"
done
# rename confounds file
for f in ${bids}/sub-${subj}/func/*confounds_timeseries.tsv; do
    newname="${f//timeseries.tsv/timeseries_ori.tsv}"
    mv "$f" "$newname"
done

## 1.2 copy anat data in individual space
mkdir ${bids}/sub-${subj}/anat
cp ${PublicDataPath_sMRI}/${subj}_V1_MR/MNINonLinear/T1w.nii.gz               ${bids}/sub-${subj}/anat/sub-${subj}_desc-preproc_T1w.nii.gz
cp ${PublicDataPath_sMRI}/${subj}_V1_MR/MNINonLinear/brainmask_fs.nii.gz      ${bids}/sub-${subj}/anat/sub-${subj}_desc-brain_mask.nii.gz
cp ${PublicDataPath_sMRI}/${subj}_V1_MR/MNINonLinear/ribbon.nii.gz            ${bids}/sub-${subj}/anat/sub-${subj}_desc-ribbon_T1w.nii.gz
cp ${PublicDataPath_sMRI}/${subj}_V1_MR/MNINonLinear/wmparc.nii.gz            ${bids}/sub-${subj}/anat/sub-${subj}_desc-wmparc_T1w.nii.gz
cp ${PublicDataPath_sMRI}/${subj}_V1_MR/MNINonLinear/aparc+aseg.nii.gz        ${bids}/sub-${subj}/anat/sub-${subj}_desc-aparcaseg_dseg.nii.gz

## 1.3 convert anat data to func space (6Asym)
module load afni
# T1w
funcTarget=${bids}/sub-${subj}/func/sub-${subj}_task-REST1_acq-AP_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz
anatOri=${PublicDataPath_sMRI}/${subj}_V1_MR/MNINonLinear/T1w.nii.gz
3dresample -master $funcTarget -input $anatOri -prefix ${bids}/sub-${subj}/anat/sub-${subj}_space-MNI152NLin6Asym_desc-preproc_T1w.nii.gz
# brainmask
anatOri=${PublicDataPath_sMRI}/${subj}_V1_MR/MNINonLinear/brainmask_fs.nii.gz
3dresample -master $funcTarget -input $anatOri -prefix ${bids}/sub-${subj}/anat/sub-${subj}_space-MNI152NLin6Asym_desc-preproc_brainmask.nii.gz
# ribbon
anatOri=${PublicDataPath_sMRI}/${subj}_V1_MR/MNINonLinear/ribbon.nii.gz
3dresample -master $funcTarget -input $anatOri -prefix ${bids}/sub-${subj}/anat/sub-${subj}_space-MNI152NLin6Asym_desc-preproc_ribbon.nii.gz
# wmparc
anatOri=${PublicDataPath_sMRI}/${subj}_V1_MR/MNINonLinear/wmparc.nii.gz
3dresample -master $funcTarget -input $anatOri -prefix ${bids}/sub-${subj}/anat/sub-${subj}_space-MNI152NLin6Asym_desc-preproc_wmparc.nii.gz
# aparc+aseg
anatOri=${PublicDataPath_sMRI}/${subj}_V1_MR/MNINonLinear/aparc+aseg.nii.gz
3dresample -master $funcTarget -input $anatOri -prefix ${bids}/sub-${subj}/anat/sub-${subj}_space-MNI152NLin6Asym_desc-preproc_aparc+aseg.nii.gz

## 1.4 construct custom confounds file
customCP=${Path}/custom_confounds_36p_noWM
mkdir -p ${customCP}
# AP_REST1
module load MATLAB/R2019a
filePath=${bids}/sub-${subj}/func/sub-${subj}_task-REST1_acq-AP_desc-confounds_timeseries_ori.tsv
savePath_bids=${bids}/sub-${subj}/func/sub-${subj}_task-REST1_acq-AP_desc-confounds_timeseries.tsv
rm -rf ${savePath_bids}
savePath_customCP=${customCP}/sub-${subj}_task-REST1_acq-AP_desc-confounds_timeseries.tsv
rm -rf ${savePath_customCP}
matlab -nodisplay -nosplash -nodesktop -r \
    "extract_columns_by_title_new2025('${filePath}', '${savePath_bids}', '${savePath_customCP}', {'csf','global_signal','csf_derivative1', 'global_signal_derivative1','csf_power2','global_signal_power2','csf_derivative1_power2','global_signal_derivative1_power2'}); exit;"
# PA_REST1
filePath=${bids}/sub-${subj}/func/sub-${subj}_task-REST1_acq-PA_desc-confounds_timeseries_ori.tsv
savePath_bids=${bids}/sub-${subj}/func/sub-${subj}_task-REST1_acq-PA_desc-confounds_timeseries.tsv
rm -rf ${savePath_bids}
savePath_customCP=${customCP}/sub-${subj}_task-REST1_acq-PA_desc-confounds_timeseries.tsv
rm -rf ${savePath_customCP}
matlab -nodisplay -nosplash -nodesktop -r \
    "extract_columns_by_title_new2025('${filePath}', '${savePath_bids}', '${savePath_customCP}', {'csf','global_signal','csf_derivative1', 'global_signal_derivative1','csf_power2','global_signal_power2','csf_derivative1_power2','global_signal_derivative1_power2'}); exit;"
# AP_REST2
filePath=${bids}/sub-${subj}/func/sub-${subj}_task-REST2_acq-AP_desc-confounds_timeseries_ori.tsv
savePath_bids=${bids}/sub-${subj}/func/sub-${subj}_task-REST2_acq-AP_desc-confounds_timeseries.tsv
rm -rf ${savePath_bids}
savePath_customCP=${customCP}/sub-${subj}_task-REST2_acq-AP_desc-confounds_timeseries.tsv
rm -rf ${savePath_customCP}
matlab -nodisplay -nosplash -nodesktop -r \
    "extract_columns_by_title_new2025('${filePath}', '${savePath_bids}', '${savePath_customCP}', {'csf','global_signal','csf_derivative1', 'global_signal_derivative1','csf_power2','global_signal_power2','csf_derivative1_power2','global_signal_derivative1_power2'}); exit;"
# PA_REST2
filePath=${bids}/sub-${subj}/func/sub-${subj}_task-REST2_acq-PA_desc-confounds_timeseries_ori.tsv
savePath_bids=${bids}/sub-${subj}/func/sub-${subj}_task-REST2_acq-PA_desc-confounds_timeseries.tsv
rm -rf ${savePath_bids}
savePath_customCP=${customCP}/sub-${subj}_task-REST2_acq-PA_desc-confounds_timeseries.tsv
rm -rf ${savePath_customCP}
matlab -nodisplay -nosplash -nodesktop -r \
    "extract_columns_by_title_new2025('${filePath}', '${savePath_bids}', '${savePath_customCP}', {'csf','global_signal','csf_derivative1', 'global_signal_derivative1','csf_power2','global_signal_power2','csf_derivative1_power2','global_signal_derivative1_power2'}); exit;"

## 1.5 get bold reference file
anatPath=${bids}/sub-${subj}/anat/sub-${subj}_space-MNI152NLin6Asym_desc-preproc_wmparc.nii.gz
# AP_REST1
filePath=${bids}/sub-${subj}/func/sub-${subj}_task-REST1_acq-AP_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz
savePath=${bids}/sub-${subj}/func/sub-${subj}_task-REST1_acq-AP_space-MNI152NLin6Asym_boldref.nii
matlab -nodisplay -nosplash -nodesktop -r \
    "get_boldref('${filePath}', '${anatPath}', '${savePath}'); exit;"
gzip ${savePath}
# PA_REST1
filePath=${bids}/sub-${subj}/func/sub-${subj}_task-REST1_acq-PA_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz
savePath=${bids}/sub-${subj}/func/sub-${subj}_task-REST1_acq-PA_space-MNI152NLin6Asym_boldref.nii
matlab -nodisplay -nosplash -nodesktop -r \
    "get_boldref('${filePath}', '${anatPath}', '${savePath}'); exit;"
gzip ${savePath}
# AP_REST2
filePath=${bids}/sub-${subj}/func/sub-${subj}_task-REST2_acq-AP_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz
savePath=${bids}/sub-${subj}/func/sub-${subj}_task-REST2_acq-AP_space-MNI152NLin6Asym_boldref.nii
matlab -nodisplay -nosplash -nodesktop -r \
    "get_boldref('${filePath}', '${anatPath}', '${savePath}'); exit;"
gzip ${savePath}
# PA_REST2
filePath=${bids}/sub-${subj}/func/sub-${subj}_task-REST2_acq-PA_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz
savePath=${bids}/sub-${subj}/func/sub-${subj}_task-REST2_acq-PA_space-MNI152NLin6Asym_boldref.nii
matlab -nodisplay -nosplash -nodesktop -r \
    "get_boldref('${filePath}', '${anatPath}', '${savePath}'); exit;"
gzip ${savePath}

## 1.6 add the BIDS version in the json file
cp /ibmgpfs/cuizaixu_lab/zhaoshaoling/MSC_data/HCPD/code_xcpd0.7.1rc5_hcpMiniPrepData/dataset_description.json ${bids}/dataset_description.json

## 1.7 copy T1w-MNI and MNI-T1wtxt file from HCP to HCPD
cp /ibmgpfs/cuizaixu_lab/zhaoshaoling/HCP_WMfMRI/hcp_rest_2/xcpd0.7.1rc5_rest2/bids/sub-100206/anat/sub-100206_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.txt           ${bids}/sub-${subj}/anat/sub-${subj}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.txt
cp /ibmgpfs/cuizaixu_lab/zhaoshaoling/HCP_WMfMRI/hcp_rest_2/xcpd0.7.1rc5_rest2/bids/sub-100206/anat/sub-100206_from-T1w_to-MNI152NLin6Asym_mode-image_xfm.txt         ${bids}/sub-${subj}/anat/sub-${subj}_from-T1w_to-MNI152NLin6Asym_mode-image_xfm.txt

## 1.8 construct deseg file
wmparc=${bids}/sub-${subj}/anat/sub-${subj}_space-MNI152NLin6Asym_desc-preproc_wmparc.nii.gz
brainmask=${bids}/sub-${subj}/anat/sub-${subj}_space-MNI152NLin6Asym_desc-preproc_brainmask.nii.gz
dseg=${bids}/sub-${subj}/anat/sub-${subj}_dseg.nii
matlab -nodisplay -nosplash -nodesktop -r "convert_wmparc2dseg('$wmparc', '$dseg', '$brainmask'); exit;"
gzip ${dseg}

## 1.9 add brain mask for each func data
cp ${bids}/sub-${subj}/anat/sub-${subj}_space-MNI152NLin6Asym_desc-preproc_brainmask.nii.gz ${bids}/sub-${subj}/func/sub-${subj}_task-REST1_acq-AP_space-MNI152NLin6Asym_desc-brain_mask.nii.gz
cp ${bids}/sub-${subj}/anat/sub-${subj}_space-MNI152NLin6Asym_desc-preproc_brainmask.nii.gz ${bids}/sub-${subj}/func/sub-${subj}_task-REST1_acq-PA_space-MNI152NLin6Asym_desc-brain_mask.nii.gz
cp ${bids}/sub-${subj}/anat/sub-${subj}_space-MNI152NLin6Asym_desc-preproc_brainmask.nii.gz ${bids}/sub-${subj}/func/sub-${subj}_task-REST2_acq-AP_space-MNI152NLin6Asym_desc-brain_mask.nii.gz
cp ${bids}/sub-${subj}/anat/sub-${subj}_space-MNI152NLin6Asym_desc-preproc_brainmask.nii.gz ${bids}/sub-${subj}/func/sub-${subj}_task-REST2_acq-PA_space-MNI152NLin6Asym_desc-brain_mask.nii.gz

###### 2. running xcpd ######
output=${Path}/step_2nd_24PcsfGlobal
mkdir -p ${output}
wd=${Path}/step_2nd_wd/sub-${subj}
mkdir -p ${wd}
unset PYTHONPATH
export SINGULARITYENV_TEMPLATEFLOW_HOME=$templateflow
singularity run --cleanenv \
        -B $Path/bids:/bids \
        -B $output:/output \
        -B $wd:/wd \
        -B $customCP:/custom_confounds \
        -B $fslic:/fslic \
        -B $templateflow:$templateflow \
        /ibmgpfs/cuizaixu_lab/xulongzhou/apps/singularity/xcpd-0.7.1rc5.simg \
        /bids /output participant \
        --participant_label ${subj} --task-id rest \
        --input-type fmriprep \
        --fs-license-file /fslic/license.txt \
        -w /wd --nthreads 4 --mem-gb 40 \
        --nuisance-regressors 24P --despike -c /custom_confounds \
        --lower-bpf=0.01 --upper-bpf=0.1 \
        --smoothing 0 \
        --motion-filter-type lp --band-stop-min 6 \
        --skip-parcellation \
        --fd-thresh -1

rm -rf $wd