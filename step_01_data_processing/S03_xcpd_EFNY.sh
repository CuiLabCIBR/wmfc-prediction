#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -p q_cn

module load singularity
subj=$1

fmriprep_Path=/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/fmriprep
xcpd_Path=/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/xcpd
wd=/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/wd/xcpd_rest_cifti
temp_dir=/ibmgpfs/cuizaixu_lab/congjing/pro_temp/${subj}
mkdir $temp_dir
fslic=/ibmgpfs/cuizaixu_lab/congjing/freesurfer_license
templateflow=/ibmgpfs/cuizaixu_lab/congjing/softwarepackages/templateflow

customCP=/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/custom_confounds

mkdir -p  $wd
rm -rf $wd/$subj
mkdir -p $wd/$subj

#Run xcpd
echo ""
echo "Running xcpd on participant: $subj"
echo ""


unset PYTHONPATH
export SINGULARITYENV_TEMPLATEFLOW_HOME=$templateflow
singularity run --cleanenv \
        -B $fmriprep_Path:/fmriprep \
        -B $xcpd_Path:/xcpd_Path \
        -B $wd/$subj:/wd \
        -B $customCP:/custom_confounds \
        -B $fslic:/fslic \
        -B $templateflow:$templateflow \
		-B /ibmgpfs/cuizaixu_lab/congjing/tmp:/tmp \
        -B $temp_dir:$HOME \
        /ibmgpfs/cuizaixu_lab/xulongzhou/apps/singularity/xcpd-0.7.1rc5.simg \
        /fmriprep /xcpd_Path participant \
        --participant_label ${subj} --task-id rest \
        --fs-license-file /fslic/license.txt \
        -w /wd --nthreads 2 --mem-gb 40 \
        --nuisance-regressors 24P --despike -c /custom_confounds \
        --lower-bpf=0.01 --upper-bpf=0.1 \
        --smoothing 0 \
        --motion-filter-type lp --band-stop-min 6 \
        --skip-parcellation \
        --fd-thresh -1

rm -rf $wd/$subj
rm -rf $temp_dir