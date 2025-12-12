#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -p q_cn
#SBATCH -e ./log/job.%j.log # Standard error
#SBATCH -o ./log/job.%j.out.txt # Standard output

module load singularity
sub_id=$1
subj=$(sed -n "${sub_id}p" subName_used.txt)
echo ${subj}
#subj=NDARINV02UVMTY7

Path=/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/data
fslic=/ibmgpfs/cuizaixu_lab/xulongzhou/tool/freesurfer
templateflow=/ibmgpfs/cuizaixu_lab/xulongzhou/tool/templateflow

# running xcpd
output=${Path}/step_2nd_24PcsfGlobal
mkdir -p ${output}
wd=${Path}/step_2nd_wd/sub-${subj}
mkdir -p ${wd}
unset PYTHONPATH
export SINGULARITYENV_TEMPLATEFLOW_HOME=$templateflow
singularity run --cleanenv \
        -B ${bids}:/bids \
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
        -w /wd --nthreads 2 --mem-gb 20 \
        --nuisance-regressors 24P --despike -c /custom_confounds \
        --lower-bpf=0.01 --upper-bpf=0.1 \
        --smoothing 0 \
        --motion-filter-type lp --band-stop-min 6 \
        --skip-parcellation \
        --fd-thresh -1

rm -rf $wd