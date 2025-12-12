#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH -p q_fat_c
#SBATCH -q high_c


module load singularity/3.7.0
subj=$1

bids=/ibmgpfs/cuizaixu_lab/liyang/BrainProject25/Tsinghua_data/BIDS
wd=/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/wd
output=/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/fmriprep
temp_dir=/ibmgpfs/cuizaixu_lab/congjing/pro_temp/${subj}
mkdir $temp_dir

mkdir -p  $wd
rm -rf $wd/$subj
mkdir -p $wd/$subj
mkdir -p $output
rm -rf $output/$subj
fs_license=/ibmgpfs/cuizaixu_lab/congjing/freesurfer_license
templateflow=/ibmgpfs/cuizaixu_lab/congjing/softwarepackages/templateflow

#Run fmriprep
echo ""
echo "Running fmriprep on participant: $subj"
echo ""

#Run fmriprep
export SINGULARITYENV_TEMPLATEFLOW_HOME=$templateflow

unset PYTHONPATH; 
singularity run --cleanenv \
    -B $wd/$subj:/wd \
    -B $bids:/bids \
    -B $output:/output \
    -B $fs_license:/fs_license \
    -B $templateflow:$templateflow \
	-B /ibmgpfs/cuizaixu_lab/congjing/tmp:/tmp \
    -B $temp_dir:$HOME \
	/ibmgpfs/cuizaixu_lab/xulongzhou/apps/singularity/fmriprep-24.0.0.simg \
	/bids /output participant -w /wd \
    --n-cpus 2 --omp-nthreads 2 --mem-mb 40G \
    --participant_label $subj \
    --fs-license-file /fs_license/license.txt \
    --ignore slicetiming \
    --skip_bids_validation \
    --task-id rest \
    --output-spaces T1w MNI152NLin6Asym:res-2 \
    --stop-on-first-crash \
    --resource-monitor \
	--fs-no-reconall \
    --notrack
	
rm -rf $wd/$subj

