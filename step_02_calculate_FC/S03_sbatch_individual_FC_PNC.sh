#!/bin/bash
#SBATCH -p q_fat
#SBATCH --ntasks=1 # Run a single serial task
#SBATCH --cpus-per-task=2

module load MATLAB/R2019a
cd /ibmgpfs/cuizaixu_lab/congjing/WM_prediction/PNC/code/3rd_calFC
i_sub=$1 
matlab -nodisplay -nosplash -r "s02_get_individualFCmatrix('$i_sub')"