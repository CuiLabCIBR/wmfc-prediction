#!/usr/bin/env bash
set -euo pipefail


SUBJECTS_FILE="/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/results/check_data/has_smooth_data_new.txt"  
PYTHON_BIN="${PYTHON_BIN:-python}"
SCRIPT="/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/code/4th_calFC/s02_calFC.py"    

PARTITION="q_fat_c"
CPUS_PER_TASK=2

LOG_DIR="/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/code/4th_calFC/logs"


[[ -f "$SCRIPT" ]] || { echo "ERROR: No Python script found: $SCRIPT"; exit 1; }
[[ -f "$SUBJECTS_FILE" ]] || { echo "ERROR: No participant list found: $SUBJECTS_FILE"; exit 1; }
mkdir -p "$LOG_DIR"


mapfile -t SUBJECTS < <(grep -Ev '^\s*($|#)' "$SUBJECTS_FILE" | tr -d '\r')
(( ${#SUBJECTS[@]} > 0 )) || { echo "ERROR: $SUBJECTS_FILE No valid participants ID"; exit 1; }

echo "Submit ${#SUBJECTS[@]} jobs to the partition: ${SUBJECTS[*]}"


SLURM_OPTS=(
  -p "$PARTITION"
  --cpus-per-task="$CPUS_PER_TASK"
)


for sub in "${SUBJECTS[@]}"; do
  jobname="fc_${sub}"

  CMD="export OMP_NUM_THREADS=$CPUS_PER_TASK MKL_NUM_THREADS=$CPUS_PER_TASK OPENBLAS_NUM_THREADS=$CPUS_PER_TASK NUMEXPR_NUM_THREADS=$CPUS_PER_TASK; \
       SUBJECT='${sub}' ${PYTHON_BIN} '${SCRIPT}'"

  echo "Submitting ${jobname} ..."
  sbatch "${SLURM_OPTS[@]}" -J "$jobname" \
    -o "${LOG_DIR}/%x_%j.out" -e "${LOG_DIR}/%x_%j.err" \
    --export=ALL,SUBJECT="${sub}" \
    --wrap "$CMD"


done

echo "All jobs have been submitted."

