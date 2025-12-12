#!/usr/bin/env bash
set -euo pipefail

SUBJECTS_FILE="/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/sublist_final_532.txt"
PYTHON_BIN="${PYTHON_BIN:-python}"
SCRIPT="/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/code/3rd_cal_FC/s02_calFC.py"

PARTITION="q_cn"
CPUS_PER_TASK=2
LOG_DIR="/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/code/3rd_cal_FC/logs"

[[ -f "$SCRIPT" ]] || { echo "ERROR: No Python script found: $SCRIPT"; exit 1; }
[[ -f "$SUBJECTS_FILE" ]] || { echo "ERROR: No participant list found: $SUBJECTS_FILE"; exit 1; }
mkdir -p "$LOG_DIR"

mapfile -t SUBJECTS < <(grep -Ev '^\s*($|#)' "$SUBJECTS_FILE" | tr -d '\r')
(( ${#SUBJECTS[@]} > 0 )) || { echo "ERROR: $SUBJECTS_FILE No valid participants ID"; exit 1; }

echo "Submit ${#SUBJECTS[@]} jobs to the partition ${PARTITION}"

SLURM_OPTS=(-p "$PARTITION" --cpus-per-task="$CPUS_PER_TASK")

JOBIDS=()
for sub in "${SUBJECTS[@]}"; do
  jobname="fc_${sub}"
  CMD="export OMP_NUM_THREADS=$CPUS_PER_TASK MKL_NUM_THREADS=$CPUS_PER_TASK OPENBLAS_NUM_THREADS=$CPUS_PER_TASK NUMEXPR_NUM_THREADS=$CPUS_PER_TASK; \
       ${PYTHON_BIN} '${SCRIPT}'"
  echo "Submitting ${jobname} ..."

  out=$(sbatch "${SLURM_OPTS[@]}" -J "$jobname" \
        -o "${LOG_DIR}/%x_%j.out" -e "${LOG_DIR}/%x_%j.err" \
        --export=ALL,SUBJECT="${sub}" \
        --wrap "$CMD")
  jid=$(awk '{print $4}' <<<"$out")
  JOBIDS+=("$jid")
done


if ((${#JOBIDS[@]} > 0)); then
  dep="afterok:$(IFS=:; echo "${JOBIDS[*]}")"
  echo "Submitting merge job with dependency: $dep"
  sbatch -p "$PARTITION" -J "merge_run_usage" \
    -o "${LOG_DIR}/%x_%j.out" -e "${LOG_DIR}/%x_%j.err" \
    --dependency="$dep" --export=ALL \
    --wrap "MERGE_RUN_USAGE=1 ${PYTHON_BIN} '${SCRIPT}'"
fi

echo "All jobs have been submitted."
