#!/bin/bash
#-----------------------------------------------------------------------
# Submit the remaining case-study arrays as QoS capacity frees.
#
# The `short` QoS caps submitted jobs per user at 100 and every array element
# counts, so three 30-51 task arrays cannot all be queued at once. This waits
# for room and submits in priority order rather than requiring someone to sit
# and retry.
#
#   quiron       the substantive Section 3.2 dataset
#   communities  highest measured shift; the stress test for the slope model
#   acs_pums     replication of the BRFSS design in a different subject area
#
# Usage:  nohup bash slurm/submit_case_queue.sh > logs/submit_queue.log 2>&1 &
#-----------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
LIMIT=100
MODEL=gi_hd_slope

# dataset:array_size
QUEUE=("quiron:51" "communities:30" "acs_pums:51")

for entry in "${QUEUE[@]}"; do
  ds="${entry%%:*}"
  n="${entry##*:}"
  while true; do
    used=$(squeue -u "$USER" -h -r -o "%i" 2>/dev/null | wc -l)
    room=$(( LIMIT - used ))
    if [ "$room" -ge "$n" ]; then
      echo "$(date +%H:%M:%S) submitting $ds ($n tasks), room $room"
      if sbatch --array=1-"$n"%15 slurm/06_case_study.sh "$ds" "$MODEL"; then
        sleep 30
        break
      fi
      echo "$(date +%H:%M:%S) sbatch refused $ds; retrying"
    else
      echo "$(date +%H:%M:%S) waiting for $ds: need $n, room $room"
    fi
    sleep 180
  done
done
echo "$(date +%H:%M:%S) all case-study arrays submitted"
