#!/usr/bin/env bash
#
# Usage: ./submit_temps.sh 300 310 320
#

if [ $# -lt 1 ]; then
  echo "Usage: $0 TEMP1 [TEMP2 ...]"
  exit 1
fi

for TEMP in "$@"; do
  OUTPREFIX="cgtriple${TEMP}K"
  echo "Submitting T=${TEMP} K → prefix=${OUTPREFIX}"
  sbatch \
    --job-name="${OUTPREFIX}" \
    --output="${OUTPREFIX}.out" \
    --error="${OUTPREFIX}.err" \
    --export=ALL,TEMP="${TEMP}",OUTPREFIX="${OUTPREFIX}" \
    submit.sh
done
