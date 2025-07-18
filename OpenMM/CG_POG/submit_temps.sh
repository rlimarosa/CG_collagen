#!/usr/bin/env bash
#
# Usage: ./submit_temps.sh 300 310 320 ...

if [ $# -lt 1 ]; then
  echo "Usage: $0 TEMP1 [TEMP2 ...]"
  exit 1
fi

systems=(8 10 12 14 16 18 20)
replicas=(1 2 3)

for TEMP in "$@"; do
  for sys in "${systems[@]}"; do
    system_name="POG${sys}"
    
    for rep in "${replicas[@]}"; do
      replica_dir="./${system_name}/replicas/replica${rep}"
      if [[ ! -d "$replica_dir" ]]; then
        echo "⚠️ Missing directory: $replica_dir, skipping"
        continue
      fi

      cd "$replica_dir" || continue
      OUTPREFIX="${system_name}_${TEMP}K"
      echo "🚀 Submitting → $system_name | replica${rep} | ${TEMP}K"

      sbatch \
        --job-name="${OUTPREFIX}" \
        --output="${OUTPREFIX}.out" \
        --error="${OUTPREFIX}.err" \
        --export=ALL,TEMP="${TEMP}",OUTPREFIX="${OUTPREFIX}" \
        submit.sh

      cd - > /dev/null
    done
  done
done
