#!/usr/bin/env bash
set -euo pipefail

SCRIPT="Rg_process_replica.py"

# ensure the script exists here
if [[ ! -f "$SCRIPT" ]]; then
  echo "Error: $SCRIPT not found in $(pwd)"
  exit 1
fi

# loop over each replica folder
for rep_dir in replicas/replica*/; do
  echo "=== Setting up $rep_dir ==="
  
  # copy the processor script into the replica folder
  cp "$SCRIPT" "$rep_dir"
  echo "→ Copied $SCRIPT to $rep_dir"

  # run it inside that folder
  (
    cd "$rep_dir"
    echo "⏳ Running $SCRIPT in $rep_dir"
    python "$SCRIPT"
    echo "✅ Finished $rep_dir"
  )
done

echo "🎉 All replicas processed."
