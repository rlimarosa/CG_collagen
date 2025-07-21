#!/usr/bin/env bash
set -euo pipefail

SCRIPT="E2E_process_replica.py"

if [[ ! -f "$SCRIPT" ]]; then
  echo "Error: $SCRIPT not found here." >&2
  exit 1
fi

for rep_dir in replicas/replica*/; do
  echo "=== Processing $rep_dir ==="
  cp "$SCRIPT" "$rep_dir"
  (
    cd "$rep_dir"
    python "$SCRIPT"
  )
done

echo "🎉 All replicas processed."
