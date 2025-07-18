#!/bin/bash

# Define system names
systems=(8 10 12 14 16 18 20)

# Loop through each system folder
for sys in "${systems[@]}"; do
    system_dir="POG${sys}"
    script_path="./${system_dir}/fraction-analysis.py"

    if [[ -f "$script_path" ]]; then
        echo "📊 Running fraction analysis in ${system_dir}"
        (cd "$system_dir" && python fraction-analysis.py)
    else
        echo "⚠️  Missing script: ${script_path}"
    fi
done
