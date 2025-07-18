#!/bin/bash

# List of systems and replicas
systems=(8 10 12 14 16 18 20)
replicas=(1 2 3)

# Loop over each system and replica
for sys in "${systems[@]}"; do
    system_name="POG${sys}"
    
    for rep in "${replicas[@]}"; do
        rep_dir="./${system_name}/replicas/replica${rep}"
        script="${rep_dir}/contacts-analysis.py"

        if [[ -f "$script" ]]; then
            echo "🔍 Running contact analysis for ${system_name} replica${rep}"
            (cd "$rep_dir" && python contacts-analysis.py)
        else
            echo "⚠️ Script not found: $script"
        fi
    done
done
