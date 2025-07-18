#!/bin/bash

# Define system and replica ranges
systems=(8 10 12 14 16 18 20)
replicas=(1 2 3)

# Define paths
structure_dir="./structures"
source_dir="./source"

# Required source files (only those placed inside replica folders)
replica_files=("setup_cgsimulation.py" "contacts-analysis.py" "submit.sh")

# Check for fraction-analysis.py for top-level
fraction_script="${source_dir}/fraction-analysis.py"
if [[ ! -f "$fraction_script" ]]; then
    echo "❌ Missing fraction-analysis.py in source/"
    exit 1
fi

# Validate replica-scope files
for file in "${replica_files[@]}"; do
    if [[ ! -f "${source_dir}/${file}" ]]; then
        echo "❌ Missing file in source: ${file}"
        exit 1
    fi
done

# Create directories and copy files
for sys in "${systems[@]}"; do
    system_name="POG${sys}"
    pdb_file="${system_name}_1trimervs3.pdb"
    pdb_path="${structure_dir}/${pdb_file}"

    if [[ ! -f "$pdb_path" ]]; then
        echo "❌ PDB not found: $pdb_path"
        continue
    fi

    # Copy fraction-analysis.py to system-level
    mkdir -p "./${system_name}/replicas"
    cp "$fraction_script" "./${system_name}/fraction-analysis.py"

    for rep in "${replicas[@]}"; do
        replica_dir="./${system_name}/replicas/replica${rep}"
        mkdir -p "$replica_dir"

        cp "$pdb_path" "$replica_dir/"
        for file in "${replica_files[@]}"; do
            cp "${source_dir}/${file}" "$replica_dir/"
        done

        echo "✅ Setup for ${system_name}/replica${rep}"
    done
done
