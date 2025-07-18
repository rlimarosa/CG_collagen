#!/bin/bash -l
#
# ───────────────── SLURM DIRECTIVES ─────────────────
#SBATCH --job-name=${OUTPREFIX}
#SBATCH --output=${OUTPREFIX}.out
#SBATCH --error=${OUTPREFIX}.err         
#SBATCH --partition=gpu             
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=1
#SBATCH --ntasks-per-node=1              
#SBATCH --time=24:00:00              

# ─────────────────── PREAMBLE ───────────────────
module --force purge
module load modules/2.3-20240529 cuda/12.3.2

# Activate virtual environment
source ~/openmm-env/bin/activate

cd $SLURM_SUBMIT_DIR

# Determine which PDB to use (should be only one in folder matching POG*.pdb)
pdb_file=$(ls POG*_1trimervs3.pdb | head -n 1)

# Run setup
python setup_cgsimulation.py \
    --pdb "$pdb_file" \
    --outprefix "${OUTPREFIX}" \
    --temp "${TEMP}"
