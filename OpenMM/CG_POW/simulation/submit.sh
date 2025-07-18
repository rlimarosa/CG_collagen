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
#SBATCH --time=24:00:00               # walltime hh:mm:ss
         
# ─────────────────── PREAMBLE ───────────────────
module --force purge
module load modules/2.3-20240529 cuda/12.3.2

# Activate your Python venv 
source ~/openmm-env/bin/activate

# Change into the directory where you submitted this job
cd $SLURM_SUBMIT_DIR

# ─────────────────── RUN THE PYTHON SCRIPT ───────────────────
# You can pass --pdb, --Tstar, and --outprefix if desired
#python run_simple-vs9.py --pdb rebuilt_structure.pdb --outprefix cgtriple300K


# ─────────────────── RUN THE PYTHON SCRIPT ───────────────────
# Expect two environment variables:
#   TEMP      → temperature in K
#   OUTPREFIX → prefix for all output files
python setup_cgsimulation_POW.py \
    --pdb POG5POWPOG6_1trimervs3.pdb \
    --outprefix ${OUTPREFIX} \
    --temp ${TEMP}