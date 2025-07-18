project_root/
│
├── structures/                         # Contains input PDBs for each system (e.g., POG14_1trimervs3.pdb)
│
├── source/                             # Template scripts copied into system/replica folders
│   ├── setup_cgsimulation.py           # Main simulation script  - FF and MD simulation parameters (used by submit.sh)
│   ├── contacts-analysis.py            # Analyzes hydrogen bond contacts per frame - HBD-HBA beads contacts in 4 angstrons
│   ├── fraction-analysis.py            # Analyzes average contact fraction across replicas 
│
├── POG14/                              # Example system folder (one for each POGx system)
│   ├── fraction-analysis.py            # Calculates melting curve from all replicas
│   └── replicas/
│       ├── replica1/
│       │   ├── setup_cgsimulation.py
│       │   ├── contacts-analysis.py
│       │   ├── POG14_1trimervs3.pdb
│       │   ├── POG14_300K.pdb          # PDB of the last frame of the simulation
│       │   └── POG14_300K.dcd
│       ├── replica2/
│       └── replica3/
│
├── setup.sh                            # Creates folders and copies template scripts into each system/replica
├── submit.sh                           # SLURM job script to run a simulation using setup_cgsimulation.py
├── submit_temps.sh                     # Submits jobs for all systems/replicas at specified temperatures
├── run_all_contacts-analysis.sh        # Runs contacts-analysis.py in every replica folder
├── run_all_fractions-analysis.sh       # Runs fraction-analysis.py in every system folder


Workflow

Setup folders for each desire system
./setup.sh

Run simulations for desired temperatures (e.g. 300 310 320)
./submit_temps.sh 300 310 320

Wait for all jobs to finish, then compute contacts:
./run_all_contacts-analysis.sh

After all contact files are generated, calculate melting curves:
./run_all_fractions-analysis.sh
