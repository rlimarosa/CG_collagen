POW Project – Simulation and Analysis of Collagen-Like Peptides with a Defect

───────────────────────────────────────────────────────────────────────────────
📁 Folder Structure
───────────────────────────────────────────────────────────────────────────────

project_root/
│
├── modelling/                          # Contains scripts for model building
│   ├── clp_builder_defect.py          # System builder defining BBW and HBW beads
│   └── cgmodel_defect.ipynb           # Jupyter notebook to tune geometry and connect atoms
│
├── simulation/                         # Main simulation folder
│   ├── setup_cgsimulation_POW.py      # Contains FF and MD setup for POW (20 ns run)
│   │                                   # BBW/HBW beads mimic backbone/hydrogen interactions
│   │                                   # Modified nonbonded sigmas:
│   │                                     - BB-HB (with W): 0.8 (vs. 0.65)
│   │                                     - HB-HB (with W): 0.85 (vs. 0.70)
│   ├── submit.sh                       # SLURM script to run a single temperature
│   └── submit_temps.sh                # Automates submission of simulations at various temperatures
│
└── README.txt                          # This file


───────────────────────────────────────────────────────────────────────────────
⚙️ POW System Overview
───────────────────────────────────────────────────────────────────────────────

• Base peptide: (POG)₁₂ → replaced 1 repeat with a W bead (BBW + HBW)
• BBW: mimics standard backbone beads (BBP, BBO, BBG)
• HBW: mimics hydrogen-cap beads (HBP, HBG)
• Custom bonds, angles, torsions defined in `setup_cgsimulation_POW.py`
• Nonbonded sigma modifications:
    - BB-HB (if W involved): 0.80 (default 0.65)
    - HB-HB (if W involved): 0.85 (default 0.70)
• No changes in WCA potential
• Simulation duration: 20 ns

