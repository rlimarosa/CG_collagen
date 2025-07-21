project_root/
│
│ ├── RG_process_replica.py                    # compute per‑frame Rg → CSV
│ ├── E2E_process_replica.py                   # compute per‑frame end‑to‑end distance → CSV
│
├── replicas/                                  # each system’s replica subfolders
│ ├── POG12/
│ │ └── replicas/
│ │ ├── replica1/
│ │ ├── replica2/
│ │ └── replica3/
│
├── Rg_aggregate_replicas.py                    # aggregate per‑replica Rg → frame‑wise mean CSV
├── E2E_aggregate_replicas.py                   # aggregate per‑replica E2E → frame‑wise mean CSV
│
├── run_all-Rg.sh                               # copy & run RG_process_replica.py in every replica folder
├── run_all-E2E.sh                              # copy & run E2E_process_replica.py in every replica folder
│
├── notebooks/                                  # Jupyter notebooks for plotting
│ ├── Rg-analysis_Plots.ipynb                   # reads *_meanRg.csv → time vs Rg plots
│ └── E2E-analysis.ipynb                        # reads meanRg/E2E CSVs → 2D free‐energy surfaces
│
└── README.md                                   # this file

Workflow

Copy and run python scripts for each replica:
   ./run_all-Rg.sh
   ./run_all-E2E.sh

Compute all the Rg and E2E doing an avg:
  python Rg_aggregate_replicas.py
  python E2E_aggregate_replicas.py

Use the notebooks to generate the plots:
  Rg-analysis_Plots.ipynb
  E2E-analysis.ipynb
  E2E-analysis.ipynb
