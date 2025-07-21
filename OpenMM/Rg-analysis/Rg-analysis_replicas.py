#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

# ─── Gather all per‑replica CSVs ─────────────────────────────────────────────
csvs = sorted(Path('replicas').glob('replica*/*_Rg.csv'))
if not csvs:
    raise FileNotFoundError("No '*_Rg.csv' files found under 'replicas/'")

# ─── Read & concatenate ──────────────────────────────────────────────────────
df = pd.concat((pd.read_csv(p) for p in csvs), ignore_index=True)

# ─── Compute mean Rg per temperature, chain, and frame ────────────────────────
mean_rg = (
    df
    .groupby(['Temperature (K)', 'Chain', 'Frame'])['Rg (Å)']
    .mean()
    .reset_index(name='Mean Rg (Å)')
    .sort_values(['Temperature (K)', 'Chain', 'Frame'])
)

# ─── Save final summary ──────────────────────────────────────────────────────
system_name = csvs[0].stem.split('_')[0]  # e.g. "POG12"
out = Path(f"{system_name}_replicas_framewise_meanRg.csv")
mean_rg.to_csv(out, index=False)
print(f"✅ Saved frame‐wise aggregated summary to {out}")
