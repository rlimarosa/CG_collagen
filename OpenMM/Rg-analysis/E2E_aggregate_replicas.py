#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

# ─── Gather all per‑replica CSVs ─────────────────────────────────────────────
csvs = sorted(Path('replicas').glob('replica*/*_e2e.csv'))
if not csvs:
    raise FileNotFoundError("No '*_e2e.csv' files found under 'replicas/'")

# ─── Read & concatenate ──────────────────────────────────────────────────────
df = pd.concat((pd.read_csv(p) for p in csvs), ignore_index=True)

# ─── Compute frame‑wise mean E2E (only on the numeric column) ───────────────
mean_e2e = (
    df
    .groupby(['Temperature (K)', 'Chain', 'Frame'], as_index=False)['E2E_Å']
    .mean()
    .rename(columns={'E2E_Å': 'Mean E2E (Å)'})
    .sort_values(['Temperature (K)', 'Chain', 'Frame'])
)

# ─── Save final summary ──────────────────────────────────────────────────────
system_name = csvs[0].stem.split('_')[0]  # e.g. "POG12"
out = Path(f"{system_name}_replicas_framewise_meanE2E.csv")
mean_e2e.to_csv(out, index=False)
print(f"✅ Saved aggregated summary to {out}")
