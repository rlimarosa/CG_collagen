#!/usr/bin/env python3
import mdtraj as md
import numpy as np
import pandas as pd
import glob, re
from pathlib import Path

# ─── Settings ────────────────────────────────────────────────────────────────
backbone_names = ['BBP','BBO','BBG']
chains         = ['A','B','C']

# ─── Identify files ──────────────────────────────────────────────────────────
cwd = Path().resolve()
replica = cwd.name  # e.g. "replica1"

# Grab anything ending in 'K.dcd' (e.g. POG12_300K.dcd)
dcd_files = sorted(glob.glob('*K.dcd'))
if not dcd_files:
    raise FileNotFoundError("No '*K.dcd' files found in this folder.")

# Extract system name (e.g. 'POG12' from 'POG12_300K.dcd')
system_name = dcd_files[0].split('_')[0]

# Find one matching PDB for alignment/selection
pdb_refs = sorted(glob.glob(f'{system_name}_*K.pdb'))
if not pdb_refs:
    raise FileNotFoundError(f"No '{system_name}_*K.pdb' found in this folder.")
ref = md.load(pdb_refs[0])

# ─── Precompute atom selections ──────────────────────────────────────────────
sel_idx = {}
for i, ch in enumerate(chains):
    sel_str = (
        f"chainid {i} and ("
        + " or ".join(f"name {b}" for b in backbone_names)
        + ")"
    )
    idx = ref.topology.select(sel_str)
    if len(idx) == 0:
        raise ValueError(f"No atoms selected for chain {ch} with '{sel_str}'")
    sel_idx[ch] = idx
    print(f"Chain {ch}: {len(idx)} atoms")

# ─── Process each trajectory ─────────────────────────────────────────────────
records = []
for dcd in dcd_files:
    m = re.search(rf'{system_name}_(\d+)K\.dcd', dcd)
    if not m:
        continue
    temp = int(m.group(1))
    top  = dcd.replace('.dcd', '.pdb')
    if not Path(top).exists():
        raise FileNotFoundError(f"Missing topology {top} for {dcd}")

    traj = md.load(dcd, top=top)
    traj.superpose(ref)

    for ch in chains:
        rg = md.compute_rg(traj.atom_slice(sel_idx[ch])) * 10.0
        for frame, val in enumerate(rg):
            records.append({
                'Replica':         replica,
                'Temperature (K)': temp,
                'Chain':           ch,
                'Frame':           frame,
                'Rg (Å)':          val
            })

# ─── Save CSV ────────────────────────────────────────────────────────────────
df = pd.DataFrame(records)
out = cwd / f"{system_name}_{replica}_Rg.csv"
df.to_csv(out, index=False)
print(f"✅ Saved data to {out}")
