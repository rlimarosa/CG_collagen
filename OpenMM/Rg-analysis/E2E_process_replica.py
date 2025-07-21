#!/usr/bin/env python3
import mdtraj as md
import numpy as np
import pandas as pd
import glob, re
from pathlib import Path

# ─── Settings ────────────────────────────────────────────────────────────────
backbone_beads = ['BBP','BBO','BBG']
chain_map      = {0:'A', 1:'B', 2:'C'}

# ─── Locate DCD & PDB ─────────────────────────────────────────────────────────
cwd       = Path().resolve()
replica   = cwd.name  # e.g. "replica1"
dcd_files = sorted(glob.glob('*K.dcd'))
if not dcd_files:
    raise FileNotFoundError("No '*K.dcd' files found in this folder.")

# extract system name (e.g. "POG12")
system_name = dcd_files[0].split('_')[0]

# find a PDB for alignment
pdb_refs = sorted(glob.glob(f'{system_name}_*K.pdb'))
if not pdb_refs:
    raise FileNotFoundError(f"No '{system_name}_*K.pdb' found in this folder.")
ref = md.load(pdb_refs[0])

# ─── Process trajectories ────────────────────────────────────────────────────
records = []
for dcd in dcd_files:
    m = re.search(rf'{system_name}_(\d+)K\.dcd', dcd)
    if not m:
        continue
    temp = int(m.group(1))
    pdb  = dcd.replace('.dcd', '.pdb')
    if not Path(pdb).exists():
        raise FileNotFoundError(f"Missing topology {pdb}")

    traj = md.load(dcd, top=pdb)
    traj.superpose(ref)
    coords = traj.xyz  # nm

    for chain_idx, chain_name in chain_map.items():
        # identify first/last residue indices for this chain
        residues = [r for r in ref.topology.residues if r.chain.index == chain_idx]
        if not residues:
            continue
        first_i = min(r.index for r in residues)
        last_i  = max(r.index for r in residues)

        # atom indices for backbone beads
        first_idxs = [a.index for a in ref.topology.atoms
                      if a.residue.index == first_i and a.name in backbone_beads]
        last_idxs  = [a.index for a in ref.topology.atoms
                      if a.residue.index == last_i  and a.name in backbone_beads]

        # compute COM → distance (Å)
        com1   = coords[:, first_idxs, :].mean(axis=1)
        com2   = coords[:, last_idxs,  :].mean(axis=1)
        e2e_nm = np.linalg.norm(com2 - com1, axis=1)
        e2e_A  = e2e_nm * 10.0

        for i in range(traj.n_frames):
            records.append({
                'Replica':         replica,
                'Temperature (K)': temp,
                'Chain':           chain_name,
                'Frame':           i,
                'Time_ps':         float(traj.time[i]),
                'Time_ns':         float(traj.time[i]) / 1000.0,
                'E2E_Å':           e2e_A[i]
            })

# ─── Save CSV ────────────────────────────────────────────────────────────────
df = pd.DataFrame(records)
out = cwd / f"{system_name}_{replica}_e2e.csv"
df.to_csv(out, index=False)
print(f"✅ Wrote {len(df)} rows to {out}")
