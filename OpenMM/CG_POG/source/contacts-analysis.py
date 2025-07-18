import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis import distances
from pathlib import Path
import matplotlib.pyplot as plt
import re
import os

# Infer system and replica name from current working directory
cwd = Path.cwd()
system = cwd.parts[-3] if 'replicas' in cwd.parts else 'UNKNOWN'
replica = cwd.parts[-1]

print(f"📁 System: {system}, Replica: {replica}")

# Setup paths
data_dir = Path('./')
out_dir  = data_dir / f"{replica}_contacts_dat"
out_dir.mkdir(exist_ok=True)

stride = 10
cutoff = 4.0   # Å

# Find PDB–DCD pairs
pdb_files = sorted(data_dir.glob('*.pdb'))
pairs = []
for pdb in pdb_files:
    dcd = pdb.with_suffix('.dcd')
    if dcd.exists():
        pairs.append((pdb, dcd))
    else:
        print(f'⚠️  no DCD for {pdb.name}, skipping')

print(f'→ will process {len(pairs)} trajectories')

# Analyze each PDB/DCD pair
for pdb, dcd in pairs:
    u = mda.Universe(pdb, dcd)
    donors = u.select_atoms('name HBG')
    acceptors = u.select_atoms('name HBP')
    records = []

    for ts in u.trajectory[::stride]:
        frame = ts.frame
        D = distances.distance_array(donors.positions, acceptors.positions, box=u.dimensions)
        found = False

        for i in range(D.shape[0]):
            for j in range(D.shape[1]):
                if (donors[i].segid != acceptors[j].segid) and (D[i, j] <= cutoff):
                    records.append((frame, donors.indices[i], acceptors.indices[j], D[i, j]))
                    found = True

        if not found:
            records.append((frame, 0, 0, 0.0))

    df = pd.DataFrame(records, columns=['Frame', 'HBG_index', 'HBP_index', 'Distance (Å)'])
    out_file = out_dir / f'{pdb.stem}_{replica}_contacts.dat'
    df.to_csv(out_file, sep='\t', index=False)
    print(f'✅ Wrote {out_file.name} ({len(df)} rows)')

# Count per frame
for contact_file in out_dir.glob(f'*_{replica}_contacts.dat'):
    df = pd.read_csv(contact_file, sep='\t')
    all_frames = sorted(df['Frame'].unique())
    real = df[df['HBG_index'] != 0]
    counts = real.groupby('Frame').size()
    counts = counts.reindex(all_frames, fill_value=0)
    count_per_frame = counts.reset_index(name='NumContacts')

    out_count = contact_file.with_name(f"{contact_file.stem}_count.dat")
    count_per_frame.to_csv(out_count, sep='\t', index=False)
    print(f'✅ Wrote {out_count.name}')
