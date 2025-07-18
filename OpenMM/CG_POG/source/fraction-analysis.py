# fraction-analysis.py

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit
import os

# ─── DETERMINE SYSTEM NAME FROM CURRENT DIRECTORY ─────────────────────────────
script_path = Path(__file__).resolve()
system_dir = script_path.parent  # This should be ./POG12
system_match = re.search(r'POG\d+', system_dir.name)

if not system_match:
    print("❌ Could not determine system name from script path.")
    exit(1)

system_name = system_match.group(0)

# ─── SET max_pairs BASED ON SYSTEM ────────────────────────────────────────────
max_pairs_dict = {
    'POG8': 22,
    'POG10': 28,
    'POG12': 34,
    'POG14': 40,
    'POG16': 46,
    'POG18': 52,
    'POG20': 58
}
max_pairs = max_pairs_dict.get(system_name)
if max_pairs is None:
    print(f"❌ No max_pairs defined for {system_name}")
    exit(1)

print(f"📦 System: {system_name}, using max_pairs = {max_pairs}")

# ─── USER SETTINGS ────────────────────────────────────────────────────────────
replicas_dir = Path('./replicas')
start_frame  = 9703
threshold    = 350

# ─── COLLECT REPLICA FRACTIONS ────────────────────────────────────────────────
fractions_by_temp = {}

for rep in sorted(replicas_dir.glob('replica*')):
    count_dir = rep / f"{rep.name}_contacts_dat"
    if not count_dir.exists():
        print(f"⚠️  missing {count_dir}, skipping")
        continue

    for f in sorted(count_dir.glob('*_count.dat')):
        df = pd.read_csv(f, sep='\t')
        df2 = df[df['Frame'] >= start_frame]
        if df2.empty:
            continue

        frac = df2['NumContacts'].mean() / max_pairs
        m = re.search(r'(\d+)K', f.stem)
        if not m:
            continue
        temp = int(m.group(1))

        fractions_by_temp.setdefault(temp, []).append(frac)

# ─── BUILD SUMMARY DATAFRAME ──────────────────────────────────────────────────
temps       = sorted(fractions_by_temp)
mean_fracts = [np.mean(fractions_by_temp[T]) for T in temps]
std_fracts  = [np.std (fractions_by_temp[T], ddof=1) for T in temps]
n_reps      = [len(fractions_by_temp[T]) for T in temps]
se_fracts   = [std_fracts[i] / np.sqrt(n_reps[i]) for i in range(len(temps))]

df_frac = pd.DataFrame({
    'Temperature (K)': temps,
    'Fraction':         mean_fracts,
    'StdDev':           std_fracts,
    'StdErr':           se_fracts
})

# ─── FIT MELTING CURVE ────────────────────────────────────────────────────────
T_data = df_frac['Temperature (K)'].values
F_data = df_frac['Fraction'].values
E_data = df_frac['StdErr'].values

extra_temps  = np.arange(threshold + 50, 501, 50)
extra_fracts = np.zeros_like(extra_temps, dtype=float)
extra_errs   = np.zeros_like(extra_temps, dtype=float)

T_all = np.concatenate([T_data, extra_temps])
F_all = np.concatenate([F_data, extra_fracts])
E_all = np.concatenate([E_data, extra_errs])

idx    = np.argsort(T_all)
T_plot = T_all[idx]
F_plot = F_all[idx]
E_plot = E_all[idx]

def melting_curve_top(T, Tm, k, top):
    return top / (1.0 + np.exp((T - Tm) / k))

p0 = [np.median(T_data), 5.0, np.max(F_data)]
popt, _ = curve_fit(melting_curve_top, T_data, F_data, p0=p0, maxfev=3000)
Tm_fit, k_fit, top_fit = popt
print(f"✅ Fitted Tm = {Tm_fit:.2f} K, k = {k_fit:.2f}, top = {top_fit:.2f}")

T_fit = np.linspace(T_plot.min(), T_plot.max(), 300)
F_fit = melting_curve_top(T_fit, Tm_fit, k_fit, top_fit)

# ─── PREPARE & SAVE RESULTS ───────────────────────────────────────────────────
results_df = df_frac.copy()
results_df['Tm_fit']   = Tm_fit
results_df['k_fit']    = k_fit
results_df['top_frac'] = top_fit

output_path = Path(f'./{system_name}_replicas_fraction_summary.csv')
results_df.to_csv(output_path, index=False)
print(f"📄 Saved summary to {output_path}")
