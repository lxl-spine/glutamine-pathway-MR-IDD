"""
SuppFig S5: Reverse MR (IDD -> expression).
51 estimable pairs; shows direction is mixed and underpowered.
"""
import sys
try:
    from pathlib import Path
    import sys, os
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # <ROOT>/pipeline
    from config.paths import ROOT, RESULTS_DIR
    OUTPUT_FIG = ROOT / "figures" / "final"
    SCIPILOT_SCRIPTS = os.environ.get("SCIPILOT_SKILL_SCRIPTS")
    if SCIPILOT_SCRIPTS:
        sys.path.insert(0, SCIPILOT_SCRIPTS)
except ImportError as exc:
    raise SystemExit(f"[figure] 无法导入 pipeline/config/paths.py: {exc}; 如使用 scipilot 样式，请设置环境变量 SCIPILOT_SKILL_SCRIPTS")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless rendering (portable; no Tk needed)
import matplotlib.pyplot as plt
from setup_style import setup_style

setup_style(journal='nature', lang='en', constrained_layout=False, use_sciplots=False)

okabe = {
    'blue':   '#0072B2',
    'vermil': '#D55E00',
    'grey':   '#999999',
    'lgrey':  '#dddddd',
}

rmr = pd.read_csv(str(RESULTS_DIR / "reverse_mr_IDD_to_expr.csv"), sep='\t')

# forward vs reverse beta scatter
fig, ax = plt.subplots(figsize=(7.2, 6.58))

same = rmr[rmr['same_direction'] == 'Y']
opp  = rmr[rmr['same_direction'] == 'N']

ax.scatter(same['eQTL_b'], same['reverse_b'],
           s=14, c=okabe['blue'], alpha=0.6, edgecolors='none', label=f'Same direction (n={len(same)})')
ax.scatter(opp['eQTL_b'], opp['reverse_b'],
           s=14, c=okabe['vermil'], alpha=0.6, edgecolors='none', label=f'Opposite (n={len(opp)})')

ax.axhline(0, c='#cccccc', lw=0.5)
ax.axvline(0, c='#cccccc', lw=0.5)

# highlight core genes
for gene, color in [('GFPT1', okabe['vermil']), ('SLC25A13', okabe['blue'])]:
    for _, row in rmr[(rmr['gene']==gene) & (rmr['tissue']=='Whole_Blood')].iterrows():
        ax.scatter(row['eQTL_b'], row['reverse_b'], s=50, c=color,
                   edgecolors='black', linewidths=0.8, zorder=5)
        ax.annotate(gene, (row['eQTL_b'], row['reverse_b']),
                    xytext=(row['eQTL_b']-0.12, row['reverse_b']-0.3),
                    fontsize=6.5, fontweight='bold', color=color, ha='left')

ax.set_xlim(-0.12, 0.30)
ax.set_ylim(-1.8, 6.5)
ax.set_xlabel(r'Forward $\beta_{MR}$ (expression -> IDD)')
ax.set_ylabel(r'Reverse $\beta_{MR}$ (IDD -> expression)')
ax.legend(frameon=False, fontsize=6, loc='upper left')

plt.subplots_adjust(left=0.18, right=0.97, top=0.95, bottom=0.15)

out = str(OUTPUT_FIG / "FigureS5_reverse_MR.png")
fig.set_size_inches(7.2, 6.58, forward=True)
fig.set_dpi(300)
fig.savefig(out, dpi=300)
fig.savefig(out.replace('.png', '.pdf'))
fig.savefig(out.replace('.png', '.svg'))
print(f'Saved: {out}')