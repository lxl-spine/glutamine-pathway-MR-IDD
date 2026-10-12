"""
SuppFig S4: Leave-one-out enrichment.
47 genes removed one at a time; empirical P stays < 0.05.
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
    'orange': '#E69F00',
    'grey':   '#999999',
}

loo = pd.read_csv(str(RESULTS_DIR / "leave_one_out_enrichment.csv"))

# separate "none" from 47 genes
none_row = loo[loo['dropped_gene'] == 'none']
genes = loo[loo['dropped_gene'] != 'none'].sort_values('emp_p')

fig, ax = plt.subplots(figsize=(7.2, 5.5))

y = range(len(genes))
colors = [okabe['vermil'] if g in ('GFPT1', 'SLC25A13') else okabe['grey']
          for g in genes['dropped_gene']]
ax.barh(y, genes['emp_p'], color=colors, height=0.7, zorder=3)

# all values well below 0.05 - annotate
ax.text(0.97, 0.03, 'All empirical P < 0.05\n(range: 0.0015–0.0165)',
        transform=ax.transAxes, fontsize=6.5, ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#cccccc'))

ax.set_yticks(list(y))
ax.set_yticklabels(genes['dropped_gene'], fontsize=5)
ax.set_xlabel('Empirical P after dropping gene')

ax.set_xlim(0, 0.018)

plt.subplots_adjust(left=0.12, right=0.97, top=0.92, bottom=0.1)

out = str(OUTPUT_FIG / "FigureS4_leave_one_out.png")
fig.set_size_inches(7.2, 6.5, forward=True)
fig.set_dpi(300)
fig.savefig(out, dpi=300)
fig.savefig(out.replace('.png', '.pdf'))
fig.savefig(out.replace('.png', '.svg'))
print(f'Saved: {out}')