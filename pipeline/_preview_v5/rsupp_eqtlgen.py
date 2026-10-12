"""
Supplementary Figure: eQTLGen vs GTEx WB replication scatter.
198 genes with eQTLGen whole-blood FDR < 0.05.
"""
import sys, os
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
from matplotlib.lines import Line2D
from setup_style import setup_style

setup_style(journal='nature', lang='en', constrained_layout=False, use_sciplots=False)

# ---- data ----
eq = pd.read_csv(str(RESULTS_DIR / "replication_eQTLGen_vs_GTExWB.csv"), sep='\t')

okabe = {
    'blue':   '#0072B2',
    'vermil': '#D55E00',
    'orange': '#E69F00',
    'sky':    '#56B4E9',
    'green':  '#009E73',
    'grey':   '#999999',
}

# ---- figure ----
fig, ax = plt.subplots(figsize=(7.2, 6.58))

# direction-consistent vs not
consistent = eq[eq['same_direction'] == 'Y']
opposite   = eq[eq['same_direction'] == 'N']

ax.scatter(consistent['GTExWB_b_MR'], consistent['eQTLGen_b_MR'],
           s=14, c=okabe['blue'], alpha=0.5, edgecolors='none', zorder=3)
ax.scatter(opposite['GTExWB_b_MR'], opposite['eQTLGen_b_MR'],
           s=14, c=okabe['grey'], alpha=0.5, edgecolors='none', zorder=3)

# 1:1 line
lim = 1.0
ax.plot([-lim, lim], [-lim, lim], '--', c='#666666', lw=0.7, zorder=2)
ax.axhline(0, c='#cccccc', lw=0.5, zorder=1)
ax.axvline(0, c='#cccccc', lw=0.5, zorder=1)

# highlight core genes
labels = {
    'GFPT1':    (0.04, 0.08),   # right + up
    'SLC25A13': (0.0, -0.12),   # directly below
}
for gene, color in [('GFPT1', okabe['vermil']), ('SLC25A13', okabe['blue'])]:
    row = eq[eq['gene'] == gene].iloc[0]
    ax.scatter(row['GTExWB_b_MR'], row['eQTLGen_b_MR'],
               s=60, c=color, edgecolors='black', linewidths=0.8, zorder=5)
    dx, dy = labels[gene]
    ha = 'left'
    ax.annotate(gene, (row['GTExWB_b_MR'], row['eQTLGen_b_MR']),
                xytext=(row['GTExWB_b_MR'] + dx, row['eQTLGen_b_MR'] + dy),
                fontsize=7, fontweight='bold', color=color, zorder=6, ha=ha)

# stats text
n_con = len(consistent)
n_tot = len(eq)
ax.text(0.03, 0.97, f'{n_con}/{n_tot} direction-consistent (83%)\nP = 6.4e-22 (binomial)',
        transform=ax.transAxes, fontsize=6.5, va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#cccccc', alpha=0.9))

ax.set_xlim(-0.55, 0.55)
ax.set_ylim(-1.15, 0.70)
ax.set_xlabel(r'GTEx WB $\beta_{MR}$')
ax.set_ylabel(r'eQTLGen WB $\beta_{MR}$')

# legend
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=okabe['blue'],
           markersize=5, label='Same direction (164)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=okabe['grey'],
           markersize=5, label='Opposite direction (34)'),
]
ax.legend(handles=legend_elements, loc='lower right', frameon=False, fontsize=6.5)

plt.subplots_adjust(left=0.16, right=0.97, top=0.97, bottom=0.14)

out_png = str(OUTPUT_FIG / "FigureS1_eQTLGen_replication.png")
out_pdf = out_png.replace('.png', '.pdf')
out_svg = out_png.replace('.png', '.svg')
fig.set_size_inches(7.2, 6.58, forward=True)
fig.set_dpi(300)
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf)
fig.savefig(out_svg)
print(f'Saved: {out_png}')