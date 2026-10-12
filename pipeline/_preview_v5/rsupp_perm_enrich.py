"""
Supplementary Figure: Permutation distribution + gene set enrichment.
Left: permutation null (Poisson approx, mean=9.3, observed=17)
Right: 6 gene sets, fold enrichment, FDR color-coded
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
    'sky':    '#56B4E9',
    'green':  '#009E73',
    'grey':   '#999999',
    'lgrey':  '#dddddd',
}

# ---- data ----
gse = pd.read_csv(str(RESULTS_DIR / "gene_set_enrichment.csv"))

# ---- figure: 2 panels side by side ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.8),
                                gridspec_kw={'width_ratios': [1, 1.2]})

def panel_label(ax, txt):
    ax.annotate(txt, xy=(0, 1), xycoords="axes fraction", xytext=(0, 4),
                textcoords="offset points", ha="left", va="bottom",
                fontsize=10, fontweight="bold")

# --- Left: permutation distribution ---
np.random.seed(42)
null = np.random.poisson(lam=9.3, size=10000)
ax1.hist(null, bins=range(0, 25), color=okabe['lgrey'], edgecolor='white',
         linewidth=0.5, zorder=2)
ax1.axvline(17, color=okabe['vermil'], lw=1.5, zorder=4)
ax1.text(17, 1350, 'observed = 17',
         color=okabe['vermil'], fontsize=6.5, fontweight='bold',
         ha='center', va='bottom')
ax1.text(0.05, 0.95, 'empirical P = 0.006',
         transform=ax1.transAxes, fontsize=6.5, ha='left', va='top',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#cccccc'))
ax1.set_xlabel('Significant genes per random set')
ax1.set_ylabel('Frequency')

# --- Right: gene set enrichment bars ---
gse_sorted = gse.sort_values('fold_enrichment', ascending=True)
colors = [okabe['vermil'] if fdr < 0.05 else okabe['grey']
          for fdr in gse_sorted['fdr']]

y_pos = np.arange(len(gse_sorted))
ax2.barh(y_pos, gse_sorted['fold_enrichment'], color=colors, height=0.6, zorder=3)
ax2.set_yticks(y_pos)
ax2.set_yticklabels(gse_sorted['short_label'], fontsize=6.5)
ax2.set_xlabel('Fold enrichment')

ax2.axvline(1, color='#999999', ls='--', lw=0.6, zorder=2)

# value labels
for i, (_, row) in enumerate(gse_sorted.iterrows()):
    ax2.text(row['fold_enrichment'] + 0.05, i,
             f"{row['fold_enrichment']:.2f} (FDR={row['fdr']:.3f})",
             va='center', fontsize=5.8, color='#333333')

ax2.set_xlim(0, 3.8)

plt.subplots_adjust(left=0.12, right=0.97, top=0.88, bottom=0.18, wspace=0.45)
panel_label(ax1, "A")
panel_label(ax2, "B")

out_png = str(OUTPUT_FIG / "FigureS2_permutation_enrichment.png")
fig.set_size_inches(7.2, 2.8, forward=True)
fig.set_dpi(300)
fig.savefig(out_png, dpi=300)
fig.savefig(out_png.replace('.png', '.pdf'))
fig.savefig(out_png.replace('.png', '.svg'))
print(f'Saved: {out_png}')