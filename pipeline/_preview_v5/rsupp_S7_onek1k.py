"""
SuppFig S7: OneK1K cell-type MR + Power distribution (proper histogram).
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

# --- OneK1K ---
ok1k = pd.read_csv(str(RESULTS_DIR / "onek1k_scMR_IDD_partial.csv"))

# --- Power ---
pow_ = pd.read_csv(str(RESULTS_DIR / "power_per_gene_all.csv"), sep='\t')

# stacked layout: A (43-row barh) tall, B (histogram) wide & short
fig = plt.figure(figsize=(7.2, 7.0))
gs = fig.add_gridspec(2, 1, height_ratios=[4.8, 2.2], hspace=0.32)
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

def panel_label(ax, txt):
    ax.annotate(txt, xy=(0, 1), xycoords="axes fraction", xytext=(0, 4),
                textcoords="offset points", ha="left", va="bottom",
                fontsize=10, fontweight="bold")

# OneK1K: bar of -log10(P) by cell type, colored by gene
ok1k_sorted = ok1k.sort_values('p_mr')
colors = [okabe['vermil'] if g == 'GFPT1' else okabe['blue'] if g == 'SLC25A13'
          else okabe['grey'] for g in ok1k_sorted['gene']]

bars = ax1.barh(range(len(ok1k_sorted)), -np.log10(ok1k_sorted['p_mr']),
                color=colors, height=0.7, zorder=3)
# emphasise nominal-significant pairs (here: SLC7A7) with a dark edge + star
sig = (ok1k_sorted['p_mr'] < 0.05).values
for b, s, pv in zip(bars, sig, ok1k_sorted['p_mr']):
    if s:
        b.set_edgecolor('#111111')
        b.set_linewidth(0.9)
        xend = b.get_width()
        ax1.annotate('*', xy=(xend, b.get_y() + b.get_height() / 2),
                     xytext=(2.5, 0), textcoords='offset points',
                     ha='left', va='center', fontsize=7, color='#111111')
ax1.axvline(-np.log10(0.05), color='#666666', ls='--', lw=0.6)
ax1.set_yticks(range(len(ok1k_sorted)))
ax1.set_yticklabels([f"{r.celltype} ({r.gene})" for r in ok1k_sorted.itertuples()],
                    fontsize=5)
ax1.set_xlabel(r'$-\log_{10}(P)$ cell-type MR')
ax1.set_xlim(0, 1.6)

# Power: histogram of detect_OR (clipped to <5), bicolor by 1.3 threshold
dOR = pow_['detect_OR'].clip(upper=5)
n, bins, patches = ax2.hist(dOR, bins=40, color=okabe['lgrey'],
                            edgecolor='white', linewidth=0.3, zorder=3)
for p, b in zip(patches, bins):
    if b <= 1.3:
        p.set_facecolor(okabe['blue'])
med = pow_['detect_OR'].median()
n_le = int((pow_['detect_OR'] <= 1.3).sum())
ax2.axvline(med, color=okabe['vermil'], lw=1.2)
ax2.text(med+0.08, ax2.get_ylim()[1]*0.88, f'median OR={med:.2f}',
         color=okabe['vermil'], fontsize=6.5)
ax2.axvline(1.3, color=okabe['blue'], lw=1, ls='--')
ax2.text(1.32, ax2.get_ylim()[1]*0.68, 'OR=1.3',
         color=okabe['blue'], fontsize=6)
ax2.set_xlabel('Detectable OR (80% power)')
ax2.set_ylabel('Pairs')
from matplotlib.patches import Patch
ax2.legend(handles=[
    Patch(facecolor=okabe['blue'], edgecolor='white', label=f'OR \u2264 1.3 ({n_le:,} pairs)'),
    Patch(facecolor=okabe['lgrey'], edgecolor='white', label=f'OR > 1.3 ({len(pow_)-n_le:,} pairs)'),
], loc='upper right', frameon=False, fontsize=6.5)

ax2.set_xlim(1.0, 5.0)

plt.subplots_adjust(left=0.14, right=0.97, top=0.96, bottom=0.07)
# B keeps the same left edge as A for a clean, aligned stack
pos = ax2.get_position()
ax2.set_position([0.14, pos.y0, 0.83, pos.height])
panel_label(ax1, "A")
panel_label(ax2, "B")

out = str(OUTPUT_FIG / "FigureS7_onek1k_power.png")
fig.set_size_inches(7.2, 7.0, forward=True)
fig.set_dpi(300)
fig.savefig(out, dpi=300)
fig.savefig(out.replace('.png', '.pdf'))
fig.savefig(out.replace('.png', '.svg'))
print(f'Saved: {out}')