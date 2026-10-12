"""
SuppFig S3: Expression validation in independent disc datasets.
Left: GSE70362 log2FC for 47 genes
Right: GSE124272 log2FC, colored by direction match with MR
Core genes GFPT1/SLC25A13 highlighted.
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

g1 = pd.read_csv(str(RESULTS_DIR / "ge_validation_GSE70362_47genes.csv"))
g2 = pd.read_csv(str(RESULTS_DIR / "ge_validation_GSE124272_47genes.csv"))

core = {'GFPT1': okabe['vermil'], 'SLC25A13': okabe['blue']}

g1s = g1.sort_values('log2FC').reset_index(drop=True)
g2s = g2.sort_values('log2FC').reset_index(drop=True)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 6.5),
                                gridspec_kw={'width_ratios': [1, 1]})

def panel_label(ax, txt):
    ax.annotate(txt, xy=(0, 1), xycoords="axes fraction", xytext=(0, 4),
                textcoords="offset points", ha="left", va="bottom",
                fontsize=10, fontweight="bold")

# --- GSE70362 ---
colors1 = []
for _, r in g1s.iterrows():
    if r['gene'] in core:
        colors1.append(core[r['gene']])
    elif r['FDR'] < 0.05:
        colors1.append(okabe['grey'])
    else:
        colors1.append(okabe['lgrey'])

ax1.barh(range(len(g1s)), g1s['log2FC'], color=colors1, height=0.7, zorder=3)
ax1.set_yticks(range(len(g1s)))
ax1.set_yticklabels(g1s['gene'], fontsize=4.8)
# bold core gene labels
for i, gene in enumerate(g1s['gene']):
    if gene in core:
        ax1.get_yticklabels()[i].set_fontweight('bold')
        ax1.get_yticklabels()[i].set_color(core[gene])
ax1.axvline(0, c='#666666', lw=0.5)
ax1.set_xlabel('log2(FC) degenerate vs control')

ax1.set_ylim(-0.5, len(g1s)-0.5)

# --- GSE124272 ---
colors2 = []
for _, r in g2s.iterrows():
    if r['gene'] in core:
        colors2.append(core[r['gene']])
    elif r['FDR'] < 0.05 and r['dir_match'] == 1:
        colors2.append(okabe['grey'])
    elif r['FDR'] < 0.05 and r['dir_match'] == 0:
        colors2.append(okabe['vermil'])
    else:
        colors2.append(okabe['lgrey'])

ax2.barh(range(len(g2s)), g2s['log2FC'], color=colors2, height=0.7, zorder=3)
ax2.set_yticks(range(len(g2s)))
ax2.set_yticklabels(g2s['gene'], fontsize=4.8)
for i, gene in enumerate(g2s['gene']):
    if gene in core:
        ax2.get_yticklabels()[i].set_fontweight('bold')
        ax2.get_yticklabels()[i].set_color(core[gene])
ax2.axvline(0, c='#666666', lw=0.5)
ax2.set_xlabel('log2(FC) degenerate vs control')

ax2.set_ylim(-0.5, len(g2s)-0.5)

plt.subplots_adjust(left=0.18, right=0.97, top=0.95, bottom=0.06, wspace=0.55)
panel_label(ax1, "A")
panel_label(ax2, "B")

out = str(OUTPUT_FIG / "FigureS3_expression_validation.png")
fig.set_size_inches(7.2, 5.5, forward=True)
fig.set_dpi(300)
fig.savefig(out, dpi=300)
fig.savefig(out.replace('.png', '.pdf'))
fig.savefig(out.replace('.png', '.svg'))
print(f'Saved: {out}')