"""
SuppFig S6: GTEx tissue panel for GFPT1 and SLC25A13.
Dot plot showing lead eQTL effect size across tissues.
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
}

gtex = pd.read_csv(str(RESULTS_DIR / "gtex_alltissues_eqtl_2genes.csv"))
lead = gtex.loc[gtex.groupby(['gene', 'tissue'])['p'].idxmin()].copy()
lead['neg_log10p'] = -np.log10(lead['p'])

# simplify tissue names
lead['tissue_simple'] = lead['tissue'].str.replace('_', ' ').str.replace('Cells - ', '')

# pivot: tissues as rows, genes as columns
gfpt1 = lead[lead['gene'] == 'GFPT1'].set_index('tissue_simple')
slc = lead[lead['gene'] == 'SLC25A13'].set_index('tissue_simple')

# common tissues
common = sorted(set(gfpt1.index) & set(slc.index))
gfpt1 = gfpt1.loc[common]
slc = slc.loc[common]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 6),
                                gridspec_kw={'width_ratios': [1, 1]})

def panel_label(ax, txt):
    ax.annotate(txt, xy=(0, 1), xycoords="axes fraction", xytext=(0, 4),
                textcoords="offset points", ha="left", va="bottom",
                fontsize=10, fontweight="bold")

y = np.arange(len(common))

# GFPT1
for i, tissue in enumerate(common):
    is_wb = tissue == 'Whole Blood'
    ax1.scatter(gfpt1['nes'].iloc[i], i,
                s=gfpt1['neg_log10p'].iloc[i]*2.5,
                c=okabe['vermil'], alpha=0.8 if is_wb else 0.6,
                edgecolors='black' if is_wb else 'none',
                linewidths=1.0 if is_wb else 0, zorder=5 if is_wb else 3)
ax1.axvline(0, c='#666666', lw=0.5)
ax1.set_yticks(y)
ax1.set_yticklabels(common, fontsize=6)
for i, tissue in enumerate(common):
    if tissue == 'Whole Blood':
        ax1.get_yticklabels()[i].set_fontweight('bold')
ax1.set_xlabel('eQTL effect size (NES)')

ax1.invert_yaxis()

# SLC25A13
for i, tissue in enumerate(common):
    is_wb = tissue == 'Whole Blood'
    ax2.scatter(slc['nes'].iloc[i], i,
                s=slc['neg_log10p'].iloc[i]*2.5,
                c=okabe['blue'], alpha=0.8 if is_wb else 0.6,
                edgecolors='black' if is_wb else 'none',
                linewidths=1.0 if is_wb else 0, zorder=5 if is_wb else 3)
ax2.axvline(0, c='#666666', lw=0.5)
ax2.set_yticks(y)
ax2.set_yticklabels([])
ax2.set_xlabel('eQTL effect size (NES)')

ax2.invert_yaxis()

plt.subplots_adjust(left=0.22, right=0.97, top=0.95, bottom=0.06, wspace=0.15)
panel_label(ax1, "A")
panel_label(ax2, "B")

out = str(OUTPUT_FIG / "FigureS6_GTEx_tissues.png")
fig.set_size_inches(7.2, 6.0, forward=True)
fig.set_dpi(300)
fig.savefig(out, dpi=300)
fig.savefig(out.replace('.png', '.pdf'))
fig.savefig(out.replace('.png', '.svg'))
print(f'Saved: {out}')