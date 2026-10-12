"""
Fig5: Replication + GSE245147 volcano + scRNA heatmap
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
from matplotlib.colors import LinearSegmentedColormap
from setup_style import setup_style

setup_style(journal='nature', lang='en', use_sciplots=False)

R = str(RESULTS_DIR)
OUT = str(OUTPUT_FIG)

cohorts = ['FinnGen\nIDD', 'UKB\nIDD', 'UKB\nback\npain', 'UKB\ndisc\nprol.']
gfpt1_p = [4.74e-9, 0.0596, 0.1797, 0.313]
slc_p   = [2.36e-6, 0.0200, 0.0477, 0.0105]

gse = pd.read_csv(os.path.join(R, "GSE245147_validation_47genes.csv"))
gse = gse[gse.in_data == True].copy()
gse['nlog10FDR'] = -np.log10(gse.FDR_ttest.clip(lower=1e-10))
gse['sig'] = gse.FDR_ttest < 0.05

sc = pd.read_csv(os.path.join(R, "scRNA_GSE153066_subtype.csv"))
sc['lfc'] = np.log2((sc.IDD+0.1)/(sc.CTL+0.1))
pivot = sc.pivot_table(index='gene', columns='celltype', values='lfc')
cell_order = ['NP_like','AF_like','Mesenchymal','Hypertrophic','Endothelial','Immune','Unclassified']
pivot = pivot.reindex(columns=[c for c in cell_order if c in pivot.columns])

fig = plt.figure(figsize=(7.2, 3.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.0, 1.15], wspace=0.28)

def panel_label(ax, txt):
    ax.annotate(txt, xy=(0, 1), xycoords="axes fraction", xytext=(0, 4),
                textcoords="offset points", ha="left", va="bottom",
                fontsize=10, fontweight="bold")

# panel a
ax = fig.add_subplot(gs[0, 0])
x = np.arange(len(cohorts))
w = 0.35
ax.bar(x - w/2, -np.log10(gfpt1_p), w, color='#D55E00', label='GFPT1-WB', edgecolor='none')
ax.bar(x + w/2, -np.log10(slc_p), w, color='#0072B2', label='SLC25A13-WB', edgecolor='none')
ax.axhline(-np.log10(0.05), color='gray', ls='--', lw=0.6)
# P=0.05 label just above the line, near left
ax.text(0.05, -np.log10(0.05)+0.2, 'P=0.05', fontsize=5.5, color='gray')
ax.set_xticks(x)
ax.set_xticklabels(cohorts, fontsize=6)
ax.set_ylabel(r'$-log_{10}(P_{MR})$', fontsize=6.5)

ax.legend(fontsize=6, frameon=False, loc='upper right')
ax.set_ylim(0, 9)
ax.set_xlim(-0.6, 3.6)
ax.tick_params(labelsize=6)
for s in ['top','right']: ax.spines[s].set_visible(False)

# panel b
ax2 = fig.add_subplot(gs[0, 1])
ax2.scatter(gse[~gse.sig].log2FC, gse[~gse.sig].nlog10FDR, s=12, color='#BDBDBD', edgecolor='none', alpha=0.7)
ax2.scatter(gse[gse.sig].log2FC, gse[gse.sig].nlog10FDR, s=18, color='#D55E00', edgecolor='none', alpha=0.9)
for _, r in gse[gse.gene.isin(['GFPT1','SLC25A13'])].iterrows():
    c = '#D55E00' if r.FDR_ttest < 0.05 else '#888888'
    ax2.scatter(r.log2FC, r.nlog10FDR, s=70, marker='D', facecolor=c,
                edgecolor='black', linewidth=0.6, zorder=5)
    ax2.annotate(r.gene, xy=(r.log2FC, r.nlog10FDR), xytext=(6, 2),
                 textcoords='offset points', fontsize=6.5, fontweight='bold', color=c)
ax2.axhline(-np.log10(0.05), color='gray', ls='--', lw=0.6)
ax2.axvline(0, color='black', lw=0.4)
ax2.set_xlabel(r'$log_{2}(FC)$ degenerate vs control', fontsize=6.5)
ax2.set_ylabel(r'$-log_{10}(FDR)$', fontsize=6.5)

ax2.tick_params(labelsize=6)
for s in ['top','right']: ax2.spines[s].set_visible(False)

# panel c
ax3 = fig.add_subplot(gs[0, 2])
vmax = np.ceil(np.nanmax(np.abs(pivot.values)))
cmap = LinearSegmentedColormap.from_list('rdbu', ['#0072B2','#f7f7f7','#D55E00'])
im = ax3.imshow(pivot.values, aspect='auto', cmap=cmap, vmin=-vmax, vmax=vmax)
ax3.set_xticks(range(len(pivot.columns)))
ax3.set_xticklabels([c.replace('_','-') for c in pivot.columns], rotation=35, ha='right', fontsize=5.2)
ax3.set_yticks(range(len(pivot.index)))
ax3.set_yticklabels(pivot.index, fontsize=6)

for i, g in enumerate(pivot.index):
    if g == 'GFPT1':
        ax3.get_yticklabels()[i].set_color('#D55E00')
        ax3.get_yticklabels()[i].set_fontweight('bold')
    elif g == 'SLC25A13':
        ax3.get_yticklabels()[i].set_color('#0072B2')
        ax3.get_yticklabels()[i].set_fontweight('bold')
cbar = fig.colorbar(im, ax=ax3, fraction=0.05, pad=0.04)
cbar.set_label(r'$log_{2}(IDD/CTL)$', fontsize=5.5)
cbar.ax.tick_params(labelsize=5)

fig.set_size_inches(7.2, 3.4, forward=True)
fig.set_dpi(300)
panel_label(ax, "A")
panel_label(ax2, "B")
panel_label(ax3, "C")
fig.savefig(os.path.join(OUT, "Figure3_crossoutcome.png"), dpi=300)
fig.savefig(os.path.join(OUT, "Figure3_crossoutcome.pdf"))
fig.savefig(os.path.join(OUT, "Figure3_crossoutcome.svg"))
print("Fig3 done")