"""
Fig2 Discovery: 47x4 heatmap (Ward) + ranked beta_MR
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
from scipy.cluster.hierarchy import linkage, dendrogram
from setup_style import setup_style

setup_style(journal='nature', lang='en', use_sciplots=False)

R = str(RESULTS_DIR)
OUT = str(OUTPUT_FIG)

gse = pd.read_csv(os.path.join(R, "GSE245147_validation_47genes.csv"))
genes47 = sorted(gse.loc[gse.in_data==True, 'gene'].tolist())

mr = pd.read_csv(os.path.join(R, "mr_results_GTEx4tissues.csv"), sep='\t')
mr47 = mr[mr.gene.isin(genes47)].copy()

tissues = ['Whole_Blood','Muscle_Skeletal','Nerve_Tibial','Adipose_Subcutaneous']
tlabels = ['Whole blood','Muscle-skeletal','Nerve-tibial','Adipose']

piv = mr47.pivot_table(index='gene', columns='tissue', values='p_MR').reindex(columns=tissues)
piv = piv.fillna(1.0)
logp = -np.log10(piv.clip(lower=1e-300))

Z = linkage(logp.values, method='ward')
order = dendrogram(Z, no_plot=True)['leaves']
logp = logp.iloc[order]

fig = plt.figure(figsize=(7.2, 4.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.1], wspace=0.22)

def panel_label(ax, txt):
    ax.annotate(txt, xy=(0, 1), xycoords="axes fraction", xytext=(0, 4),
                textcoords="offset points", ha="left", va="bottom",
                fontsize=10, fontweight="bold")

# ---- panel a: heatmap ----
ax = fig.add_subplot(gs[0, 0])
vmax = 4.0
cmap = LinearSegmentedColormap.from_list('bg', ['#f7f7f7','#c6dbef','#6baed6','#08519c'])
im = ax.imshow(logp.values, aspect='auto', cmap=cmap, vmin=0, vmax=vmax)

hits = mr47[mr47.FDR_global < 0.05]
star_colors = {'GFPT1': '#D55E00', 'SLC25A13': '#0072B2'}
for _, r in hits.iterrows():
    gi = list(logp.index).index(r.gene)
    ti = tissues.index(r.tissue)
    ax.plot(ti, gi, marker='*', color=star_colors[r.gene], markersize=10, zorder=5,
            markeredgecolor='white', markeredgewidth=0.5)

ax.set_xticks(range(4))
ax.set_xticklabels(tlabels, rotation=25, ha='right', fontsize=6.5)
ax.set_yticks(range(len(logp.index)))
ax.set_yticklabels(logp.index, fontsize=5.0)

cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
cbar.set_label(r'$-log_{10}(P_{MR})$', fontsize=6.5)
cbar.ax.tick_params(labelsize=5.5)
cbar.set_ticks([0,1,2,3,4])

# ---- panel b: ranked beta ----
ax2 = fig.add_subplot(gs[0, 1])
d = mr47[['gene','tissue','b_MR','p_MR','FDR_global']].copy()
d['sig'] = d.FDR_global < 0.05
d = d.sort_values('b_MR').reset_index(drop=True)
ypos = np.arange(len(d))
hit_colors = {'GFPT1': '#D55E00', 'SLC25A13': '#0072B2'}
colors = ['#BDBDBD'] * len(d)
for i, row in d.iterrows():
    if row.sig:
        colors[i] = hit_colors[row.gene]
ax2.barh(ypos, d.b_MR, color=colors, height=0.8, edgecolor='none')
ax2.axvline(0, color='black', lw=0.5)
ax2.set_yticks([])
ax2.set_xlabel(r'$\beta_{MR}$ (SD IDD per SD eQTL)', fontsize=6.5)
ax2.set_ylabel('170 gene–tissue pairs', fontsize=6.5)

# label hits with arrows so staggered text stays clearly attached to its bar
hit_rows = d[d.sig].sort_values('b_MR')
prev_y = -999
for _, r in hit_rows.iterrows():
    orig_y = list(d.b_MR).index(r.b_MR)
    y = orig_y
    if abs(y - prev_y) < 3:
        y += 3
    prev_y = y
    short_tissue = r.tissue.replace('_Subcutaneous','').replace('_Skeletal','').replace('_Tibial','').replace('_',' ')
    label = f"{r.gene} · {short_tissue}"
    ax2.annotate(label, xy=(r.b_MR, orig_y), xytext=(r.b_MR + 0.02, y),
                 fontsize=5.8, va='center', color=hit_colors[r.gene], fontweight='bold',
                 arrowprops=dict(arrowstyle='-', lw=0.5, color='#999999',
                                 shrinkA=0, shrinkB=0))

# legend: global gene-level color coding shared with Figures 3-5
from matplotlib.lines import Line2D
legend_handles = [
    Line2D([0],[0], marker='s', color='w', markerfacecolor='#0072B2',
           markeredgecolor='none', markersize=6, label='SLC25A13 (Tier 1)'),
    Line2D([0],[0], marker='s', color='w', markerfacecolor='#D55E00',
           markeredgecolor='none', markersize=6, label='GFPT1 (Tier 2)'),
    Line2D([0],[0], marker='s', color='w', markerfacecolor='#BDBDBD',
           markeredgecolor='none', markersize=6, label='Other pairs'),
]
ax2.legend(handles=legend_handles, loc='lower right', fontsize=5.5,
           frameon=True, framealpha=0.9, edgecolor='#cccccc')

ax2.set_xlim(-0.45, 0.65)
ax2.tick_params(labelsize=6)
for s in ['top','right']: ax2.spines[s].set_visible(False)

fig.set_size_inches(7.2, 4.4, forward=True)
fig.set_dpi(300)
panel_label(ax, "A")
panel_label(ax2, "B")
fig.savefig(os.path.join(OUT, "Figure2_discovery.png"), dpi=300)
fig.savefig(os.path.join(OUT, "Figure2_discovery.pdf"))
fig.savefig(os.path.join(OUT, "Figure2_discovery.svg"))
print("Fig2 done")