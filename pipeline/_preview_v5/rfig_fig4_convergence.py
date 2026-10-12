"""
Fig3 Convergence: SMR-HEIDI quadrant + coloc PP.H4
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

setup_style(journal='nature', lang='en', use_sciplots=False)

R = str(RESULTS_DIR)
OUT = str(OUTPUT_FIG)

smr = pd.read_csv(os.path.join(R, "smr_glut_17genes_4tissues.csv"), sep='\t')
coloc = pd.read_csv(os.path.join(R, "coloc_results.csv"), sep='\t')
smr = smr.dropna(subset=['p_HEIDI']).copy()
smr['x'] = -np.log10(smr.p_SMR.clip(lower=1e-300))
smr['y'] = -np.log10(smr.p_HEIDI.clip(lower=1e-300))

tiss_colors = {
    'Whole_Blood': '#56B4E9',
    'Muscle_Skeletal': '#E69F00',
    'Nerve_Tibial': '#009E73',
    'Adipose_Subcutaneous': '#CC79A7',
}

core = smr[
    ((smr.gene=='GFPT1') & (smr.tissue=='Whole_Blood')) |
    ((smr.gene=='GFPT1') & (smr.tissue=='Nerve_Tibial')) |
    ((smr.gene=='SLC25A13') & (smr.tissue=='Whole_Blood'))
].copy()

fig = plt.figure(figsize=(7.2, 3.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.22)

def panel_label(ax, txt):
    ax.annotate(txt, xy=(0, 1), xycoords="axes fraction", xytext=(0, 4),
                textcoords="offset points", ha="left", va="bottom",
                fontsize=10, fontweight="bold")
ax = fig.add_subplot(gs[0, 0])

ax.axhspan(0, 1.3, xmin=0.15, xmax=1.0, color='#D5E8D4', alpha=0.45, zorder=0, edgecolor='none')
ax.axhspan(1.3, 4.5, xmin=0.15, xmax=1.0, color='#F8CECC', alpha=0.35, zorder=0, edgecolor='none')
ax.axvline(1.3, color='gray', ls='--', lw=0.6, zorder=1)
ax.axhline(1.3, color='gray', ls='--', lw=0.6, zorder=1)

for tiss, sub in smr.groupby('tissue'):
    ax.scatter(sub.x, sub.y, s=12, color=tiss_colors[tiss],
               alpha=0.45, edgecolor='none', zorder=2)

for _, r in core.iterrows():
    c = '#0072B2' if r.gene=='SLC25A13' else '#D55E00'
    ax.scatter(r.x, r.y, s=100, color=c, edgecolor='black', linewidth=0.8, zorder=5)

# labels - positioned clear of dashed lines, data points and quadrant text (no white boxes)
# SLC25A13-WB at (5.53, 0.36) -> label at (4.2, 1.0), left of the 'colocalized' text (x~6.9-8.0)
ax.annotate('SLC25A13 (WB)', xy=(5.53, 0.36), xytext=(4.2, 1.0),
            fontsize=6.5, fontweight='bold', color='#0072B2',
            arrowprops=dict(arrowstyle='-', lw=0.5, color='#0072B2',
                            shrinkA=14, shrinkB=3))
# GFPT1-WB at (7.03, 1.77) -> label directly above, clear of the fail text
ax.annotate('GFPT1 (WB)', xy=(7.03, 1.77), xytext=(7.35, 2.5),
            fontsize=6.5, fontweight='bold', color='#D55E00',
            arrowprops=dict(arrowstyle='-', lw=0.5, color='#D55E00',
                            shrinkA=14, shrinkB=3))
# GFPT1-Nerve at (5.18, 1.68), a HEIDI-fail signal -> label high in the fail quadrant,
# clear of the GFPT1-Adipose point (3.85, 3.59), the GFPT1 (WB) label and the linkage text
ax.annotate('GFPT1 (Nerve)', xy=(5.18, 1.68), xytext=(5.3, 2.8),
            fontsize=6.0, color='#D55E00',
            arrowprops=dict(arrowstyle='-', lw=0.5, color='#D55E00',
                            shrinkA=14, shrinkB=3))

ax.set_xlabel(r'SMR $-log_{10}(P)$', fontsize=7)
ax.set_ylabel(r'HEIDI $-log_{10}(P)$', fontsize=7)
ax.set_xlim(-0.3, 8.5)
ax.set_ylim(0, 4.5)

# quadrant labels tucked into the lower corners (kept short, never crowd gene labels)
ax.text(0.15, 0.9, 'colocalized', fontsize=5, ha='left',
        color='#2d6a2d', style='italic')
ax.text(0.15, 3.3, 'linkage', fontsize=5, ha='left',
        color='#a33', style='italic')

legend_elems = [Line2D([0],[0], marker='o', color='w', markerfacecolor=c,
                       markeredgecolor='none',
                       markersize=6, alpha=0.7, label=t.replace('_',' '))
               for t,c in tiss_colors.items()]
ax.legend(handles=legend_elems, fontsize=5.5, frameon=False, ncol=4,
          loc='upper center', bbox_to_anchor=(0.5, -0.22), handletextpad=0.4,
          columnspacing=1.0)

# panel b: coloc.abf vs coloc.susie dumbbell
ax2 = fig.add_subplot(gs[0, 1])
show = coloc[coloc.gene.isin(['GFPT1','SLC25A13'])].copy()
show = show.sort_values(['gene','tissue'])
show['label'] = show.gene + ' · ' + show.tissue.replace('WholeBlood','WB').replace('_Subcutaneous','').replace('_Skeletal','').replace('_Tibial','')

# PP.H4 under coloc.susie (multi-signal model); from Supplementary Table S4 (coloc.susie runs)
SUSIE = {
    ('GFPT1','Adipose'): 0.000002,
    ('GFPT1','Muscle'): 0.002300,
    ('GFPT1','WholeBlood'): 0.219000,
    ('SLC25A13','WholeBlood'): 0.904000,
}
susie_vals = [SUSIE.get((r.gene, r.tissue)) for _, r in show.iterrows()]

ypos = np.arange(len(show))
colors = ['#D55E00' if g=='GFPT1' else '#0072B2' for g in show.gene]

for i, (_, r) in enumerate(show.iterrows()):
    abf = r['PP.H4']
    sus = susie_vals[i]
    if sus is not None:
        ax2.plot([abf, sus], [i, i], color=colors[i], lw=0.9, alpha=0.75, zorder=2)
        ax2.scatter(abf, i, s=26, color='white', edgecolor=colors[i], linewidth=1.1, zorder=3)
        ax2.scatter(sus, i, s=26, color=colors[i], edgecolor=colors[i], linewidth=0.8, zorder=4)
    else:
        ax2.scatter(abf, i, s=26, color=colors[i], edgecolor=colors[i], linewidth=0.8, zorder=4)

ax2.axvline(0.75, color='black', ls='--', lw=0.7)

ax2.set_yticks(ypos)
ax2.set_yticklabels(show.label, fontsize=5.8)
ax2.set_xlim(0, 1.0)
ax2.set_ylim(-0.8, 7.8)
ax2.set_xlabel('Coloc posterior (PP.H4)', fontsize=6.5)

ax2.invert_yaxis()

for i, (_, r) in enumerate(show.iterrows()):
    abf = r['PP.H4']
    sus = susie_vals[i]
    if abf > 0.1:
        ax2.text(abf, i+0.42, f"{abf:.2f}", va='bottom', ha='center',
                 fontsize=5.3, color='#333333')
    if sus is not None and sus > 0.1:
        ax2.text(sus, i-0.42, f"{sus:.2f}", va='top', ha='center',
                 fontsize=5.3, color=colors[i])

legend_elems2 = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor='white',
           markeredgecolor='#333333', markersize=5, label='coloc.abf'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#333333',
           markersize=5, label='coloc.susie'),
]
ax2.legend(handles=legend_elems2, fontsize=5.5, frameon=True,
           facecolor='white', edgecolor='none', framealpha=0.9,
           loc='upper left', bbox_to_anchor=(0.02, 1.02), ncol=2)

# threshold line only (0.75 text removed to avoid crossing the dashed line)

fig.set_size_inches(7.2, 3.5, forward=True)
fig.set_dpi(300)
panel_label(ax, "A")
panel_label(ax2, "B")
fig.savefig(os.path.join(OUT, "Figure4_convergence.png"), dpi=300)
fig.savefig(os.path.join(OUT, "Figure4_convergence.pdf"))
fig.savefig(os.path.join(OUT, "Figure4_convergence.svg"))
print("Fig4 done")