"""
Fig4 Locus: SLC25A13 chr7 — GTEx WB eQTL (top) vs FinnGen IDD GWAS marginal+conditional (bottom)

v9 rewrite (2026-10-10):
  - fixed 7.2in canvas @300dpi (2160px wide), NO bbox_inches='tight'
    -> consistent physical scaling across all figures in the docx
  - unified type scale (panel label 10 / title 9 / axis label 8 / tick 7 / legend 7 / annotation 7)
  - position fixes: rs10631575 label back into the data area; legend moved off the
    lead diamond; SLC25A13 gene label placed above the truncated eQTL wall;
    secondary SNP annotation arrowed; truncation cap annotated
  - panel labels a/b at the conventional upper-left outside position
"""
import sys, os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # <ROOT>/pipeline
from config.paths import ROOT, RESULTS_DIR

OUTPUT_FIG = ROOT / "figures" / "final"
SCIPILOT_SCRIPTS = os.environ.get("SCIPILOT_SKILL_SCRIPTS")
if SCIPILOT_SCRIPTS:
    sys.path.insert(0, SCIPILOT_SCRIPTS)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless rendering (no Tk/Tcl needed; portable)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from setup_style import setup_style

setup_style(journal="nature", lang="en", constrained_layout=False, use_sciplots=False)
# unified type scale (final display size ~6.2 in -> all labels stay >= 6 pt)
plt.rcParams.update(
    {
        "font.size": 7.5,
        "axes.labelsize": 8,
        "axes.titlesize": 8.5,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "lines.linewidth": 1.0,
        "axes.linewidth": 0.6,
    }
)

R = str(RESULTS_DIR)
OUT = str(OUTPUT_FIG)

# ---------------- data ----------------
eqtl = pd.read_csv(os.path.join(R, "cis_snps_full_Whole_Blood.csv"), sep="\t")
eq = eqtl[eqtl.gene == "SLC25A13"].copy().rename(columns={"SNP_bp": "bp", "p": "p_eqtl"})
gwas = pd.read_csv(os.path.join(R, "conditional_analysis_slc25a13_region.csv"))
gwas = gwas.rename(columns={"BP": "bp", "p": "p_gwas", "p_cond": "p_gwas_cond"})

lo, hi = 95.4e6, 96.1e6
eq = eq[(eq.bp >= lo) & (eq.bp <= hi)].copy()
gwas = gwas[(gwas.bp >= lo) & (gwas.bp <= hi)].copy()

# eQTL -log10(P); truncate at 10^-12 and annotate the cap (avoids an unexplained flat wall)
CAP = 12.0
eq["nlogp"] = np.where(eq.p_eqtl > 10**-CAP, -np.log10(eq.p_eqtl), CAP)
eq["trunc"] = eq.p_eqtl <= 10**-CAP

lead_bp = 95811083  # rs4727335 (GWAS lead / conditional panel lead)

def ld_color(bp):
    d = abs(bp - lead_bp)
    if d < 50000:
        return "#08519C"
    if d < 100000:
        return "#3182BD"
    if d < 200000:
        return "#9ECAE1"
    return "#C8C8C8"

eq["c"] = [ld_color(bp) for bp in eq.bp]
gwas["c"] = [ld_color(bp) for bp in gwas.bp]

# ---------------- figure ----------------
fig, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(7.2, 4.7), sharex=True,
    gridspec_kw={"height_ratios": [1, 1.35], "hspace": 0.12},
)

fig.subplots_adjust(left=0.085, right=0.985, top=0.94, bottom=0.105, hspace=0.16)

def panel_label(ax, txt):
    ax.annotate(
        txt, xy=(0, 1), xycoords="axes fraction", xytext=(-4, 4),
        textcoords="offset points", ha="left", va="bottom",
        fontsize=10, fontweight="bold",
    )

# shared region highlight (lead LD block)
for ax in (ax1, ax2):
    ax.axvspan(95.70, 95.85, color="#f0f0f0", alpha=0.55, zorder=0)

# ---- panel a: GTEx WB eQTL ----
ax1.scatter(eq.bp / 1e6, eq.nlogp, c=eq.c, s=6, edgecolor="none", alpha=0.85, zorder=2)
ax1.axhline(CAP, color="#888888", lw=0.5, ls=(0, (3, 3)), alpha=0.6, zorder=1)
ax1.set_ylim(0, 13.0)
ax1.set_ylabel(r"GTEx WB eQTL" "\n" r"$-\log_{10}(P)$", fontsize=8)
# gene label above the capped wall (no white box; sits in the empty top margin)
ax1.text(
    0.5, 0.96, "SLC25A13", transform=ax1.transAxes, ha="center", va="top",
    fontsize=8, fontstyle="italic", color="#333333", zorder=8,
)
ax1.text(
    0.985, 0.98, "P capped at 10$^{-12}$", transform=ax1.transAxes,
    ha="right", va="top", fontsize=6.5, style="italic", color="#555555", zorder=8,
)

lead_eq = eqtl[eqtl.gene == "SLC25A13"].copy()
lead_eq["bp"] = lead_eq.SNP_bp
lead_eq["p_eqtl"] = lead_eq.p
lead_eq = lead_eq[(lead_eq.bp >= lo) & (lead_eq.bp <= hi)]
lead_eq = lead_eq[lead_eq.SNP == "rs10631575"]
if len(lead_eq):
    x0, y0 = lead_eq.bp.values[0] / 1e6, np.where(
        lead_eq.p_eqtl.values[0] > 10**-CAP,
        -np.log10(lead_eq.p_eqtl.values[0]), CAP,
    )
    ax1.scatter(x0, y0, s=44, facecolor="none", edgecolor="#111111", linewidth=1.0, zorder=5)
    ax1.annotate(
        "rs10631575", xy=(x0, y0), xytext=(10, 4), textcoords="offset points",
        fontsize=7, fontweight="bold", zorder=6,
        arrowprops=dict(arrowstyle="-", color="#111111", lw=0.6),
    )

# ---- panel b: FinnGen IDD GWAS marginal + conditional ----
ax2.scatter(
    gwas.bp / 1e6, -np.log10(gwas.p_gwas.clip(lower=1e-300)),
    c=gwas.c, s=6, edgecolor="none", alpha=0.85, zorder=2,
)
ax2.scatter(
    gwas.bp / 1e6, -np.log10(gwas.p_gwas_cond.clip(lower=1e-300)),
    facecolor="none", edgecolor="#0072B2", s=9, linewidth=0.45, alpha=0.5, zorder=3,
)
ax2.set_ylim(0, 9.5)
ax2.set_ylabel(r"FinnGen IDD GWAS" "\n" r"$-\log_{10}(P)$", fontsize=8)
ax2.set_xlabel("Chromosome 7 position (Mb)", fontsize=8)
ax2.set_xticks(np.arange(95.4, 96.11, 0.1))

lead_gw = gwas[gwas.snp == "rs4727335"]
if len(lead_gw):
    x1, y1 = lead_gw.bp.values[0] / 1e6, -np.log10(lead_gw.p_gwas.values[0])
    ax2.scatter(x1, y1, s=46, marker="D", facecolor="none", edgecolor="#111111", linewidth=1.1, zorder=6)
    ax2.annotate(
        "rs4727335 (lead)", xy=(x1, y1), xytext=(-18, 12), textcoords="offset points",
        ha="right", fontsize=7, fontweight="bold", zorder=7,
        arrowprops=dict(arrowstyle="-", color="#111111", lw=0.6),
    )

sec_gw = gwas[gwas.snp == "rs201828368"]
if len(sec_gw):
    x2, y2 = sec_gw.bp.values[0] / 1e6, -np.log10(sec_gw.p_gwas.values[0])
    ax2.scatter(x2, y2, s=38, marker="D", facecolor="#E69F00", edgecolor="#111111", linewidth=0.7, zorder=6)
    ax2.annotate(
        "rs201828368 (secondary)", xy=(x2, y2), xytext=(12, 14), textcoords="offset points",
        fontsize=7, fontweight="bold", color="#8a5a00", zorder=7,
        arrowprops=dict(arrowstyle="-", color="#8a5a00", lw=0.6),
    )

# legend: upper right (y > 6 zone is free of points; avoids covering the low-P cloud)
ax2.legend(
    handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#08519C", markersize=5.5,
               markeredgewidth=0, label="marginal association"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor="#0072B2",
               markersize=5.5, markeredgewidth=0.9, label="conditional on lead (rs4727335)"),
    ],
    loc="upper right", frameon=True, framealpha=0.92, edgecolor="#dddddd",
    borderaxespad=0.6, handletextpad=0.4, fontsize=7,
)

panel_label(ax1, "A")
panel_label(ax2, "B")

for ax in (ax1, ax2):
    ax.tick_params(labelsize=7)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)

# ---------------- export (fixed canvas, no bbox_inches='tight') ----------------
# two guards keep the raster at exactly 7.2in x 300dpi = 2160px wide:
#  1) use_sciplots=False above  (SciencePlots mis-scales savefig under matplotlib>=3.11)
#  2) lock inches+dpi right before savefig (also neutralizes the 3.11 dpi quirk)
fig.set_size_inches(7.2, 4.7, forward=True)
fig.set_dpi(300)
fig.savefig(os.path.join(OUT, "Figure5_locus.png"), dpi=300)
fig.savefig(os.path.join(OUT, "Figure5_locus.pdf"))
fig.savefig(os.path.join(OUT, "Figure5_locus.svg"))
print("Fig5 done ->", os.path.join(OUT, "Figure5_locus.png"))
