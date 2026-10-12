# -*- coding: utf-8 -*-
"""23_gene_set_enrichment.py
从 mr_results_GTEx4tissues.csv 现算 6 个谷氨酰胺通路基因集的富集统计，
输出 results/gene_set_enrichment.csv（Figure 3f 的正式数据源，可复现）。
富集定义与 17_enrichment_permutation.py 完全一致：
  - 显著基因 = 层内 MR p<0.05 的基因数（gene 级 min p）
  - 期望 = bg_rate × n_genes_available，bg_rate = (全部基因 min p_MR < 0.05).mean()
  - 同时输出全局 FDR<0.05 计数供核对。
v2 (2026-09-23): 增加超几何富集 p 值 (hypergeom, 背景=全部 MR 基因)
  + BH-FDR 校正 (6 个基因集)，使 Figure 3f 能渲染显著性星号。
"""
import pandas as pd
import numpy as np
from pathlib import Path
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P
from scipy.stats import hypergeom

MR_FILE = P.RESULTS_DIR / "mr_results_GTEx4tissues.csv"
OUT_FILE = P.RESULTS_DIR / "gene_set_enrichment.csv"

mr = pd.read_csv(MR_FILE, sep="\t")
mr.columns = [c.strip() for c in mr.columns]

GENE_SETS = {
    "KEGG_D-Glutamine_and_D-glutamate_metabolism": ["GLS", "GLS2", "GLUD1", "GLUD2", "GOT1", "GOT2", "GPT", "GPT2", "GLUL"],
    "KEGG_Alanine_aspartate_glutamate_metabolism": ["GOT1", "GOT2", "GPT", "GPT2", "GAD1", "GAD2", "GLUD1", "GLUD2", "ASNS", "ASNSD1", "ASL", "ASS1", "CAD", "CTPS1", "CTPS2", "GLS", "GLS2", "GLUL", "PPAT", "NFS1", "CPS1", "SLC25A12", "SLC25A13", "SLC1A5", "SLC38A3", "SLC38A4"],
    "Reactome_Glutamate_neurotransmitter_release": ["SLC1A5", "SLC7A5", "SLC38A1", "SLC38A2", "SLC38A3", "SLC38A4", "SLC38A7", "SLC38A9", "SLC7A6", "SLC7A7", "SLC7A11", "GLS", "GLUD1", "GOT1", "GOT2", "GLUL"],
    "GO_Glutamine_metabolic_process": ["GLS", "GLS2", "GLUL", "GLUD1", "GLUD2", "GPT", "GPT2", "GOT1", "GOT2", "ASNS", "ASNSD1", "CAD", "CTPS1", "CTPS2", "GFPT1", "GFPT2"],
    "GO_Glutamate_metabolic_process": ["GLS", "GLS2", "GLUD1", "GLUD2", "GAD1", "GAD2", "GPT", "GPT2", "GOT1", "GOT2", "ASNS", "ASNSD1", "GLUL", "GCLM", "GCLC", "GSR", "GSS"],
    "custom_transport_mitochondria": ["SLC1A5", "SLC38A1", "SLC38A2", "SLC38A3", "SLC38A4", "SLC38A5", "SLC38A7", "SLC38A9", "SLC25A12", "SLC25A13", "SLC7A5", "SLC7A6", "SLC7A7", "SLC7A11"],
}
SHORT = {
    "KEGG_D-Glutamine_and_D-glutamate_metabolism": "KEGG Gln/Glu metabolism",
    "KEGG_Alanine_aspartate_glutamate_metabolism": "KEGG Ala/Asp/Glu metab.",
    "Reactome_Glutamate_neurotransmitter_release": "Reactome Glu release",
    "GO_Glutamine_metabolic_process": "GO Glutamine proc.",
    "GO_Glutamate_metabolic_process": "GO Glutamate proc.",
    "custom_transport_mitochondria": "Custom: transport/mito",
}
SOURCE = {
    "KEGG_D-Glutamine_and_D-glutamate_metabolism": "KEGG",
    "KEGG_Alanine_aspartate_glutamate_metabolism": "KEGG",
    "Reactome_Glutamate_neurotransmitter_release": "Reactome",
    "GO_Glutamine_metabolic_process": "GO",
    "GO_Glutamate_metabolic_process": "GO",
    "custom_transport_mitochondria": "Custom",
}

all_genes = set(mr.gene.unique())
bg_rate = (mr.groupby("gene")["p_MR"].min() < 0.05).mean()
N_all = len(all_genes)
K_bg = int((mr.groupby("gene")["p_MR"].min() < 0.05).sum())

rows = []
for name, gs in GENE_SETS.items():
    gs_avail = [g for g in gs if g in all_genes]
    gs_sub = mr[mr.gene.isin(gs_avail)]
    n_sig = int(gs_sub.groupby("gene")["p_MR"].min().lt(0.05).sum())   # layer-level p<0.05
    n_sig_fdr = int(gs_sub.groupby("gene")["FDR_global"].min().lt(0.05).sum())  # global FDR<0.05
    expect = bg_rate * len(gs_avail)
    fold = n_sig / expect if expect > 0 else np.nan
    # 超几何富集 p（右尾：观测到 ≥ k 个显著基因的概率）
    p_hyper = float(hypergeom.sf(n_sig - 1, N_all, K_bg, len(gs_avail))) if len(gs_avail) > 0 else np.nan
    rows.append({
        "gene_set": name,
        "short_label": SHORT[name],
        "source": SOURCE[name],
        "n_genes": len(gs_avail),
        "n_sig_genes": n_sig,
        "significant_genes_globalFDR": n_sig_fdr,
        "expected": round(expect, 2),
        "fold_enrichment": round(fold, 2),
        "p_value": p_hyper,
    })

out = pd.DataFrame(rows)
# BH-FDR across the 6 gene sets
m = len(out)
p_vals = out["p_value"].astype(float).values
order = np.argsort(p_vals)
q = np.empty(m)
q[order] = np.minimum.accumulate((p_vals[order] * m / np.arange(1, m + 1))[::-1])[::-1]
out["fdr"] = q
out.to_csv(OUT_FILE, index=False, encoding="utf-8-sig")
print(out[["short_label", "n_genes", "n_sig_genes", "expected",
           "fold_enrichment", "p_value", "fdr"]].to_string(index=False))
print(f"\nsaved {OUT_FILE}")
