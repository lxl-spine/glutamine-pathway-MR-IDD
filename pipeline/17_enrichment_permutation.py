# -*- coding: utf-8 -*-
"""富集置换检验 v2：10000 次置换 + 工具变量强度匹配 + 全局 FDR
改进（2026-09-22）：
 1) 置换次数 1000 -> 10000（更稳的经验 p）
 2) 零假设基因与 47 目标基因匹配"lead eQTL 工具强度"分布（按最显著组织 p_eQTL 分位数分层抽样）
 3) 新增全局 BH-FDR（76,953 行合并）下 47 基因显著计数
"""
import pandas as pd, numpy as np, json
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

mr = pd.read_csv(P.RESULTS_DIR / "mr_results_GTEx4tissues.csv", sep="\t")
mr.columns = [c.strip() for c in mr.columns]

genes = ["ARG1","ARG2","ASL","ASNS","ASNSD1","ASS1","CAD","CPS1","CTPS1","GAD1","GCLC","GCLM",
         "GFPT1","GFPT2","GLS","GLUD1","GLUL","GMPS","GOT1","GOT2","GPT","GPT2","GSR","GSS",
         "IDH1","IDH2","NFS1","OGDH","PPAT","SLC1A5","SLC25A12","SLC25A13","SLC38A1","SLC38A2",
         "SLC38A3","SLC38A4","SLC38A7","SLC38A9","SLC7A11","SLC7A5","SLC7A6","SLC7A7","SUCLA2","SUCLG1",
         "CTPS2","SLC38A5","GLUD2"]
all_genes = mr.gene.unique()
bg_genes = np.array([g for g in all_genes if g not in genes])
print(f"背景基因: {len(bg_genes)}（剔除 47 目标）")

# ---- 工具变量强度：每基因最显著组织的 p_eQTL ----
def lead_strength(g):
    gg = mr[mr.gene == g]
    return gg["p_eQTL"].min()

strength = pd.Series({g: lead_strength(g) for g in bg_genes})
target_strength = pd.Series({g: lead_strength(g) for g in genes})
# 分位数分层（10 层）
qs = strength.quantile(np.linspace(0, 1, 11)).values
def layer(p):
    for i in range(10):
        if p <= qs[i + 1] + 1e-300:
            return i
    return 9
strength_layer = strength.apply(layer)
target_layer = target_strength.apply(layer)
layer_counts = target_layer.value_counts().sort_index()
print("目标基因按工具强度分层:", dict(layer_counts))
bg_by_layer = {i: bg_genes[strength_layer.values == i] for i in range(10)}
for i in range(10):
    if layer_counts.get(i, 0) > 0 and len(bg_by_layer[i]) == 0:
        print(f"  !! 层 {i} 无背景基因，抽样将失败")

# ---- 观察值 ----
sub = mr[mr.gene.isin(genes)]
obs_genes = (sub.groupby("gene")["p_MR"].min() < 0.05).sum()
obs_pairs_total = (sub["p_MR"] < 0.05).sum()
# 全局 FDR（76,953 行合并 BH——修复：accumulate min 须在排序空间）
p_all = mr["p_MR"].values
m = len(p_all)
order = np.argsort(p_all)
q_sorted = p_all[order] * m / np.arange(1, m + 1)
q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
q_global = np.empty(m)
q_global[order] = q_sorted
mr["FDR_global"] = q_global
sub = mr[mr.gene.isin(genes)]
obs_genes_global = (sub.groupby("gene")["FDR_global"].min() < 0.05).sum()
obs_pairs_global = (sub["FDR_global"] < 0.05).sum()
print(f"观察(层内 p<0.05): 显著基因 {obs_genes}/47, gene×tissue {obs_pairs_total}")
print(f"观察(全局 FDR<0.05): 显著基因 {obs_genes_global}/47, gene×tissue {obs_pairs_global}")

# ---- 10000 次特征匹配置换 ----
rng = np.random.default_rng(42)
N_PERM = 10000
perm_genes = np.zeros(N_PERM); perm_pairs = np.zeros(N_PERM)
for i in range(N_PERM):
    chosen = []
    for layer_i, n in layer_counts.items():
        pool = bg_by_layer[layer_i]
        chosen.extend(rng.choice(pool, size=int(n), replace=False))
    chosen = np.array(chosen)
    gg = mr[mr.gene.isin(chosen)]
    perm_genes[i] = (gg.groupby("gene")["p_MR"].min() < 0.05).sum()
    perm_pairs[i] = (gg["p_MR"] < 0.05).sum()
    if (i + 1) % 2000 == 0:
        print(f"  置换 {i+1}/{N_PERM} ... 当前 emp_p_genes={(perm_genes[:i+1] >= obs_genes).mean():.4f}")

p_emp_genes = (perm_genes >= obs_genes).mean()
p_emp_pairs = (perm_pairs >= obs_pairs_total).mean()
print(f"\n置换(特征匹配, N={N_PERM}): emp_p_genes={p_emp_genes} (perm均值 {perm_genes.mean():.2f}, 观察 {obs_genes})")
print(f"                    emp_p_pairs={p_emp_pairs} (perm均值 {perm_pairs.mean():.2f}, 观察 {obs_pairs_total})")

# ---- 多基因集（沿用） ----
gene_sets = {
    "KEGG_D-Glutamine_and_D-glutamate_metabolism": ["GLS","GLS2","GLUD1","GLUD2","GOT1","GOT2","GPT","GPT2","GLUL"],
    "KEGG_Alanine_aspartate_glutamate_metabolism": ["GOT1","GOT2","GPT","GPT2","GAD1","GAD2","GLUD1","GLUD2","ASNS","ASNSD1","ASL","ASS1","CAD","CTPS1","CTPS2","GLS","GLS2","GLUL","PPAT","NFS1","CPS1","SLC25A12","SLC25A13","SLC1A5","SLC38A3","SLC38A4"],
    "Reactome_Glutamate_neurotransmitter_release": ["SLC1A5","SLC7A5","SLC38A1","SLC38A2","SLC38A3","SLC38A4","SLC38A7","SLC38A9","SLC7A6","SLC7A7","SLC7A11","GLS","GLUD1","GOT1","GOT2","GLUL"],
    "GO_Glutamine_metabolic_process": ["GLS","GLS2","GLUL","GLUD1","GLUD2","GPT","GPT2","GOT1","GOT2","ASNS","ASNSD1","CAD","CTPS1","CTPS2","GFPT1","GFPT2"],
    "GO_Glutamate_metabolic_process": ["GLS","GLS2","GLUD1","GLUD2","GAD1","GAD2","GPT","GPT2","GOT1","GOT2","ASNS","ASNSD1","GLUL","GCLM","GCLC","GSR","GSS"],
    "自定义_谷氨酰胺转运与线粒体": ["SLC1A5","SLC38A1","SLC38A2","SLC38A3","SLC38A4","SLC38A5","SLC38A7","SLC38A9","SLC25A12","SLC25A13","SLC7A5","SLC7A6","SLC7A7","SLC7A11"],
}
bg_rate = (mr.groupby("gene")["p_MR"].min() < 0.05).mean()
set_results = []
for name, gs in gene_sets.items():
    gs_avail = [g for g in gs if g in all_genes]
    gs_sig = mr[(mr.gene.isin(gs_avail)) & (mr["p_MR"] < 0.05)].gene.nunique()
    expect = bg_rate * len(gs_avail)
    set_results.append({"基因集": name, "基因数": len(gs_avail), "显著基因": int(gs_sig),
                        "富集倍数": round(gs_sig / expect, 2) if expect > 0 else np.nan,
                        "期望": round(expect, 2)})

json.dump({"n_perm": N_PERM, "feature_matched": True,
           "emp_p_genes": float(p_emp_genes), "emp_p_pairs": float(p_emp_pairs),
           "obs_genes": int(obs_genes), "obs_pairs": int(obs_pairs_total),
           "perm_mean_genes": float(perm_genes.mean()), "perm_mean_pairs": float(perm_pairs.mean()),
           "obs_genes_globalFDR": int(obs_genes_global), "obs_pairs_globalFDR": int(obs_pairs_global),
           "gene_sets": set_results},
          open(P.GEO_WORK / "enrichment_permutation_results.json", "w", encoding="utf-8"), ensure_ascii=False, indent=2)
mr.to_csv(P.RESULTS_DIR / "mr_results_GTEx4tissues.csv", sep="\t", index=False, encoding="utf-8-sig")
print("\n已存 enrichment_permutation_results.json（含全局 FDR 与特征匹配）")
