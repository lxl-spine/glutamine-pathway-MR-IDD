# -*- coding: utf-8 -*-
"""GSE153066 亚群定位（marker 粗分群）：核心基因在各细胞类型中的表达（CTL vs IDD）"""
# ============================================================
# 步骤 20 · 单细胞亚群定位（marker 粗分群，探索性/方法学对照）
# ------------------------------------------------------------
# 功能: 每细胞 marker 得分粗分群 + 各亚群核心基因 CTL vs IDD
# 注意: 该粗分群信号（SLC25A13 Immune p=0.004）经正式 Seurat
#       聚类（步骤 19）样本级检验后不成立——仅供方法学对照
# 输入: GSE153066 counts
# 输出: results/scRNA_GSE153066_subtype.csv
# ============================================================
import gzip
import pandas as pd, numpy as np
from scipy import stats
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

MARKERS = {
    "NP_like": ["COL2A1", "ACAN", "SOX9"],          # 髓核样（软骨）
    "Hypertrophic": ["COL10A1", "IBSP", "SPP1"],     # 肥大软骨样
    "AF_like": ["COL1A1", "COL1A2"],                 # 纤维环样
    "Immune": ["PTPRC", "CD3D", "CD68"],             # 免疫
    "Endothelial": ["PECAM1", "VWF", "CDH5"],        # 内皮
    "Mesenchymal": ["PDGFRA", "THY1", "ENG"],        # 间充质/祖细胞
}
CORE = ["GFPT1","SLC25A13","SLC25A12","NFS1","GLS","GLUD1","SLC38A1","SLC7A5","ARG1","SLC7A11"]
need = set()
for m in MARKERS.values(): need.update(m)
need.update(CORE)
need.update(["COL2A1","ACAN","SOX9"])

with gzip.open(P.GEO_WORK / "GSE153066_AllSample.counts.tsv.gz", "rt") as f:
    header = f.readline()
    cells = header.strip().split("\t")[1:]
    rows = {}
    colsum = np.zeros(len(cells), dtype=np.float64)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        g = parts[0]
        vals = np.array([int(x) for x in parts[1:]], dtype=np.float64)
        colsum += vals
        if g in need:
            rows[g] = vals
print(f"读取 {len(rows)} 基因 x {len(cells)} 细胞")

# 每细胞 CPM 归一化
lib = colsum / 1e6
norm = {g: v / lib for g, v in rows.items()}
df = pd.DataFrame(norm)  # 细胞 x 基因

# marker 粗分群：细胞分配给 marker 得分最高的类别（平均 marker CPM，需 > 阈值）
score = pd.DataFrame({k: df[list(v)].mean(axis=1) for k, v in MARKERS.items()})
thresh = 1.0
label = score.idxmax(axis=1)
label[score.max(axis=1) < thresh] = "Unclassified"
print("\n细胞类型分布:", label.value_counts().to_dict())

# 核心基因在各亚群 CTL vs IDD（按样本聚合 CPM 后检验，避免细胞数偏差）
samples = sorted(set(c.split("_")[0] for c in cells))
sample_of = {c: c.split("_")[0] for c in cells}
df["cell"] = cells; df["sample"] = [sample_of[c] for c in cells]; df["label"] = label.values

results = []
for ct in label.unique():
    sub = df[df.label == ct]
    if len(sub) < 30: continue
    ctrl = [s for s in samples if s.startswith("CTL")]
    idd = [s for s in samples if s.startswith("IDD")]
    for g in CORE:
        if g not in sub.columns: continue
        a = sub[sub["sample"].isin(ctrl)][g].mean() if len(sub[sub["sample"].isin(ctrl)]) else np.nan
        b = sub[sub["sample"].isin(idd)][g].mean() if len(sub[sub["sample"].isin(idd)]) else np.nan
        # 按样本聚合再做检验（每样本平均）
        ps = {s: sub[sub["sample"] == s][g].mean() for s in samples if len(sub[sub["sample"] == s]) > 5}
        if len(ps) >= 8:
            ac = [ps[s] for s in ctrl if s in ps]; bi = [ps[s] for s in idd if s in ps]
            if len(ac) >= 3 and len(bi) >= 3:
                t, p = stats.ttest_ind(bi, ac, equal_var=False)
            else: p = np.nan
        else: p = np.nan
        results.append({"celltype": ct, "gene": g, "CTL": a, "IDD": b, "p": p})
res = pd.DataFrame(results)
# P0-7 修复：FDR 用 mask 赋值，避免 dropna 长度不匹配导致错位/崩溃
mask = res["p"].notna()
if mask.sum() > 0:
    res.loc[mask, "FDR"] = stats.false_discovery_control(res.loc[mask, "p"].values, method="bh")
res = res.sort_values("p")
print("\n=== 核心基因各亚群 CTL vs IDD（伪 bulk 均值）===")
print(res.round(3).to_string(index=False))
res.to_csv(P.RESULTS_DIR / "scRNA_GSE153066_subtype.csv", index=False, encoding="utf-8-sig")
print("\n已存 results\\scRNA_GSE153066_subtype.csv")
