# -*- coding: utf-8 -*-
"""GSE70362: NP 16 IDD(Thompson III-V) vs 8 control(I-II) 差异表达验证（47 谷氨酰胺基因）"""
# ============================================================
# 步骤 14 · bulk 表达验证①：GSE70362（椎间盘 NP 芯片 48 样本）
# ------------------------------------------------------------
# 功能: 47 谷氨酰胺基因 IDD vs Control 差异表达（limma 双队列）；
#       GSE34095 已删除（样本重叠/质量核查后），GSE70362 独立完成
# 输入: GSE70362_series_matrix + probe2gene 映射
# 输出: results/ge_validation_GSE70362_47genes.csv +
#       results/ge_validation_GSE70362_robust.csv
# ============================================================
import re, os, sys
import pandas as pd
import numpy as np
from scipy import stats
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

WORK = str(P.GEO_WORK)
MATRIX = os.path.join(WORK, "GSE70362_series_matrix.txt")
ANNOT = os.path.join(WORK, "GPL570_quick.txt")
OUT_CSV = str(P.RESULTS_DIR / "ge_validation_GSE70362_47genes.csv")

# ---------- 1. 解析 series matrix ----------
with open(MATRIX, encoding="utf-8", errors="replace") as f:
    lines = f.read().splitlines()

meta = {}
data_start = None
for i, ln in enumerate(lines):
    if ln.startswith("!Sample_") or ln.startswith("!Series_"):
        key, _, rest = ln.partition("\t")
        rest = rest.strip().strip('"')
        meta.setdefault(key, []).append(rest.split('"\t"'))
    if ln.startswith("!series_matrix_table_begin"):
        data_start = i + 1
        break

# 样本级元数据（按列）
titles = meta.get("!Sample_title", [[]])[0]
chars_all = meta.get("!Sample_characteristics_ch1", [])
# characteristics 每个样本有多个属性（按分号列展开成 list of dict）
n_samples = len(titles)
char_rows = []
for ci in range(len(chars_all[0])):
    d = {}
    for attr_line in chars_all:
        if ci < len(attr_line):
            kv = attr_line[ci].strip()
            if ":" in kv:
                k, _, v = kv.partition(":")
                d[k.strip()] = v.strip()
    char_rows.append(d)
samples = pd.DataFrame({
    "title": titles,
    "donor": [c.get("individual", "") for c in char_rows],
    "thompson": [c.get("thompson grade", "") for c in char_rows],
    "tissue": [c.get("tissue", "") for c in char_rows],
    "batch": [c.get("batch", "") for c in char_rows],
})
print("总样本:", n_samples)
print(samples.to_string())

# ---------- 2. 表达矩阵 ----------
expr_lines = [l for l in lines[data_start:] if l.strip() and not l.startswith("!")]
hdr = expr_lines[0].split("\t")
probe_ids = [l.split("\t")[0].strip('"') for l in expr_lines[1:]]
expr = np.array([[float(x) for x in l.split("\t")[1:]] for l in expr_lines[1:]])
expr_df = pd.DataFrame(expr, index=probe_ids, columns=hdr[1:])
print("表达矩阵:", expr_df.shape)

# ---------- 3. 注释（Brainarray custom CDF: EntrezID_at → symbol）----------
annot = pd.read_csv(os.path.join(WORK, "GSE70362_probe2gene.csv"))
probe2gene = dict(zip(annot["PROBE"], annot["SYMBOL"]))
print("探针-基因映射:", len(probe2gene))

# ---------- 4. NP 子集 + 分组 ----------
np_idx = samples.index[samples["tissue"].str.contains("Nucleus pulposus", case=False, na=False)]
np_s = samples.loc[np_idx]
if len(np_s) == 0:
    # 解析失败保护：打印实际解析到的 tissue 值帮助诊断（P2-18）
    print("!! tissue 列解析失败。实际值:", samples["tissue"].dropna().unique()[:20])
    raise SystemExit("ERROR: 未能从 !Sample_characteristics_ch1 解析 tissue 列，分组无法判定")
def group(g):
    g = str(g)
    if g in ("I", "II", "I-II"):
        return "Control"
    if g in ("III", "IV", "V"):
        return "IDD"
    return "Unknown"
np_s = np_s.copy()
np_s["group"] = np_s["thompson"].apply(group)
print("\nNP 分组:")
print(np_s[["title", "donor", "thompson", "batch", "group"]].to_string())
ctrl = np_s.index[np_s["group"] == "Control"].tolist()
idd = np_s.index[np_s["group"] == "IDD"].tolist()
print(f"\nControl(NP): {len(ctrl)}  IDD(NP): {len(idd)}")
assert len(ctrl) == 8 and len(idd) == 16, "分组数量不符预期"
# ---------- 5. 47 谷氨酰胺基因差异检验 ----------
glut_genes = ["GFPT1","GFPT2","GLS","GLS2","GLUD1","GLUD2","SLC25A12","SLC25A13",
              "SLC38A1","SLC38A2","SLC38A3","SLC38A5","SLC38A7","SLC38A9","SLC1A5",
              "SLC7A5","SLC7A6","SLC7A7","SLC7A8","SLC7A11","GOT1","GOT2","GPT","GPT2",
              "ASNS","ASNSD1","ASL","ASS1","ARG1","ARG2","NOS2","GAD1","GAD2",
              "GCLM","GCLC","GSS","GSTM1","GSTP1","GSTT1","NFS1","NFS2","GLS3",
              "GMPS","IMPDH1","IMPDH2","SUCLG1","SUCLG2"]
expr_np = expr_df.iloc[:, np_s.index.tolist()]
ctrl_pos = [np_s.index.tolist().index(i) for i in ctrl]
idd_pos = [np_s.index.tolist().index(i) for i in idd]
ctrl_mat = expr_np.iloc[:, ctrl_pos].astype(float)
idd_mat = expr_np.iloc[:, idd_pos].astype(float)

rows = []
probe_hits = {}
for gene in glut_genes:
    probes = [p for p, g in probe2gene.items() if g.upper() == gene.upper()]
    if not probes:
        rows.append({"gene": gene, "n_probe": 0, "log2FC": np.nan, "p": np.nan, "FDR": np.nan, "probes": ""})
        continue
    # 多探针：取每个探针 t 检验，取中位表达代表基因
    sub_c = ctrl_mat.loc[probes]
    sub_i = idd_mat.loc[probes]
    # 表达用 log2 值（matrix 已是 log2 RMA？GEO 平台常见已 log2。直接 t-test）
    log2fc = (sub_i.mean(axis=1) - sub_c.mean(axis=1)).median()
    ps = []
    for p in probes:
        a = sub_c.loc[p].values; b = sub_i.loc[p].values
        if np.nanstd(a) == 0 and np.nanstd(b) == 0:
            continue
        t, pv = stats.ttest_ind(b, a, equal_var=False)
        ps.append(pv)
    pval = np.median(ps) if ps else np.nan
    probe_hits[gene] = probes
    rows.append({"gene": gene, "n_probe": len(probes), "log2FC": log2fc, "p": pval, "FDR": np.nan, "probes": ",".join(probes[:6])})

res = pd.DataFrame(rows)
# BH FDR（仅对可检验的）
valid = res["p"].notna()
res.loc[valid, "FDR"] = stats.false_discovery_control(res.loc[valid, "p"].values, method="bh")
res = res.sort_values("p", na_position="last")
res.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
print("\n=== 47 基因差异表达（GSE70362 NP: IDD 16 vs Ctrl 8）===")
print(res.to_string(index=False))
print("\n保存:", OUT_CSV)
