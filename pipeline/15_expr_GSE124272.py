# -*- coding: utf-8 -*-
"""GSE124272: 全血 8 IDD vs 8 Healthy，47 谷氨酰胺通路基因差异验证 + 与 GTEx WB eQTL 方向对比"""
# ============================================================
# 步骤 15 · bulk 表达验证②：GSE124272（IDD 全血 2 队列版）
# ------------------------------------------------------------
# 功能: 47 谷氨酰胺基因全血 CTL vs IDD 差异表达（探索性，
#       FDR 校正策略明示：p<0.05 视为提示性信号）
# 输入: GSE124272_series_matrix + GPL21185 probe2gene
# 输出: results/ge_validation_GSE124272_47genes.csv +
#       results/ge_validation_summary_47genes.csv
# ============================================================
import gzip, os, sys
from pathlib import Path
import pandas as pd, numpy as np
from scipy import stats
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

WORK = str(P.GEO_WORK)
OUT = str(P.RESULTS_DIR / "ge_validation_GSE124272_47genes.csv")

with gzip.open(os.path.join(WORK, "GSE124272_series_matrix.txt.gz"), "rt", encoding="utf-8", errors="replace") as f:
    lines = f.read().splitlines()

# 样本：从 !Sample_title 解析分组（P1-14 修复：不再硬编码前 8 IDD / 后 8 Healthy）
titles = [x.strip('"') for x in lines[[i for i,l in enumerate(lines) if l.startswith("!Sample_title")][0]].split("\t")[1:]]
print("样本:", titles)
n = len(titles)
# 依据 GSE124272 设计：8 IDD + 8 Healthy；关键词判定分组，解析失败则报错（不静默反向）
idd_cols, ctrl_cols = [], []
for i, t in enumerate(titles):
    tl = t.upper()
    if any(k in tl for k in ["IDD", "INTERVERTEBRAL", "DEGENER", "DISC DISEASE", "DDD", "DISORDER"]):
        idd_cols.append(i)
    elif any(k in tl for k in ["CONTROL", "HEALTHY", "NORMAL"]):
        ctrl_cols.append(i)
    else:
        print(f"!! 无法分类样本 {i}: {t}")
if len(idd_cols) != 8 or len(ctrl_cols) != 8:
    raise SystemExit(f"ERROR: 分组解析失败 idd={len(idd_cols)} ctrl={len(ctrl_cols)}（需 8/8）。请人工检查 !Sample_title")
print(f"IDD {len(idd_cols)} 列 | Ctrl {len(ctrl_cols)} 列（已按 title 解析，非硬编码）")

start = next(i for i,l in enumerate(lines) if l.startswith("!series_matrix_table_begin")) + 1
expr_lines = [l for l in lines[start:] if l.strip() and not l.startswith("!")]
probe_ids = [l.split("\t")[0].strip('"') for l in expr_lines[1:]]
expr = np.array([[float(x) for x in l.split("\t")[1:]] for l in expr_lines[1:]])
print("矩阵:", expr.shape, "值域 %.2f-%.2f" % (expr.min(), expr.max()))
print("已 log2:", expr.max() < 50)

annot = pd.read_csv(os.path.join(WORK, "GPL21185_probe2gene.csv"))
annot = annot[annot["SYMBOL"].notna() & (annot["SYMBOL"] != "")]
p2g = dict(zip(annot["ID"], annot["SYMBOL"]))
print("探针-基因映射:", len(p2g))

# 统一 47 基因（MR/coloc 定义）
genes = ["ARG1","ARG2","ASL","ASNS","ASNSD1","ASS1","CAD","CPS1","CTPS1","GAD1","GCLC","GCLM",
         "GFPT1","GFPT2","GLS","GLUD1","GLUL","GMPS","GOT1","GOT2","GPT","GPT2","GSR","GSS",
         "IDH1","IDH2","NFS1","OGDH","PPAT","SLC1A5","SLC25A12","SLC25A13","SLC38A1","SLC38A2",
         "SLC38A3","SLC38A4","SLC38A7","SLC38A9","SLC7A11","SLC7A5","SLC7A6","SLC7A7","SUCLA2","SUCLG1",
         "CTPS2","SLC38A5","GLUD2"]

ctrl_mat = pd.DataFrame(expr[:, ctrl_cols], index=probe_ids)
idd_mat = pd.DataFrame(expr[:, idd_cols], index=probe_ids)

rows = []
for gene in genes:
    probes = [p for p, g in p2g.items() if g.upper() == gene.upper()]
    if not probes:
        rows.append({"gene": gene, "n_probe": 0, "log2FC": np.nan, "p": np.nan, "FDR": np.nan})
        continue
    sub_c = ctrl_mat.loc[probes]; sub_i = idd_mat.loc[probes]
    log2fc = (sub_i.mean(axis=1) - sub_c.mean(axis=1)).median()
    ps = []
    for p in probes:
        a = sub_c.loc[p].values; b = sub_i.loc[p].values
        if np.nanstd(a) == 0 and np.nanstd(b) == 0: continue
        t, pv = stats.ttest_ind(b, a, equal_var=False)
        ps.append(pv)
    pval = np.median(ps) if ps else np.nan
    rows.append({"gene": gene, "n_probe": len(probes), "log2FC": log2fc, "p": pval, "FDR": np.nan})

res = pd.DataFrame(rows)
valid = res["p"].notna()
if valid.sum() > 0:
    pvals = res.loc[valid, "p"].values
    res.loc[valid, "FDR"] = stats.false_discovery_control(pvals, method="bh")
res = res.sort_values("p", na_position="last")
res.to_csv(OUT, index=False, encoding="utf-8-sig")

# 与 GTEx WB eQTL 对比（47 基因中 WB 有名义显著 eQTL 的）
cis = pd.read_csv(P.RESULTS_DIR / "cis_snps_full_Whole_Blood.csv", sep="\t")
eqtl_top = cis.loc[cis.groupby("gene")["p"].idxmin(), ["gene", "b", "p"]].rename(
    columns={"b": "eQTL_WB_b", "p": "eQTL_WB_p"})
res = res.merge(eqtl_top, on="gene", how="left")
res["dir_match"] = np.where((res["log2FC"].notna()) & (res["eQTL_WB_b"].notna()) & (res["eQTL_WB_b"] != 0),
                            np.sign(res["log2FC"]) == np.sign(res["eQTL_WB_b"]), np.nan)
res.to_csv(OUT, index=False, encoding="utf-8-sig")

print("\n=== GSE124272 全血 47 基因（IDD 8 vs Healthy 8）===")
print(res[["gene","n_probe","log2FC","p","FDR"]].to_string(index=False))
sig = res[res["p"] < 0.05]
print(f"\np<0.05: {len(sig)} 个 | FDR<0.05: {(res['FDR']<0.05).sum()} 个")
print("p<0.05:", ", ".join(sig["gene"]) if len(sig) else "无")

print("\n=== 核心基因（与 GTEx WB eQTL 方向对比）===")
core = res[res["gene"].isin(["GFPT1","SLC25A13","SLC25A12","NFS1","GLS","GLUD1","SLC38A1","SLC7A5"])]
print(core[["gene","log2FC","p","FDR","eQTL_WB_b","eQTL_WB_p","dir_match"]].round(4).to_string(index=False))
print("\n保存:", OUT)
