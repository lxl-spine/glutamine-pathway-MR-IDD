# -*- coding: utf-8 -*-
"""双向 MR：IDD GWAS 独立信号 → 47 基因表达（反向因果检查）
# ============================================================
# 步骤 12 · 双向 MR（SLC25A13/GFPT1 等关键基因方向性）
# ------------------------------------------------------------
# 功能: 基因表达 → IDD（正向，已由主分析覆盖）+ IDD → 表达
#       （反向可行性：FinnGen IDD 自身 cis 工具强度检查 R4）
# 输入: mr_results + FinnGen + GTEx eQTL
# 输出: results/reverse_mr_IDD_to_expr.csv +
#       results/reverse_mr_overlap_hits.csv
# ============================================================
步骤：显著 SNP + 位置 → plink clump（1000G LD）→ 独立信号 → 与 47 基因 cis 区域重叠 → 反向 Wald ratio
"""
import os, subprocess, pandas as pd, numpy as np, math
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

WORK = str(P.SMR_WORK)
OUT = str(P.RESULTS_DIR / "reverse_mr_IDD_to_expr.csv")
PLINK = str(P.PLINK)
BIM = str(P.PANEL1000G) + ".bim"
BFILE = str(P.PANEL1000G)

gwas = pd.read_csv(os.path.join(WORK, "gwas_summary.txt"), sep="\t", dtype={"snp": str})
bim = pd.read_csv(BIM, sep=r"\s+", header=None, names=["chr", "snp", "cm", "bp", "A1", "A2"])
bim = bim.drop_duplicates("snp")

# 显著 SNP + 位置
sig = gwas[gwas["p"] < 5e-6].copy()
sig = sig.merge(bim[["snp", "chr", "bp"]], on="snp", how="left")
print("p<5e-6 且有位置:", len(sig))

# plink clump 用标准 assoc 格式（CHR SNP BP A1 F_A F_U A2 CHISQ P OR，带头）
sig = sig.dropna(subset=["chr", "bp"]).copy()
sig["chr"] = sig["chr"].astype(int); sig["bp"] = sig["bp"].astype(int)
sig = sig.merge(bim[["snp", "A1", "A2"]], on="snp", how="left")
sig["F_A"] = 0.5; sig["F_U"] = 0.5
sig["CHISQ"] = (sig["beta"] / sig["se"]) ** 2
sig["OR"] = np.exp(sig["beta"])
assoc = sig[["chr", "snp", "bp", "A1", "F_A", "F_U", "A2", "CHISQ", "p", "OR"]].dropna(subset=["A1", "A2"])
assoc.columns = ["CHR", "SNP", "BP", "A1", "F_A", "F_U", "A2", "CHISQ", "P", "OR"]
assoc.to_csv(os.path.join(WORK, "idd_sig.assoc"), sep="\t", index=False)
cmd = [PLINK, "--bfile", BFILE, "--clump", os.path.join(WORK, "idd_sig.assoc"),
       "--clump-p1", "5e-6", "--clump-p2", "5e-6", "--clump-r2", "0.01", "--clump-kb", "1000",
       "--out", os.path.join(WORK, "idd_clump"), "--allow-extra-chr", "--silent"]
print("运行 clump...")
subprocess.run(cmd, check=False, capture_output=True)
cl = pd.read_csv(os.path.join(WORK, "idd_clump.clumped"), sep=r"\s+")
ren = {"TOTAL": "TOT"}
cl = cl.rename(columns=ren)
print("独立信号数:", len(cl))
print(cl[["CHR", "BP", "SNP", "P", "TOT", "SP2"]].head(30).to_string(index=False))

# 47 基因 cis 区域
gene_cis = pd.DataFrame({
    "gene": ["GFPT1","GFPT2","GLS","GLS2","GLUD1","GLUD2","SLC25A12","SLC25A13",
             "SLC38A1","SLC38A2","SLC38A3","SLC38A5","SLC38A7","SLC38A9","SLC1A5",
             "SLC7A5","SLC7A6","SLC7A7","SLC7A8","SLC7A11","GOT1","GOT2","GPT","GPT2",
             "ASNS","ASNSD1","ASL","ASS1","ARG1","ARG2","NOS2","GAD1","GAD2",
             "GCLM","GCLC","GSS","GSTM1","GSTP1","GSTT1","NFS1","NFS2","GLS3",
             "GMPS","IMPDH1","IMPDH2","SUCLG1","SUCLG2"],
})
# 基因位置从 cis_snps 文件提取（已有 4 组织 cis_snps_full_*.csv 含 chr/bp 范围）
import glob
cis_files = glob.glob(str(P.RESULTS_DIR / "cis_snps_full_*.csv"))
pos = {}
for f in cis_files:
    d = pd.read_csv(f, sep="\t")
    for g, c, b in zip(d["gene"], d["chr"], d["probe_bp"]):
        pos.setdefault(g, []).append((c, b))
gene_pos = {}
for g, lst in pos.items():
    chrs = [c for c, b in lst]
    c = max(set(chrs), key=chrs.count)
    bps = [b for cc, b in lst if cc == c]
    gene_pos[g] = (c, min(bps), max(bps))
print("\n基因位置获取:", len(gene_pos))

# 重叠：独立信号在基因 ±2Mb
rows = []
for _, r in cl.iterrows():
    for g, (c, gstart, gend) in gene_pos.items():
        if r["CHR"] == c and (gstart - 2_000_000) <= r["BP"] <= (gend + 2_000_000):
            rows.append({"clump_SNP": r["SNP"], "chr": c, "bp": r["BP"], "GWAS_p": r["P"],
                         "gene": g, "gene_tss_low": gstart, "gene_tss_high": gend})
ov = pd.DataFrame(rows)
print("\n与 47 基因 cis 区域重叠的独立信号:", len(ov))
if len(ov):
    print(ov.to_string(index=False))
    ov.to_csv(P.RESULTS_DIR / "reverse_mr_overlap_hits.csv", index=False, encoding="utf-8-sig")

# ============================================================
# P0-4 修复：真正计算反向 Wald ratio（IDD → 表达）并输出
# reverse_mr_IDD_to_expr.csv（原代码只做 clump+重叠，21 号断链）
# 暴露 = IDD GWAS（b/se），结局 = GTEx eQTL（b/se，按组织）
# θ_rev = b_eqtl / b_gwas；se_rev 用完整 delta method
# ============================================================
if len(ov):
    # GWAS 效应（对齐用）
    gwas_sig = sig.set_index("snp")
    # eQTL 面板（4 组织 cis_snps_full）
    rev_rows = []
    for f in sorted(glob.glob(str(P.RESULTS_DIR / "cis_snps_full_*.csv"))):
        tissue = os.path.basename(f).replace("cis_snps_full_", "").replace(".csv", "")
        d = pd.read_csv(f, sep="\t", dtype={"SNP": str})
        for _, o in ov.iterrows():
            g = o["gene"]
            sub = d[(d["gene"] == g) & (d["SNP"] == o["clump_SNP"])]
            if sub.empty:
                continue
            sub = sub.iloc[0]
            gg = gwas_sig.loc[o["clump_SNP"]]
            ea1, ea2 = sub["A1"], sub["A2"]
            ga1, ga2 = gg["a1"], gg["a2"]
            if {ga1, ga2} != {ea1, ea2}:
                continue
            b_e = float(sub["b"])
            if ga1 == ea2:
                b_e = -b_e
            se_e = float(sub["se"])
            b_g = float(gg["beta"]); se_g = float(gg["se"])
            if b_g == 0:
                continue
            # 反向 Wald ratio + 完整 delta method SE
            b_rev = b_e / b_g
            se_rev = math.sqrt((se_e * se_e) / (b_g * b_g) + (b_e * b_e * se_g * se_g) / (b_g ** 4))
            z_rev = b_rev / se_rev
            p_rev = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z_rev) / math.sqrt(2.0))))
            same = (b_e > 0) == (b_g > 0)
            rev_rows.append({
                "clump_SNP": o["clump_SNP"], "gene": g, "tissue": tissue,
                "IDD_beta": b_g, "IDD_p": o["GWAS_p"],
                "eQTL_b": b_e, "eQTL_se": se_e, "eQTL_p": float(sub["p"]),
                "reverse_b": b_rev, "reverse_se": se_rev, "reverse_p": p_rev,
                "same_direction": "Y" if same else "N",
                "eQTL_sig": "Y" if float(sub["p"]) < 5e-5 else "N",
            })
    rev = pd.DataFrame(rev_rows)
    if len(rev):
        rev.to_csv(OUT, sep="\t", index=False, encoding="utf-8-sig")
        print(f"\n反向 MR 结果: {len(rev)} 条 -> {OUT}")
        print(rev.to_string(index=False))
    else:
        print("\n反向 MR: 无有效 SNP×基因×组织重叠（输出空文件占位）")
        pd.DataFrame(columns=["clump_SNP", "gene", "tissue", "IDD_beta", "IDD_p",
                              "eQTL_b", "eQTL_se", "eQTL_p", "reverse_b", "reverse_se",
                              "reverse_p", "same_direction", "eQTL_sig"]).to_csv(OUT, sep="\t", index=False)
else:
    print("无重叠信号，跳过反向 MR")
print("DONE 12")
