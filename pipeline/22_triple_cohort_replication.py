# -*- coding: utf-8 -*-
"""步骤 22 · 三队列复制（FinnGen IDD 主 + UKB IDD + UKB 背痛 118k 扩展）
表型差异声明：ukb-b-9838 为"背痛"（症状表型），非结构性 IDD，仅作扩展复制。
"""
import pandas as pd, numpy as np, json, math
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

BASE = P.RESULTS_DIR
mr = pd.read_csv(BASE / "mr_results_GTEx4tissues.csv", sep="\t")
mr.columns = [c.strip() for c in mr.columns]

# ---- 1. 核心基因四组织 FinnGen MR ----
core = mr[mr.gene.isin(["GFPT1", "SLC25A13"])][["gene", "tissue", "SNP", "eA1", "eA2",
    "b_eQTL", "se_eQTL", "p_eQTL", "b_GWAS", "se_GWAS", "p_GWAS", "b_MR", "se_MR", "p_MR", "FDR_MR"]].copy()
core["F_eQTL"] = (core["b_eQTL"] / core["se_eQTL"]) ** 2

# ---- 2. UKB IDD (ukb-b-19807) ----
ukb19807 = pd.read_csv(BASE / "ukb_replication_47genes.csv")
ukb19807.columns = [c.strip() for c in ukb19807.columns]
u19807 = ukb19807[ukb19807.gene.isin(["GFPT1", "SLC25A13"])][["gene", "tissue", "SNP", "b_MR_ukb", "se_MR_ukb", "p_MR_ukb"]].copy()
# SLC25A13 补 rs4727335 proxy（UKB IDD）
d1 = json.load(open(P.GEO_WORK / "ukb_19807_proxy_rs4727335.json", encoding="utf-8"))[0]
b_e = 0.240163; se_e = 0.030035
u19807 = pd.concat([u19807, pd.DataFrame([{"gene": "SLC25A13", "tissue": "Whole_Blood",
    "SNP": "rs4727335(proxy)", "b_MR_ukb": d1["beta"] / b_e, "se_MR_ukb": d1["se"] / b_e,
    "p_MR_ukb": d1["p"]}])], ignore_index=True)

# ---- 3. UKB 背痛 (ukb-b-9838, 118,471 cases) ----
d2 = json.load(open(P.GEO_WORK / "ukb_b9838_core_snps.json", encoding="utf-8"))
# eQTL 效应映射（GTEx WB；GFPT1 用各组织 lead 的 WB eQTL 近似——如实标注组织）
eqtl_map = {
    ("GFPT1", "Adipose_Subcutaneous", "rs6758228"): (-0.151209, 0.022926),
    ("GFPT1", "Muscle_Skeletal", "rs7585658"): (0.182524, 0.028061),
    ("GFPT1", "Nerve_Tibial", "rs145160777"): (-0.205058, 0.028040),
    ("GFPT1", "Whole_Blood", "rs7370877"): (0.163872, 0.018249),
    ("SLC25A13", "Muscle_Skeletal", "rs10250928"): (0.220449, 0.032000),
    ("SLC25A13", "Nerve_Tibial", "rs1987584"): (-0.104206, 0.032705),
    ("SLC25A13", "Whole_Blood", "rs4727335"): (0.240163, 0.030035),
    ("SLC25A13", "Whole_Blood", "rs10631575"): (0.256325, 0.029768),
}
rows = []
for a in d2:
    rsid = a["rsid"]; b_g = a["beta"]; se_g = a["se"]; p_g = a["p"]; ea = a["ea"]
    # 匹配基因/组织：通过 rsid 反查
    matched = None
    for (g, t, snp), (b_e, se_e) in eqtl_map.items():
        if snp == rsid:
            matched = (g, t, snp, b_e, se_e); break
    if not matched:
        continue
    g, t, snp, b_e, se_e = matched
    # 等位对齐：ukb ea 与 eQTL A1 方向——用已知 A1
    # eQTL 各 SNP A1：rs6758228 A1=C; rs7585658 A1=G; rs145160777 A1=?; rs7370877 A1=G;
    # rs10250928 A1=?; rs1987584 A1=A; rs4727335 A1=A; rs10631575 A1=C
    a1_map = {"rs6758228": "C", "rs7585658": "G", "rs145160777": "T", "rs7370877": "G",
              "rs10250928": "A", "rs1987584": "A", "rs4727335": "A", "rs10631575": "C"}
    align = 1 if ea == a1_map.get(snp) else -1
    b_g_al = b_g * align
    b_mr = b_g_al / b_e
    se_mr = se_g / abs(b_e)
    z = b_mr / se_mr
    p_mr = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    rows.append({"gene": g, "tissue": t, "SNP": snp, "b_MR_ukb_bp": b_mr, "se_MR_ukb_bp": se_mr,
                 "p_MR_ukb_bp": p_mr, "dir_ukb_bp": "+" if b_mr > 0 else "-"})
ukb9838 = pd.DataFrame(rows)

# ---- 4. 组装三队列表 ----
out = core.merge(u19807[["gene", "tissue", "SNP", "b_MR_ukb", "se_MR_ukb", "p_MR_ukb"]],
                 on=["gene", "tissue"], how="left", suffixes=("", "_ukb"))
out = out.merge(ukb9838[["gene", "tissue", "b_MR_ukb_bp", "se_MR_ukb_bp", "p_MR_ukb_bp", "dir_ukb_bp"]],
                on=["gene", "tissue"], how="left")
# 方向一致性
def dir_of(b):
    return "+" if b > 0 else "-"
out["dir_fin"] = out["b_MR"].apply(dir_of)
out["dir_19807"] = np.where(out["b_MR_ukb"].notna(), out["b_MR_ukb"].apply(dir_of), np.nan)
out["dir_9838"] = out.get("dir_ukb_bp", np.nan)
out["triple_dir"] = out.apply(lambda r: "一致(3)" if r["dir_fin"] == r["dir_19807"] == r["dir_9838"] else
                              ("一致(2)" if (r["dir_fin"] == r["dir_19807"]) or (r["dir_fin"] == r["dir_9838"]) else "NA"), axis=1)
out.to_csv(BASE / "triple_cohort_replication.csv", index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(out[["gene", "tissue", "SNP", "F_eQTL", "b_MR", "p_MR", "FDR_MR", "b_MR_ukb", "p_MR_ukb",
           "b_MR_ukb_bp", "p_MR_ukb_bp", "triple_dir"]].round(4).to_string(index=False))
print("\n已存 triple_cohort_replication.csv")
