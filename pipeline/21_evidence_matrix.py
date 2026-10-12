# -*- coding: utf-8 -*-
"""证据矩阵 → 一致性矩阵 v2.0：加入 UKB 复制、条件分析、Steiger、置换检验"""
# ============================================================
# 步骤 21 · 证据矩阵汇总 v2（六层证据一致性矩阵）
# ------------------------------------------------------------
# 功能: 每基因×组织合并 FinnGen MR / eQTLGen / SMR / coloc.abf /
#       双向 MR / UKB 复制六层证据，输出方向+显著性矩阵
# 输入: 上述全部结果 CSV
# 输出: results/evidence_matrix_v2.csv
# ============================================================
import pandas as pd, numpy as np, os
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

# ============================================================
# P2-17 修复：evidence_matrix_glut_17genes.csv（基础矩阵）若缺失则
# 自动从各步骤结果构建（原代码依赖无生成脚本的手工产物，管线断链）
# ============================================================
BASE_MATRIX = P.RESULTS_DIR / "evidence_matrix_glut_17genes.csv"
if not os.path.exists(BASE_MATRIX):
    print("构建 evidence_matrix_glut_17genes.csv（从步骤 03/07/06/08 结果）...", flush=True)
    mr = pd.read_csv(P.RESULTS_DIR / "mr_results_GTEx4tissues.csv", sep="\t")
    base = mr[["tissue", "gene", "b_MR", "FDR_MR", "p_MR", "b_GWAS", "p_GWAS"]].copy()
    # eQTLGen p（按基因合并，17 基因在 eQTLGen 有记录的）
    eqtlgen = pd.read_csv(P.RESULTS_DIR / "mr_results_eQTLGen_blood.csv", sep="\t")
    eqtlgen_top = eqtlgen.sort_values("p_MR").drop_duplicates("gene")[["gene", "p_MR"]].rename(columns={"p_MR": "eQTLGen_p_MR"})
    base = base.merge(eqtlgen_top, on="gene", how="left")
    # coloc.abf PP.H4
    coloc = pd.read_csv(P.RESULTS_DIR / "coloc_results.csv", sep="\t")
    coloc_ren = {"Adipose": "Adipose_Subcutaneous", "Muscle": "Muscle_Skeletal",
                 "Nerve": "Nerve_Tibial", "WholeBlood": "Whole_Blood"}
    coloc["tissue"] = coloc["tissue"].map(coloc_ren).fillna(coloc["tissue"])
    base = base.merge(coloc[["tissue", "gene", "PP.H4", "PP.H3"]], on=["tissue", "gene"], how="left")
    # SMR-HEIDI
    smr = pd.read_csv(P.RESULTS_DIR / "smr_glut_17genes_4tissues.csv", sep="\t")
    base = base.merge(smr[["tissue", "gene", "b_SMR", "p_SMR", "p_HEIDI"]], on=["tissue", "gene"], how="left")
    # 17 基因主证据子集
    genes17 = ["GFPT1","SLC25A13","ARG1","ASL","GAD1","GCLM","GLUD1","GMPS","GSS",
               "NFS1","SLC25A12","SLC38A1","SLC38A3","SLC7A5","SLC7A6","SLC7A7","SUCLG1"]
    base = base[base.gene.isin(genes17)].reset_index(drop=True)
    base.to_csv(BASE_MATRIX, sep="\t", index=False)
    print(f"evidence_matrix_glut_17genes.csv: {len(base)} 行", flush=True)

em = pd.read_csv(BASE_MATRIX, sep="\t")
ukb = pd.read_csv(P.RESULTS_DIR / "ukb_replication_47genes.csv")
ukb["gene_tissue"] = ukb["gene"] + "_" + ukb["tissue"]
em["gene_tissue"] = em["gene"] + "_" + em["tissue"]
# UKB 复制合并（方向一致 + UKB p）
em = em.merge(ukb[["gene_tissue","dir_consistent","p","p_MR_ukb"]].rename(
    columns={"p": "UKB_p", "p_MR_ukb": "UKB_MR_p"}), on="gene_tissue", how="left")
# 双向 MR
rev = pd.read_csv(P.RESULTS_DIR / "reverse_mr_IDD_to_expr.csv", sep="\t")
rev["gene_tissue"] = rev["gene"] + "_" + rev["tissue"]
em = em.merge(rev[["gene_tissue","clump_SNP","IDD_p","eQTL_p","same_direction"]].rename(
    columns={"clump_SNP": "bidir_SNP", "IDD_p": "bidir_IDD_p", "eQTL_p": "bidir_eQTL_p",
             "same_direction": "bidir_same_dir"}), on="gene_tissue", how="left")

em.to_csv(P.RESULTS_DIR / "evidence_matrix_v2.csv", index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250)
show = em[em.gene.isin(["SLC25A13","GFPT1","SLC25A12","NFS1","GLS","GLUD1","SLC38A1","SLC7A5"])]
print(show[["gene","tissue","b_MR","FDR_MR","eQTLGen_p_MR","PP.H4","b_SMR","p_SMR","p_HEIDI",
            "dir_consistent","UKB_p","bidir_SNP","bidir_IDD_p","bidir_same_dir"]].round(4).to_string(index=False))
print("\n已存 evidence_matrix_v2.csv")
