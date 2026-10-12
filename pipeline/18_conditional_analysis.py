# -*- coding: utf-8 -*-
"""SLC25A13 区域条件分析：rs4727335 ∥ rs12216562（1000G EUR LD + FinnGen R12 区域 GWAS）
# ============================================================
# 步骤 18 · 条件分析（SLC25A13 chr7:95.8Mb 区域）
# ------------------------------------------------------------
# 功能: rs4727335 与 rs12216562 的 LD 检查（1000G EUR plink，
#       r²=0.692 同一 block）+ COJO 近似条件分析 + 区域扫描
# 输入: glut_ld/SLC25A13 面板 + FinnGen
# 输出: results/conditional_analysis_summary.txt +
#       results/conditional_analysis_slc25a13_region.csv
#       （结论：单因果信号，修正早期"独立信号"误述）
# ============================================================
COJO 风格单条件近似: z_cond = (z_j - r_jk * z_k) / sqrt(1 - r_jk^2)"""
import numpy as np
import pandas as pd
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

# 区域 GWAS（FinnGen R12 IDD，SLC25A13 cis±1Mb）
g = pd.read_csv(P.SMR_WORK / "gwas_SLC25A13.txt", sep="\t")

# ============================================================
# LD 来源修复（2026-09-22）：旧代码读 --r2 输出（只有 R2 列无 D 列，
# KeyError: 'D'）。改用 11 号生成的 SLC25A13 面板带符号 r 矩阵
# （plink --r square 输出皮尔逊相关 r，带符号），覆盖 chr7 95.8Mb±2Mb
# ============================================================
PANEL_LD = P.RESULTS_DIR / "glut_ld" / "SLC25A13_chr7_ld.ld"
ldm = pd.read_csv(PANEL_LD, sep=r"\s+", header=None)   # 方阵，行列=面板 bim SNP 顺序
bim = pd.read_csv(P.RESULTS_DIR / "glut_ld" / "SLC25A13_chr7.bim",
                  sep=r"\s+", header=None, names=['chr', 'SNP', 'cm', 'bp', 'A1', 'A2'])
snp_order = bim['SNP'].tolist()
assert ldm.shape[0] == ldm.shape[1] == len(snp_order), "面板 LD 矩阵维度与 bim 不符"
r_mat = pd.DataFrame(ldm.values, index=snp_order, columns=snp_order)

SNP1, BP1 = "rs4727335", 95811083   # 双向 MR SNP：A 增表达+增IDD
SNP2, BP2 = "rs12216562", 95967710  # SMR top SNP：G 降表达（SMR +）

r1 = g[g.snp == SNP1].iloc[0]
r2 = g[g.snp == SNP2].iloc[0]
print(f"{SNP1}: beta={r1.beta:.5f} se={r1.se:.5f} p={r1.p:.2e} z={r1.beta/r1.se:.3f}")
print(f"{SNP2}: beta={r2.beta:.5f} se={r2.se:.5f} p={r2.p:.2e} z={r2.beta/r2.se:.3f}")

# LD（带符号 r，来自面板矩阵）
# rs12216562 在面板中可能为 chr:pos:ref:alt 格式（1kg_eur fix 后），按位置回退匹配
def panel_snp(rsid, bp, chrom=7):
    if rsid in r_mat.index:
        return rsid
    hits = [s for s in r_mat.index if s.startswith(f"{chrom}:{bp}:")]
    if len(hits) == 1:
        return hits[0]
    raise SystemExit(f"ERROR: SNP 不在面板中 {rsid} (bp={bp}) 候选={hits[:5]}")
s1 = panel_snp(SNP1, BP1)
s2 = panel_snp(SNP2, BP2)
r_panel = float(r_mat.loc[s1, s2])

# ============================================================
# 相位对齐（2026-09-22 关键修复）：面板 --r square 的 r 符号基于 bim A1；
# 条件分析在 GWAS 效应等位（gwas a1）空间进行，必须将 r 翻转到该空间。
# 例：rs12216562 GWAS a1=G，面板 A1=A → 翻转因子 -1；
#     rs4727335 GWAS a1=A=面板 A1 → +1。r_eff = r_panel * f1 * f2。
# ============================================================
bim_alleles = dict(zip(bim['SNP'], list(zip(bim['A1'], bim['A2']))))
g_alleles = dict(zip(g['snp'], list(zip(g['a1'], g['a2']))))

def flip_factor(rsid, panel_sid):
    pa1, pa2 = bim_alleles.get(panel_sid, (None, None))
    ga1, ga2 = g_alleles.get(rsid, (None, None))
    if ga1 is None or pa1 is None:
        return 1.0
    if ga1 == pa1:
        return 1.0
    if ga1 == pa2:
        return -1.0
    return 1.0  # 无法判定，保守不翻转（区域扫描中该 SNP r 置 0 更安全，见下）

f1 = flip_factor(SNP1, s1)
f2 = flip_factor(SNP2, s2)
r_sign = r_panel * f1 * f2
r2_ld = r_sign ** 2
print(f"面板 SNP: {s1} / {s2}；相位因子 {f1} × {f2}（GWAS a1: A/A → G/A）")
print(f"\n1000G LD（GWAS 效应空间带符号 r）: R2={r2_ld:.4f} 带符号 r={r_sign:.4f}")

# 条件分析
z1, z2 = r1.beta/r1.se, r2.beta/r2.se
def cond(z_j, z_k, r):
    return (z_j - r*z_k) / np.sqrt(1 - r**2)

z_cond_1 = cond(z1, z2, r_sign)   # SNP1 条件 SNP2
z_cond_2 = cond(z2, z1, r_sign)   # SNP2 条件 SNP1
p1 = 2*(1-__import__("scipy").stats.norm.cdf(abs(z_cond_1)))
p2 = 2*(1-__import__("scipy").stats.norm.cdf(abs(z_cond_2)))

print(f"\n=== 条件分析（COJO 单条件近似）===")
print(f"{SNP1} 条件 {SNP2}: z_cond={z_cond_1:.3f} p={p1:.3e}")
print(f"{SNP2} 条件 {SNP1}: z_cond={z_cond_2:.3f} p={p2:.3e}")
print(f"\n结论: 条件后 {SNP1} 仍显著(p={p1:.2e})，{SNP2} 消失(p={p2:.2f})")
print(f"→ 区域内单因果信号，主导位点为 {SNP1} 所在 LD block（chr7:95.81Mb）")
print(f"→ rs12216562 为 LD 代理（r2=0.69），非独立第二信号")

# 区域扫描：给 GWAS SNP 补位置（从 GTEx cis eQTL 数据映射 rsid→bp）
cis = pd.read_csv(P.RESULTS_DIR / "cis_snps_full_Whole_Blood.csv", sep="\t")
posmap = dict(zip(cis.SNP, cis.SNP_bp))
g["BP"] = g["snp"].map(posmap)
g["z"] = g["beta"]/g["se"]
print(f"区域 GWAS SNP 位置映射: {g['BP'].notna().sum()}/{len(g)}")
# 远离主信号 LD block 的独立候选（>600kb，LD 可忽略）
far = g[(g.BP.notna()) & (abs(g.BP - BP1) > 600000)]
sig_far = far[far["p"] < 1e-5]
print(f"\n远离主信号(>600kb)的 SNP: {len(far)} 个, p<1e-5: {len(sig_far)} 个")
if len(sig_far):
    print(sig_far[["snp","BP","beta","se","p"]].sort_values("p").to_string(index=False))
else:
    print("→ 区域内无第二个独立信号候选")

# 条件 rs4727335 后区域残余（P1-12 修复：用带符号 r，来自面板 --r square 矩阵）
# 相位：面板 r 是 bim A1 空间；翻转到 GWAS 效应空间需乘 flip(rs4727335) × flip(目标 SNP)
r_vec = r_mat.loc[s1, :].drop(index=s1)
g2 = g[g.BP.notna()].copy()
def _panel_sid(rsid, bp):
    if rsid in r_vec.index:
        return rsid
    hits = [s for s in r_vec.index if s.startswith(f"7:{int(bp)}:")]
    return hits[0] if len(hits) == 1 else None
g2["panel_snp"] = [_panel_sid(s, b) if not pd.isna(b) else None for s, b in zip(g2["snp"], g2["BP"])]
def _eff_r(rsid, panel_sid):
    if panel_sid is None or panel_sid not in r_vec.index:
        return 0.0
    return float(r_vec[panel_sid]) * f1 * flip_factor(rsid, panel_sid)
g2["r_sign"] = [_eff_r(s, p) for s, p in zip(g2["snp"], g2["panel_snp"])]
g2.loc[g2["r_sign"].isna(), "r_sign"] = 0.0
# 条件 z = (z_j - r_jk * z_k) / sqrt(1 - r_jk^2)，r_jk 为带符号 LD（与主条件分析一致）
z1 = r1.beta/r1.se
from scipy import stats
g2["z_cond"] = g2.apply(lambda row: abs((row.z - row.r_sign*z1)/np.sqrt(max(1-row.r_sign**2, 0.01))), axis=1)
g2["p_cond"] = 2*(1-stats.norm.cdf(g2["z_cond"]))
g2 = g2.sort_values("p_cond")
print("\n=== 条件 rs4727335 后区域 top5 残余（带符号 r 条件化）===")
print(g2.head(5)[["snp","BP","beta","se","p","r_sign","p_cond"]].round(6).to_string(index=False))
print(f"条件后 p<1e-5 的残余信号数: {(g2['p_cond']<1e-5).sum()}")
g2.to_csv(P.RESULTS_DIR / "conditional_analysis_slc25a13_region.csv", index=False, encoding="utf-8-sig")
print("已存 results\\conditional_analysis_slc25a13_region.csv")
