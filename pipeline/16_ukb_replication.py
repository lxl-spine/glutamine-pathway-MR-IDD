# -*- coding: utf-8 -*-
"""UKB (ukb-b-19807) 复制分析：方向一致性 + Wald ratio + 核心基因复制表"""
# ============================================================
# 步骤 16 · UKB 独立复制（ukb-b-19807, 1,045 例/461,965 对照）
# ------------------------------------------------------------
# 功能: OpenGWAS API（JWT 认证）查询 47 基因 top SNP 的 UKB 关联，
#       Wald ratio 方向一致性复制（GFPT1 WB p=0.06、
#       rs4727335 p=0.02 方向一致——功效受限表述）
# 输入: mr_results top SNP + OpenGWAS token（环境变量 OPENGWAS_TOKEN）
# 输出: results/ukb_replication_47genes.csv
# ============================================================
import json
import pandas as pd, numpy as np
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
import os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

mr = pd.read_csv(P.RESULTS_DIR / "mr_results_GTEx4tissues.csv", sep="\t")
genes = ["ARG1","ARG2","ASL","ASNS","ASNSD1","ASS1","CAD","CPS1","CTPS1","GAD1","GCLC","GCLM",
         "GFPT1","GFPT2","GLS","GLUD1","GLUL","GMPS","GOT1","GOT2","GPT","GPT2","GSR","GSS",
         "IDH1","IDH2","NFS1","OGDH","PPAT","SLC1A5","SLC25A12","SLC25A13","SLC38A1","SLC38A2",
         "SLC38A3","SLC38A4","SLC38A7","SLC38A9","SLC7A11","SLC7A5","SLC7A6","SLC7A7","SUCLA2","SUCLG1",
         "CTPS2","SLC38A5","GLUD2"]
sub = mr[mr.gene.isin(genes)].copy()

# ============================================================
# P2-17 修复：ukb_19807_associations.json 缺失时自动从 OpenGWAS API
# 拉取（ukb-b-19807），不再依赖手工放置的缓存文件
# ============================================================
UKB_JSON = str(P.GEO_WORK / "ukb_19807_associations.json")
import os, json
if os.path.exists(UKB_JSON):
    ukb = json.load(open(UKB_JSON))
else:
    print("ukb_19807_associations.json 缺失，从 OpenGWAS API 拉取...", flush=True)
    import urllib.request
    # token 来源：环境变量 OPENGWAS_TOKEN → <GEO_WORK>/opengwas_token.txt（统一走 paths.get_opengwas_token()）
    token = P.get_opengwas_token()
    if not token:
        raise SystemExit("ERROR: 需要 OPENGWAS_TOKEN 环境变量（或 <GEO_WORK>/opengwas_token.txt）")
    ukb_id = "ukb-b-19807"
    # 47 基因 top SNP（去重）
    top_snps = sub.sort_values("p_MR").drop_duplicates("gene")["SNP"].tolist()
    ukb = {}
    API = "https://api.opengwas.io/api/associations"
    # 分批（每批 ≤64）
    for i in range(0, len(top_snps), 60):
        batch = top_snps[i:i+60]
        body = json.dumps({"variant": batch}).encode()
        req = urllib.request.Request(API, data=body, method="POST",
              headers={"Content-Type": "application/json", "Authorization": f"Bearer {token}"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = json.loads(resp.read().decode())
        for a in data:
            if a.get("variant"):
                ukb[a["variant"]] = {"beta": a.get("beta"), "se": a.get("se"),
                                     "p": a.get("p"), "ea": a.get("ea"), "nea": a.get("nea")}
        print(f"  batch {i//60+1}: {len(data)} variants", flush=True)
    with open(UKB_JSON, 'w') as f:
        json.dump(ukb, f)
    print(f"API 拉取完成 {len(ukb)} variants -> {UKB_JSON}", flush=True)
ukb_df = pd.DataFrame(ukb).T.reset_index().rename(columns={"index": "SNP"})

# 合并（对齐等位基因：FinnGen MR 的 b_GWAS 是 eA1 方向；UKB beta 是 ea 方向）
m = sub.merge(ukb_df, on="SNP", how="inner", suffixes=("", "_ukb"))
print(f"UKB 有记录的 MR 行: {len(m)} / {len(sub)} | 基因: {m.gene.nunique()} | 组织×基因: {m.groupby(['gene','tissue']).ngroups}")

# 等位基因对齐：flip UKB beta 使其对齐 eA1
def align(row):
    b = row["beta"]
    if row["ea"] == row["eA1"] and row["nea"] == row["eA2"]:
        return b
    if row["ea"] == row["eA2"] and row["nea"] == row["eA1"]:
        return -b
    return np.nan  # 等位基因不匹配（palindrome 或多等位）

m["b_ukb_align"] = m.apply(align, axis=1)
m["dir_consistent"] = np.sign(m["b_ukb_align"]) == np.sign(m["b_GWAS"])
m["ukb_nominal"] = m["p"] < 0.05

# Wald ratio 复制：b_ukb/b_eQTL
from scipy import stats
m["b_MR_ukb"] = m["b_ukb_align"] / m["b_eQTL"]
m["se_MR_ukb"] = abs(m["se"] / m["b_eQTL"])
z = m["b_MR_ukb"] / m["se_MR_ukb"]
m["p_MR_ukb"] = [2*(1-stats.norm.cdf(abs(x))) if np.isfinite(x) else np.nan for x in z]

print("\n=== 17 候选基因 UKB 复制（Wald ratio 方向）===")
cand = ["GFPT1","SLC25A13","SLC25A12","NFS1","SLC7A7","SLC38A1","SLC38A3","GAD1","ASL","GCLM",
        "ARG1","GLUD1","GSS","GMPS","SUCLG1","SLC7A5","SLC7A6"]
mc = m[m.gene.isin(cand)].copy()
mc = mc.sort_values("p")
show = mc[["tissue","gene","SNP","b_MR","p_MR","FDR_MR","b_ukb_align","p","dir_consistent","p_MR_ukb"]].round(4)
print(show.to_string(index=False))

# 汇总统计
print("\n=== UKB 复制汇总 ===")
print(f"UKB 有记录: {len(m)} 行 | 方向一致: {m['dir_consistent'].sum()} | UKB 名义 p<0.05: {m['ukb_nominal'].sum()}")
print(f"FinnGen 名义显著 (p_MR<0.05) 且有 UKB 记录: {len(m[m['p_MR']<0.05])}")
fg_sig = m[m['p_MR']<0.05]
if len(fg_sig):
    print(f"  其中方向一致: {fg_sig['dir_consistent'].sum()}/{len(fg_sig)}")
    print(f"  其中 UKB p<0.05: {fg_sig['ukb_nominal'].sum()}")

# 核心三基因
print("\n=== 核心基因 UKB 复制明细 ===")
core = m[m.gene.isin(["SLC25A13","GFPT1","SLC25A12","NFS1"])]
print(core[["tissue","gene","SNP","b_eQTL","p_eQTL","b_GWAS","p_GWAS","b_ukb_align","p","dir_consistent","p_MR_ukb"]].round(5).to_string(index=False))

m.to_csv(P.RESULTS_DIR / "ukb_replication_47genes.csv", index=False, encoding="utf-8-sig")
print("\n保存: results\\ukb_replication_47genes.csv")
