# -*- coding: utf-8 -*-
"""druggability + 竞争情报汇总（2026-09-22 检索）"""
# ============================================================
# 步骤 13 · 可成药性注释（不 claim 靶点，仅讨论素材）
# ------------------------------------------------------------
# 功能: DGIdb/Open Targets 对 47 谷氨酰胺基因的可成药性标注
# 输入: 候选基因清单
# 输出: results/druggability_summary.csv
# ============================================================
import pandas as pd
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

rows = [
    # gene, 可药性, 已知药物/工具化合物, 孟德尔病/功能证据, 竞争占用状态
    ["GFPT1", "高（可药靶点）",
     "IOX1（DGIdb，分子对接+酶活抑制+斑马鱼拯救）；DON/Azaserine（谷氨酰胺类似物不可逆抑制剂）；anticapsin（Ki=9.5µM）",
     "先天性肌无力综合征（CMS）；已糖胺通路限速酶",
     "已占用×5：Commun Biol 2026（跨种族+斑马鱼IOX1，Tian JCB）；JOR Spine 2025（TWAS+MR，PP.H4=0.81）；Spinal Stenosis MR 2026（eQTLGen OR=1.42）；Front Mol Neurosci 2026（DLK1+GFPT1 靶点）；EHR GWAS（GFPT1+NFU1 位点）"],
    ["SLC25A13", "低（无小分子药物；代谢管理）",
     "无直接药物（citrin 缺乏症：MCT 补充、低糖饮食等管理）",
     "Citrin 缺乏症（常隐，东亚高频，明确孟德尔病）；线粒体天冬氨酸-谷氨酸载体",
     "IDD/DDD 领域无人占用（本次检索确认）"],
    ["SLC25A12", "低（同族无药）",
     "无直接药物",
     "Aralar1（SLC25A13 同源载体）",
     "无 IDD 占用；但 SMR 显著+多 IV 复现（Muscle）"],
    ["NFS1", "低-中",
     "无直接药物；铁硫簇合成限速酶",
     "中山 GLS1-NFS1 轴（FRBM 2025 PMID 39710108）",
     "机制桥（谷氨酰胺→铁硫簇）"],
]
df = pd.DataFrame(rows, columns=["gene", "可药性", "已知药物/工具化合物", "孟德尔病/功能证据", "竞争占用状态"])
df.to_csv(P.RESULTS_DIR / "druggability_summary.csv", index=False, encoding="utf-8-sig")
print(df.to_string(index=False))
print("\n已保存: results\\druggability_summary.csv")
