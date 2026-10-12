# -*- coding: utf-8 -*-
"""
统一路径配置（可移植版）—— pipeline 全部脚本的路径唯一入口。

设计约定（详见 pipeline/config/REFACTOR_CONVENTIONS.md）：
  1. 项目根目录由本文件位置自动推导（本文件位于 <ROOT>/pipeline/config/paths.py，
     上溯两级即项目根），因此整个仓库换机器、换目录均无需修改任何脚本内的项目路径。
  2. 机器专属的外部工具/数据目录集中在本文件下方“外部工具与机器专属路径”区：
     发布版不内置任何开发机路径——RSCRIPT / PLINK / SMR_EXE / R_LIBS / SMR_WORK /
     GEO_WORK 均须通过环境变量提供（环境变量优先；未设置时导入即报错或打印警告）。
  3. 环境变量一览：
        REP_MR_IDD_RSCRIPT  Rscript 可执行文件路径（R 步骤 06/09/19 需要）
        REP_MR_IDD_PLINK    PLINK 可执行文件路径（08/09/10/11/12/18 需要）
        REP_MR_IDD_SMR_EXE  SMR 可执行文件路径（默认位于 data/ 下）
        REP_MR_IDD_SMR_WORK SMR/条件分析中间工作目录（必须设置）
        REP_MR_IDD_GEO      GEO 原始数据与 UKB/OpenGWAS 缓存目录（必须设置）
        REP_MR_IDD_RLIBS    R 包库目录（R 脚本 .libPaths 使用；可选）
        OPENGWAS_TOKEN      OpenGWAS 访问 token（16 使用；也可置于 <GEO_WORK>/opengwas_token.txt）
        REP_MR_IDD_RESULTS  结果输出目录（默认 <ROOT>/results；重跑验证时重定向到临时目录）
  4. 可执行文件（RSCRIPT/PLINK/SMR_EXE）在导入时做存在性检查，缺失立即报错并提示如何配置。
"""
import os
from pathlib import Path

# ================= 项目根（自动推导，换电脑/换目录无需修改） =================
ROOT = Path(__file__).resolve().parents[2]
PIPELINE_DIR = ROOT / "pipeline"
CONFIG_DIR = PIPELINE_DIR / "config"
DATA_DIR = ROOT / "data"
RESULTS_DIR = Path(os.environ.get("REP_MR_IDD_RESULTS", str(ROOT / "results")))

# ================= 项目内数据文件（随仓库整体迁移，无需改动） =================
FINNGEN_MA = DATA_DIR / "finngen_R12_M13_INTERVERTEB.ma"
FINNGEN_RSID_MA = DATA_DIR / "finngen_R12_M13_INTERVERTEB.rsid.ma"
FINNGEN_SMR_MA = DATA_DIR / "finngen_R12_M13_INTERVERTEB.smr.ma"
GTEX_DIR = DATA_DIR / "GTEx_v8_westlake"
EQTLGEN_LEAD_CSV = DATA_DIR / "lead_eqtl_eQTLGen.csv"
EQTLGEN_LEAD_RDS = DATA_DIR / "lead_eqtl_eQTLGen.rds"
PANEL1000G = DATA_DIR / "1000G_plink" / "1kg_eur"          # plink 面板前缀（.bed/.bim/.fam）
GLUT_GENES = DATA_DIR / "glut_genes.txt"                    # 47 谷氨酰胺基因清单（01b/11/12 等）
GLUTAMINE_GENES = CONFIG_DIR / "glutamine_genes.txt"        # 47 基因（17/24）
GENE_SETS_CSV = CONFIG_DIR / "gene_sets.csv"                # 6 组 gene set（23）
SMR_EXE = DATA_DIR / "smr-1.3.1-win-x86_64" / "smr-1.3.1-win.exe"


def besd_prefix(tissue):
    """GTEx v8 BESD 前缀：<GTEX_DIR>/<tissue>/<tissue>（.besd/.epi/.esi）。"""
    return GTEX_DIR / tissue / tissue


# ================= 外部工具与机器专属路径 =================
def _env(key, default, what):
    v = os.environ.get(key)
    if v and v.strip():
        return Path(v.strip())
    if default:
        return Path(default)
    return None


def _require(p, key, what):
    if p is None or not p.exists():
        raise RuntimeError(
            f"[paths.py] 找不到 {what}（未配置默认值）\n"
            f"  请设置环境变量 {key}（指向正确路径）。"
        )
    return p


RSCRIPT = _require(_env("REP_MR_IDD_RSCRIPT", "", "Rscript"),
                   "REP_MR_IDD_RSCRIPT", "Rscript 可执行文件")
# PLINK 无合理默认路径（机器专属）：必须通过环境变量 REP_MR_IDD_PLINK 提供
PLINK = _require(_env("REP_MR_IDD_PLINK", "", "plink"),
                 "REP_MR_IDD_PLINK", "PLINK 可执行文件")
SMR_EXE = _require(_env("REP_MR_IDD_SMR_EXE", str(SMR_EXE), "smr"),
                   "REP_MR_IDD_SMR_EXE", "SMR 可执行文件")
R_LIBS = _env("REP_MR_IDD_RLIBS", "", "R 包库")


def _work_dir(key, what):
    v = os.environ.get(key, "").strip()
    if not v:
        print(f"[paths.py] 提示：{what} 未设置（环境变量 {key}）；如需该功能请配置后重跑。")
        return None
    return Path(v)


SMR_WORK = _work_dir("REP_MR_IDD_SMR_WORK", "SMR 工作目录")
GEO_WORK = _work_dir("REP_MR_IDD_GEO", "GEO 工作目录")

for _p, _k, _w in ((R_LIBS, "REP_MR_IDD_RLIBS", "R 包库"),):
    if _p is not None and not _p.exists():
        print(f"[paths.py] 警告：{_w} 不存在: {_p}（如需使用请设置环境变量 {_k}）")


def get_opengwas_token():
    """OpenGWAS token：环境变量 OPENGWAS_TOKEN → <GEO_WORK>/opengwas_token.txt。"""
    t = os.environ.get("OPENGWAS_TOKEN", "").strip()
    if not t:
        tf = GEO_WORK / "opengwas_token.txt"
        if tf.exists():
            t = tf.read_text(encoding="utf-8", errors="ignore").strip()
    return t
