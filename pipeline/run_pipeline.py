# -*- coding: utf-8 -*-
"""rep_MR_IDD 管线重跑驱动：按依赖顺序执行 01→24（修复后）
用法: python run_pipeline.py [step_start] [step_end]  （可断点续跑，如 run_pipeline.py 12 24）
"""
import os, sys, subprocess, time, json
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

BASE = str(P.ROOT)
PIPE = str(P.PIPELINE_DIR)
RES = str(P.RESULTS_DIR)
DATA = str(P.DATA_DIR)
RSCRIPT = str(P.RSCRIPT)
PLINK = str(P.PLINK)
SMR = str(P.SMR_EXE)

os.makedirs(RES, exist_ok=True)

def run(cmd, tag, cwd=None, timeout=None, env=None):
    t0 = time.time()
    print(f"\n=== [{tag}] {time.strftime('%H:%M:%S')} START: {cmd}", flush=True)
    try:
        r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True,
                           encoding="utf-8", errors="replace", timeout=timeout)
        out = (r.stdout or "") + ("\n[STDERR]\n" + r.stderr if r.stderr else "")
    except subprocess.TimeoutExpired:
        print(f"=== [{tag}] TIMEOUT after {timeout}s", flush=True)
        return False
    ok = (r.returncode == 0) if hasattr(r, 'returncode') else False
    # 打印末尾 30 行
    tail = out.strip().splitlines()[-30:] if out.strip() else ["(no output)"]
    for l in tail:
        print(f"  {l}", flush=True)
    print(f"=== [{tag}] END rc={getattr(r,'returncode','?')} ({time.time()-t0:.0f}s)", flush=True)
    return ok

def py(script, tag, args=None, timeout=None):
    cmd = [sys.executable, os.path.join(PIPE, script)] + (args or [])
    return run(cmd, tag, cwd=PIPE, timeout=timeout)

def r(script, tag, timeout=None):
    # R 侧 bootstrap 依赖环境变量定位项目根/results/geo/smr 目录，显式导出保持一致
    env = dict(os.environ)
    env["REP_MR_IDD_ROOT"] = BASE
    env["REP_MR_IDD_RESULTS"] = RES
    env["REP_MR_IDD_GEO"] = str(P.GEO_WORK)
    env["REP_MR_IDD_SMR_WORK"] = str(P.SMR_WORK)
    env["REP_MR_IDD_RLIBS"] = str(P.R_LIBS)
    return run([RSCRIPT, os.path.join(PIPE, script)], tag, cwd=PIPE, timeout=timeout, env=env)

# 47 基因清单
GENES = str(P.GLUT_GENES)
TISSUES = ['Adipose_Subcutaneous', 'Muscle_Skeletal', 'Nerve_Tibial', 'Whole_Blood']
SHORT = {'Adipose_Subcutaneous': 'Adipose', 'Muscle_Skeletal': 'Muscle',
         'Nerve_Tibial': 'Nerve', 'Whole_Blood': 'WholeBlood'}

def step01():
    oks = []
    for t in TISSUES:
        pre = str(P.besd_prefix(t))
        out = os.path.join(RES, f"lead_eqtl_{SHORT[t]}.csv")
        if os.path.exists(out) and os.path.getsize(out) > 1e5:
            print(f"  {t} lead 已存在，跳过", flush=True); oks.append(True); continue
        oks.append(py("01_data_prep_lead_eqtl.py", f"01|{t}", [pre, out], timeout=1800))
    return all(oks)

def step01b():
    oks = []
    for t in TISSUES:
        pre = str(P.besd_prefix(t))
        out = os.path.join(RES, f"cis_snps_full_{t}.csv")
        if os.path.exists(out) and os.path.getsize(out) > 1e5:
            print(f"  {t} cis 已存在，跳过", flush=True); oks.append(True); continue
        oks.append(py("01b_cis_snps_extract.py", f"01b|{t}", [pre, GENES, out], timeout=1800))
    return all(oks)

def step02(): return py("02_power_per_gene.py", "02 功效", timeout=900)
def step03(): return py("03_mr_main.py", "03 MR主分析", timeout=900)
def step04(): return py("04_tissue_stratify.py", "04 组织分层", timeout=600)
def step05(): return py("05_coloc_prep.py", "05 coloc输入", timeout=900)
def step06(): return r("06_coloc_abf.R", "06 coloc.abf", timeout=7200)
def step07(): return py("07_eqtlgen_replication.py", "07 eQTLGen复制", timeout=900)
def step11(): return py("11_glut_panels.py", "11 面板+LD", timeout=3600)
def step08(): return py("08_smr_heidi.py", "08 SMR-HEIDI", timeout=7200)
def step09(): return r("09_coloc_susie.R", "09 coloc.susie", timeout=10800)
def step10(): return py("10_multi_iv_mr.py", "10 多IV敏感性", timeout=1800)
def step12(): return py("12_bidirectional_mr.py", "12 双向MR", timeout=1800)
def step13(): return py("13_druggability.py", "13 可成药性", timeout=300)
def step14(): return py("14_expr_GSE70362.py", "14 GSE70362", timeout=900)
def step15(): return py("15_expr_GSE124272.py", "15 GSE124272", timeout=900)
def step16(): return py("16_ukb_replication.py", "16 UKB复制", timeout=1800)
def step17(): return py("17_enrichment_permutation.py", "17 富集置换", timeout=1800)
def step18(): return py("18_conditional_analysis.py", "18 条件分析", timeout=600)
def step19(): return r("19_scRNA_seurat.R", "19 scRNA Seurat", timeout=7200)
def step20(): return py("20_scRNA_subtype.py", "20 scRNA subtype", timeout=3600)
def step21(): return py("21_evidence_matrix.py", "21 证据矩阵", timeout=900)
def step22(): return py("22_triple_cohort_replication.py", "22 三队列复制", timeout=600)
def step23(): return py("23_gene_set_enrichment.py", "23 gene set富集", timeout=600)
def step24(): return py("24_loo_tissue_permutation.py", "24 LOO+组织分层", timeout=7200)

STEPS = [
    ("01", step01), ("01b", step01b), ("02", step02), ("03", step03), ("04", step04),
    ("05", step05), ("06", step06), ("07", step07),
    ("11", step11), ("08", step08), ("09", step09), ("10", step10),
    ("12", step12), ("13", step13), ("14", step14), ("15", step15),
    ("16", step16), ("17", step17), ("18", step18), ("19", step19),
    ("20", step20), ("21", step21), ("22", step22), ("23", step23), ("24", step24),
]

def main():
    start, end = 0, len(STEPS)
    if len(sys.argv) >= 2:
        start = max(0, next((i for i, (s, _) in enumerate(STEPS) if s == sys.argv[1]), 0))
    if len(sys.argv) >= 3:
        end = min(len(STEPS), next((i for i, (s, _) in enumerate(STEPS) if s == sys.argv[2]), len(STEPS)) + 1)
    log = os.path.join(RES, "pipeline_run_log.jsonl")
    results = {}
    for name, fn in STEPS[start:end]:
        if name in results and results[name]:
            continue
        ok = fn()
        results[name] = ok
        with open(log, 'a') as f:
            f.write(json.dumps({"step": name, "ok": ok, "time": time.strftime('%Y-%m-%d %H:%M:%S')}) + "\n")
        if not ok:
            print(f"!!! [{name}] 失败——管线终止（可用 run_pipeline.py {name} 24 续跑）", flush=True)
            break
    print(f"\n===== 管线结束 ok={ {k:v for k,v in results.items()} } =====", flush=True)

if __name__ == '__main__':
    main()
