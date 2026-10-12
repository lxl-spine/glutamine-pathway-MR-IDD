# ============================================================
# 步骤 03 · MR 主分析：单 IV Wald ratio（手算，独立于 TwoSampleMR）
# ------------------------------------------------------------
# 功能: GTEx 四组织 lead cis-eQTL × FinnGen R12 IDD，Wald ratio +
#       完整 delta method SE；每组织层内 BH-FDR（R9）
# 输入: results/lead_eqtl_{组织}.csv + data/finngen_R12_*.rsid.ma
# 输出: results/mr_results_GTEx4tissues.csv（76,953 基因-组织对）
# ============================================================
# 04_mr_main.py — 单 IV Wald ratio MR 主分析（手算，独立于 TwoSampleMR）
# 暴露: GTEx v8 四组织 lead cis-eQTL (BESD 解析输出)
# 结局: FinnGen R12 IDD (rsID 版 .ma)
# 输出: results/mr_results_GTEx4tissues.csv（每组织×每基因 Wald ratio + 层内 BH-FDR）
import sys, math, csv
from collections import defaultdict

# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

TISSUES = ['Adipose', 'Muscle', 'Nerve', 'WholeBlood']
TISSUE_LABEL = {'Adipose': 'Adipose_Subcutaneous', 'Muscle': 'Muscle_Skeletal',
                'Nerve': 'Nerve_Tibial', 'WholeBlood': 'Whole_Blood'}

# 1. 读 FinnGen GWAS (rsID 版)
print("loading GWAS...", flush=True)
gwas = {}
with open(str(P.FINNGEN_RSID_MA), 'rt') as f:
    for line in f:
        c = line.rstrip('\n').split('\t')
        rs = c[0]
        if not rs.startswith('rs'):
            continue
        gwas[rs] = (c[1], c[2], float(c[4]), float(c[5]), float(c[6]))  # A1,A2,b,se,p
print(f"GWAS loaded: {len(gwas)} rsid SNPs", flush=True)

def bh_fdr(ps):
    m = len(ps)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: ps[i])
    ranked = [0.0] * m
    prev = 1.0
    for k in range(m - 1, -1, -1):
        idx = order[k]
        val = ps[idx] * m / (k + 1)
        prev = min(prev, val)
        ranked[idx] = prev
    return [min(x, 1.0) for x in ranked]

all_rows = []
for tissue in TISSUES:
    lead = str(P.RESULTS_DIR / f"lead_eqtl_{tissue}.csv")
    rows = []
    with open(lead, 'rt') as f:
        rd = csv.DictReader(f, delimiter='\t')
        for r in rd:
            rows.append(r)
    print(f"{tissue}: {len(rows)} genes", flush=True)
    mr_rows = []
    n_miss = 0
    n_flip = 0
    for r in rows:
        rs = r['SNP']
        g = gwas.get(rs)
        if g is None:
            n_miss += 1
            continue
        ga1, ga2, gb, gse, gp = g
        ea1, ea2 = r['A1'], r['A2']
        b_e = float(r['b'])
        se_e = float(r['se'])
        # 等位基因对检查（P1-9 修复）：{ga1,ga2} 必须与 {ea1,ea2} 完全匹配
        if {ga1, ga2} != {ea1, ea2}:
            n_miss += 1
            continue
        if ga1 == ea2:
            b_e = -b_e
            n_flip += 1
        # Wald ratio + 完整 delta method SE
        # se(MR) = sqrt( se_Y^2/beta_X^2 + beta_Y^2*se_X^2/beta_X^4 )
        # 一阶近似 se_Y/|beta_X| 遗漏第二项，SE 系统性低估（中位 23%，最大 2.5 倍）
        b_mr = gb / b_e
        se_mr = math.sqrt((gse * gse) / (b_e * b_e) + (gb * gb * se_e * se_e) / (b_e ** 4))
        z = b_mr / se_mr
        if abs(z) > 8.5:
            # 渐近近似（P1-8 修复：双尾需乘 2）: p ≈ 2*exp(-z^2/2) / (|z| * sqrt(2*pi))
            p_mr = 2.0 * math.exp(-z * z / 2.0) / (abs(z) * math.sqrt(2.0 * math.pi))
        else:
            p_mr = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z) / math.sqrt(2.0))))
        p_e = float(r['p'])
        if p_e == 0:
            p_e = float(r['p'])
        mr_rows.append({
            'tissue': TISSUE_LABEL[tissue], 'probeID': r['probeID'], 'gene': r['gene'],
            'chr': r['chr'], 'probe_bp': r['probe_bp'],
            'SNP': rs, 'SNP_bp': r['SNP_bp'],
            'eA1': ea1, 'eA2': ea2, 'gA1': ga1, 'gA2': ga2,
            'b_eQTL': b_e, 'se_eQTL': se_e, 'p_eQTL': p_e,
            'b_GWAS': gb, 'se_GWAS': gse, 'p_GWAS': gp,
            'b_MR': b_mr, 'se_MR': se_mr, 'p_MR': p_mr,
            'p_MR_log10': math.log10(max(p_mr, 1e-320)),
            'flipped': 1 if b_e != float(r['b']) else 0
        })
    # 层内 BH-FDR (组织内全部基因, R9)
    ps = [x['p_MR'] for x in mr_rows]
    fdr = bh_fdr(ps)
    for x, f in zip(mr_rows, fdr):
        x['FDR_MR'] = f
    print(f"  {len(mr_rows)} MR rows, miss {n_miss}, flip {n_flip}", flush=True)
    all_rows.extend(mr_rows)

if not all_rows:
    raise SystemExit("ERROR: MR 主分析无任何有效行——请检查 01 输出与 GWAS 是否匹配")

out = str(P.RESULTS_DIR / "mr_results_GTEx4tissues.csv")
with open(out, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(all_rows)
print(f"DONE {len(all_rows)} rows -> {out}")
