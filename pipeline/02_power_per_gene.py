# ============================================================
# 步骤 02 · 功效门槛预计算（§1.6 逐基因×组织，实际数据驱动）
# ------------------------------------------------------------
# 功能: 用实际 beta/se（完整 delta method）计算每基因×组织可检测
#       最小 OR；OR>1.30 记为"功效不足"（R5 分层前置判定）
# 输入: results/lead_eqtl_{Adipose,Muscle,Nerve,WholeBlood}.csv
#       （步骤 01 输出，每基因一行）+ FinnGen R12 .ma
# 输出: results/power_per_gene_all.csv
# 依赖: 01（本脚本不依赖 03，修正原依赖颠倒问题）
# ============================================================
# step1_power_per_gene.py — v4.0 R5 逐基因×组织功效门槛
# 单 IV Wald ratio，完整 delta method SE（与 03 主分析一致）
# se_MR = sqrt(se_Y^2/beta_X^2 + beta_Y^2*se_X^2/beta_X^4)
# 可检测 OR = exp( (z_alpha + z_beta) * se_MR )，z 对应 alpha=5e-5 双侧、power=0.80
import csv, math, sys
from math import sqrt, exp, log

# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

ALPHA = 5e-5
POWER = 0.80
OR_CUT = 1.30
TISSUES = {'Adipose': 'Adipose_Subcutaneous', 'Muscle': 'Muscle_Skeletal',
           'Nerve': 'Nerve_Tibial', 'WholeBlood': 'Whole_Blood'}

# z_alpha = qnorm(1-ALPHA/2) 双侧; z_beta = qnorm(POWER)
def qnorm(p):
    # Acklam 有理逼近
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p <= 0 or p >= 1:
        return float('nan')
    if p < plow:
        q = sqrt(-2 * log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    elif p <= phigh:
        q = p - 0.5
        r = q * q
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    else:
        q = sqrt(-2 * log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)

z_alpha = qnorm(1 - ALPHA / 2)
z_beta = qnorm(POWER)
K = z_alpha + z_beta
print(f"z_alpha={z_alpha:.4f} z_beta={z_beta:.4f} K={K:.4f}")

# 1. 读 FinnGen GWAS (rsID 版) — 与 03 主分析完全同源
gwas = {}
with open(str(P.FINNGEN_RSID_MA), 'rt') as f:
    for line in f:
        c = line.rstrip('\n').split('\t')
        if not c[0].startswith('rs'):
            continue
        gwas[c[0]] = (c[1], c[2], float(c[4]), float(c[5]), float(c[6]))  # A1,A2,b,se,p
print(f"GWAS loaded: {len(gwas)} rsid SNPs", flush=True)

# 2. 逐组织读 lead eQTL（01 输出），等位基因对齐同 03
out_rows = []
n_miss = 0; n_flip = 0
for t_short, t_full in TISSUES.items():
    lead = str(P.RESULTS_DIR / f"lead_eqtl_{t_short}.csv")
    with open(lead, 'rt') as f:
        rd = csv.DictReader(f, delimiter='\t')
        for r in rd:
            g = gwas.get(r['SNP'])
            if g is None:
                n_miss += 1
                continue
            ga1, ga2, gb, gse, gp = g
            ea1, ea2 = r['A1'], r['A2']
            b_e = float(r['b'])
            if ga1 == ea1:
                pass
            elif ga1 == ea2:
                b_e = -b_e
                n_flip += 1
            else:
                n_miss += 1
                continue
            se_e = float(r['se'])
            if b_e == 0 or se_e <= 0 or gse <= 0:
                continue
            # 完整 delta method SE（与 03 完全一致）
            se_mr = math.sqrt((gse * gse) / (b_e * b_e) + (gb * gb * se_e * se_e) / (b_e ** 4))
            detect_OR = exp(K * se_mr)
            power_ok = detect_OR <= OR_CUT
            out_rows.append({
                'tissue': t_full, 'gene': r['gene'], 'probeID': r['probeID'],
                'SNP': r['SNP'], 'b_eQTL': b_e, 'se_eQTL': se_e, 'p_eQTL': r['p'],
                'se_GWAS': gse, 'p_MR': gp,
                'detect_OR': detect_OR, 'power_ok_OR1.30': power_ok
            })
print(f"miss={n_miss} flip={n_flip} rows={len(out_rows)}", flush=True)

if not out_rows:
    raise SystemExit("ERROR: 无有效输入。请先运行步骤 01 生成 lead_eqtl_{组织}.csv")

with open(str(P.RESULTS_DIR / "power_per_gene_all.csv"), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(out_rows)

from collections import Counter
for t in ['Adipose_Subcutaneous', 'Muscle_Skeletal', 'Nerve_Tibial', 'Whole_Blood']:
    sub = [x for x in out_rows if x['tissue'] == t]
    if not sub:
        continue
    ok = sum(1 for x in sub if x['power_ok_OR1.30'])
    print(f"{t}: n={len(sub)} power_ok={ok} ({ok/len(sub)*100:.1f}%)")
