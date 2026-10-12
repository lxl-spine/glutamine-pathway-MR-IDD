# ============================================================
# 步骤 07 · eQTLGen 全血扩展复制（跨平台一致性）
# ------------------------------------------------------------
# 功能: eQTLGen lead（N=31,684）Z→beta/se，Wald ratio MR + 层内 BH；
#       与 GTEx Whole_Blood 方向一致性对比
# 输入: results/lead_eqtl_eQTLGen.csv
# 输出: results/mr_results_eQTLGen_blood.csv +
#       results/replication_eQTLGen_vs_GTExWB.csv（一致率 81.3%）
# ============================================================
# 08_eqtlgen_replication.py — eQTLGen 全血扩展复制（R6）
# Zscore -> beta/se (Zhu 2016 转换), MAF 用 FinnGen af_alt 近似（登记）
# Wald ratio MR + 层内 BH-FDR, 与 GTEx Whole_Blood 方向对照
import csv, math, bisect, sys

# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

N_EQTLGEN = 31684

# 加载 GWAS
print("loading GWAS...", flush=True)
gwas = {}
with open(str(P.FINNGEN_RSID_MA)) as f:
    for line in f:
        c = line.rstrip('\n').split('\t')
        if not c[0].startswith('rs'):
            continue
        gwas[c[0]] = (c[1], c[2], float(c[3]), float(c[4]), float(c[5]), float(c[6]))
print(f"GWAS loaded {len(gwas)}", flush=True)

# eQTLGen lead（外部输入数据，位于 data\，非管线产物）
rows = []
with open(str(P.EQTLGEN_LEAD_CSV)) as f:
    rd = csv.DictReader(f, delimiter='\t')
    for r in rd:
        rs = r['SNP']
        gd = gwas.get(rs)
        if gd is None:
            continue
        ga1, ga2, gfreq, gb, gse, gp = gd
        z = float(r['Zscore'])
        n = int(r['NrSamples'])
        maf = gfreq
        if not (0 < maf < 1):
            continue
        # Z -> beta/se（标准化线性模型近似）
        var_e = 2 * maf * (1 - maf) * (n + z * z)
        b_e = z / math.sqrt(var_e)
        se_e = 1.0 / math.sqrt(var_e)
        # 等位基因对齐
        ea = r['AssessedAllele']
        if ga1 == ea:
            pass
        elif ga2 == ea:
            b_e = -b_e
        else:
            continue
        # Wald ratio + 完整 delta method SE（与 04_mr_main.py 一致）
        b_mr = gb / b_e
        se_mr = math.sqrt((gse * gse) / (b_e * b_e) + (gb * gb * se_e * se_e) / (b_e ** 4))
        zz = b_mr / se_mr
        if abs(zz) > 8.5:
            # P1-8 修复：双尾渐近需乘 2
            p_mr = 2.0 * math.exp(-zz * zz / 2.0) / (abs(zz) * math.sqrt(2.0 * math.pi))
        else:
            p_mr = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(zz) / math.sqrt(2.0))))
        rows.append({
            'gene': r['GeneSymbol'], 'probeID': r['Gene'], 'SNP': rs,
            'b_eQTL': b_e, 'se_eQTL': se_e, 'z_eQTL': z,
            'b_GWAS': gb, 'se_GWAS': gse, 'p_GWAS': gp,
            'b_MR': b_mr, 'se_MR': se_mr, 'p_MR': p_mr,
            'p_MR_log10': math.log10(max(p_mr, 1e-320)),
        })
print(f"MR rows: {len(rows)}", flush=True)

# 层内 BH-FDR（P1-11 修复：用索引排序，不用 p 值做 key——重复 p 不再互相覆盖）
m = len(rows)
order = sorted(range(m), key=lambda i: rows[i]['p_MR'])
ranked = [0.0] * m
prev = 1.0
for k in range(m - 1, -1, -1):
    idx = order[k]
    val = rows[idx]['p_MR'] * m / (k + 1)
    prev = min(prev, val)
    ranked[idx] = min(prev, 1.0)
for x, q in zip(rows, ranked):
    x['FDR_MR'] = q

with open(str(P.RESULTS_DIR / "mr_results_eQTLGen_blood.csv"), 'w', newline='', encoding='utf-8') as f:
    # rows 中已含 FDR_MR，勿重复追加（此前 bug 产生 FDR_MR/FDR_MR.1 双列）
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(rows)

sig = [x for x in rows if x['FDR_MR'] < 0.05]
print(f"eQTLGen blood FDR<0.05: {len(sig)}", flush=True)

# 与 GTEx Whole_Blood 对照
gtex_wb = {}
with open(str(P.RESULTS_DIR / "mr_results_GTEx4tissues.csv")) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['tissue'] == 'Whole_Blood':
            gtex_wb[r['gene']] = r
agree, disagree = 0, 0
with open(str(P.RESULTS_DIR / "replication_eQTLGen_vs_GTExWB.csv"), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['gene', 'eQTLGen_b_MR', 'eQTLGen_p_MR', 'eQTLGen_FDR', 'GTExWB_b_MR', 'GTExWB_FDR', 'same_direction'])
    for x in sig:
        g = gtex_wb.get(x['gene'])
        if g is None:
            continue
        b1, b2 = x['b_MR'], float(g['b_MR'])
        same = (b1 > 0) == (b2 > 0)
        agree += same; disagree += (not same)
        w.writerow([x['gene'], f"{b1:.4f}", x['p_MR'], f"{x['FDR_MR']:.4f}",
                    f"{b2:.4f}", g['FDR_MR'], 'Y' if same else 'N'])
print(f"GTExWB overlap: {agree+disagree}, same direction: {agree} ({agree/(agree+disagree)*100:.1f}%)" if (agree+disagree) else "GTExWB overlap: 0（无交集，跳过方向一致性）", flush=True)
print("DONE")
