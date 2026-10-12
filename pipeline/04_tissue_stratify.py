# ============================================================
# 步骤 04 · 组织分层（R5 前置锁定规则 + R1 lead SNP 共享矩阵）
# ------------------------------------------------------------
# 功能: 共享候选（≥2 组织 FDR<0.05 同向且功效足够）/ 组织特异 /
#       血特异 / 方向冲突 四类；lead SNP 跨组织共享矩阵（R1）
# 输入: mr_results_GTEx4tissues.csv + power_per_gene_all.csv
# 输出: results/tier_shared.csv / tier_tissue_specific.csv /
#       tier_blood_specific.csv / tier_conflict.csv
# ============================================================
# 05_tissue_stratify.py — v4.0 R1/R5 组织分层（分析前锁定规则）
# 共享候选: ≥2 组织 FDR<0.05 且方向一致，且全部功效足够
# 组织特异候选: 仅 1 组织 FDR<0.05（功效足够）
# 血液特异: 仅 Whole_Blood 显著
# 同时输出 lead SNP 共享矩阵（R1: 共享 lead SNP 的组织对 = 统计非独立）
import csv, sys
from collections import defaultdict

# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

TISSUES = ['Adipose_Subcutaneous', 'Muscle_Skeletal', 'Nerve_Tibial', 'Whole_Blood']

# 读 MR 结果
# P0-1 修复：01 已保证每基因×组织唯一 lead SNP 行；此处显式校验，
# 若出现重复 gene+tissue（非法输入）则按 p_MR 最小保留，而非字典覆盖"最后一行"
mr = defaultdict(dict)   # gene -> tissue -> row
dup = 0
with open(str(P.RESULTS_DIR / "mr_results_GTEx4tissues.csv")) as f:
    rd = csv.DictReader(f, delimiter='\t')
    for r in rd:
        old = mr[r['gene']].get(r['tissue'])
        if old is None:
            mr[r['gene']][r['tissue']] = r
        else:
            dup += 1
            # 保留 p_MR 更小的行（与 01 lead 定义一致：最显著 SNP）
            if float(r['p_MR']) < float(old['p_MR']):
                mr[r['gene']][r['tissue']] = r
if dup:
    print(f"WARN: 检测到 {dup} 个重复 gene×tissue 行（已按 p_MR 最小保留）")

# 读功效门槛
pw = {}
with open(str(P.RESULTS_DIR / "power_per_gene_all.csv")) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        pw[(r['tissue'], r['gene'])] = r['power_ok_OR1.30'] == 'True'

genes = sorted(mr.keys())
print(f"{len(genes)} genes with MR results")

shared, conflict, tissue_spec, blood_spec = [], [], [], []
lead_matrix = {}   # gene -> {tissue: rsid}

for g in genes:
    sig = {}
    for t in TISSUES:
        r = mr[g].get(t)
        if r is None:
            continue
        fdr = float(r['FDR_MR'])
        pw_ok = pw.get((t, g), False)
        if fdr < 0.05 and pw_ok:
            sig[t] = r
    if len(sig) == 0:
        continue
    # lead SNP 记录
    snps = {t: r['SNP'] for t, r in sig.items()}
    lead_matrix[g] = snps
    # 方向
    signs = {t: float(r['b_MR']) > 0 for t, r in sig.items()}
    uniq_signs = set(signs.values())
    if len(sig) >= 2 and len(uniq_signs) == 1:
        shared.append((g, sig, snps))
    elif len(sig) >= 2 and len(uniq_signs) > 1:
        # R5 方向冲突类别（此前被静默丢弃）
        conflict.append((g, sig, snps))
    elif len(sig) == 1:
        t0 = list(sig.keys())[0]
        if t0 == 'Whole_Blood':
            blood_spec.append((g, sig, snps))
        else:
            tissue_spec.append((g, sig, snps))

print(f"shared (>=2 tissues, same direction): {len(shared)}")
print(f"conflict (>=2 tissues, opposite directions): {len(conflict)}")
print(f"tissue-specific: {len(tissue_spec)}")
print(f"blood-specific: {len(blood_spec)}")

# 输出分层表
def dump(fname, items):
    with open(str(P.RESULTS_DIR / fname), 'w', newline='', encoding='utf-8') as f:
        f.write("gene\ttissues\tn_tissues\tdirection\tlead_SNP_per_tissue\tb_MR_per_tissue\tFDR_per_tissue\n")
        for g, sig, snps in items:
            ts = sorted(sig.keys())
            bstr = ';'.join(f"{t}={float(sig[t]['b_MR']):.3f}" for t in ts)
            fstr = ';'.join(f"{t}={float(sig[t]['FDR_MR']):.4f}" for t in ts)
            sstr = ';'.join(f"{t}={snps[t]}" for t in ts)
            signs = {t: '+' if float(sig[t]['b_MR']) > 0 else '-' for t in ts}
            dstr = ''.join(signs[t] for t in ts)
            f.write(f"{g}\t{','.join(ts)}\t{len(ts)}\t{dstr}\t{sstr}\t{bstr}\t{fstr}\n")

dump('tier_shared.csv', shared)
dump('tier_conflict.csv', conflict)
dump('tier_tissue_specific.csv', tissue_spec)
dump('tier_blood_specific.csv', blood_spec)

# 共享 lead SNP 矩阵（R1）
print("\n--- shared candidates: lead SNP 共享检查 ---")
dup_snp = 0
for g, sig, snps in shared:
    ts = sorted(sig.keys())
    rsids = [snps[t] for t in ts]
    if len(set(rsids)) < len(rsids):
        dup_snp += 1
        print(f"  {g}: {' vs '.join(ts)} share lead SNP {set(rsids)}")
print(f"shared genes with duplicate lead SNP: {dup_snp}/{len(shared)}")
print("DONE")
