# ============================================================
# 步骤 05 · coloc 输入准备：cis±2Mb 全 SNP + FinnGen 合并
# ------------------------------------------------------------
# 功能: 从 BESD 提取候选基因 cis 区域全量 SNP，与 FinnGen GWAS
#       等位基因对齐 → coloc 输入文件（同源数据，R3）
# 输入: results/cis_snps_full_{组织}.csv + FinnGen .ma
# 输出: results/coloc_input/{组织}_{基因}.txt（1,770 对）
# ============================================================
# 06_coloc_prep.py — 合并 eQTL cis SNP + FinnGen GWAS → coloc 输入
# 输出: results/coloc_input/{tissue}_{gene}.txt  (每基因×组织)
import csv, os, math, sys
from collections import defaultdict

# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

TISSUES = ['Adipose', 'Muscle', 'Nerve', 'WholeBlood']
TISSUE_FULL = {'Adipose': 'Adipose_Subcutaneous', 'Muscle': 'Muscle_Skeletal',
               'Nerve': 'Nerve_Tibial', 'WholeBlood': 'Whole_Blood'}
OUTDIR = P.RESULTS_DIR / "coloc_input"
OUTDIR.mkdir(parents=True, exist_ok=True)

# GTEx v8 组织样本量（公开元数据）
N_EQTL = {'Adipose': 581, 'Muscle': 706, 'Nerve': 532, 'WholeBlood': 670}
NCASE, NCONT = 51683, 353224

# 加载 GWAS (rsID -> A1,A2,freq,b,se,p)
print("loading GWAS...", flush=True)
gwas = {}
with open(str(P.FINNGEN_RSID_MA)) as f:
    for line in f:
        c = line.rstrip('\n').split('\t')
        if not c[0].startswith('rs'):
            continue
        gwas[c[0]] = (c[1], c[2], float(c[3]), float(c[4]), float(c[5]), float(c[6]))
print(f"GWAS loaded {len(gwas)}", flush=True)

total = 0
for tissue in TISSUES:
    # P0-3 修复：统一读 01b 输出的 cis_snps_full_{组织}.csv（原代码读短名 cis_snps_{tissue}.csv）
    inp = str(P.RESULTS_DIR / f"cis_snps_full_{TISSUE_FULL[tissue]}.csv")
    full = TISSUE_FULL[tissue]
    if not os.path.exists(inp):
        print(f"!! 缺少 {inp}——请先运行步骤 01b (cis_snps_extract)")
        continue
    by_gene = defaultdict(list)
    with open(inp) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            by_gene[r['gene']].append(r)
    print(f"{tissue}: {len(by_gene)} genes", flush=True)
    for g, rows in by_gene.items():
        out_rows = []
        for r in rows:
            gd = gwas.get(r['SNP'])
            if gd is None:
                continue
            ga1, ga2, gfreq, gb, gse, gp = gd
            ea1, ea2 = r['A1'], r['A2']
            b_e = float(r['b'])
            # P1-9 修复：等位基因对完全匹配检查
            if {ga1, ga2} != {ea1, ea2}:
                continue
            if ga1 == ea2:
                b_e = -b_e
            out_rows.append((r['SNP'], b_e, float(r['se']), gb, gse, gfreq))
        if len(out_rows) < 50:
            continue
        with open(str(OUTDIR / f"{full}_{g}.txt"), 'w') as o:
            o.write("snp\tb_eqtl\tse_eqtl\tb_gwas\tse_gwas\tmaf\n")
            for x in out_rows:
                o.write('\t'.join(map(str, x)) + '\n')
        total += 1
print(f"DONE {total} gene-tissue coloc input files -> {OUTDIR}")

# P2-17 修复：生成 09 (coloc.susie) 依赖的 todo_full.txt（gene×tissue 清单）
# 09 读取 results/coloc_susie/todo_full.txt，每行 "tissue gene"（空格分隔）
TODO_DIR = P.RESULTS_DIR / "coloc_susie"
TODO_DIR.mkdir(parents=True, exist_ok=True)
todo_path = str(TODO_DIR / "todo_full.txt")
with open(todo_path, 'w') as f:
    for fn in sorted(os.listdir(OUTDIR)):
        if not fn.endswith('.txt'):
            continue
        # 文件名: {TissueFull}_{Gene}.txt（TissueFull 含下划线，如 Adipose_Subcutaneous）
        stem = fn[:-4]
        # 逆匹配组织名（从已知 full 名集合）
        matched = None
        for tf in ['Adipose_Subcutaneous', 'Muscle_Skeletal', 'Nerve_Tibial', 'Whole_Blood']:
            if stem.startswith(tf + '_'):
                matched = (tf, stem[len(tf) + 1:])
                break
        if matched:
            f.write(f"{matched[0]} {matched[1]}\n")
print(f"todo_full.txt 已生成: {todo_path}")
