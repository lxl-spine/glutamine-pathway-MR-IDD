# ============================================================
# 步骤 08 · SMR-HEIDI 多效性判别（R8）
# ------------------------------------------------------------
# 功能: GTEx allpairs 区域全 SNP × IDD GWAS → SMR + HEIDI 检验；
#       P_HEIDI<0.05 降级为"关联但方向不确定"
# 输入: smr.exe + GTEx 区域 eQTL + 1000G EUR LD + FinnGen
# 输出: results/smr_glut_17genes_4tissues.csv
# ============================================================
# run_smr_batch_regional2.py — 区域 GWAS 子集 + 区域面板 SMR-HEIDI（17 基因 × 4 组织）
import pandas as pd, subprocess, os, glob, sys

# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

SMR = str(P.SMR_EXE)
OUT = str(P.SMR_WORK)
GWAS = os.path.join(OUT, "gwas_summary.txt")
GLD = str(P.RESULTS_DIR / "glut_ld")
TISSUES = ['Adipose_Subcutaneous', 'Muscle_Skeletal', 'Nerve_Tibial', 'Whole_Blood']
GENES = ['GFPT1', 'SLC25A13', 'ARG1', 'ASL', 'GAD1', 'GCLM', 'GLUD1', 'GMPS', 'GSS',
         'NFS1', 'SLC25A12', 'SLC38A1', 'SLC38A3', 'SLC7A5', 'SLC7A6', 'SLC7A7', 'SUCLG1']

pan_map = {}
for f in os.listdir(GLD):
    if f.endswith(".bim") and not f.endswith("_qc.bim"):
        g = f.replace(".bim", "")
        parts = g.rsplit("_chr", 1)
        if len(parts) == 2:
            pan_map[parts[0]] = g
PAN_DIR = str(P.SMR_WORK / "panels")

# ============================================================
# P2-17 修复：管线自足输入准备（原代码依赖手工放置的 gwas_summary /
# {t}_besd / panels，均无生成步骤）
# 1) gwas_summary.txt ← FinnGen R12 .ma（a1,a2,beta,se,p 列名）
# ============================================================
def prepare_inputs():
    os.makedirs(PAN_DIR, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    # 1. gwas_summary.txt
    if not os.path.exists(GWAS):
        print("生成 gwas_summary.txt（从 FinnGen .ma）...", flush=True)
        rows = []
        with open(str(P.FINNGEN_RSID_MA)) as f:
            for line in f:
                c = line.rstrip('\n').split('\t')
                if not c[0].startswith('rs'):
                    continue
                rows.append((c[0], c[1], c[2], c[4], c[5], c[6]))  # snp,a1,a2,beta,se,p
        pd.DataFrame(rows, columns=['snp', 'a1', 'a2', 'beta', 'se', 'p']).to_csv(GWAS, sep='\t', index=False)
        print(f"gwas_summary.txt: {len(rows)} SNPs", flush=True)
    # 2. {t}_besd 组织 BESD 前缀（从 data/GTEx_v8_westlake 复制）
    for t in TISSUES:
        src = str(P.besd_prefix(t))
        dst = os.path.join(OUT, f"{t}_besd")
        for ext in ['.besd', '.epi', '.esi']:
            if os.path.exists(src + ext) and not os.path.exists(dst + ext):
                import shutil
                shutil.copyfile(src + ext, dst + ext)
    # 3. SMR 面板（从 results/glut_ld 复制到 SMR_WORK/panels）
    for f in os.listdir(GLD):
        if f.endswith(('.bed', '.bim', '.fam')) and not '_qc.' in f:
            src = os.path.join(GLD, f)
            dst = os.path.join(PAN_DIR, f)
            if not os.path.exists(dst):
                import shutil
                shutil.copyfile(src, dst)
    print("SMR 输入准备完成", flush=True)

prepare_inputs()

def run(cmd, tag):
    r = subprocess.run(cmd, capture_output=True, text=True, encoding='gbk', errors='replace')
    if r.returncode != 0:
        print(f"  ERR {tag}: {(r.stderr or r.stdout or '')[-200:]}", flush=True)
    return r

# 1) 按基因提取区域 GWAS 子集（面板 SNP 交集）
gwas_full = pd.read_csv(GWAS, sep='\t', dtype={'snp': str})
gwas_full = gwas_full.set_index('snp')
for gene in GENES:
    pan = pan_map[gene]
    gfile = os.path.join(OUT, f"gwas_{gene}.txt")
    if os.path.exists(gfile):
        continue
    bim = pd.read_csv(os.path.join(GLD, f"{pan}.bim"), sep=r'\s+', header=None, names=['chr', 'SNP', 'cm', 'bp', 'a1', 'a2'])
    snps = set(bim['SNP'])
    sub = gwas_full.loc[gwas_full.index.isin(snps)]
    sub.to_csv(gfile, sep='\t')
    print(f"gwas_{gene}: {len(sub)} SNPs", flush=True)

# 2) 逐基因 × 组织 SMR
done = 0
for t in TISSUES:
    besd_pre = os.path.join(OUT, f"{t}_besd")
    for gene in GENES:
        if gene not in pan_map:
            continue
        out = os.path.join(OUT, f"{t}_{gene}_smr")
        if os.path.exists(out + ".smr") and os.path.getsize(out + ".smr") > 100:
            done += 1
            continue
        pan = os.path.join(PAN_DIR, pan_map[gene])
        r = run([SMR, '--bfile', pan, '--gwas-summary', os.path.join(OUT, f"gwas_{gene}.txt"),
                 '--beqtl-summary', besd_pre, '--smr',
                 '--peqtl-smr', '0.05', '--peqtl-heidi', '0.05',
                 '--disable-freq-ck', '--out', out, '--thread-num', '2'],
                f"{t}_{gene}")
        if os.path.exists(out + ".smr") and os.path.getsize(out + ".smr") > 100:
            done += 1
            print(f"  OK {t}_{gene}", flush=True)
print(f"ALL DONE ({done})", flush=True)

# ============================================================
# P2 修复（2026-09-22）：聚合 SMR-HEIDI 结果 → results/smr_glut_17genes_4tissues.csv
# （21 号证据矩阵依赖此文件；原代码只产出 SMR_WORK/{t}_{gene}_smr.smr 未聚合）
# ============================================================
agg = []
for t in TISSUES:
    for gene in GENES:
        if gene not in pan_map:
            continue
        out = os.path.join(OUT, f"{t}_{gene}_smr")
        if not os.path.exists(out + ".smr"):
            continue
        try:
            df = pd.read_csv(out + ".smr", sep=r"\s+")
        except Exception as e:
            print(f"  PARSE-ERR {t}_{gene}: {e}", flush=True)
            continue
        if df.empty:
            continue
        top = df.iloc[0]
        agg.append({
            'tissue': t, 'gene': gene,
            'topSNP': top['topSNP'], 'A1': top['A1'], 'A2': top['A2'],
            'b_eQTL': top['b_eQTL'], 'p_eQTL': top['p_eQTL'],
            'b_GWAS': top['b_GWAS'], 'p_GWAS': top['p_GWAS'],
            'b_SMR': top['b_SMR'], 'se_SMR': top['se_SMR'], 'p_SMR': top['p_SMR'],
            'p_HEIDI': top['p_HEIDI'], 'nsnp_HEIDI': top['nsnp_HEIDI'],
            'HEIDI_pass': top['p_HEIDI'] >= 0.05,
        })
agg_df = pd.DataFrame(agg)
agg_out = str(P.RESULTS_DIR / "smr_glut_17genes_4tissues.csv")
agg_df.to_csv(agg_out, sep="\t", index=False)
print(f"SMR 聚合: {len(agg_df)} 行 -> {agg_out}", flush=True)
