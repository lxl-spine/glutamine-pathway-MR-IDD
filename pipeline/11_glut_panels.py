# ============================================================
# 步骤 11 · 谷氨酰胺通路 47 基因区域面板 + LD 矩阵
# ------------------------------------------------------------
# 功能: 从 1000G EUR plink 面板提取 47 谷氨酰胺基因 ±2Mb 区域
#       面板 + 平方 LD 矩阵（供条件分析/多 IV/区域扫描）
# 输入: data/1000G_plink/1kg_eur
# 输出: results/glut_ld/{基因}_chr{chr}.{bed,bim,fam,ld}
# ============================================================
# build_glut_panels.py — 从 1kg_eur 提取 47 谷氨酰胺基因区域面板 + LD 矩阵
import pandas as pd, subprocess, os, glob
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

PLINK = str(P.PLINK)
BFILE = str(P.PANEL1000G)
GLD = str(P.RESULTS_DIR / "glut_ld")
os.makedirs(GLD, exist_ok=True)

# 基因位置：从四组织 cis_snps_full 合并提取（每基因 probe chr/bp）
# P3 修复：drop_duplicates('gene') 只保留首个组织位置可能选错（不同组织 probe_bp
# 近似相同但 chr 可能不同）——改为按基因取众数 chr + 中位 bp
frames = []
for t in ['Adipose_Subcutaneous', 'Muscle_Skeletal', 'Nerve_Tibial', 'Whole_Blood']:
    p = P.RESULTS_DIR / f"cis_snps_full_{t}.csv"
    if os.path.exists(p):
        frames.append(pd.read_csv(p, sep='\t', usecols=['gene', 'chr', 'probe_bp']))
cis = pd.concat(frames, ignore_index=True)
loc_list = {}
for g, grp in cis.groupby('gene'):
    c = grp['chr'].mode().iloc[0]
    loc_list[g] = (int(c), int(grp.loc[grp['chr'] == c, 'probe_bp'].median()))
loc = pd.DataFrame(loc_list, index=['chr', 'probe_bp']).T
print(f"{len(loc)} genes located", flush=True)

W = 2000000  # ±2Mb
for gene, row in loc.iterrows():
    c, bp = int(row['chr']), int(row['probe_bp'])
    pan = f"{gene}_chr{c}"
    bed = os.path.join(GLD, pan)
    if os.path.exists(bed + ".bed"):
        continue
    lo, hi = max(1, bp - W), bp + W
    r = subprocess.run([PLINK, "--bfile", BFILE, "--chr", str(c),
                        "--from-bp", str(lo), "--to-bp", str(hi),
                        "--make-bed", "--out", bed],
                       capture_output=True, text=True, encoding='gbk', errors='replace')
    if r.returncode != 0:
        print(f"  {gene} panel ERR: {(r.stderr or '')[-200:]}", flush=True)
    else:
        print(f"{gene} panel {os.path.exists(bed+'.bed')}", flush=True)

# LD 矩阵（每基因）
for gene, row in loc.iterrows():
    c = int(row['chr'])
    pan = f"{gene}_chr{c}"
    bed = os.path.join(GLD, pan)
    if not os.path.exists(bed + ".bed"):
        continue
    r = subprocess.run([PLINK, "--bfile", bed, "--r", "square",
                        "--out", os.path.join(GLD, f"{pan}_ld")],
                       capture_output=True, text=True, encoding='gbk', errors='replace')
    ok = os.path.exists(os.path.join(GLD, f"{pan}_ld.ld"))
    if not ok:
        print(f"  {gene} LD ERR: {(r.stderr or '')[-200:]}", flush=True)

# P2-17 修复：生成 09 (coloc.susie) 依赖的 glut_panels.txt（gene -> panel 名）
# 09 读取 results/coloc_susie/glut_panels.txt（两列: gene, panel，无表头）
WORK09 = str(P.RESULTS_DIR / "coloc_susie")
os.makedirs(WORK09, exist_ok=True)
with open(os.path.join(WORK09, "glut_panels.txt"), 'w') as f:
    for gene, row in loc.iterrows():
        c = int(row['chr'])
        f.write(f"{gene}\t{gene}_chr{c}\n")
print("glut_panels.txt 已生成:", os.path.join(WORK09, "glut_panels.txt"))
print("ALL PANELS DONE", flush=True)
