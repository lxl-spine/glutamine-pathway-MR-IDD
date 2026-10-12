# ============================================================
# 步骤 10 · 多 IV MR 敏感性矩阵（IVW/Egger/加权中位数/Steiger）
# ------------------------------------------------------------
# 功能: 对 8 个双证据基因 + GFPT1/SLC25A13，clump 阈值
#       r²∈{0.001,0.01,0.1} 多 IV 分析；Steiger 方向过滤
# 输入: glut_ld 面板（plink LD）+ mr_results
# 输出: results/multi_iv_mr_sensitivity.csv（全部 Steiger_dir_ok=True）
# ============================================================
# run_multi_iv_mr.py — 多 IV MR 敏感性分析（IVW/Egger/加权中位数/Steiger）
# 对象：8 个双证据基因×组织 + GFPT1/SLC25A13（WB 重跑验证）
import pandas as pd, numpy as np, os
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

GLD = str(P.RESULTS_DIR / "glut_ld")
RES = str(P.RESULTS_DIR)

# 面板名（非 _qc）
def pan_name(gene):
    for f in os.listdir(GLD):
        if f.endswith(".bim") and not f.endswith("_qc.bim") and f.startswith(gene + "_chr"):
            return f[:-4]
    return None

def load_panel(gene):
    pan = pan_name(gene)
    bim = pd.read_csv(os.path.join(GLD, pan + ".bim"), sep=r'\s+', header=None,
                      names=['chr','SNP','cm','bp','A1','A2']).drop_duplicates('SNP')
    ld = pd.read_csv(os.path.join(GLD, pan + "_ld.ld"), sep='\t', header=None)
    return bim, ld.values

def clump(snps, ld, r2_thresh=0.01, window=1e6):
    """snps: DataFrame(SNP, bp, p) 按 p 升序；贪心 clump"""
    kept = []
    for _, s in snps.iterrows():
        drop = False
        for k in kept:
            if abs(s['bp'] - k['bp']) < window:
                i = snps.index[snps['SNP'] == k['SNP']][0]
                j = snps.index[snps['SNP'] == s['SNP']][0]
                if i < len(ld) and j < len(ld) and ld[i, j] ** 2 > r2_thresh:
                    drop = True
                    break
        if not drop:
            kept.append(s)
    return pd.DataFrame(kept)

def ivw(b, se):
    w = 1 / se ** 2
    b_ivw = np.sum(w * b) / np.sum(w)
    se_ivw = np.sqrt(1 / np.sum(w))
    Q = np.sum(w * (b - b_ivw) ** 2)
    tau2 = max(0, (Q - (len(b) - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w)))
    se_rand = np.sqrt(1 / np.sum(w) * (1 + tau2 * np.sum(w))) if tau2 > 0 else se_ivw
    z = b_ivw / se_rand
    return b_ivw, se_rand, 2 * (1 - norm_cdf(abs(z)))

def norm_cdf(x):
    return 0.5 * (1 + erf(x / np.sqrt(2)))

def erf(x):
    # Abramowitz-Stegun 近似
    sign = 1 if x >= 0 else -1
    x = abs(x)
    t = 1 / (1 + 0.3275911 * x)
    y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * np.exp(-x * x)
    return sign * y

def egger(b_x, b_y, se_y):
    """标准 MR-Egger（P0-6 修复）：b_Y ~ a + g*b_X 加权回归，权重 1/se_Y^2
    斜率 g = 因果效应；截距 a = 定向水平多效性（检验其是否为 0）"""
    w = 1 / se_y ** 2
    X = np.column_stack([np.ones(len(b_x)), b_x])
    W = np.diag(w)
    XtWX = X.T @ W @ X
    XtWy = X.T @ W @ b_y
    try:
        beta = np.linalg.solve(XtWX, XtWy)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan, np.nan
    resid = b_y - X @ beta
    sigma2 = np.sum(w * resid ** 2) / (len(b_x) - 2)
    cov = sigma2 * np.linalg.inv(XtWX)
    intercept = beta[0]; slope = beta[1]
    se_i = np.sqrt(cov[0, 0]); se_s = np.sqrt(cov[1, 1])
    p_i = 2 * (1 - norm_cdf(abs(intercept / se_i)))
    p_s = 2 * (1 - norm_cdf(abs(slope / se_s)))
    return slope, se_s, p_s, (intercept, p_i)

def wm(b, se):
    # 加权中位数（Bowden 2016）
    w = 1 / se ** 2
    ord = np.argsort(b)
    b_s = b[ord]; w_s = w[ord]
    sw = np.cumsum(w_s) / w_s.sum()
    k = np.searchsorted(sw, 0.5)
    return b_s[k]

def steiger(b_eqtl, se_eqtl, N_eqtl, b_gwas, se_gwas, N_gwas):
    # R2 ≈ Z^2 / (Z^2 + N)
    z_e = b_eqtl / se_eqtl; z_g = b_gwas / se_gwas
    r2e = z_e ** 2 / (z_e ** 2 + N_eqtl)
    r2g = z_g ** 2 / (z_g ** 2 + N_gwas)
    return r2e.sum(), r2g.sum()

# 分析对象（8 双证据 + 2 核心）
TARGETS = [
    ('Whole_Blood', 'GFPT1'), ('Whole_Blood', 'SLC25A13'),
    ('Muscle_Skeletal', 'SLC25A12'), ('Nerve_Tibial', 'NFS1'),
    ('Nerve_Tibial', 'SLC38A1'), ('Whole_Blood', 'SLC7A7'),
    ('Adipose_Subcutaneous', 'SLC38A3'), ('Whole_Blood', 'GAD1'),
    ('Whole_Blood', 'ASL'), ('Muscle_Skeletal', 'GCLM'),
]
N_EQTL = {'Adipose_Subcutaneous': 581, 'Muscle_Skeletal': 706, 'Nerve_Tibial': 532, 'Whole_Blood': 670}
N_GWAS = 404907

# 预读 MR 结果（含 GWAS 效应）
mr = pd.read_csv(P.RESULTS_DIR / "mr_results_GTEx4tissues.csv", sep='\t')

rows = []
for tissue, gene in TARGETS:
    cis = pd.read_csv(P.RESULTS_DIR / f"cis_snps_full_{tissue}.csv", sep='\t')
    sub = cis[(cis['gene'] == gene)].copy()
    # cis-eQTL 显著 SNP（分级阈值）
    thresh_used = None
    for thr in [5e-8, 5e-6, 5e-5, 1e-4]:
        iv_cand = sub[sub['p'] < thr].copy()
        if len(iv_cand) >= 2:
            thresh_used = f'{thr:.0e}'
            break
    if thresh_used is None:
        rows.append({'tissue': tissue, 'gene': gene, 'note': 'IV<2 (p<1e-4)'})
        continue
    bim, ld = load_panel(gene)
    iv_cand = iv_cand.merge(bim[['SNP', 'bp']], on='SNP', how='left').dropna(subset=['bp'])
    # LD clump（r²<0.001 标准；不足则放宽 0.01）
    idx_map = {s: i for i, s in enumerate(bim['SNP'])}
    ld_small = np.zeros((len(iv_cand), len(iv_cand)))
    for a in range(len(iv_cand)):
        for b_ in range(len(iv_cand)):
            ia, ib = idx_map[iv_cand['SNP'].iloc[a]], idx_map[iv_cand['SNP'].iloc[b_]]
            ld_small[a, b_] = ld[ia, ib]
    def do_clump(r2_cut):
        ivs = []
        used = set()
        order = iv_cand.sort_values('p').index.tolist()
        for i in order:
            if i in used: continue
            ivs.append(i)
            for j in order:
                if j in used or i == j: continue
                if abs(iv_cand.loc[i, 'bp'] - iv_cand.loc[j, 'bp']) < 1e6 and ld_small[iv_cand.index.get_loc(i), iv_cand.index.get_loc(j)] ** 2 > r2_cut:
                    used.add(j)
        return ivs
    r2_used = None
    ivs = []
    # P0-6 修复：clump 阈值固定 {0.001, 0.01, 0.1}，禁止放宽到 r²=0.8
    #（r²=0.8 会把高度连锁 SNP 当独立 IV，违反 MR 独立性假设）
    for r2c in [0.001, 0.01, 0.1]:
        ivs = do_clump(r2c)
        if len(ivs) >= 2:
            r2_used = r2c
            break
    if r2_used is None or len(ivs) < 2:
        rows.append({'tissue': tissue, 'gene': gene, 'note': 'IV<2 even r2<0.1 (single LD block)'})
        continue
    iv_df = iv_cand.loc[ivs]
    # GWAS 效应（从 gwas_summary 全量提取；对齐用 cis 的 A1 与 GWAS 的 a1）
    gwas = pd.read_csv(P.SMR_WORK / "gwas_summary.txt", sep='\t', dtype={'snp': str})
    gwas_eff = gwas[gwas['snp'].isin(iv_df['SNP'])][['snp', 'a1', 'a2', 'beta', 'se', 'p']].rename(
        columns={'a1': 'gwas_a1', 'a2': 'gwas_a2', 'beta': 'gwas_beta', 'se': 'gwas_se', 'p': 'gwas_p'})
    iv_df = iv_df.merge(gwas_eff, left_on='SNP', right_on='snp', how='left').dropna(subset=['gwas_beta', 'gwas_se'])
    if len(iv_df) < 2:
        rows.append({'tissue': tissue, 'gene': gene, 'note': 'IV<2 after GWAS merge'})
        continue
    # GWAS 效应对齐到面板等位基因方向（bim 的 A1/A2 为基准）
    pan_alleles = bim.set_index('SNP')[['A1', 'A2']].rename(columns={'A1': 'panA1', 'A2': 'panA2'})
    iv_df = iv_df.merge(pan_alleles, on='SNP', how='left')
    align = np.where(iv_df['gwas_a1'] == iv_df['panA1'], 1,
            np.where(iv_df['gwas_a1'] == iv_df['panA2'], -1, 0))
    if (align == 0).any():
        rows.append({'tissue': tissue, 'gene': gene, 'note': f'{int((align==0).sum())} IV unaligned to panel'})
        continue
    # P0-6 修复：eQTL 效应也必须对齐到面板方向（原代码漏掉此项，
    # 若 eQTL A1 ≠ 面板 A1 则 Wald ratio 符号/量值错误）
    align_e = np.where(iv_df['A1'] == iv_df['panA1'], 1,
              np.where(iv_df['A1'] == iv_df['panA2'], -1, 0))
    if (align_e == 0).any():
        rows.append({'tissue': tissue, 'gene': gene, 'note': f'{int((align_e==0).sum())} IV eQTL unaligned to panel'})
        continue
    n_iv = len(iv_df)
    b_g = iv_df['gwas_beta'].values * align      # 面板方向 GWAS 效应
    b_e = iv_df['b'].values * align_e            # 面板方向 eQTL 效应（P0-6 修复）
    se_g = iv_df['gwas_se'].values
    se_e = iv_df['se'].values
    theta = b_g / b_e                             # Wald ratio（MR 因果效应）
    se_theta = se_g / np.abs(b_e)                 # delta 一阶近似
    b_ivw, se_ivw, p_ivw = ivw(theta, se_theta)
    if n_iv >= 3:
        # P0-6 修复：标准 MR-Egger —— b_gwas ~ a + g*b_eqtl（b_X 是 eQTL 效应，非 se）
        egg = egger(b_e, b_g, se_g)
        egger_b, egger_p, egger_i, egger_ip = egg[0], egg[2], egg[3][0], egg[3][1]
    else:
        egger_b = egger_p = egger_i = egger_ip = np.nan
    b_wm = wm(theta, se_theta)
    r2e, r2g = steiger(b_e, se_e, N_EQTL[tissue], b_g, se_g, N_GWAS)
    snps = ','.join(iv_df['SNP'].tolist())
    rows.append({'tissue': tissue, 'gene': gene, 'n_IV': n_iv, 'thresh': thresh_used, 'clump_r2': r2_used,
                 'IVW_b': b_ivw, 'IVW_se': se_ivw, 'IVW_p': p_ivw,
                 'Egger_b': egger_b, 'Egger_p': egger_p, 'Egger_intercept': egger_i,
                 'Egger_intercept_p': egger_ip,
                 'WM_b': b_wm, 'R2_expr': r2e, 'R2_gwas': r2g,
                 'Steiger_dir_ok': r2e > r2g, 'SNPs': snps,
                 'Wald_ratio_top': theta[0], 'MR_single_p': iv_df['gwas_p'].iloc[0]})

out = pd.DataFrame(rows)
out.to_csv(P.RESULTS_DIR / "multi_iv_mr_sensitivity.csv", sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20)
print(out.to_string(index=False))
