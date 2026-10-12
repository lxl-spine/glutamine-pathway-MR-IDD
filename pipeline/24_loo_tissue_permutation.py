# -*- coding: utf-8 -*-
"""Leave-one-out enrichment + tissue-stratified MR enrichment.

Same convention as 17_enrichment_permutation.py:
  - target set  : 47 glutamine-pathway genes
  - outcome     : gene-level min p_MR < 0.05 (layer-internal), observed = 17
  - permutation : instrument-strength matched (10 quantile layers of lead p_eQTL),
                  background = genome-wide genes not in the 47-gene set
  - leave-one-out: re-run enrichment dropping each of the 47 genes (N=2000)
  - tissue-stratified: per tissue, enrich among genes with cis-instruments in
                  that tissue, strength layers computed within tissue (N=1000)
Outputs:
  results/leave_one_out_enrichment.csv
  results/tissue_enrichment_MR_permutation.csv
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
# ---- 统一路径入口（pipeline/config/paths.py；换电脑无需改本文件）----
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P

MR_FILE = P.RESULTS_DIR / "mr_results_GTEx4tissues.csv"
GENES_FILE = P.GLUTAMINE_GENES
OUT_LOO = P.RESULTS_DIR / "leave_one_out_enrichment.csv"
OUT_TIS = P.RESULTS_DIR / "tissue_enrichment_MR_permutation.csv"

rng = np.random.default_rng(42)

mr = pd.read_csv(MR_FILE, sep="\t")
mr.columns = [c.strip() for c in mr.columns]
mr["p_MR"] = mr["p_MR"].astype(float)
mr["p_eQTL"] = mr["p_eQTL"].astype(float)

genes47 = [g.strip() for g in GENES_FILE.read_text(encoding="utf-8").splitlines() if g.strip()]
genes47 = [g for g in genes47 if g in set(mr["gene"])]
print(f"target genes: {len(genes47)}")

all_genes = mr["gene"].unique()
bg_genes = np.array([g for g in all_genes if g not in set(genes47)])
print(f"background genes: {len(bg_genes)}")

# ---- per-gene instruments: lead eQTL strength (most significant tissue p_eQTL) ----
def lead_strength(sub):
    return sub.groupby("gene")["p_eQTL"].min()

strength = lead_strength(mr)                     # all genes, cross-tissue lead
qs = strength.quantile(np.linspace(0, 1, 11)).values

def layer_of(p):
    for i in range(10):
        if p <= qs[i + 1] + 1e-300:
            return i
    return 9

layer_all = strength.apply(layer_of)

# per-gene significance (gene-level min p_MR < 0.05)
minp = mr.groupby("gene")["p_MR"].min()
sig_all = (minp < 0.05).astype(int)
sig_val = sig_all.to_dict()          # gene -> 0/1
layer_val = layer_all.to_dict()      # gene -> 0..9

# target layer distribution
tl = layer_all.reindex(genes47)
tl_counts = tl.value_counts().sort_index()

# background pools per layer
bg_by_layer = {i: np.array([g for g in bg_genes if layer_val.get(g) == i]) for i in range(10)}

def perm_emp(obs, target_genes, layer_counts, n_perm):
    """Strength-matched permutation: empirical P(perm >= obs)."""
    cnt = np.zeros(n_perm)
    for i in range(n_perm):
        chosen = []
        for li, n in layer_counts.items():
            pool = bg_by_layer.get(li)
            if pool is None or len(pool) == 0:
                # no background in layer -> draw from nearest non-empty layer
                for li2 in range(10):
                    if len(bg_by_layer[li2]) > 0:
                        pool = bg_by_layer[li2]; break
            chosen.extend(rng.choice(pool, size=int(n), replace=False))
        cnt[i] = sum(sig_val.get(g, 0) for g in chosen)
    return float((cnt >= obs).mean()), float(cnt.mean())

# ---------------------------------------------------------------- LOO ----
print("\n== Leave-one-out (N=2000) ==")
loo_rows = []
# full-set baseline (N=10000)
obs_full = int(sum(sig_val.get(g, 0) for g in genes47))
emp_full, mean_full = perm_emp(obs_full, genes47, tl_counts, 10000)
loo_rows.append({"dropped_gene": "none", "n_targets": len(genes47),
                 "obs_sig": obs_full, "emp_p": emp_full, "perm_mean": mean_full})
print(f"  full set: obs={obs_full} emp_p={emp_full:.4f} perm_mean={mean_full:.2f}")

for g in genes47:
    tg = [x for x in genes47 if x != g]
    tlc = pd.Series({li: sum(1 for x in tg if layer_val.get(x) == li)
                     for li in range(10)}, dtype=int)
    tlc = tlc[tlc > 0]
    obs = int(sum(sig_val.get(x, 0) for x in tg))
    emp, mean_perm = perm_emp(obs, tg, tlc, 2000)
    loo_rows.append({"dropped_gene": g, "n_targets": len(tg),
                     "obs_sig": obs, "emp_p": emp, "perm_mean": mean_perm})
    print(f"  drop {g:>12s}: obs={obs} emp_p={emp:.4f} perm_mean={mean_perm:.2f}")

pd.DataFrame(loo_rows).to_csv(OUT_LOO, index=False, encoding="utf-8-sig")
print(f"  -> {OUT_LOO}")

# ------------------------------------------------- tissue-stratified ----
print("\n== Tissue-stratified (N=1000) ==")
tis_rows = []
for t in mr["tissue"].unique():
    sub = mr[mr["tissue"] == t]
    tg = [g for g in genes47 if g in set(sub["gene"])]
    bg_t = np.array([g for g in sub["gene"].unique() if g not in set(genes47)])
    if len(tg) < 3 or len(bg_t) < 50:
        print(f"  {t}: insufficient ({len(tg)} target / {len(bg_t)} bg) - skip")
        continue
    # within-tissue strength layering
    st_t = sub.groupby("gene")["p_eQTL"].min()
    qs_t = st_t.quantile(np.linspace(0, 1, 11)).values
    def layer_t(p):
        for i in range(10):
            if p <= qs_t[i + 1] + 1e-300:
                return i
        return 9
    lv_t = {g: layer_t(p) for g, p in st_t.items()}
    # within-tissue significance: gene has p_MR < 0.05 in this tissue
    sig_t = (sub.groupby("gene")["p_MR"].min() < 0.05).to_dict()
    obs = int(sum(sig_t.get(g, 0) for g in tg))
    bg_rate = float((sub.groupby("gene")["p_MR"].min() < 0.05).mean())
    expect = bg_rate * len(tg)
    fold = obs / expect if expect > 0 else np.nan
    # strength-matched permutation within tissue
    tlc = pd.Series({li: sum(1 for g in tg if lv_t.get(g) == li) for li in range(10)}, dtype=int)
    tlc = tlc[tlc > 0]
    bg_by_layer_t = {i: np.array([g for g in bg_t if lv_t.get(g) == i]) for i in range(10)}
    cnt = np.zeros(1000)
    for i in range(1000):
        chosen = []
        for li, n in tlc.items():
            pool = bg_by_layer_t.get(li)
            if pool is None or len(pool) == 0:
                for li2 in range(10):
                    if len(bg_by_layer_t[li2]) > 0:
                        pool = bg_by_layer_t[li2]; break
            chosen.extend(rng.choice(pool, size=int(n), replace=False))
        cnt[i] = sum(sig_t.get(g, 0) for g in chosen)
    emp = float((cnt >= obs).mean())
    tis_rows.append({"tissue": t, "n_targets": len(tg), "obs_sig": obs,
                     "bg_sig_rate": round(bg_rate, 4), "expected": round(expect, 2),
                     "fold_enrichment": round(fold, 2), "emp_p": emp,
                     "perm_mean": round(float(cnt.mean()), 2)})
    print(f"  {t}: n={len(tg)} obs={obs} expect={expect:.2f} fold={fold:.2f} emp_p={emp:.4f}")

pd.DataFrame(tis_rows).to_csv(OUT_TIS, index=False, encoding="utf-8-sig")
print(f"  -> {OUT_TIS}")

summary = {"loo_full_obs": obs_full, "loo_full_emp_p": emp_full,
           "loo_always_sig": int(all(r["emp_p"] < 0.05 for r in loo_rows[1:])),
           "tissue_folds": {r["tissue"]: r["fold_enrichment"] for r in tis_rows},
           "tissue_emp": {r["tissue"]: r["emp_p"] for r in tis_rows}}
json.dump(summary, open(P.RESULTS_DIR / "loo_tissue_summary.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=2)
print("\nsummary:", json.dumps(summary, ensure_ascii=False, indent=2))
