# Mendelian randomization across the glutamine pathway implicates SLC25A13 in intervertebral disc degeneration (IDD)

Genome-wide FDR-controlled cis-eQTL Mendelian randomization (MR) screen across pre-specified glutamine-pathway genes for intervertebral disc degeneration.

## What this repository contains

Reproducible analysis pipeline (`pipeline/`, numbered scripts `01` → `24`, run order = script number) covering:

- Cis-eQTL instrument preparation (GTEx v8, four tissues) and MR estimation (Wald ratio, allele-alignment checks)
- Genome-wide screen (76,953 gene–tissue pairs) with within-layer Benjamini–Hochberg global FDR control
- Bayesian colocalization (coloc.abf, coloc.susie), SMR/HEIDI, conditional analysis
- Independent consistency: eQTLGen, UK Biobank (direction-consistent, same-LD proxy), triple-cohort replication
- Pathway enrichment with permutation (10,000 iterations, instrument-strength-matched), leave-one-out and tissue-stratified stability
- Expression validation (GSE70362, GSE124272, GSE245147) and single-cell atlas (GSE153066)
- Figure scripts for the main and supplementary figures

Key result: **SLC25A13–whole blood** carries the strongest convergent genetic evidence (robust colocalization under both models, passing HEIDI filter, direction-consistent UK Biobank estimates); **GFPT1** is a literature-supported secondary signal.

## Getting started

Full setup instructions (dependencies, environment variables, external data download list, run order, resume-from-step) are in **[`pipeline/README_pipeline.md`](pipeline/README_pipeline.md)**.

Quick start:

```bash
# 1. set machine-specific paths via environment variables (see pipeline/README_pipeline.md)
# 2. download external data listed in pipeline/README_pipeline.md (GTEx v8 BESD, FinnGen R12, 1000G EUR, GEO datasets)
pip install pandas numpy scipy requests
cd pipeline
python run_pipeline.py            # run 01 → 24
python run_pipeline.py 12 24      # resume from step 12
```

All paths are resolved through `pipeline/config/paths.py`; no hard-coded machine paths. Results are written to `results/` (redirectable via `REP_MR_IDD_RESULTS`).

## Data availability

All input data are publicly available (GTEx v8, FinnGen R12, 1000 Genomes Phase 3 EUR, eQTLGen, UK Biobank OpenGWAS, GEO datasets). See `pipeline/README_pipeline.md` for the download checklist. Analysis results can be fully reproduced with this pipeline.

## License

MIT — see [LICENSE](LICENSE).
