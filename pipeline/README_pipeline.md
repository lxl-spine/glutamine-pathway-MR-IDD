# rep_MR_IDD — 完整分析管线（pipeline/ 编号脚本）

> 定位：`pipeline/` 是**从原始数据到论文结果表的完整可复现管线**。
> 每个脚本 = 一个可独立运行的步骤；编号即执行顺序（01 → 24）。
> **可移植性**：全部脚本不内嵌机器专属绝对路径，统一入口为 `pipeline/config/paths.py`。
> 换电脑只需：① 把整个仓库放到任意目录；② 按下方"外部依赖与配置"设置环境变量；
> ③ 补齐外部数据（见下载清单）。所有脚本均已在开发机跑通并产出全部核心结果。

---

## 目录结构

```
glutamine-pathway-MR-IDD/        # 仓库根（GitHub 仓库名）
├── pipeline/                  # 全部分析脚本（编号即执行顺序）
│   ├── config/
│   │   ├── paths.py            # ★ 唯一路径入口（换电脑只设环境变量，不改文件）
│   │   ├── glutamine_genes.txt # 47 谷氨酰胺基因清单
│   │   └── gene_sets.csv      # 6 组通路基因集
│   ├── 01_…py … 24_…py        # 步骤脚本
│   ├── run_pipeline.py         # 管线驱动（断点续跑）
│   └── README_pipeline.md
├── data/                       # 外部数据下载位置（见下表；不随仓库提交）
└── results/                    # 全部输出（可用 REP_MR_IDD_RESULTS 重定向）
```

---

## 管线总览（数据流）

```
原始数据                                  中间产物                       最终产物
─────────                                ─────────                     ─────────
GTEx v8 BESD ×4 ──01──▶ lead eQTL ──┐
FinnGen R12 .ma ──▶ (02 功效) ──────┤
                                    ├─03──▶ MR 主表(76,953对) ──04──▶ 组织分层
                                    │                                  │
BESD cis 区域 ──05──▶ coloc 输入 ──06──▶ coloc.abf(184对,2对PP.H4≥0.75)    │
                                    ├─07──▶ eQTLGen 复制(82.8%一致)    │
1000G EUR LD ──11──▶ 面板+LD ──08/09──▶ SMR-HEIDI + coloc.susie        │
                                    ├─10──▶ 多IV敏感性(Steiger全OK)    │
                                    ├─12──▶ 双向MR(方向性)             │
47 谷氨酰胺基因 ──13──▶ 可成药性 ──┐                                   │
                                  ├─17──▶ 富集置换(特征匹配10000次) ───┤
                                  ├─16──▶ UKB 复制(方向一致) ─────────┤
                                  ├─18──▶ 条件分析(单信号修正) ────────┤
GSE70362 ──14──▶ 表达验证① ──────┤                                   │
GSE124272 ──15──▶ 表达验证② ─────┤                                   │
GSE153066 ──19/20──▶ scRNA ──────┘                                   │
                                                                    ▼
                                                    21──▶ 证据矩阵 v2（论文核心表）
                                                    22──▶ 三队列复制
                                                    23──▶ gene set 富集（Figure 3f）
                                                    24──▶ LOO + 组织分层置换
```

---

## 步骤清单（执行顺序 = 编号）

| # | 脚本 | 功能 | 关键产物 | 运行环境 |
|---|---|---|---|---|
| 01 | `01_data_prep_lead_eqtl.py` | GTEx v8 四组织 BESD → **每基因 lead cis-eQTL**（cis±2Mb 内 p 最小 SNP，`argmin(p)`） | `results/lead_eqtl_{4组织}.csv`（每基因一行） | Python |
| 01b | `01b_cis_snps_extract.py` | 同一 BESD → 目标基因 cis±2Mb **全量 SNP**（coloc/SMR/多IV 输入） | `results/cis_snps_full_{4组织}.csv` | Python |
| 02 | `02_power_per_gene.py` | 逐基因×组织功效门槛（读 01 lead + GWAS，**完整 delta method**，不依赖 03） | `results/power_per_gene_all.csv` | Python |
| 03 | `03_mr_main.py` | **MR 主分析**：单 IV Wald ratio（等位基因对检查、双尾渐近 p）+ 层内 BH-FDR | `results/mr_results_GTEx4tissues.csv` | Python |
| 04 | `04_tissue_stratify.py` | 组织分层（R5）+ lead SNP 共享矩阵（R1，显式防重复） | `results/tier_{shared,tissue_specific,blood_specific,conflict}.csv` | Python |
| 05 | `05_coloc_prep.py` | cis_snps_full（继承 01b 的 cis±2Mb 窗口）× FinnGen 对齐合并 → coloc 输入（`coloc_input/*.txt`） | `results/coloc_input/*.txt` | Python |
| 06 | `06_coloc_abf.R` | coloc.abf（PP.H4≥0.75；断点续跑合并旧结果，无 res 残留） | `results/coloc_results.csv` | R |
| 07 | `07_eqtlgen_replication.py` | eQTLGen 全血复制（双尾 p + 索引排序 FDR）+ GTEx WB 一致性 | `results/replication_eQTLGen_vs_GTExWB.csv` | Python |
| 08 | `08_smr_heidi.py` | SMR-HEIDI（17 基因 × 4 组织；**自足生成 gwas_summary/{t}_besd/panels**） | `results/smr_glut_17genes_4tissues.csv` | Python + smr.exe |
| 09 | `09_coloc_susie.R` | coloc.susie 精细共定位（LD 表头防御） | `results/coloc_susie/` | R |
| 10 | `10_multi_iv_mr.py` | 多 IV 敏感性：IVW/**标准 MR-Egger**/加权中位数/Steiger（**eQTL 对齐面板、clump 上限 r²=0.1**） | `results/multi_iv_mr_sensitivity.csv` | Python |
| 11 | `11_glut_panels.py` | 47 基因 ±2Mb 面板 + LD 矩阵 + **glut_panels.txt**（09 依赖） | `results/glut_ld/*` | Python + plink |
| 12 | `12_bidirectional_mr.py` | 双向 MR：clump → 重叠 → **真正反向 Wald ratio**（P0-4 补全） | `results/reverse_mr_IDD_to_expr.csv` + overlap | Python |
| 13 | `13_druggability.py` | DGIdb/Open Targets 可成药性注释 | `results/druggability_summary.csv` | Python |
| 14 | `14_expr_GSE70362.py` | bulk 表达验证①（NP 芯片 48 样本；解析失败显式报错） | `results/ge_validation_GSE70362_*.csv` | Python |
| 15 | `15_expr_GSE124272.py` | bulk 表达验证②（**解析 !Sample_title 分组**，不再硬编码列） | `results/ge_validation_GSE124272_*.csv` | Python |
| 16 | `16_ukb_replication.py` | UKB 复制（ukb-b-19807；JSON 缺失时**自动调 OpenGWAS API**） | `results/ukb_replication_47genes.csv` | Python + API |
| 17 | `17_enrichment_permutation.py` | 富集置换（**10000 次，工具强度 10 分位匹配**；**背景剔除 47 目标基因**）+ 多基因集；含全局 BH-FDR | `<GEO_WORK>/enrichment_permutation_results.json`；回写 `results/mr_results_GTEx4tissues.csv`（加 FDR_global 列） | Python |
| 18 | `18_conditional_analysis.py` | 条件分析（**带符号 LD r 条件化**，LD r²=0.69 单信号） | `results/conditional_analysis_*.csv` | Python + plink 面板 |
| 19 | `19_scRNA_seurat.R` | 单细胞正式 Seurat 注释（**列名 g 前缀核对断言**） | `results/scRNA_GSE153066_seurat_subtype.csv`；`<GEO_WORK>/GSE153066_seurat.rds` | R + Seurat |
| 20 | `20_scRNA_subtype.py` | 亚群定位（marker 粗分群；**FDR mask 赋值**；方法学对照） | `results/scRNA_GSE153066_subtype.csv` | Python |
| 21 | `21_evidence_matrix.py` | **证据矩阵 v2**（**自足构建基础矩阵**，六层证据） | `results/evidence_matrix_v2.csv` | Python |
| 22 | `22_triple_cohort_replication.py` | 三队列复制（FinnGen IDD 主 + UKB IDD ukb-b-19807 + UKB 背痛 ukb-b-9838）方向一致性 | `results/triple_cohort_replication.csv` | Python |
| 23 | `23_gene_set_enrichment.py` | 6 组通路基因集富集（hypergeom + BH-FDR；与 17 同富集定义）——Figure 3f 正式数据源 | `results/gene_set_enrichment.csv` | Python |
| 24 | `24_loo_tissue_permutation.py` | Leave-one-out 富集（N=2000/基因，全集基线 N=10000）+ 组织分层置换（N=1000，种子 42） | `results/leave_one_out_enrichment.csv`、`results/tissue_enrichment_MR_permutation.csv`、`results/loo_tissue_summary.json` | Python |

---

## 环境安装

- **Python 3.13+**：`pandas`、`numpy`、`scipy`、`requests`。
  ```
  pip install pandas numpy scipy requests
  ```
- **R 4.6.x** + 包：`Seurat`（5.5.x）、`coloc`、`susieR`、`data.table`、`Matrix`。
  R 包库目录通过 `REP_MR_IDD_RLIBS` 指定（可选；不设则用 R 默认库位置）。
- **PLINK v1.9**：放置任意路径，把路径写入 `REP_MR_IDD_PLINK`（或改 paths.py 默认值）。
- **SMR v1.3.1 (Windows)**：解压到 `data/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe`
  （即 `data/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe`；可用 `REP_MR_IDD_SMR_EXE` 覆盖）。

## 外部数据下载清单（01 之前需就位）

| 数据 | 放置位置 | 用途 | 来源 |
|---|---|---|---|
| GTEx v8 BESD ×4 组织（Adipose_Subcutaneous / Muscle_Skeletal / Nerve_Tibial / Whole_Blood） | `data/GTEx_v8_westlake/{组织}/{组织}.besd/.epi/.esi` | 01/03/05/08 | GTEx Portal v8（Westlake 处理版 BESD） |
| FinnGen R12 IDD（M13 椎间盘退变）GWAS：`.ma`、`.rsid.ma`、`.smr.ma` 三个切片 | `data/finngen_R12_M13_INTERVERTEB*.ma` | 03/05/08/12/18 | FinnGen R12 公开汇总数据 |
| eQTLGen lead | `data/lead_eqtl_eQTLGen.rds/.csv` | 07 | eQTLGen 联盟公开数据 |
| 1000G EUR plink 面板（`.bed/.bim/.fam`） | `data/1000G_plink/1kg_eur*` | 08/09/10/11/18 | 1000 Genomes EUR 面板 |
| 47 基因清单 | `data/glut_genes.txt`、`pipeline/config/glutamine_genes.txt` | 01b/11/12/17/24 | 随仓库 |
| GSE153066 计数矩阵（85,205 细胞） | `<GEO_WORK>/GSE153066_AllSample.counts.tsv.gz` | 19/20 | GEO GSE153066 |
| UKB proxy JSON（ukb-b-19807 rs4727335；ukb-b-9838 core SNPs） | `<GEO_WORK>/ukb_19807_proxy_rs4727335.json`、`ukb_b9838_core_snps.json` | 22 | OpenGWAS 导出缓存 |
| SLC25A13 区域 GWAS 切片 | `<SMR_WORK>/gwas_SLC25A13.txt` | 18 | FinnGen R12 cis±1Mb 切片 |
| OpenGWAS token | 环境变量 `OPENGWAS_TOKEN` 或 `<GEO_WORK>/opengwas_token.txt` | 16 | OpenGWAS 官网申请 |

> `<GEO_WORK>` / `<SMR_WORK>` 为外部工作目录，**发布版不内置任何本机默认路径**，
> 必须通过环境变量指定（见下表），可指向任意位置。

## 外部依赖与配置（paths.py + 环境变量）

`pipeline/config/paths.py` 自动推导项目根（本文件上溯两级），项目内路径**无需配置**。
机器专属外部路径按下表配置（环境变量优先；不设则用 paths.py 内默认值）：

| 环境变量 | 含义 | 说明 |
|---|---|---|
| `REP_MR_IDD_RESULTS` | 结果输出目录（重跑验证时可重定向到临时目录） | 默认 `<仓库根>/results` |
| `REP_MR_IDD_GEO` | GEO 原始数据 / UKB 缓存 / enrichment JSON 目录 | **必设** |
| `REP_MR_IDD_SMR_WORK` | SMR 与条件分析中间工作目录 | **必设**（08/10/12/18 需要） |
| `REP_MR_IDD_RSCRIPT` | Rscript 可执行文件路径 | 必设，或确保 `Rscript` 在系统 PATH |
| `REP_MR_IDD_PLINK` | PLINK v1.9 可执行文件路径 | **必设**（08/09/10/11/12/18 需要） |
| `REP_MR_IDD_SMR_EXE` | SMR 可执行文件路径 | 默认 `data/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe` |
| `REP_MR_IDD_RLIBS` | R 包库目录（R 脚本 `.libPaths`） | 可选 |
| `OPENGWAS_TOKEN` | OpenGWAS 访问 token（16 使用） | — |

R 脚本（06/09/19）顶部有独立样板块：优先读 `REP_MR_IDD_ROOT` / `REP_MR_IDD_RESULTS` /
`REP_MR_IDD_GEO` 等环境变量；经 `run_pipeline.py` 调用时自动注入；单独手动跑 R 脚本时
请先 `$env:REP_MR_IDD_ROOT='<项目根>'`（或在 pipeline/ 目录下运行，脚本自动以 `..` 为根）。

## 执行顺序与断点续跑

- **顺序**：01 → 01b → 02 → 03 → 04 → 05 → 06 → 07 → 11 → 08 → 09 → 10 → 12 → 13 → 14 → 15 → 16 → 17 → 18 → 19 → 20 → 21 → 22 → 23 → 24。
  依赖关系：02 读 01；03/04 读 01；05 读 01b；06 读 05 的 `coloc_input/*.txt`（全量）；08 读 01b/11（面板+自足输入）；09/10/18 读 11（面板/LD）；10 读 01b+05；12 读 05+11；16-18 读 03；19-20 读原始单细胞数据；21 读全部；22 读 03+16+UKB JSON；23/24 读 03。
- **一键/断点续跑**：
  ```
  python run_pipeline.py            # 从 01 跑到 24
  python run_pipeline.py 12 24      # 只跑 12 → 24（断点续跑）
  python run_pipeline.py 23 23      # 只跑 23
  ```
  大文件步骤（01/01b/05/08/11）已有结果时跳过重算；任一步失败即停，修好用 `run_pipeline.py <失败步> 24` 续跑。
- **幂等性**：每个脚本可重复运行（输出覆盖）；17 会把加了 `FDR_global` 列的 MR 主表原样回写，属设计内行为。

## 修订记录（管线开发中的关键修正摘要）

- 01：lead eQTL 选择器为每基因最小 p SNP（`argmin(p)`，19,269 行实测）；04 显式防重复。
- 02：功效计算独立读取 01 lead + GWAS（完整 delta method），不依赖 03 输出。
- 01b：新增 cis±2Mb 全量 SNP 提取脚本（coloc/SMR/多 IV 输入）。
- 12：补全反向 Wald ratio（原仅 clump + 重叠）。
- 06：修复断点续跑覆盖旧结果与 `exists(res)` 残留。
- 10：eQTL 方向对齐面板、clump 上限固定 r²=0.1、MR-Egger 采用标准 `b_Y ~ b_X`（权重 1/se_Y²）。
- 20：修复 FDR 长度不匹配（mask 赋值）。
- 03/07：大 z 渐近 p 双尾乘 2；03/05/07 增加等位基因对 {ga1,ga2}=={ea1,ea2} 配对检查。
- 07：FDR 改用索引排序（重复 p 不再覆盖）。
- 18：区域扫描改为带符号 r 条件化。
- 17：置换背景剔除 47 目标基因。
- 15：分组改为解析 !Sample_title（不再硬编码样本序）。
- 19：增加 Seurat 列名 `g` 前缀核对断言。
- 08/16/21：补齐中间文件生成脚本；14 解析失败显式报错。

## 已知限定（诚实声明）

1. **SLC25A13 免疫亚群信号**：步骤 20（粗分群）发现 p=0.004 上调，经步骤 19（正式 Seurat 样本级检验）**不成立**——最终以 19 为准，20 仅作方法学对照，不进正文。
2. **UKB 复制**：ukb-b-19807 仅 1,045 例，功效极低——表述为"方向一致外部观察"，不宣称强复制。
3. **条件分析**：rs4727335 与 rs12216562 LD r²=0.69，为同一 LD block 单因果信号（修正早期"独立信号"误述）。
4. **反向 MR**：FinnGen IDD 无强 cis 工具（仅 92 对 p<5e-8）——放弃并如实报告。

---

## 图表脚本（投稿图，位于 `pipeline/_preview_v5/`）

> 定位：主图 + 补充图共 11 个脚本，全部**可移植**（统一走 `config/paths.py`；
> 图样式依赖 `scipilot-figure-skill`，运行时需设置
> `SCIPILOT_SKILL_SCRIPTS=<skill>/scripts`）。输出 `figures/final/{png,pdf,svg}`。
> 图注正文：`Manuscript_body_IDD_MR.docx` / `Supplementary_Materials_IDD_MR.docx`
> 中已嵌入对应渲染版本；重跑图后须重新替换嵌入图（见"替换"）。

| 脚本 | 输出（figures/final/） | 内容 |
|---|---|---|
| `rfig_fig2_discovery.py` | `Figure1_discovery` | 主图 1：47 基因 × 4 组织热图 + 170 对 MR 效应分布（3 个全局 FDR 命中高亮） |
| `rfig_fig3_crossoutcome.py` | `Figure2_crossoutcome` | 主图 2：跨结局 −log10P / GSE245147 差异表达 / 单细胞热图 |
| `rfig_fig4_convergence.py` | `Figure3_convergence` | 主图 3：SMR–HEIDI 象限 + coloc.abf/susie 双模型 |
| `rfig_fig5_locus.py` | `Figure4_locus` | 主图 4：SLC25A13 chr7 位点双轨（eQTL / GWAS + conditional） |
| `rsupp_eqtlgen.py` | `FigureS1_eQTLGen_replication` | 补充图 S1：eQTLGen 独立复制 |
| `rsupp_perm_enrich.py` | `FigureS2_permutation_enrichment` | 补充图 S2：置换 + 通路富集 |
| `rsupp_S3_expr.py` | `FigureS3_expression_validation` | 补充图 S3：bulk 表达验证 |
| `rsupp_S4_loo.py` | `FigureS4_leave_one_out` | 补充图 S4：leave-one-out 富集稳定性 |
| `rsupp_S5_revMR.py` | `FigureS5_reverse_MR` | 补充图 S5：反向 MR |
| `rsupp_S6_tissues.py` | `FigureS6_GTEx_tissues` | 补充图 S6：GTEx 四组织 eQTL 全谱 |
| `rsupp_S7_onek1k.py` | `FigureS7_onek1k_power` | 补充图 S7：OneK1K 细胞型 MR + 功效直方图 |

重跑全部图（PowerShell）：
```powershell
$env:SCIPILOT_SKILL_SCRIPTS = "C:\...\scipilot-figure-skill\scripts"
$env:REP_MR_IDD_PLINK     = "<plink 或 smr exe 路径>"   # 仅 rfig_fig5_locus 用
cd pipeline\_preview_v5
python rfig_fig2_discovery.py; python rfig_fig3_crossoutcome.py
python rfig_fig4_convergence.py; python rfig_fig5_locus.py
python rsupp_eqtlgen.py; python rsupp_perm_enrich.py; python rsupp_S3_expr.py
python rsupp_S4_loo.py; python rsupp_S5_revMR.py; python rsupp_S6_tissues.py
python rsupp_S7_onek1k.py
```
替换 docx 内嵌图：`_scratch\_audit_extract\move_figs3.py`（归档重命名）+ `replace_docx_images.py`（替换正文/SI 内嵌）。

> 历史图脚本（五图版 workflow 等）已归档到 `pipeline/_archive/`，不参与当前投稿。
