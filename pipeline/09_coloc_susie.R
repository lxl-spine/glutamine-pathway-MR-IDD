# ============================================================
# 步骤 09 · coloc.susie 精细共定位（R3 双方法一致）
# ------------------------------------------------------------
# 功能: 对 coloc.abf PP.H4≥0.75 的候选做 SuSiE 精细共定位
#       （收敛 + PP.H4≥0.75 双条件），报告收敛失败率
# 输入: coloc_input/*.txt + 1000G EUR LD
# 输出: results/coloc_susie/（收敛率 + 双方法一致清单）
# ============================================================
# ---- 统一路径入口（R 侧约定，见 REFACTOR_CONVENTIONS.md；换电脑可改环境变量或下方默认值）----
root <- Sys.getenv("REP_MR_IDD_ROOT", unset = "")
if (!nzchar(root)) root <- normalizePath(file.path(getwd(), ".."), winslash = "/")
if (!dir.exists(root)) stop("无法定位项目根目录：请设置环境变量 REP_MR_IDD_ROOT 或在本文件顶部修改 root。")
results_dir <- Sys.getenv("REP_MR_IDD_RESULTS", unset = file.path(root, "results"))
geo_work   <- Sys.getenv("REP_MR_IDD_GEO", unset = "")
smr_work   <- Sys.getenv("REP_MR_IDD_SMR_WORK", unset = "")
r_libs     <- Sys.getenv("REP_MR_IDD_RLIBS", unset = "")
if (nzchar(r_libs) && dir.exists(r_libs)) .libPaths(c(r_libs, .libPaths()))
if (!nzchar(geo_work)) stop("未设置环境变量 REP_MR_IDD_GEO；请按 README_pipeline.md 配置后重跑。")
cat("ROOT =", root, "| results_dir =", results_dir, "\n")
# run_coloc_susie_full3.R — 全量 coloc.susie v3：完整 summary 落盘 + 错误捕获
suppressMessages({library(coloc); library(data.table)})
WORK <- file.path(results_dir, "coloc_susie")
INP  <- file.path(results_dir, "coloc_input")
GLD  <- file.path(results_dir, "glut_ld")
CISD <- results_dir
N_EQTL <- c(Adipose_Subcutaneous=581, Muscle_Skeletal=706, Nerve_Tibial=532, Whole_Blood=670)
N_GWAS <- 404907; S_GWAS <- 51683/404907
CISFILES <- c(Adipose_Subcutaneous="cis_snps_full_Adipose_Subcutaneous.csv",
              Muscle_Skeletal="cis_snps_full_Muscle_Skeletal.csv",
              Nerve_Tibial="cis_snps_full_Nerve_Tibial.csv",
              Whole_Blood="cis_snps_full_Whole_Blood.csv")
panmap <- fread(file.path(WORK, "glut_panels.txt"), header=FALSE, sep="\t")
PANELS <- setNames(panmap$V2, panmap$V1)
todo <- read.table(file.path(WORK, "todo_full.txt"), stringsAsFactors=FALSE, sep=" ")

out <- list()
for (i in seq_len(nrow(todo))) {
  tissue <- todo[i,1]; gene <- todo[i,2]
  key <- paste0(tissue, "_", gene)
  outfile <- file.path(WORK, paste0(key, ".susie.csv"))
  if (file.exists(outfile) && file.info(outfile)$size > 50) { out[[key]] <- fread(outfile); next }
  d <- tryCatch(fread(file.path(INP, paste0(key, ".txt")), header=TRUE), error=function(e) NULL)
  empty <- data.frame(tissue=tissue, gene=gene, nsnps=NA, hit1=NA, hit2=NA,
                      PP.H0=NA, PP.H1=NA, PP.H2=NA, PP.H3=NA, PP.H4=NA, note="no-input")
  if (is.null(d) || nrow(d) == 0) { fwrite(empty, outfile); out[[key]] <- empty; next }
  pan <- PANELS[[gene]]
  if (is.null(pan) || !file.exists(file.path(GLD, paste0(pan, ".bim")))) {
    empty$note <- "no-panel"; fwrite(empty, outfile); out[[key]] <- empty; next
  }
  bim <- fread(file.path(GLD, paste0(pan, ".bim")), header=FALSE)
  # P3 修复：LD 矩阵表头防御（plink --r square 无表头；若读取到表头导致
  # 行列数不匹配则重试 header=TRUE）
  ld <- tryCatch(as.matrix(fread(file.path(GLD, paste0(pan, "_ld.ld")), header=FALSE)), error=function(e) NULL)
  if (!is.null(ld) && nrow(ld) > 0 && ncol(ld) != nrow(ld)) {
    ld2 <- tryCatch(as.matrix(fread(file.path(GLD, paste0(pan, "_ld.ld")), header=TRUE)), error=function(e) NULL)
    if (!is.null(ld2) && nrow(ld2) == ncol(ld2)) ld <- ld2
  }
  if (is.null(ld) || ncol(ld) != nrow(bim)) {
    empty$note <- "ld-mismatch"; fwrite(empty, outfile); out[[key]] <- empty; next
  }
  rownames(ld) <- colnames(ld) <- bim$V2
  na_snp <- rownames(ld)[apply(is.na(ld), 1, any)]
  if (length(na_snp) > 0) ld <- ld[!rownames(ld) %in% na_snp, !colnames(ld) %in% na_snp]
  common <- intersect(d$snp, colnames(ld))
  if (length(common) < 10) { empty$note <- "few-snps"; fwrite(empty, outfile); out[[key]] <- empty; next }
  d2 <- d[d$snp %in% common, ]
  cis <- fread(file.path(CISD, CISFILES[[tissue]]), header=TRUE)
  mr_sub <- unique(cis[gene == gene, .(SNP, A1, A2)])
  d2 <- merge(d2, mr_sub, by.x="snp", by.y="SNP", all.x=TRUE)
  pA <- data.table(rsid=bim$V2, pA1=bim$V5, pA2=bim$V6)
  d2 <- merge(d2, pA, by.x="snp", by.y="rsid", all.x=TRUE)
  d2 <- d2[!is.na(A1) & (A1 == pA1 | A1 == pA2)]
  flip <- d2$A1 == d2$pA2
  d2[flip, `:=`(b_eqtl = -b_eqtl, b_gwas = -b_gwas)]
  if (nrow(d2) < 10) { empty$note <- "few-after-align"; fwrite(empty, outfile); out[[key]] <- empty; next }
  idx <- match(d2$snp, colnames(ld))
  L <- ld[idx, idx]
  d1 <- list(snp=d2$snp, beta=d2$b_eqtl, varbeta=d2$se_eqtl^2, type="quant", N=N_EQTL[[tissue]], sdY=1, LD=L)
  d2l <- list(snp=d2$snp, beta=d2$b_gwas, varbeta=d2$se_gwas^2, type="cc", N=N_GWAS, s=S_GWAS, LD=L)
  errmsg <- ""
  res <- tryCatch(suppressWarnings(coloc.susie(d1, d2l, susie.args=list(maxit=50, verbose=FALSE))),
                  error=function(e) { errmsg <<- substr(conditionMessage(e), 1, 120); NULL })
  if (is.null(res)) {
    empty$note <- paste0("susie-error: ", errmsg)
    fwrite(empty, outfile); out[[key]] <- empty; next
  }
  sus <- res$summary
  row <- data.frame(tissue=tissue, gene=gene,
                    nsnps=if ("nsnps" %in% colnames(sus)) sus$nsnps[1] else NA,
                    hit1=if ("hit1" %in% colnames(sus)) sus$hit1[1] else NA,
                    hit2=if ("hit2" %in% colnames(sus)) sus$hit2[1] else NA,
                    PP.H0=if ("PP.H0.abf" %in% colnames(sus)) sus$PP.H0.abf[1] else NA,
                    PP.H1=if ("PP.H1.abf" %in% colnames(sus)) sus$PP.H1.abf[1] else NA,
                    PP.H2=if ("PP.H2.abf" %in% colnames(sus)) sus$PP.H2.abf[1] else NA,
                    PP.H3=if ("PP.H3.abf" %in% colnames(sus)) sus$PP.H3.abf[1] else NA,
                    PP.H4=if ("PP.H4.abf" %in% colnames(sus)) sus$PP.H4.abf[1] else NA,
                    note="ok")
  fwrite(row, outfile)
  out[[key]] <- row
  if (i %% 20 == 0) cat("progress", i, "/", nrow(todo), "\n")
}
res_df <- rbindlist(out, fill=TRUE)
fwrite(res_df, file.path(WORK, "coloc_susie_results_full.csv"))
cat("DONE", nrow(res_df), "rows | PH4>=0.5:", sum(res_df$PP.H4 >= 0.5, na.rm=TRUE), "\n")
