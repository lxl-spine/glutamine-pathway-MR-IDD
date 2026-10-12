# ============================================================
# 步骤 06 · coloc.abf 共定位（R3，阈值 PP.H4≥0.75 统一）
# ------------------------------------------------------------
# 功能: 每基因×组织 coloc.abf，eQTL N 用 GTEx v8 公开样本量，
#       GWAS cc N=404,907 s=0.1276（FinnGen R12）
# 输入: results/coloc_input/*.txt
# 输出: results/coloc_results.csv（1,770 对全部完成，117 对 PP.H4≥0.75）
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
cat("ROOT =", root, "| results_dir =", results_dir, "\n")

suppressMessages({library(coloc); library(data.table)})

INP <- file.path(results_dir, "coloc_input")
OUT <- file.path(results_dir, "coloc_results.csv")
N_EQTL <- c(Adipose=581, Muscle=706, Nerve=532, WholeBlood=670)
NCASE <- 51683; NCONT <- 353224

files <- list.files(INP, pattern="\\.txt$")
cat("files:", length(files), "\n")

# 已完成的断点
done <- character(0)
out <- list()
if (file.exists(OUT)) {
  prev <- fread(OUT, sep="\t")
  if (nrow(prev) > 0) {
    done <- paste0(prev$tissue, "_", prev$gene)
    # P0-5 修复：把已完成的 prev 并入 out，避免部分保存时覆盖旧结果
    out <- lapply(seq_len(nrow(prev)), function(i) prev[i, ])
  }
  cat("resume from", length(done), "done\n")
}
for (i in seq_along(files)) {
  fn <- files[i]
  tissue <- ifelse(grepl("^Adipose_", fn), "Adipose",
            ifelse(grepl("^Muscle_", fn), "Muscle",
            ifelse(grepl("^Nerve_", fn), "Nerve",
            ifelse(grepl("^Whole_", fn), "WholeBlood", NA))))
  if (is.na(tissue)) next
  gene <- sub("^[^_]+_[^_]+_", "", sub("\\.txt$", "", fn))
  key <- paste0(tissue, "_", gene)
  if (key %in% done) next
  dt <- tryCatch(fread(file.path(INP, fn)), error=function(e) NULL)
  if (is.null(dt) || nrow(dt) < 50) next
  dt <- dt[is.finite(se_eqtl) & se_eqtl>0 & is.finite(se_gwas) & se_gwas>0 & maf>0 & maf<1, ]
  if (nrow(dt) < 50) next
  D1 <- list(type="quant", snp=dt$snp, beta=dt$b_eqtl, varbeta=dt$se_eqtl^2, MAF=dt$maf, N=N_EQTL[[tissue]])
  D2 <- list(type="cc", snp=dt$snp, beta=dt$b_gwas, varbeta=dt$se_gwas^2, MAF=dt$maf, N=NCASE+NCONT, s=NCASE/(NCASE+NCONT))
  # P0-5 修复：每次循环前清空 res，防止上一轮残留被 exists("res") 误用
  rm(list="res"); res <- NULL
  fit <- tryCatch(suppressWarnings(suppressMessages(
                    invisible(capture.output(res <- coloc.abf(D1, D2))))),
                  error=function(e) NULL)
  if (!is.null(fit) && exists("res") && !is.null(res)) {
    out[[length(out)+1]] <- data.frame(tissue=tissue, gene=gene, nsnp=nrow(dt),
                                       PP.H4=res$summary[["PP.H4.abf"]],
                                       PP.H3=res$summary[["PP.H3.abf"]], stringsAsFactors=FALSE)
  }
  if (i %% 200 == 0) {
    cat("  ", i, "/", length(files), "\n")
    if (length(out) > 0) {
      tmp <- do.call(rbind, out)
      fwrite(tmp, OUT, sep="\t")
      cat("    partial saved", nrow(tmp), "\n")
    }
  }
}
if (length(out) > 0) {
  tmp <- do.call(rbind, out)
  fwrite(tmp, OUT, sep="\t")
  cat("DONE", nrow(tmp), "coloc results ->", OUT, "\n")
} else {
  cat("no new results\n")
}
