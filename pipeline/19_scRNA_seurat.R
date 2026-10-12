# ============================================================
# 步骤 19 · 单细胞 GSE153066 正式 Seurat 注释（85,205 细胞）
# ------------------------------------------------------------
# 功能: QC→Normalize→2000 HVG→PCA(30)→UMAP→FindClusters(24)→
#       marker 注释 6 类；核心基因各亚群 CTL vs IDD（样本级 n=16）
# 输入: GSE153066_AllSample.counts.tsv.gz
# 输出: geo_work/GSE153066_seurat.rds +
#       results/scRNA_GSE153066_seurat_subtype.csv
#       （核心基因各亚群无显著差异——粗分群信号经校正不成立）
# ============================================================
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
options(future.globals.maxSize = 8 * 1024^3)
suppressMessages({
  library(data.table); library(Matrix); library(Seurat)
})

cat("=== GSE153066 正式 Seurat 流程 ===\n")
t0 <- Sys.time()

# 复现种子：固定 Seurat 随机步骤（RunUMAP / FindClusters）的随机数流，
# 使整套单细胞流程在换电脑后可逐位复现（2026-10-10 遗留事项 #4 新增）
set.seed(42)

# 1. 读 counts（data.table fread，gzip 自动识别）
cat("读取计数矩阵...\n")
mat <- fread(file.path(geo_work, "GSE153066_AllSample.counts.tsv.gz"), header = TRUE, sep = "\t")
genes <- mat[[1]]
mat[, 1 := NULL]
cells <- colnames(mat)
mat <- as.matrix(mat)
rownames(mat) <- genes
cat(sprintf("矩阵: %d 基因 x %d 细胞\n", nrow(mat), ncol(mat)))

# 2. 创建 Seurat 对象
so <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200)
cat(sprintf("Seurat 对象: %d 细胞\n", ncol(so)))
rm(mat); gc()

# 3. QC：线粒体比例
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
cat("QC 过滤前:", ncol(so), "细胞\n")
print(summary(so$nFeature_RNA))
# 保守过滤（85k 细胞，髓核 RNA 低含量，放宽阈值）
so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
cat("QC 过滤后:", ncol(so), "细胞\n")

# 4. 归一化 + 高变基因
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)

# 5. Scale + PCA
so <- ScaleData(so, vars.to.regress = c("percent.mt"))
so <- RunPCA(so, npcs = 30, verbose = FALSE)

# 6. UMAP + 聚类
so <- RunUMAP(so, dims = 1:20, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:20)
so <- FindClusters(so, resolution = 0.8, verbose = FALSE)
cat("聚类数:", length(unique(so$seurat_clusters)), "\n")

# 7. 样本信息
so$sample <- sapply(strsplit(colnames(so), "_"), `[`, 1)
so$group <- ifelse(grepl("^CTL", so$sample), "CTL", "IDD")
cat("样本分布:", table(so$sample), "\n")

# 8. marker 注释（已知 NP 亚群 marker，取各 cluster 平均表达）
marker_genes <- c("COL2A1","ACAN","SOX9","COL1A1","COL1A2","COL10A1","IBSP","SPP1",
                  "PTPRC","CD3D","CD68","PECAM1","VWF","PDGFRA","THY1","PRG4","KRT19")
avail <- intersect(marker_genes, rownames(so))
avail
avail <- avail[avail %in% rownames(so)]
cat("可用 marker:", paste(avail, collapse=", "), "\n")
avg <- AverageExpression(so, features = avail, group.by = "seurat_clusters")$RNA
avg <- as.data.frame(avg)
# P1-15 防御：核对列名是否为 "g<cluster>" 形式（该版本 Seurat 返回 g 前缀，
# 已实测有效；若版本变化导致列名不匹配，此处显式报错而非静默 NA）
if (!all(grepl("^g[0-9]+$", colnames(avg)))) {
  stop("AverageExpression 列名非 g 前缀形式: ", paste(head(colnames(avg), 10), collapse=", "))
}
cat("AverageExpression 列名核对通过（g 前缀）:", length(colnames(avg)), "clusters\n")
# 手动注释规则：基于 marker 平均表达（AUCell 太慢，用 max 归一化打分）
score <- t(apply(avg, 2, function(x) (x - min(x)) / (max(x) - min(x) + 1e-9)))
score <- as.data.frame(score)
# 简化注释逻辑
annot <- sapply(rownames(score), function(cl) {
  s <- score[cl, , drop = FALSE]
  cols <- colnames(s)
  vals <- as.numeric(s[1, ])
  names(vals) <- cols
  np <- mean(vals[c("COL2A1","ACAN","SOX9")], na.rm=TRUE)
  af <- mean(vals[c("COL1A1","COL1A2")], na.rm=TRUE)
  hyp <- mean(vals[c("COL10A1","IBSP","SPP1")], na.rm=TRUE)
  imm <- mean(vals[c("PTPRC","CD3D","CD68")], na.rm=TRUE)
  endo <- mean(vals[c("PECAM1","VWF")], na.rm=TRUE)
  mes <- mean(vals[c("PDGFRA","THY1")], na.rm=TRUE)
  lab <- c("NP_like"=np, "AF_like"=af, "Hypertrophic"=hyp, "Immune"=imm,
           "Endothelial"=endo, "Mesenchymal"=mes)
  if (max(lab) < 0.3) return("Unclassified")
  names(lab)[which.max(lab)]
})
so$celltype <- unname(annot[paste0("g", as.character(so$seurat_clusters))])
cat("\n细胞类型分布:\n"); print(table(so$celltype))

# 9. 核心基因在各亚群 CTL vs IDD（伪 bulk：每样本每亚群平均表达）
core <- c("GFPT1","GFPT2","SLC25A13","SLC25A12","NFS1","GLS","GLUD1","GLUL","SLC38A1","SLC38A3",
          "SLC7A5","SLC7A6","SLC7A7","SLC7A11","GSS","GMPS","SLC1A5","GOT2","ASNS","ARG1","GAD1","GCLM","SUCLG1")
core <- intersect(core, rownames(so))
Idents(so) <- "celltype"
res_list <- list()
for (ct in unique(so$celltype)) {
  sub <- subset(so, idents = ct)
  if (ncol(sub) < 30) next
  expr <- as.matrix(GetAssayData(sub, assay = "RNA", layer = "data"))[core, , drop = FALSE]
  samples <- unique(sub$sample)
  ctrl <- samples[grepl("^CTL", samples)]; idd <- samples[grepl("^IDD", samples)]
  for (g in core) {
    ps <- sapply(samples, function(s) mean(expr[g, sub$sample == s]))
    if (length(ctrl) >= 3 & length(idd) >= 3) {
      t <- t.test(ps[idd], ps[ctrl])
      pv <- t$p.value
    } else pv <- NA
    res_list[[length(res_list)+1]] <- data.frame(celltype=ct, gene=g,
      CTL=mean(ps[ctrl], na.rm=TRUE), IDD=mean(ps[idd], na.rm=TRUE), p=pv)
  }
}
res <- do.call(rbind, res_list)
res$FDR <- p.adjust(res$p, method = "BH")
res <- res[order(res$p), ]
cat("\n=== 核心基因各亚群 CTL vs IDD（Seurat 正式注释）===\n")
print(head(res, 30), row.names = FALSE)
write.csv(res, file.path(results_dir, "scRNA_GSE153066_seurat_subtype.csv"), row.names = FALSE)
cat("\n已存 results/scRNA_GSE153066_seurat_subtype.csv\n")

# 10. 保存 Seurat 对象（供后续使用）
saveRDS(so, file.path(geo_work, "GSE153066_seurat.rds"))
cat("已存 GSE153066_seurat.rds（geo_work）\n")
cat(sprintf("\n总耗时: %.1f 分钟\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))
