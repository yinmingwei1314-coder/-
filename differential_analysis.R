##############################################################################
# 差异表达分析脚本 - 多GEO数据集
# 筛选条件: |logFC| > 1 (即 logFC > 1 或 logFC < -1), adj.P.Val < 0.05
##############################################################################

# ============== 加载所需包 ==============
required_packages <- c("GEOquery", "limma", "DESeq2", "dplyr", "ggplot2",
                        "ggrepel", "pheatmap", "EnhancedVolcano")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("GEOquery", "limma", "DESeq2", "EnhancedVolcano")) {
      BiocManager::install(pkg, ask = FALSE)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# ============== 全局参数 ==============
LOGFC_CUTOFF <- 1        # |logFC| > 1
PVALUE_CUTOFF <- 0.05    # adj.P.Val < 0.05

dir.create("results", showWarnings = FALSE)

##############################################################################
# 1. GSE162729 - RNA-seq数据 (DESeq2分析)
##############################################################################
analyze_GSE162729 <- function() {
  cat("\n========== 分析 GSE162729 (RNA-seq) ==========\n")

  # 读取原始计数矩阵
  counts <- read.delim("GSE162729/GSE162729_raw_counts_GRCh38.p13_NCBI.tsv",
                        row.names = 1, check.names = FALSE)

  # 读取series matrix获取样本分组信息
  gse <- getGEO(filename = "GSE162729/GSE162729_series_matrix.txt",
                 getGPL = FALSE)
  pdata <- pData(gse)

  # 提取分组信息 (根据实际数据调整列名)
  # 通常在 characteristics_ch1 或 source_name_ch1 中
  group_col <- grep("disease|condition|group|treatment|status",
                     colnames(pdata), value = TRUE, ignore.case = TRUE)

  if (length(group_col) == 0) {
    group_col <- grep("characteristics_ch1", colnames(pdata),
                       value = TRUE)[1]
  }

  cat("使用分组列:", group_col[1], "\n")
  groups <- factor(pdata[[group_col[1]]])
  cat("分组水平:", levels(groups), "\n")

  # 确保样本顺序一致
  common_samples <- intersect(colnames(counts), rownames(pdata))
  counts <- counts[, common_samples]
  groups <- groups[match(common_samples, rownames(pdata))]

  # 过滤低表达基因
  keep <- rowSums(counts >= 10) >= min(table(groups))
  counts <- counts[keep, ]
  cat("过滤后基因数:", nrow(counts), "\n")

  # DESeq2分析
  coldata <- data.frame(condition = groups, row.names = common_samples)
  dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                 colData = coldata,
                                 design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, alpha = PVALUE_CUTOFF)
  res_df <- as.data.frame(res) %>%
    mutate(gene = rownames(res)) %>%
    arrange(padj)

  # 筛选差异基因: |logFC| > 1 且 adj.P < 0.05
  deg <- res_df %>%
    filter(abs(log2FoldChange) > LOGFC_CUTOFF & padj < PVALUE_CUTOFF)

  up_genes <- deg %>% filter(log2FoldChange > LOGFC_CUTOFF)
  down_genes <- deg %>% filter(log2FoldChange < -LOGFC_CUTOFF)

  cat("上调基因数:", nrow(up_genes), "\n")
  cat("下调基因数:", nrow(down_genes), "\n")
  cat("总差异基因:", nrow(deg), "\n")

  # 保存结果
  write.csv(res_df, "results/GSE162729_all_results.csv", row.names = FALSE)
  write.csv(deg, "results/GSE162729_DEGs_logFC1.csv", row.names = FALSE)

  # 火山图
  plot_volcano(res_df, "log2FoldChange", "padj", "GSE162729")

  # 热图 (top 50 DEGs)
  if (nrow(deg) > 0) {
    top_genes <- head(deg$gene, min(50, nrow(deg)))
    norm_counts <- counts(dds, normalized = TRUE)
    plot_heatmap(norm_counts[top_genes, ], coldata, "GSE162729")
  }

  return(list(all = res_df, deg = deg))
}

##############################################################################
# 2. GSE28623 - 芯片数据 (limma分析)
##############################################################################
analyze_GSE28623 <- function() {
  cat("\n========== 分析 GSE28623 (Microarray) ==========\n")
  analyze_microarray("GSE28623/GSE28623_series_matrix.txt",
                      "GSE28623/GPL4133-12599.txt",
                      "GSE28623")
}

##############################################################################
# 3. GSE34608 - 芯片数据 (limma分析)
##############################################################################
analyze_GSE34608 <- function() {
  cat("\n========== 分析 GSE34608 (Microarray) ==========\n")
  analyze_microarray("GSE34608/GSE34608-GPL6480_series_matrix.txt",
                      "GSE34608/GPL6480-9577.txt",
                      "GSE34608")
}

##############################################################################
# 4. GSE83456 - 芯片数据 (limma分析)
##############################################################################
analyze_GSE83456 <- function() {
  cat("\n========== 分析 GSE83456 (Microarray) ==========\n")
  analyze_microarray("GSE83456/GSE83456_series_matrix.txt",
                      "GSE83456/GPL10558-50081.txt",
                      "GSE83456")
}

##############################################################################
# 通用芯片数据分析函数 (limma)
##############################################################################
analyze_microarray <- function(matrix_file, platform_file, dataset_name) {

  # 读取GEO数据
  gse <- getGEO(filename = matrix_file, getGPL = FALSE)
  expr <- exprs(gse)
  pdata <- pData(gse)

  cat("样本数:", ncol(expr), "\n")
  cat("探针数:", nrow(expr), "\n")

  # 提取分组信息
  group_col <- grep("disease|condition|group|treatment|status|source",
                     colnames(pdata), value = TRUE, ignore.case = TRUE)

  if (length(group_col) == 0) {
    group_col <- grep("characteristics_ch1", colnames(pdata),
                       value = TRUE)[1]
  }

  cat("使用分组列:", group_col[1], "\n")
  groups <- factor(pdata[[group_col[1]]])
  cat("分组水平:", levels(groups), "\n")

  # 检查数据是否需要log2转换
  if (max(expr, na.rm = TRUE) > 50) {
    cat("对表达数据进行log2转换\n")
    expr[expr <= 0] <- NA
    expr <- log2(expr)
  }

  # 去除NA行
  expr <- expr[complete.cases(expr), ]

  # limma差异分析
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)

  fit <- lmFit(expr, design)

  # 构建对比 (第二组 vs 第一组)
  contrast_name <- paste0(levels(groups)[2], "-", levels(groups)[1])
  contrast_matrix <- makeContrasts(contrasts = contrast_name, levels = design)

  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  # 获取所有结果
  results <- topTable(fit2, number = Inf, sort.by = "P")
  results$probe_id <- rownames(results)

  # 读取平台注释文件进行探针到基因的映射
  platform <- tryCatch({
    read.delim(platform_file, comment.char = "#", check.names = FALSE)
  }, error = function(e) {
    cat("警告: 无法读取平台文件, 使用探针ID\n")
    NULL
  })

  if (!is.null(platform)) {
    # 查找基因符号列
    gene_col <- grep("gene.?symbol|symbol|gene_symbol|GENE_SYMBOL",
                      colnames(platform), value = TRUE, ignore.case = TRUE)[1]
    id_col <- grep("^ID$|^ID_REF$|probe", colnames(platform),
                    value = TRUE, ignore.case = TRUE)[1]

    if (!is.na(gene_col) && !is.na(id_col)) {
      gene_map <- platform[, c(id_col, gene_col)]
      colnames(gene_map) <- c("probe_id", "gene_symbol")
      gene_map <- gene_map[gene_map$gene_symbol != "" &
                             !is.na(gene_map$gene_symbol), ]
      results <- merge(results, gene_map, by = "probe_id", all.x = TRUE)
    }
  }

  # 筛选差异基因: |logFC| > 1 且 adj.P.Val < 0.05
  deg <- results %>%
    filter(abs(logFC) > LOGFC_CUTOFF & adj.P.Val < PVALUE_CUTOFF)

  up_genes <- deg %>% filter(logFC > LOGFC_CUTOFF)
  down_genes <- deg %>% filter(logFC < -LOGFC_CUTOFF)

  cat("上调基因数:", nrow(up_genes), "\n")
  cat("下调基因数:", nrow(down_genes), "\n")
  cat("总差异基因:", nrow(deg), "\n")

  # 保存结果
  write.csv(results,
            paste0("results/", dataset_name, "_all_results.csv"),
            row.names = FALSE)
  write.csv(deg,
            paste0("results/", dataset_name, "_DEGs_logFC1.csv"),
            row.names = FALSE)

  # 火山图
  plot_volcano(results, "logFC", "adj.P.Val", dataset_name)

  # 热图 (top 50 DEGs)
  if (nrow(deg) > 0) {
    top_probes <- head(deg$probe_id, min(50, nrow(deg)))
    coldata <- data.frame(condition = groups, row.names = colnames(expr))
    plot_heatmap(expr[top_probes, ], coldata, dataset_name)
  }

  return(list(all = results, deg = deg))
}

##############################################################################
# 可视化函数
##############################################################################

# 火山图
plot_volcano <- function(df, fc_col, pval_col, dataset_name) {

  df$significance <- "Not Significant"
  df$significance[df[[fc_col]] > LOGFC_CUTOFF &
                    df[[pval_col]] < PVALUE_CUTOFF] <- "Up-regulated"
  df$significance[df[[fc_col]] < -LOGFC_CUTOFF &
                    df[[pval_col]] < PVALUE_CUTOFF] <- "Down-regulated"
  df$significance <- factor(df$significance,
                             levels = c("Down-regulated",
                                        "Not Significant",
                                        "Up-regulated"))

  p <- ggplot(df, aes(x = .data[[fc_col]],
                       y = -log10(.data[[pval_col]]),
                       color = significance)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Down-regulated" = "#2166AC",
                                   "Not Significant" = "grey60",
                                   "Up-regulated" = "#B2182B")) +
    geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF),
               linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(PVALUE_CUTOFF),
               linetype = "dashed", color = "grey40") +
    labs(title = paste0(dataset_name, " - Volcano Plot"),
         subtitle = paste0("|logFC| > ", LOGFC_CUTOFF,
                           ", adj.P.Val < ", PVALUE_CUTOFF),
         x = "log2 Fold Change",
         y = "-log10(adjusted P-value)") +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(paste0("results/", dataset_name, "_volcano.pdf"),
         p, width = 8, height = 6)
  ggsave(paste0("results/", dataset_name, "_volcano.png"),
         p, width = 8, height = 6, dpi = 300)
  cat("火山图已保存:", dataset_name, "\n")
}

# 热图
plot_heatmap <- function(expr_matrix, coldata, dataset_name) {

  # 标准化用于热图展示
  scaled_data <- t(scale(t(expr_matrix)))
  scaled_data[scaled_data > 2] <- 2
  scaled_data[scaled_data < -2] <- -2

  annotation_col <- data.frame(Group = coldata$condition)
  rownames(annotation_col) <- rownames(coldata)

  pdf(paste0("results/", dataset_name, "_heatmap.pdf"),
      width = 10, height = 12)
  pheatmap(scaled_data,
           annotation_col = annotation_col,
           show_rownames = (nrow(scaled_data) <= 50),
           show_colnames = FALSE,
           clustering_method = "ward.D2",
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           main = paste0(dataset_name, " - Top DEGs Heatmap"))
  dev.off()
  cat("热图已保存:", dataset_name, "\n")
}

##############################################################################
# 运行所有分析
##############################################################################
cat("=====================================================\n")
cat("差异表达分析 - 筛选条件: |logFC| > 1, adj.P < 0.05\n")
cat("=====================================================\n\n")

results_list <- list()

# 逐一分析各数据集
tryCatch({
  results_list[["GSE162729"]] <- analyze_GSE162729()
}, error = function(e) cat("GSE162729 分析出错:", e$message, "\n"))

tryCatch({
  results_list[["GSE28623"]] <- analyze_GSE28623()
}, error = function(e) cat("GSE28623 分析出错:", e$message, "\n"))

tryCatch({
  results_list[["GSE34608"]] <- analyze_GSE34608()
}, error = function(e) cat("GSE34608 分析出错:", e$message, "\n"))

tryCatch({
  results_list[["GSE83456"]] <- analyze_GSE83456()
}, error = function(e) cat("GSE83456 分析出错:", e$message, "\n"))

##############################################################################
# 汇总统计
##############################################################################
cat("\n\n=====================================================\n")
cat("汇总统计\n")
cat("=====================================================\n")

summary_df <- data.frame(
  Dataset = character(),
  Total_Genes = integer(),
  DEGs = integer(),
  Up_regulated = integer(),
  Down_regulated = integer(),
  stringsAsFactors = FALSE
)

for (name in names(results_list)) {
  res <- results_list[[name]]
  if (!is.null(res)) {
    fc_col <- ifelse("log2FoldChange" %in% colnames(res$all),
                      "log2FoldChange", "logFC")
    pval_col <- ifelse("padj" %in% colnames(res$all), "padj", "adj.P.Val")

    summary_df <- rbind(summary_df, data.frame(
      Dataset = name,
      Total_Genes = nrow(res$all),
      DEGs = nrow(res$deg),
      Up_regulated = sum(res$deg[[fc_col]] > 0),
      Down_regulated = sum(res$deg[[fc_col]] < 0)
    ))
  }
}

print(summary_df)
write.csv(summary_df, "results/summary_statistics.csv", row.names = FALSE)

cat("\n分析完成! 结果保存在 results/ 目录下\n")
