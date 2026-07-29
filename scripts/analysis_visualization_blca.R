#!/usr/bin/env Rscript
# TCGA-BLCA visualization and statistics starter script.
# Run after Python pipeline has produced files in data/processed.

`%||%` <- function(x, y) if (is.null(x)) y else x
args_file <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_file, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE) else "scripts/analysis_visualization_blca.R"
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
processed_dir <- file.path(project_root, "data", "processed")
figure_dir <- file.path(project_root, "results", "figures")
table_dir <- file.path(project_root, "results", "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install it before running this script.", pkg), call. = FALSE)
  }
}

need("data.table")
need("ggplot2")
need("survival")

clinical <- data.table::fread(file.path(processed_dir, "blca_clinical_survival.tsv"))
sample_map <- data.table::fread(file.path(processed_dir, "blca_rna_sample_map.tsv"))

# Clinical overview
p_stage <- ggplot2::ggplot(
  clinical[!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != ""],
  ggplot2::aes(x = ajcc_pathologic_stage)
) +
  ggplot2::geom_bar(fill = "#3478b8", width = 0.72) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
  ggplot2::labs(x = "AJCC pathologic stage", y = "Patients", title = "TCGA-BLCA stage distribution")
ggplot2::ggsave(file.path(figure_dir, "clinical_stage_distribution.pdf"), p_stage, width = 7, height = 4)

# Basic survival by stage, without fancy dependencies.
surv_data <- clinical[!is.na(os_time_days) & !is.na(os_event) & !is.na(ajcc_pathologic_stage)]
surv_data[, stage_group := fifelse(grepl("Stage I|Stage II", ajcc_pathologic_stage), "Stage I-II", "Stage III-IV/Other")]
fit <- survival::survfit(survival::Surv(os_time_days, os_event) ~ stage_group, data = surv_data)
pdf(file.path(figure_dir, "overall_survival_by_stage.pdf"), width = 6, height = 5)
plot(fit, col = c("#3478b8", "#c43c39"), lwd = 2, xlab = "Days", ylab = "Overall survival")
legend("bottomleft", legend = names(fit$strata), col = c("#3478b8", "#c43c39"), lwd = 2, bty = "n")
dev.off()

# Differential expression starter:
# For publication-grade DE, install DESeq2 and run the block below on count matrix.
# It compares Primary Tumor vs Solid Tissue Normal when normal samples are available.
if (requireNamespace("DESeq2", quietly = TRUE)) {
  counts <- data.table::fread(file.path(processed_dir, "blca_rna_counts_matrix.tsv.gz"))
  gene_id <- counts[[1]]
  counts <- as.data.frame(counts[, -1])
  rownames(counts) <- gene_id
  coldata <- as.data.frame(sample_map[match(colnames(counts), matrix_column)])
  rownames(coldata) <- coldata$matrix_column
  coldata$condition <- ifelse(coldata$sample_type == "Solid Tissue Normal", "Normal", "Tumor")
  keep <- coldata$condition %in% c("Tumor", "Normal")
  dds <- DESeq2::DESeqDataSetFromMatrix(round(as.matrix(counts[, keep])), coldata[keep, ], design = ~ condition)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("condition", "Tumor", "Normal"))
  res_dt <- data.table::as.data.table(res, keep.rownames = "gene_id")
  data.table::fwrite(res_dt, file.path(table_dir, "DESeq2_tumor_vs_normal.tsv"), sep = "\t")
}

# Drug-response module starter:
# Recommended packages for downstream prediction: oncoPredict or pRRophetic.
# Use blca_rna_tpm_matrix.tsv.gz or transformed count data as the expression input,
# then validate key targets with mutation/CNV/RPPA layers before claiming mechanism.
