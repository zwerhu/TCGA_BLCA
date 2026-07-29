#!/usr/bin/env Rscript
# TCGA-BLCA exploratory differential expression screen.
#
# Uses Wilcoxon tests on log2(TPM + 1) tumor vs normal samples. This is useful
# for fast exploration and visualization; for manuscript-grade RNA-seq DE,
# rerun with DESeq2/edgeR on raw counts.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE) else file.path(getwd(), "scripts", "differential_expression_screen.R")
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

processed_dir <- file.path(project_root, "data", "processed")
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE)
palette_main <- c("#2F6F95", "#D9603D", "#4B8F5E", "#8E5EA2", "#D5A021", "#6A7FDB")

save_png <- function(file, expr, width = 1800, height = 1200, res = 180) {
  png(filename = file.path(figure_dir, file), width = width, height = height, res = res)
  old <- par(no.readonly = TRUE)
  on.exit({
    par(old)
    dev.off()
  }, add = TRUE)
  force(expr)
}

sample_map <- read.delim(file.path(processed_dir, "blca_rna_sample_map.tsv"), check.names = FALSE)
expr <- read.delim(gzfile(file.path(processed_dir, "blca_log2_tpm_gene_symbol_all_samples.tsv.gz")),
                   row.names = 1, check.names = FALSE)
expr <- as.matrix(expr)
mode(expr) <- "numeric"

tumor_cols <- sample_map$matrix_column[sample_map$sample_type == "Primary Tumor"]
normal_cols <- sample_map$matrix_column[sample_map$sample_type == "Solid Tissue Normal"]
tumor_cols <- intersect(tumor_cols, colnames(expr))
normal_cols <- intersect(normal_cols, colnames(expr))
if (length(tumor_cols) < 10 || length(normal_cols) < 3) {
  stop("Not enough tumor/normal samples for differential expression screening.")
}

gene_mean <- rowMeans(expr, na.rm = TRUE)
gene_sd <- apply(expr, 1, sd, na.rm = TRUE)
eligible <- which(is.finite(gene_mean) & is.finite(gene_sd) & gene_mean > 0.1 & gene_sd > 0.05)

message("Tumor samples: ", length(tumor_cols), "; normal samples: ", length(normal_cols))
message("Testing genes: ", length(eligible))

pvals <- rep(NA_real_, length(eligible))
tumor_mean <- rowMeans(expr[eligible, tumor_cols, drop = FALSE], na.rm = TRUE)
normal_mean <- rowMeans(expr[eligible, normal_cols, drop = FALSE], na.rm = TRUE)
log2fc <- tumor_mean - normal_mean

for (i in seq_along(eligible)) {
  idx <- eligible[i]
  pvals[i] <- tryCatch(
    wilcox.test(expr[idx, tumor_cols], expr[idx, normal_cols], exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  if (i %% 5000 == 0) message("  finished ", i, "/", length(eligible))
}

res <- data.frame(
  gene = rownames(expr)[eligible],
  tumor_mean_log2tpm = tumor_mean,
  normal_mean_log2tpm = normal_mean,
  log2fc_tumor_vs_normal = log2fc,
  pvalue = pvals,
  FDR = p.adjust(pvals, method = "BH"),
  mean_log2tpm = gene_mean[eligible],
  sd_log2tpm = gene_sd[eligible],
  stringsAsFactors = FALSE
)
res <- res[order(res$FDR, -abs(res$log2fc_tumor_vs_normal)), ]
write.table(res, file.path(table_dir, "wilcoxon_de_tumor_vs_normal.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(head(res, 200), file.path(table_dir, "wilcoxon_de_top200.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

sig <- res$FDR < 0.05 & abs(res$log2fc_tumor_vs_normal) >= 1
de_summary <- data.frame(
  tumor_samples = length(tumor_cols),
  normal_samples = length(normal_cols),
  genes_tested = nrow(res),
  significant_abs_log2fc_1_fdr_0_05 = sum(sig, na.rm = TRUE),
  upregulated = sum(sig & res$log2fc_tumor_vs_normal > 0, na.rm = TRUE),
  downregulated = sum(sig & res$log2fc_tumor_vs_normal < 0, na.rm = TRUE)
)
write.table(de_summary, file.path(table_dir, "wilcoxon_de_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

save_png("de_wilcoxon_volcano.png", {
  x <- res$log2fc_tumor_vs_normal
  y <- -log10(pmax(res$FDR, .Machine$double.xmin))
  plot(x, y, pch = 16, cex = 0.45, col = adjustcolor("gray35", alpha.f = 0.35),
       xlab = "Tumor - normal mean log2(TPM + 1)", ylab = "-log10(FDR)",
       main = "Exploratory differential expression")
  points(x[sig & x > 0], y[sig & x > 0], pch = 16, cex = 0.55, col = adjustcolor(palette_main[2], alpha.f = 0.75))
  points(x[sig & x < 0], y[sig & x < 0], pch = 16, cex = 0.55, col = adjustcolor(palette_main[1], alpha.f = 0.75))
  abline(v = c(-1, 1), col = "gray60", lty = 2)
  abline(h = -log10(0.05), col = "gray60", lty = 2)
  top <- head(res[is.finite(res$FDR), ], 12)
  text(top$log2fc_tumor_vs_normal, -log10(pmax(top$FDR, .Machine$double.xmin)), labels = top$gene,
       cex = 0.65, pos = ifelse(top$log2fc_tumor_vs_normal >= 0, 4, 2), offset = 0.25)
})

top_up <- head(res[res$log2fc_tumor_vs_normal > 0 & is.finite(res$FDR), ], 25)
top_down <- head(res[res$log2fc_tumor_vs_normal < 0 & is.finite(res$FDR), ], 25)
heat_genes <- unique(c(top_up$gene, top_down$gene))
heat_expr <- expr[heat_genes, c(normal_cols, tumor_cols), drop = FALSE]
heat_expr <- heat_expr[, c(seq_len(length(normal_cols)), length(normal_cols) + seq_len(min(120, length(tumor_cols)))), drop = FALSE]
heat_scaled <- t(scale(t(heat_expr)))
heat_scaled[!is.finite(heat_scaled)] <- 0

save_png("de_top50_heatmap.png", width = 2100, height = 1700, {
  par(mar = c(4, 10, 4, 1))
  heatmap(heat_scaled, Rowv = NA, Colv = NA, scale = "none", labCol = FALSE,
          col = colorRampPalette(c("#2F6F95", "white", "#D9603D"))(80),
          main = "Top differential genes: normal plus representative tumors")
})

save_png("de_top_log2fc_bar.png", width = 2000, height = 1300, {
  top_bar <- rbind(head(res[order(res$log2fc_tumor_vs_normal), ], 15),
                   head(res[order(-res$log2fc_tumor_vs_normal), ], 15))
  vals <- top_bar$log2fc_tumor_vs_normal
  names(vals) <- make.unique(top_bar$gene)
  cols <- ifelse(vals > 0, palette_main[2], palette_main[1])
  pad <- diff(range(vals, na.rm = TRUE)) * 0.16
  par(mar = c(5, 9, 4, 1))
  bp <- barplot(vals, horiz = TRUE, las = 1, col = cols, border = NA,
                xlim = c(min(vals) - pad, max(vals) + pad),
                main = "Largest tumor-normal expression differences",
                xlab = "Tumor - normal mean log2(TPM + 1)")
  abline(v = 0, col = "gray40")
  text(vals, bp, labels = sprintf("%.2f", vals), pos = ifelse(vals > 0, 4, 2), cex = 0.65)
})

html <- c(
  "<!doctype html>",
  "<html><head><meta charset=\"utf-8\"><title>TCGA-BLCA differential expression</title>",
  "<style>body{font-family:Arial,sans-serif;margin:24px;color:#222}h1{font-size:24px}section{margin:24px 0 36px}img{max-width:100%;border:1px solid #ddd}code{background:#f3f3f3;padding:2px 4px}</style>",
  "</head><body>",
  "<h1>TCGA-BLCA differential expression</h1>",
  "<p>Exploratory Wilcoxon screen on log2(TPM + 1). Use DESeq2 on counts for final manuscript statistics.</p>",
  "<section><h2>Volcano</h2><img src=\"figures/de_wilcoxon_volcano.png\" alt=\"Wilcoxon differential expression volcano\"></section>",
  "<section><h2>Top genes heatmap</h2><img src=\"figures/de_top50_heatmap.png\" alt=\"Top differential genes heatmap\"></section>",
  "<section><h2>Largest expression differences</h2><img src=\"figures/de_top_log2fc_bar.png\" alt=\"Top expression differences\"></section>",
  "</body></html>"
)
writeLines(html, file.path(project_root, "results", "differential_expression_report.html"))

message("Significant genes (FDR < 0.05 and abs(log2FC) >= 1): ", sum(sig, na.rm = TRUE))
message("Upregulated: ", de_summary$upregulated, "; downregulated: ", de_summary$downregulated)
message("Differential expression report complete: ", file.path(project_root, "results", "differential_expression_report.html"))
