#!/usr/bin/env Rscript
# Plot mutation/CNV/RPPA support analyses for fixed TCGA risk groups.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE) else file.path(getwd(), "scripts", "plot_multiomics_support.R")
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
palette_main <- c("#2F6F95", "#D9603D", "#4B8F5E", "#8E5EA2", "#D5A021")

save_png <- function(file, expr, width = 1800, height = 1200, res = 180) {
  png(filename = file.path(figure_dir, file), width = width, height = height, res = res)
  old <- par(no.readonly = TRUE)
  on.exit({
    par(old)
    dev.off()
  }, add = TRUE)
  force(expr)
}

mut <- read.delim(file.path(table_dir, "tcga_fixed_risk_mutation_burden_by_patient.tsv"), check.names = FALSE)
mutfreq <- read.delim(file.path(table_dir, "tcga_fixed_risk_top_mutated_gene_frequencies.tsv"), check.names = FALSE)
cnv <- read.delim(file.path(table_dir, "tcga_fixed_risk_cnv_burden_by_patient.tsv"), check.names = FALSE)
rppa <- read.delim(file.path(table_dir, "tcga_fixed_risk_rppa_differential.tsv"), check.names = FALSE)

save_png("tcga_fixed_risk_mutation_support.png", width = 2100, height = 1000, {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  boxplot(nonsilent_mutations ~ lasso_risk_group, data = mut, col = palette_main[1:2], border = "gray30",
          main = "Nonsilent mutation burden", xlab = "Risk group", ylab = "Nonsilent mutations")
  stripchart(nonsilent_mutations ~ lasso_risk_group, data = mut, vertical = TRUE, method = "jitter",
             pch = 16, col = adjustcolor("gray25", alpha.f = 0.35), add = TRUE)
  top <- head(mutfreq[order(-abs(mutfreq$frequency_diff_high_minus_low)), ], 15)
  vals <- rbind(top$low_frequency, top$high_frequency)
  colnames(vals) <- top$gene
  barplot(vals, beside = TRUE, las = 2, col = palette_main[1:2], border = NA,
          main = "Top mutation-frequency differences", ylab = "Mutated fraction")
  legend("topright", legend = c("Low risk", "High risk"), fill = palette_main[1:2], bty = "n")
})

save_png("tcga_fixed_risk_cnv_support.png", width = 2100, height = 1000, {
  par(mfrow = c(1, 3), mar = c(5, 5, 4, 1))
  boxplot(cnv_segment_count ~ lasso_risk_group, data = cnv, col = palette_main[1:2], border = "gray30",
          main = "CNV segment count", xlab = "Risk group", ylab = "Segments")
  boxplot(mean_abs_segment_mean ~ lasso_risk_group, data = cnv, col = palette_main[1:2], border = "gray30",
          main = "Mean abs CNV segment mean", xlab = "Risk group", ylab = "Mean abs segment mean")
  boxplot(total_altered_length_abs_0_2 / 1e9 ~ lasso_risk_group, data = cnv, col = palette_main[1:2], border = "gray30",
          main = "Altered length |seg.mean| >= 0.2", xlab = "Risk group", ylab = "Gb")
})

save_png("tcga_fixed_risk_rppa_support.png", {
  top <- head(rppa[order(rppa$FDR), ], 25)
  vals <- top$diff_high_minus_low
  names(vals) <- top$protein_id
  vals <- sort(vals)
  cols <- ifelse(vals > 0, palette_main[2], palette_main[1])
  par(mar = c(5, 9, 4, 1))
  pad <- diff(range(vals)) * 0.18
  bp <- barplot(vals, horiz = TRUE, las = 1, col = cols, border = NA,
                xlim = c(min(vals) - pad, max(vals) + pad),
                main = "Top RPPA differences by fixed risk group",
                xlab = "Mean RPPA difference, high - low")
  abline(v = 0, col = "gray40")
  text(vals, bp, labels = sprintf("FDR %.1e", top$FDR[match(names(vals), top$protein_id)]),
       pos = ifelse(vals > 0, 4, 2), cex = 0.6)
})

html <- c(
  "<!doctype html>",
  "<html><head><meta charset=\"utf-8\"><title>Multi-omics support</title>",
  "<style>body{font-family:Arial,sans-serif;margin:24px;color:#222}h1{font-size:24px}section{margin:24px 0 36px}img{max-width:100%;border:1px solid #ddd}</style>",
  "</head><body>",
  "<h1>TCGA fixed-risk multi-omics support</h1>",
  "<section><h2>Mutation support</h2><img src=\"figures/tcga_fixed_risk_mutation_support.png\"></section>",
  "<section><h2>CNV support</h2><img src=\"figures/tcga_fixed_risk_cnv_support.png\"></section>",
  "<section><h2>RPPA support</h2><img src=\"figures/tcga_fixed_risk_rppa_support.png\"></section>",
  "</body></html>"
)
writeLines(html, file.path(project_root, "results", "multiomics_support_report.html"))
message("Multi-omics support report complete: ", file.path(project_root, "results", "multiomics_support_report.html"))
