#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- if (!is.na(script_arg)) sub("^--file=", "", script_arg) else "scripts/plot_candidate_scRNA_geo_comparison.R"
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
fig_dir <- file.path(root, "results", "figures")
table_dir <- file.path(root, "results", "tables")
report_path <- file.path(root, "results", "candidate_geo_single_cell_comparison_report.html")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

celltype <- read.delim(file.path(table_dir, "candidate_geo_reduced_model_by_inferred_celltype.tsv"), check.names = FALSE)
sample <- read.delim(file.path(table_dir, "candidate_geo_reduced_model_by_sample.tsv"), check.names = FALSE)
status <- read.delim(file.path(table_dir, "candidate_geo_dataset_status_and_summary.tsv"), check.names = FALSE)
gene <- read.delim(file.path(table_dir, "candidate_geo_model_gene_by_inferred_celltype.tsv"), check.names = FALSE)

score_cols <- c("median_score_z", "mean_score_z", "high_fraction", "n_units")
for (cc in score_cols) {
  if (cc %in% names(celltype)) celltype[[cc]] <- as.numeric(celltype[[cc]])
  if (cc %in% names(sample)) sample[[cc]] <- as.numeric(sample[[cc]])
}
gene$mean_log_norm_expr <- as.numeric(gene$mean_log_norm_expr)
gene$pct_positive <- as.numeric(gene$pct_positive)

dataset_order <- c("GSE222315", "GSE293189", "GSE171351", "GSE269877")
celltype_order <- c(
  "Endothelial", "Fibroblast", "Smooth_muscle_pericyte",
  "Epithelial_basal", "Epithelial_luminal",
  "Myeloid", "B_Plasma", "T_NK", "Proliferating", "Unknown"
)
model_genes <- c("EMP1", "AHNAK", "TNFRSF14", "CLEC2D", "GSDMB")

save_png <- function(filename, expr, width = 1800, height = 1200) {
  path <- file.path(fig_dir, filename)
  png(path, width = width, height = height, res = 180)
  par(family = "sans")
  force(expr)
  dev.off()
  path
}

p1 <- save_png("candidate_geo_celltype_score_heatmap.png", {
  dat <- subset(celltype, inferred_cell_type %in% celltype_order)
  mat <- matrix(NA_real_, nrow = length(celltype_order), ncol = length(dataset_order),
                dimnames = list(celltype_order, dataset_order))
  for (i in seq_len(nrow(dat))) {
    mat[dat$inferred_cell_type[i], dat$dataset[i]] <- dat$median_score_z[i]
  }
  mat <- mat[, colSums(!is.na(mat)) > 0, drop = FALSE]
  zlim <- max(abs(mat), na.rm = TRUE)
  pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
  br <- seq(-zlim, zlim, length.out = 102)
  par(mar = c(6, 11, 3.5, 2))
  image(
    x = seq_len(ncol(mat)), y = seq_len(nrow(mat)), z = t(mat[nrow(mat):1, , drop = FALSE]),
    col = pal, breaks = br, axes = FALSE, xlab = "", ylab = "",
    main = "Median reduced-model score by inferred cell type"
  )
  axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2)
  axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 1)
  for (i in seq_len(nrow(mat))) {
    for (j in seq_len(ncol(mat))) {
      val <- mat[nrow(mat) - i + 1, j]
      if (is.finite(val)) text(j, i, sprintf("%.2f", val), cex = 0.72)
    }
  }
})

p2 <- save_png("candidate_geo_celltype_score_dotplot.png", width = 2100, height = 1300, {
  dat <- subset(celltype, inferred_cell_type %in% celltype_order)
  dat$dataset <- factor(dat$dataset, levels = dataset_order)
  dat$inferred_cell_type <- factor(dat$inferred_cell_type, levels = rev(celltype_order))
  dat <- dat[!is.na(dat$dataset), ]
  zlim <- max(abs(dat$median_score_z), na.rm = TRUE)
  pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
  idx <- as.integer(cut(dat$median_score_z, breaks = seq(-zlim, zlim, length.out = 102), include.lowest = TRUE))
  sizes <- sqrt(pmax(dat$n_units, 1)) / sqrt(max(dat$n_units, na.rm = TRUE)) * 2.8
  par(mar = c(6, 11, 3.5, 10), xpd = NA)
  plot(
    NA, xlim = c(0.5, length(dataset_order) - 0.5), ylim = c(0.5, length(celltype_order) + 0.5),
    axes = FALSE, xlab = "", ylab = "",
    main = "Cell-type localization of the reduced-model score"
  )
  used_ds <- levels(droplevels(dat$dataset))
  axis(1, at = seq_along(used_ds), labels = used_ds, las = 2)
  axis(2, at = seq_along(levels(dat$inferred_cell_type)), labels = levels(dat$inferred_cell_type), las = 1)
  abline(h = seq_along(levels(dat$inferred_cell_type)), v = seq_along(used_ds), col = "grey90")
  x <- match(as.character(dat$dataset), used_ds)
  y <- match(as.character(dat$inferred_cell_type), levels(dat$inferred_cell_type))
  symbols(x, y, circles = sizes, inches = 0.16, bg = pal[idx], fg = "grey35", add = TRUE)
  legend("right", inset = c(-0.23, 0), legend = c("1k", "10k", "50k"),
         pt.cex = sqrt(c(1000, 10000, 50000)) / sqrt(max(dat$n_units, na.rm = TRUE)) * 2.8,
         pch = 21, pt.bg = "grey85", bty = "n", title = "Units", cex = 0.75)
})

p3 <- save_png("candidate_geo_sample_score_by_condition.png", width = 1900, height = 1300, {
  dat <- subset(sample, dataset %in% c("GSE222315", "GSE293189"))
  dat$condition <- factor(dat$condition, levels = c("Tumor", "CIS", "Tumor_variant_or_margin", "Normal", "NAT", "Unknown"))
  dat$dataset <- factor(dat$dataset, levels = c("GSE222315", "GSE293189"))
  par(mar = c(7, 5, 3.5, 2))
  x <- as.numeric(dat$dataset) + (as.numeric(dat$condition) - 3) * 0.08
  cols <- c(Tumor = "#B2182B", CIS = "#D6604D", Tumor_variant_or_margin = "#F4A582",
            Normal = "#2166AC", NAT = "#4393C3", Unknown = "grey60")
  plot(
    x, dat$median_score_z, pch = 21, bg = cols[as.character(dat$condition)], col = "grey30",
    cex = sqrt(dat$n_units) / sqrt(max(dat$n_units, na.rm = TRUE)) * 2.5,
    xaxt = "n", xlab = "", ylab = "Sample median reduced-model score (z)",
    main = "Sample-level score distribution in scRNA cohorts"
  )
  axis(1, at = seq_along(levels(dat$dataset)), labels = levels(dat$dataset), las = 2)
  abline(h = 0, lty = 2, col = "grey50")
  legend("topright", legend = names(cols), pt.bg = cols, pch = 21, bty = "n", cex = 0.72)
})

p4 <- save_png("candidate_geo_model_gene_dotplot.png", width = 2200, height = 1500, {
  dat <- subset(gene, dataset %in% c("GSE222315", "GSE293189") & gene %in% model_genes &
                  inferred_cell_type %in% c("Endothelial", "Fibroblast", "Epithelial_basal", "Epithelial_luminal", "T_NK", "Myeloid"))
  dat$key <- paste(dat$dataset, dat$inferred_cell_type, sep = " | ")
  key_order <- unique(dat$key[order(dat$dataset, match(dat$inferred_cell_type, celltype_order))])
  dat$key <- factor(dat$key, levels = rev(key_order))
  dat$gene <- factor(dat$gene, levels = model_genes)
  z <- ave(dat$mean_log_norm_expr, dat$gene, FUN = function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  })
  zlim <- max(abs(z), na.rm = TRUE)
  pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
  idx <- as.integer(cut(z, breaks = seq(-zlim, zlim, length.out = 102), include.lowest = TRUE))
  sizes <- sqrt(pmax(dat$pct_positive, 1)) / 10 * 2.2
  par(mar = c(5, 17, 3.5, 8), xpd = NA)
  plot(
    NA, xlim = c(0.5, length(model_genes) + 0.5), ylim = c(0.5, length(levels(dat$key)) + 0.5),
    axes = FALSE, xlab = "", ylab = "",
    main = "Reduced-model gene expression in major inferred cell types"
  )
  axis(1, at = seq_along(model_genes), labels = model_genes)
  axis(2, at = seq_along(levels(dat$key)), labels = levels(dat$key), las = 1)
  abline(h = seq_along(levels(dat$key)), v = seq_along(model_genes), col = "grey90")
  symbols(
    match(as.character(dat$gene), model_genes), match(as.character(dat$key), levels(dat$key)),
    circles = sizes, inches = 0.15, bg = pal[idx], fg = "grey35", add = TRUE
  )
  legend("right", inset = c(-0.16, 0), legend = c("25%", "50%", "75%"),
         pt.cex = sqrt(c(25, 50, 75)) / 10 * 2.2, pch = 21, pt.bg = "grey85",
         bty = "n", title = "Detected", cex = 0.75)
})

rank_rows <- data.frame()
for (ds in unique(celltype$dataset)) {
  sub <- subset(celltype, dataset == ds)
  high_types <- c("Endothelial", "Fibroblast", "Epithelial_basal", "Smooth_muscle_pericyte")
  low_types <- c("T_NK", "B_Plasma", "Myeloid")
  high <- weighted.mean(sub$median_score_z[sub$inferred_cell_type %in% high_types],
                        sub$n_units[sub$inferred_cell_type %in% high_types], na.rm = TRUE)
  low <- weighted.mean(sub$median_score_z[sub$inferred_cell_type %in% low_types],
                       sub$n_units[sub$inferred_cell_type %in% low_types], na.rm = TRUE)
  rank_rows <- rbind(rank_rows, data.frame(dataset = ds, stromal_basal_vs_immune_delta = high - low))
}
rank_rows <- rank_rows[order(-rank_rows$stromal_basal_vs_immune_delta), ]
write.table(rank_rows, file.path(table_dir, "candidate_geo_dataset_support_ranking.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

top_celltype <- celltype[order(celltype$dataset, -celltype$median_score_z), ]
top_celltype <- do.call(rbind, by(top_celltype, top_celltype$dataset, head, n = 3))
top_html <- paste(apply(top_celltype[, c("dataset", "assay_unit", "inferred_cell_type", "n_units", "median_score_z", "high_fraction")], 1, function(x) {
  paste0("<tr><td>", paste(x, collapse = "</td><td>"), "</td></tr>")
}), collapse = "\n")

rank_html <- paste(apply(rank_rows, 1, function(x) {
  paste0("<tr><td>", x[["dataset"]], "</td><td>", sprintf("%.3f", as.numeric(x[["stromal_basal_vs_immune_delta"]])), "</td></tr>")
}), collapse = "\n")

html <- c(
  "<!doctype html>",
  "<html><head><meta charset='utf-8'>",
  "<title>Candidate GEO single-cell comparison</title>",
  "<style>body{font-family:Arial,sans-serif;margin:28px;color:#222;max-width:1120px} h1{font-size:26px} h2{font-size:18px;margin-top:28px} img{max-width:100%;border:1px solid #ddd;margin:8px 0 18px 0} table{border-collapse:collapse;margin:10px 0 20px 0} td,th{border-bottom:1px solid #ddd;padding:6px 8px;text-align:left} code{background:#f3f3f3;padding:2px 4px}.note{color:#555}</style>",
  "</head><body>",
  "<h1>Candidate GEO single-cell/spatial validation of the reduced TCGA-BLCA model</h1>",
  "<p class='note'>Datasets tested: GSE222315, GSE293189, GSE171351, and GSE269877. Cell-type labels in GSE222315/GSE293189 are coarse marker-based inference; GSE171351 is Visium spatial spots, not single cells.</p>",
  "<h2>Best-supported dataset ranking</h2>",
  "<table><tr><th>Dataset</th><th>Stromal/basal vs immune delta</th></tr>",
  rank_html,
  "</table>",
  "<h2>Top score-localized cell/spot types</h2>",
  "<table><tr><th>Dataset</th><th>Assay</th><th>Inferred type</th><th>Units</th><th>Median score z</th><th>High fraction</th></tr>",
  top_html,
  "</table>",
  "<h2>Cross-dataset heatmap</h2>",
  "<img src='figures/candidate_geo_celltype_score_heatmap.png'>",
  "<h2>Cell-type localization</h2>",
  "<img src='figures/candidate_geo_celltype_score_dotplot.png'>",
  "<h2>Sample-level score</h2>",
  "<img src='figures/candidate_geo_sample_score_by_condition.png'>",
  "<h2>Model gene expression</h2>",
  "<img src='figures/candidate_geo_model_gene_dotplot.png'>",
  "<h2>Output tables</h2>",
  "<ul>",
  "<li><code>results/tables/candidate_geo_dataset_status_and_summary.tsv</code></li>",
  "<li><code>results/tables/candidate_geo_reduced_model_by_inferred_celltype.tsv</code></li>",
  "<li><code>results/tables/candidate_geo_reduced_model_by_sample.tsv</code></li>",
  "<li><code>results/tables/candidate_geo_dataset_support_ranking.tsv</code></li>",
  "</ul>",
  "<p class='note'>GSE269877 was not downloaded in this lightweight run because supplementary RDS files are approximately 13GB and 2.4GB compressed. It should be handled as a separate long-running analysis if needed.</p>",
  "</body></html>"
)
writeLines(html, report_path, useBytes = TRUE)

cat("Figures:\n")
cat(paste(c(p1, p2, p3, p4), collapse = "\n"), "\n")
cat("Report:\n", report_path, "\n", sep = "")
