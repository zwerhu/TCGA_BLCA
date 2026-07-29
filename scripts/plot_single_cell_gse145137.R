#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- if (!is.na(script_arg)) sub("^--file=", "", script_arg) else "scripts/plot_single_cell_gse145137.R"
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
fig_dir <- file.path(root, "results", "figures")
report_path <- file.path(root, "results", "single_cell_gse145137_report.html")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pca_file <- file.path(root, "data", "external", "single_cell", "GSE145137", "GSE145137_primary_pca_coordinates.tsv")
meta_file <- file.path(root, "data", "external", "single_cell", "GSE145137", "GSE145137_primary_cell_metadata_with_reduced_model_scores.tsv")
expr_file <- file.path(root, "data", "external", "single_cell", "GSE145137", "GSE145137_primary_reduced_model_gene_expression_long.tsv")
summary_file <- file.path(root, "results", "tables", "single_cell_gse145137_risk_by_celltype_summary.tsv")
gene_summary_file <- file.path(root, "results", "tables", "single_cell_gse145137_model_gene_celltype_summary.tsv")
contrib_file <- file.path(root, "data", "external", "single_cell", "GSE145137", "GSE145137_primary_reduced_model_gene_contributions_long.tsv")

pca <- read.delim(pca_file, check.names = FALSE)
meta <- read.delim(meta_file, check.names = FALSE)
expr <- read.delim(expr_file, check.names = FALSE)
risk_summary <- read.delim(summary_file, check.names = FALSE)
gene_summary <- read.delim(gene_summary_file, check.names = FALSE)
contrib <- read.delim(contrib_file, check.names = FALSE)

model_genes <- c("EMP1", "AHNAK", "TNFRSF14", "CLEC2D", "GSDMB")
cell_order <- risk_summary$cell_type
cell_order <- unique(cell_order)
expr$cell_type <- factor(expr$cell_type, levels = rev(cell_order))
gene_summary$cell_type <- factor(gene_summary$cell_type, levels = rev(cell_order))
gene_summary$gene <- factor(gene_summary$gene, levels = model_genes)
pca$cell_type <- factor(pca$cell_type, levels = cell_order)
meta$cell_type <- factor(meta$cell_type, levels = cell_order)

save_png <- function(filename, expr, width = 1800, height = 1300) {
  path <- file.path(fig_dir, filename)
  png(path, width = width, height = height, res = 180)
  par(family = "sans")
  force(expr)
  dev.off()
  path
}

nice_cols <- grDevices::hcl.colors(length(cell_order), "Dark 3")
names(nice_cols) <- cell_order

p1 <- save_png("single_cell_gse145137_pca_celltype.png", {
  par(mar = c(4.5, 4.5, 3.2, 2), xpd = NA)
  plot(
    pca$PC1, pca$PC2,
    col = adjustcolor(nice_cols[as.character(pca$cell_type)], alpha.f = 0.75),
    pch = 16, cex = 0.7,
    xlab = "PC1", ylab = "PC2",
    main = "GSE145137 primary BLCA scRNA: cell-type structure"
  )
  centers <- aggregate(cbind(PC1, PC2) ~ cell_type, data = pca, FUN = median)
  text(centers$PC1, centers$PC2, labels = centers$cell_type, cex = 0.7, font = 2, pos = 4)
})

p2 <- save_png("single_cell_gse145137_pca_reduced_score.png", {
  z <- pca$reduced_model_score_z
  pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
  z_clip <- pmax(pmin(z, quantile(z, 0.98, na.rm = TRUE)), quantile(z, 0.02, na.rm = TRUE))
  idx <- as.integer(cut(z_clip, breaks = 101, include.lowest = TRUE))
  par(mar = c(4.5, 4.5, 3.2, 5))
  plot(
    pca$PC1, pca$PC2,
    col = adjustcolor(pal[idx], alpha.f = 0.78),
    pch = 16, cex = 0.7,
    xlab = "PC1", ylab = "PC2",
    main = "Reduced-model score projected onto primary BLCA single cells"
  )
  usr <- par("usr")
  yseq <- seq(usr[3], usr[4], length.out = 101)
  xleft <- usr[2] + 0.02 * diff(usr[1:2])
  rect(xleft, yseq[-101], xleft + 0.025 * diff(usr[1:2]), yseq[-1], col = pal, border = NA, xpd = NA)
  text(xleft + 0.045 * diff(usr[1:2]), usr[4], "High", adj = c(0, 1), cex = 0.75, xpd = NA)
  text(xleft + 0.045 * diff(usr[1:2]), usr[3], "Low", adj = c(0, 0), cex = 0.75, xpd = NA)
})

p3 <- save_png("single_cell_gse145137_risk_by_celltype.png", height = 1400, {
  meta$cell_type_plot <- factor(as.character(meta$cell_type), levels = rev(cell_order))
  par(mar = c(5, 11, 3.2, 2))
  boxplot(
    reduced_model_score_z ~ cell_type_plot, data = meta,
    horizontal = TRUE, las = 1, outline = FALSE,
    col = "#D9E4F5", border = "#596579",
    xlab = "Reduced-model score, z-scaled across cells",
    ylab = "", main = "Single-cell reduced-model score by annotated cell type"
  )
  stripchart(
    reduced_model_score_z ~ cell_type_plot, data = meta,
    vertical = FALSE, method = "jitter", pch = 16,
    col = adjustcolor("#34568B", alpha.f = 0.22),
    cex = 0.45, add = TRUE
  )
  abline(v = 0, lty = 2, col = "grey45")
})

p4 <- save_png("single_cell_gse145137_model_gene_dotplot.png", width = 2100, height = 1300, {
  df <- gene_summary[gene_summary$gene %in% model_genes, ]
  df$x <- match(as.character(df$gene), model_genes)
  df$y <- match(as.character(df$cell_type), levels(df$cell_type))
  zlim <- max(abs(df$mean_z_by_gene), na.rm = TRUE)
  pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
  idx <- as.integer(cut(df$mean_z_by_gene, breaks = seq(-zlim, zlim, length.out = 102), include.lowest = TRUE))
  sizes <- sqrt(pmax(df$pct_positive, 1)) / sqrt(100) * 2.5
  par(mar = c(5, 11, 3.2, 10), xpd = NA)
  plot(
    NA, xlim = c(0.5, length(model_genes) + 0.5), ylim = c(0.5, length(levels(df$cell_type)) + 0.5),
    xaxt = "n", yaxt = "n", xlab = "", ylab = "",
    main = "Reduced-model gene expression across cell types"
  )
  axis(1, at = seq_along(model_genes), labels = model_genes)
  axis(2, at = seq_along(levels(df$cell_type)), labels = levels(df$cell_type), las = 1)
  abline(h = seq_along(levels(df$cell_type)), v = seq_along(model_genes), col = "grey90")
  symbols(df$x, df$y, circles = sizes, inches = 0.18, bg = pal[idx], fg = "grey35", add = TRUE)
  legend("right", inset = c(-0.24, 0), legend = c("25%", "50%", "75%"),
         pt.cex = sqrt(c(25, 50, 75)) / sqrt(100) * 2.5, pch = 21, pt.bg = "grey85",
         bty = "n", title = "Detected", cex = 0.75)
})

p5 <- save_png("single_cell_gse145137_model_gene_contribution_heatmap.png", height = 1300, {
  cdat <- merge(contrib, meta[, c("cell", "cell_type")], by = "cell")
  agg <- aggregate(contribution ~ cell_type + gene, data = cdat, FUN = mean)
  mat <- matrix(0, nrow = length(cell_order), ncol = length(model_genes),
                dimnames = list(cell_order, model_genes))
  for (i in seq_len(nrow(agg))) {
    mat[as.character(agg$cell_type[i]), as.character(agg$gene[i])] <- agg$contribution[i]
  }
  mat <- mat[rev(cell_order), model_genes, drop = FALSE]
  zlim <- max(abs(mat), na.rm = TRUE)
  pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
  br <- seq(-zlim, zlim, length.out = 102)
  par(mar = c(5, 11, 3.2, 5))
  image(
    x = seq_len(ncol(mat)), y = seq_len(nrow(mat)), z = t(mat),
    col = pal, breaks = br, axes = FALSE,
    xlab = "", ylab = "",
    main = "Average gene-level contribution to reduced-model score"
  )
  axis(1, at = seq_len(ncol(mat)), labels = colnames(mat))
  axis(2, at = seq_len(nrow(mat)), labels = rownames(mat), las = 1)
  for (i in seq_len(nrow(mat))) {
    for (j in seq_len(ncol(mat))) {
      text(j, i, sprintf("%.3f", mat[i, j]), cex = 0.62)
    }
  }
})

p6 <- save_png("single_cell_gse145137_celltype_composition.png", height = 1200, {
  counts <- risk_summary$n_cells
  names(counts) <- risk_summary$cell_type
  par(mar = c(5, 11, 3.2, 2))
  barplot(
    rev(counts), horiz = TRUE, las = 1, col = "#A9C5E8", border = NA,
    xlab = "Number of cells", main = "Annotated cell-type composition"
  )
})

html <- c(
  "<!doctype html>",
  "<html><head><meta charset='utf-8'>",
  "<title>GSE145137 single-cell reduced-model report</title>",
  "<style>body{font-family:Arial,sans-serif;margin:28px;color:#222;max-width:1100px} h1{font-size:26px} h2{font-size:18px;margin-top:28px} img{max-width:100%;border:1px solid #ddd;margin:8px 0 18px 0} table{border-collapse:collapse} td,th{border-bottom:1px solid #ddd;padding:6px 8px;text-align:left} code{background:#f3f3f3;padding:2px 4px}.note{color:#555}</style>",
  "</head><body>",
  "<h1>GSE145137 primary BLCA single-cell analysis</h1>",
  "<p class='note'>Dataset: GSE145137, primary sample BC159-T#3. Expression: processed log2TPM matrix. Cell labels: author-provided cell-type annotation. Score: sum of TCGA reduced-model coefficients multiplied by gene-wise z-scored expression across cells.</p>",
  "<h2>Key takeaways</h2>",
  "<ul>",
  "<li>All five reduced-model genes are present in the primary single-cell matrix.</li>",
  "<li>Higher single-cell model scores are concentrated in basal tumor cells, endothelial cells, and fibroblasts.</li>",
  "<li>T cells show lower scores, consistent with negative coefficients for immune-associated genes such as TNFRSF14 and CLEC2D.</li>",
  "<li>This supports interpreting the bulk TCGA risk phenotype as a basal/EMT-stromal and vascular microenvironment signal rather than a purely tumor-cell intrinsic survival score.</li>",
  "</ul>",
  "<h2>Cell-type structure</h2>",
  "<img src='figures/single_cell_gse145137_pca_celltype.png'>",
  "<h2>Reduced-model score on PCA</h2>",
  "<img src='figures/single_cell_gse145137_pca_reduced_score.png'>",
  "<h2>Score by cell type</h2>",
  "<img src='figures/single_cell_gse145137_risk_by_celltype.png'>",
  "<h2>Model gene expression</h2>",
  "<img src='figures/single_cell_gse145137_model_gene_dotplot.png'>",
  "<h2>Gene-level contribution</h2>",
  "<img src='figures/single_cell_gse145137_model_gene_contribution_heatmap.png'>",
  "<h2>Cell-type composition</h2>",
  "<img src='figures/single_cell_gse145137_celltype_composition.png'>",
  "<h2>Output tables</h2>",
  "<ul>",
  "<li><code>results/tables/single_cell_gse145137_risk_by_celltype_summary.tsv</code></li>",
  "<li><code>results/tables/single_cell_gse145137_model_gene_celltype_summary.tsv</code></li>",
  "<li><code>data/external/single_cell/GSE145137/GSE145137_primary_cell_metadata_with_reduced_model_scores.tsv</code></li>",
  "</ul>",
  "<p class='note'>Caution: GSE145137 contains one primary tumor sample; it is best used for cellular localization and biological interpretation, not as an independent prognostic validation cohort.</p>",
  "</body></html>"
)
writeLines(html, report_path, useBytes = TRUE)

cat("Figures:\n")
cat(paste(c(p1, p2, p3, p4, p5, p6), collapse = "\n"), "\n")
cat("Report:\n", report_path, "\n", sep = "")
