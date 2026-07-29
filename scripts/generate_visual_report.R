#!/usr/bin/env Rscript
# TCGA-BLCA visual report using base R plus survival.

project_root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = FALSE)
if (is.na(project_root) || !dir.exists(project_root)) {
  project_root <- normalizePath(".", mustWork = FALSE)
}

processed_dir <- file.path(project_root, "data", "processed")
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE)
palette_main <- c("#2F6F95", "#D9603D", "#4B8F5E", "#8E5EA2", "#D5A021", "#6A7FDB", "#C44E52")

clean_stage <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "Not Reported"] <- NA
  x
}

save_png <- function(file, expr, width = 1800, height = 1200, res = 180) {
  png(filename = file.path(figure_dir, file), width = width, height = height, res = res)
  old <- par(no.readonly = TRUE)
  on.exit({
    par(old)
    dev.off()
  }, add = TRUE)
  force(expr)
}

barplot_labeled <- function(values, main, ylab, col = palette_main[1], las = 2) {
  bp <- barplot(values, col = col, border = NA, las = las, main = main, ylab = ylab)
  text(bp, values, labels = values, pos = 3, cex = 0.8)
}

clinical <- read.delim(file.path(processed_dir, "blca_clinical_survival.tsv"), check.names = FALSE)
sample_map <- read.delim(file.path(processed_dir, "blca_rna_sample_map.tsv"), check.names = FALSE)
survival_sheet <- read.delim(file.path(processed_dir, "blca_survival_expression_sample_sheet.tsv"), check.names = FALSE)
mutation_counts <- read.delim(file.path(table_dir, "mutation_counts_by_patient.tsv"), check.names = FALSE)
fc <- read.delim(file.path(table_dir, "exploratory_log2fc_tumor_vs_normal.tsv"), check.names = FALSE)

clinical$ajcc_pathologic_stage <- clean_stage(clinical$ajcc_pathologic_stage)
stage_upper <- toupper(clinical$ajcc_pathologic_stage)
clinical$stage_group <- ifelse(grepl("^STAGE III|^STAGE IV", stage_upper), "Stage III-IV",
  ifelse(grepl("^STAGE I($|[ABC ]|\\s)|^STAGE II($|[ABC ]|\\s)", stage_upper), "Stage I-II", NA)
)

# 1. Clinical stage distribution
save_png("clinical_stage_distribution.png", {
  par(mar = c(7, 5, 4, 1))
  stage_tab <- table(clinical$ajcc_pathologic_stage, useNA = "no")
  stage_tab <- stage_tab[order(names(stage_tab))]
  barplot_labeled(stage_tab, "Clinical stage distribution", "Patients", col = palette_main[1], las = 2)
})

# 2. Sample-type distribution
save_png("rna_sample_type_distribution.png", {
  par(mar = c(7, 5, 4, 1))
  sample_tab <- table(sample_map$sample_type)
  barplot_labeled(sample_tab, "RNA-seq sample types", "Samples", col = palette_main[2], las = 2)
})

# 3. Overall survival by stage group
save_png("overall_survival_by_stage.png", {
  library(survival)
  dat <- clinical[!is.na(clinical$os_time_days) & !is.na(clinical$os_event) & !is.na(clinical$stage_group), ]
  fit <- survfit(Surv(os_time_days, os_event) ~ stage_group, data = dat)
  par(mar = c(5, 5, 4, 2))
  plot(fit, col = palette_main[1:2], lwd = 2, xlab = "Days", ylab = "Overall survival probability",
       main = "Overall survival by stage group", mark.time = TRUE)
  grid(col = "gray85")
  legend_labels <- names(fit$strata)
  if (length(legend_labels) == 0) legend_labels <- sort(unique(dat$stage_group))
  legend_labels <- sub("^stage_group=", "", legend_labels)
  legend("bottomleft", legend = legend_labels, col = palette_main[seq_along(legend_labels)], lwd = 2, bty = "n")
  p <- tryCatch(survdiff(Surv(os_time_days, os_event) ~ stage_group, data = dat), error = function(e) NULL)
  if (!is.null(p)) {
    pval <- pchisq(p$chisq, length(p$n) - 1, lower.tail = FALSE)
    mtext(sprintf("Log-rank p = %.3g", pval), side = 3, line = 0.2, adj = 1)
  }
})

# 4. Overall survival by sex
save_png("overall_survival_by_gender.png", {
  library(survival)
  dat <- clinical[!is.na(clinical$os_time_days) & !is.na(clinical$os_event) & !is.na(clinical$gender), ]
  fit <- survfit(Surv(os_time_days, os_event) ~ gender, data = dat)
  par(mar = c(5, 5, 4, 2))
  plot(fit, col = palette_main[3:4], lwd = 2, xlab = "Days", ylab = "Overall survival probability",
       main = "Overall survival by gender", mark.time = TRUE)
  grid(col = "gray85")
  legend_labels <- names(fit$strata)
  if (length(legend_labels) == 0) legend_labels <- sort(unique(dat$gender))
  legend_labels <- sub("^gender=", "", legend_labels)
  legend("bottomleft", legend = legend_labels, col = palette_main[2 + seq_along(legend_labels)], lwd = 2, bty = "n")
})

# 5. Mutation burden distribution
save_png("mutation_burden_distribution.png", {
  par(mar = c(5, 5, 4, 1))
  hist(mutation_counts$nonsilent_mutations, breaks = 35, col = palette_main[5], border = "white",
       main = "Nonsilent mutation burden", xlab = "Nonsilent mutations per patient", ylab = "Patients")
  abline(v = median(mutation_counts$nonsilent_mutations, na.rm = TRUE), col = palette_main[2], lwd = 2)
  legend("topright", legend = sprintf("Median = %.0f", median(mutation_counts$nonsilent_mutations, na.rm = TRUE)),
         col = palette_main[2], lwd = 2, bty = "n")
})

# 6. Top mutated genes
maf_path <- file.path(processed_dir, "blca_somatic_mutations.maf.gz")
maf <- read.delim(gzfile(maf_path), check.names = FALSE, comment.char = "", quote = "")
nonsilent_exclude <- c("Silent", "Intron", "IGR", "3'UTR", "5'UTR", "RNA", "Targeted_Region", "Flank")
maf_ns <- maf[!(maf$Variant_Classification %in% nonsilent_exclude), ]
top_gene_counts <- sort(table(unique(maf_ns[, c("Patient_Barcode", "Hugo_Symbol")])$Hugo_Symbol), decreasing = TRUE)
top_gene_counts <- head(top_gene_counts[names(top_gene_counts) != ""], 20)
write.table(data.frame(gene = names(top_gene_counts), patients = as.integer(top_gene_counts)),
            file.path(table_dir, "top_mutated_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
save_png("top_mutated_genes.png", {
  par(mar = c(5, 9, 4, 1))
  vals <- rev(top_gene_counts)
  bp <- barplot(vals, horiz = TRUE, las = 1, col = palette_main[6], border = NA,
                xlim = c(0, max(vals, na.rm = TRUE) * 1.18),
                main = "Top nonsilently mutated genes", xlab = "Mutated patients")
  text(vals, bp, labels = vals, pos = 4, cex = 0.75)
})

# 7. Exploratory fold-change bar plot
protein_coding <- fc[fc$gene_type == "protein_coding" & is.finite(fc$exploratory_log2fc_tumor_vs_normal), ]
top_up <- head(protein_coding[order(-protein_coding$exploratory_log2fc_tumor_vs_normal), ], 15)
top_down <- head(protein_coding[order(protein_coding$exploratory_log2fc_tumor_vs_normal), ], 15)
fc_plot <- rbind(top_down, top_up)
fc_plot$label <- make.unique(fc_plot$gene_symbol)
save_png("exploratory_top_log2fc_genes.png", width = 2000, height = 1300, {
  par(mar = c(5, 9, 4, 1))
  vals <- fc_plot$exploratory_log2fc_tumor_vs_normal
  names(vals) <- fc_plot$label
  cols <- ifelse(vals > 0, palette_main[2], palette_main[1])
  label_pad <- diff(range(vals, na.rm = TRUE)) * 0.16
  bp <- barplot(vals, horiz = TRUE, las = 1, col = cols, border = NA,
                xlim = c(min(vals, na.rm = TRUE) - label_pad, max(vals, na.rm = TRUE) + label_pad),
                main = "Exploratory tumor vs normal log2TPM difference", xlab = "Tumor - normal mean log2(TPM + 1)")
  abline(v = 0, col = "gray40", lwd = 1)
  text(vals, bp, labels = sprintf("%.2f", vals), pos = ifelse(vals > 0, 4, 2), cex = 0.65)
})

# 8. RNA PCA, using the most variable protein-coding genes.
tpm <- read.delim(gzfile(file.path(processed_dir, "blca_rna_tpm_matrix.tsv.gz")), row.names = 1, check.names = FALSE)
annotation <- read.delim(file.path(processed_dir, "blca_gene_annotation.tsv"), check.names = FALSE)
common_gene_ids <- intersect(rownames(tpm), annotation$gene_id)
gene_type <- annotation$gene_type[match(common_gene_ids, annotation$gene_id)]
pc_genes <- common_gene_ids[gene_type == "protein_coding"]
expr <- log2(as.matrix(tpm[pc_genes, , drop = FALSE]) + 1)
vars <- apply(expr, 1, var, na.rm = TRUE)
keep <- names(sort(vars, decreasing = TRUE))[seq_len(min(3000, length(vars)))]
pca <- prcomp(t(expr[keep, , drop = FALSE]), center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  sample_type = sample_map$sample_type[match(rownames(pca$x), sample_map$matrix_column)]
)
write.table(pca_df, file.path(table_dir, "rna_pca_coordinates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
save_png("rna_pca_tumor_normal.png", {
  par(mar = c(5, 5, 4, 1))
  cols <- ifelse(pca_df$sample_type == "Solid Tissue Normal", palette_main[2], palette_main[1])
  pch <- ifelse(pca_df$sample_type == "Solid Tissue Normal", 17, 16)
  plot(pca_df$PC1, pca_df$PC2, pch = pch, col = adjustcolor(cols, alpha.f = 0.85),
       xlab = sprintf("PC1 (%.1f%%)", summary(pca)$importance[2, 1] * 100),
       ylab = sprintf("PC2 (%.1f%%)", summary(pca)$importance[2, 2] * 100),
       main = "RNA-seq PCA: tumor vs normal")
  grid(col = "gray88")
  legend("topright", legend = c("Primary Tumor", "Solid Tissue Normal"),
         col = c(palette_main[1], palette_main[2]), pch = c(16, 17), bty = "n")
})

# 9. RPPA top-variable heatmap
rppa <- read.delim(gzfile(file.path(processed_dir, "blca_rppa_matrix.tsv.gz")), row.names = 1, check.names = FALSE)
rppa_num <- as.matrix(rppa)
mode(rppa_num) <- "numeric"
rppa_vars <- apply(rppa_num, 1, var, na.rm = TRUE)
rppa_keep <- names(sort(rppa_vars, decreasing = TRUE))[seq_len(min(35, length(rppa_vars)))]
rppa_scaled <- t(scale(t(rppa_num[rppa_keep, , drop = FALSE])))
rppa_scaled[!is.finite(rppa_scaled)] <- 0
save_png("rppa_top_variable_heatmap.png", width = 1800, height = 1600, {
  par(mar = c(3, 11, 4, 1))
  heatmap(rppa_scaled, Rowv = NA, Colv = NA, scale = "none", labCol = FALSE,
          col = colorRampPalette(c("#2F6F95", "white", "#D9603D"))(80),
          main = "RPPA top-variable proteins")
})

# 10. CNV segment counts by chromosome
cnv <- read.delim(gzfile(file.path(processed_dir, "blca_masked_cnv_segments.tsv.gz")), check.names = FALSE)
chrom_col <- if ("Chromosome" %in% names(cnv)) "Chromosome" else names(cnv)[grep("chrom", names(cnv), ignore.case = TRUE)[1]]
cnv_chr <- cnv[[chrom_col]]
cnv_chr <- gsub("^chr", "", cnv_chr, ignore.case = TRUE)
chr_levels <- c(as.character(1:22), "X", "Y")
cnv_tab <- table(factor(cnv_chr, levels = chr_levels))
save_png("cnv_segments_by_chromosome.png", {
  par(mar = c(5, 5, 4, 1))
  barplot_labeled(cnv_tab, "CNV segment counts by chromosome", "Segments", col = palette_main[4], las = 1)
})

# 11. Multi-omics availability
availability <- data.frame(
  omics = c("Clinical OS", "RNA-seq", "Mutation MAF", "CNV segments", "RPPA"),
  samples_or_patients = c(
    sum(!is.na(clinical$os_time_days)),
    sum(sample_map$sample_type == "Primary Tumor"),
    length(unique(maf$Patient_Barcode)),
    length(unique(cnv$patient_barcode)),
    ncol(rppa)
  )
)
write.table(availability, file.path(table_dir, "multiomics_availability.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
save_png("multiomics_availability.png", {
  par(mar = c(7, 5, 4, 1))
  vals <- availability$samples_or_patients
  names(vals) <- availability$omics
  barplot_labeled(vals, "Multi-omics sample availability", "Samples / patients", col = palette_main[7], las = 2)
})

# Simple HTML gallery for quick viewing.
figs <- c(
  "clinical_stage_distribution.png",
  "rna_sample_type_distribution.png",
  "overall_survival_by_stage.png",
  "overall_survival_by_gender.png",
  "rna_pca_tumor_normal.png",
  "exploratory_top_log2fc_genes.png",
  "mutation_burden_distribution.png",
  "top_mutated_genes.png",
  "cnv_segments_by_chromosome.png",
  "rppa_top_variable_heatmap.png",
  "multiomics_availability.png"
)
html <- c(
  "<!doctype html>",
  "<html><head><meta charset=\"utf-8\"><title>TCGA-BLCA visual report</title>",
  "<style>body{font-family:Arial,sans-serif;margin:24px;color:#222}h1{font-size:24px}section{margin:24px 0 36px}img{max-width:100%;border:1px solid #ddd}code{background:#f3f3f3;padding:2px 4px}</style>",
  "</head><body>",
  "<h1>TCGA-BLCA visual report</h1>",
  sprintf("<p>Generated from <code>%s</code>.</p>", project_root),
  paste(sprintf("<section><h2>%s</h2><img src=\"figures/%s\" alt=\"%s\"></section>", tools::file_path_sans_ext(figs), figs, tools::file_path_sans_ext(figs)), collapse = "\n"),
  "</body></html>"
)
writeLines(html, file.path(project_root, "results", "visual_report.html"))

message("Visual report complete: ", file.path(project_root, "results", "visual_report.html"))
