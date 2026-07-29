#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- if (!is.na(script_arg)) sub("^--file=", "", script_arg) else "scripts/extract_imvigor210_core_data.R"
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
pkg_data <- file.path(root, "data", "external", "immunotherapy", "IMvigor210", "package_extract", "IMvigor210CoreBiologies", "data")
out_dir <- file.path(root, "data", "external", "immunotherapy", "IMvigor210", "processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genes_needed <- unique(c(
  "EMP1", "AHNAK", "TNFRSF14", "CLEC2D", "GSDMB",
  "CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "IFNG", "CXCL9", "CXCL10",
  "CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "IDO1",
  "TGFB1", "TGFBI", "ACTA2", "COL1A1", "COL1A2", "COL3A1", "FAP", "PDGFRA", "PDGFRB",
  "VIM", "FN1", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "EPCAM", "KRT5", "KRT14", "KRT20",
  "PECAM1", "VWF", "KDR", "CLDN5", "LYZ", "AIF1", "CD68", "C1QA", "C1QB", "C1QC",
  "MS4A1", "CD79A", "JCHAIN", "MKI67", "TOP2A"
))

load(file.path(pkg_data, "cds.RData"))

counts <- get("counts", envir = slot(cds, "assayData"))
pheno <- slot(slot(cds, "phenoData"), "data")
feature <- slot(slot(cds, "featureData"), "data")

feature$feature_id <- rownames(feature)
pheno$sample_id <- rownames(pheno)

symbols <- as.character(feature$symbol)
keep <- which(symbols %in% genes_needed)
counts_sub <- counts[keep, , drop = FALSE]
symbols_sub <- symbols[keep]

if (any(duplicated(symbols_sub))) {
  unique_symbols <- unique(symbols_sub)
  collapsed <- matrix(0, nrow = length(unique_symbols), ncol = ncol(counts_sub))
  rownames(collapsed) <- unique_symbols
  colnames(collapsed) <- colnames(counts_sub)
  for (g in unique_symbols) {
    idx <- which(symbols_sub == g)
    collapsed[g, ] <- colSums(counts_sub[idx, , drop = FALSE], na.rm = TRUE)
  }
  counts_sub <- collapsed
} else {
  rownames(counts_sub) <- symbols_sub
}

write.table(pheno, file.path(out_dir, "IMvigor210_clinical.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(feature, file.path(out_dir, "IMvigor210_feature_annotation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(counts_sub, file.path(out_dir, "IMvigor210_selected_gene_counts.tsv"), sep = "\t", quote = FALSE, col.names = NA)

manifest <- data.frame(
  item = c("clinical", "feature_annotation", "selected_gene_counts"),
  file = c(
    file.path(out_dir, "IMvigor210_clinical.tsv"),
    file.path(out_dir, "IMvigor210_feature_annotation.tsv"),
    file.path(out_dir, "IMvigor210_selected_gene_counts.tsv")
  ),
  stringsAsFactors = FALSE
)
write.table(manifest, file.path(out_dir, "IMvigor210_extraction_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("clinical", dim(pheno), "\n")
cat("selected_gene_counts", dim(counts_sub), "\n")
cat("genes", paste(rownames(counts_sub), collapse = ","), "\n")
