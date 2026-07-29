#!/usr/bin/env Rscript
# Simplified GDSC2 ridge drug-response prediction for TCGA-BLCA.
#
# This is not a replacement for oncoPredict; it uses the same downloaded GDSC2
# training matrices but implements a compact ridge model with base R. Lower
# predicted response is interpreted as greater predicted sensitivity.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE) else file.path(getwd(), "scripts", "drug_sensitivity_gdsc2_ridge.R")
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

drug_dir <- file.path(project_root, "data", "external", "drug", "GDSC2")
processed_dir <- file.path(project_root, "data", "processed")
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
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

bh <- function(p) {
  p.adjust(p, method = "BH")
}

compare_groups <- function(values, groups) {
  high <- values[groups == "High"]
  low <- values[groups == "Low"]
  diff <- mean(high, na.rm = TRUE) - mean(low, na.rm = TRUE)
  se <- sqrt(var(high, na.rm = TRUE) / sum(!is.na(high)) + var(low, na.rm = TRUE) / sum(!is.na(low)))
  z <- diff / max(se, 1e-12)
  c(high_mean = mean(high, na.rm = TRUE), low_mean = mean(low, na.rm = TRUE),
    diff_high_minus_low = diff, z = z, pvalue = 2 * pnorm(-abs(z)))
}

message("Reading GDSC2 training data")
gdsc_expr <- readRDS(file.path(drug_dir, "GDSC2_Expr (RMA Normalized and Log Transformed).rds"))
gdsc_res <- readRDS(file.path(drug_dir, "GDSC2_Res.rds"))

message("Reading TCGA expression")
tcga_expr <- read.delim(gzfile(file.path(processed_dir, "blca_drug_prediction_expression_input.tsv.gz")),
                        row.names = 1, check.names = FALSE)
tcga_risk <- read.delim(file.path(table_dir, "tcga_lasso_risk_scores.tsv"), check.names = FALSE)
tcga_expr <- as.matrix(tcga_expr[, tcga_risk$matrix_column, drop = FALSE])
mode(tcga_expr) <- "numeric"

common <- intersect(rownames(gdsc_expr), rownames(tcga_expr))
message("Common genes: ", length(common))
gdsc_x <- t(gdsc_expr[common, , drop = FALSE])
tcga_x <- t(tcga_expr[common, , drop = FALSE])
gene_mean <- colMeans(gdsc_x, na.rm = TRUE)
gene_sd <- apply(gdsc_x, 2, sd, na.rm = TRUE)
keep <- is.finite(gene_mean) & is.finite(gene_sd) & gene_sd > 1e-8
gdsc_x <- scale(gdsc_x[, keep, drop = FALSE], center = gene_mean[keep], scale = gene_sd[keep])
tcga_x <- scale(tcga_x[, keep, drop = FALSE], center = gene_mean[keep], scale = gene_sd[keep])
gdsc_x[!is.finite(gdsc_x)] <- 0
tcga_x[!is.finite(tcga_x)] <- 0
p <- ncol(gdsc_x)

message("Computing kernels")
K_train <- tcrossprod(gdsc_x) / p
K_test <- tcga_x %*% t(gdsc_x) / p
lambda <- 10
pred <- matrix(NA_real_, nrow = nrow(tcga_x), ncol = ncol(gdsc_res),
               dimnames = list(rownames(tcga_x), colnames(gdsc_res)))

for (j in seq_len(ncol(gdsc_res))) {
  y <- gdsc_res[, j]
  obs <- is.finite(y)
  if (sum(obs) < 50) next
  y_mean <- mean(y[obs])
  y_centered <- y[obs] - y_mean
  K_obs <- K_train[obs, obs, drop = FALSE]
  alpha <- solve(K_obs + diag(lambda, sum(obs)), y_centered)
  pred[, j] <- as.vector(K_test[, obs, drop = FALSE] %*% alpha + y_mean)
  if (j %% 50 == 0) message("  drugs predicted: ", j, "/", ncol(gdsc_res))
}

write.table(t(pred), gzfile(file.path(drug_dir, "tcga_gdsc2_ridge_predicted_response.tsv.gz")),
            sep = "\t", quote = FALSE, col.names = NA)

groups <- tcga_risk$lasso_risk_group[match(rownames(pred), tcga_risk$matrix_column)]
rows <- list()
for (drug in colnames(pred)) {
  stats <- compare_groups(pred[, drug], groups)
  rows[[length(rows) + 1]] <- data.frame(drug = drug, t(stats), row.names = NULL)
}
diff <- do.call(rbind, rows)
diff$FDR <- bh(diff$pvalue)
diff <- diff[order(diff$pvalue), ]
diff$interpretation <- ifelse(diff$diff_high_minus_low < 0, "High_risk_more_sensitive", "Low_risk_more_sensitive")
write.table(diff, file.path(table_dir, "tcga_fixed_risk_gdsc2_ridge_drug_diff.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

save_png("tcga_fixed_risk_gdsc2_drug_sensitivity.png", width = 2100, height = 1300, {
  top_high <- head(diff[diff$diff_high_minus_low < 0, ], 15)
  top_low <- head(diff[diff$diff_high_minus_low > 0, ], 15)
  dat <- rbind(top_high, top_low)
  dat <- dat[order(dat$diff_high_minus_low), ]
  vals <- dat$diff_high_minus_low
  names(vals) <- gsub("_\\d+$", "", dat$drug)
  cols <- ifelse(vals < 0, palette_main[2], palette_main[1])
  par(mar = c(5, 12, 4, 1))
  pad <- diff(range(vals, na.rm = TRUE)) * 0.2
  bp <- barplot(vals, horiz = TRUE, las = 1, col = cols, border = NA,
                xlim = c(min(vals, na.rm = TRUE) - pad, max(vals, na.rm = TRUE) + pad),
                main = "GDSC2 ridge-predicted drug response by fixed risk group",
                xlab = "Predicted response difference, high - low")
  abline(v = 0, col = "gray40")
  legend("bottomright", legend = c("High-risk lower response", "Low-risk lower response"),
         fill = c(palette_main[2], palette_main[1]), bty = "n")
})

html <- c(
  "<!doctype html>",
  "<html><head><meta charset=\"utf-8\"><title>GDSC2 drug prediction</title>",
  "<style>body{font-family:Arial,sans-serif;margin:24px;color:#222}h1{font-size:24px}section{margin:24px 0 36px}img{max-width:100%;border:1px solid #ddd}code{background:#f3f3f3;padding:2px 4px}</style>",
  "</head><body>",
  "<h1>GDSC2 ridge-predicted drug sensitivity</h1>",
  "<p>Lower predicted response indicates greater predicted sensitivity. This is a simplified ridge predictor; rerun with oncoPredict for publication-grade drug analysis.</p>",
  "<section><h2>Top predicted differences</h2><img src=\"figures/tcga_fixed_risk_gdsc2_drug_sensitivity.png\"></section>",
  "</body></html>"
)
writeLines(html, file.path(project_root, "results", "drug_sensitivity_report.html"))
message("Drug prediction report complete: ", file.path(project_root, "results", "drug_sensitivity_report.html"))
