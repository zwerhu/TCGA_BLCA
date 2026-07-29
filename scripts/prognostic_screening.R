#!/usr/bin/env Rscript
# TCGA-BLCA exploratory prognosis screening using base R plus survival.
#
# This is a fast first-pass screen, not a final publication model. Use the
# outputs to choose directions, then confirm with LASSO-Cox, multivariate Cox,
# and an external validation cohort.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE) else file.path(getwd(), "scripts", "prognostic_screening.R")
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

processed_dir <- file.path(project_root, "data", "processed")
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

library(survival)
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

survival_sheet <- read.delim(file.path(processed_dir, "blca_survival_expression_sample_sheet.tsv"), check.names = FALSE)
expr <- read.delim(gzfile(file.path(processed_dir, "blca_log2_tpm_gene_symbol_tumor.tsv.gz")),
                   row.names = 1, check.names = FALSE)

survival_sheet <- survival_sheet[!is.na(survival_sheet$os_time_days) & !is.na(survival_sheet$os_event), ]
survival_sheet <- survival_sheet[survival_sheet$matrix_column %in% colnames(expr), ]
survival_sheet <- survival_sheet[match(intersect(survival_sheet$matrix_column, colnames(expr)), survival_sheet$matrix_column), ]
expr <- as.matrix(expr[, survival_sheet$matrix_column, drop = FALSE])
mode(expr) <- "numeric"

gene_mean <- rowMeans(expr, na.rm = TRUE)
gene_sd <- apply(expr, 1, sd, na.rm = TRUE)
gene_mad <- apply(expr, 1, mad, na.rm = TRUE)
eligible <- which(is.finite(gene_mean) & is.finite(gene_sd) & gene_mean > 0.2 & gene_sd > 0.1)
eligible <- eligible[order(gene_mad[eligible], decreasing = TRUE)]
eligible <- head(eligible, min(12000, length(eligible)))

y_time <- survival_sheet$os_time_days
y_event <- as.integer(survival_sheet$os_event)

cox_one <- function(gene_index) {
  x <- as.numeric(expr[gene_index, ])
  if (sd(x, na.rm = TRUE) <= 0 || any(!is.finite(x))) return(NULL)
  fit <- tryCatch(coxph(Surv(y_time, y_event) ~ x, control = coxph.control(iter.max = 30)),
                  error = function(e) NULL, warning = function(w) suppressWarnings(coxph(Surv(y_time, y_event) ~ x)))
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  data.frame(
    gene = rownames(expr)[gene_index],
    beta = unname(s$coefficients[1, "coef"]),
    HR = unname(s$coefficients[1, "exp(coef)"]),
    se = unname(s$coefficients[1, "se(coef)"]),
    z = unname(s$coefficients[1, "z"]),
    pvalue = unname(s$coefficients[1, "Pr(>|z|)"]),
    lower95 = unname(s$conf.int[1, "lower .95"]),
    upper95 = unname(s$conf.int[1, "upper .95"]),
    mean_log2_tpm = gene_mean[gene_index],
    sd_log2_tpm = gene_sd[gene_index],
    stringsAsFactors = FALSE
  )
}

message("Testing genes: ", length(eligible))
pieces <- vector("list", length(eligible))
for (i in seq_along(eligible)) {
  pieces[[i]] <- cox_one(eligible[i])
  if (i %% 1000 == 0) message("  finished ", i, "/", length(eligible))
}
cox_res <- do.call(rbind, pieces)
cox_res$FDR <- p.adjust(cox_res$pvalue, method = "BH")
cox_res$log2HR <- cox_res$beta / log(2)
cox_res <- cox_res[order(cox_res$pvalue), ]
write.table(cox_res, file.path(table_dir, "univariate_cox_expression_screen.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(head(cox_res, 100), file.path(table_dir, "univariate_cox_top100.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Cox volcano plot
save_png("prognosis_univariate_cox_volcano.png", {
  plot(cox_res$log2HR, -log10(pmax(cox_res$FDR, .Machine$double.xmin)),
       pch = 16, cex = 0.45, col = adjustcolor("gray35", alpha.f = 0.45),
       xlab = "log2(HR)", ylab = "-log10(FDR)", main = "Univariate Cox expression screen")
  sig <- cox_res$FDR < 0.05
  points(cox_res$log2HR[sig], -log10(pmax(cox_res$FDR[sig], .Machine$double.xmin)),
         pch = 16, cex = 0.55, col = adjustcolor(palette_main[2], alpha.f = 0.75))
  abline(v = 0, col = "gray70")
  abline(h = -log10(0.05), col = palette_main[1], lty = 2)
  top <- head(cox_res[is.finite(cox_res$log2HR), ], 10)
  text(top$log2HR, -log10(pmax(top$FDR, .Machine$double.xmin)), labels = top$gene,
       cex = 0.7, pos = ifelse(top$log2HR >= 0, 4, 2), offset = 0.25)
})

km_gene_plot <- function(gene, panel_title = gene) {
  x <- as.numeric(expr[gene, ])
  group <- ifelse(x >= median(x, na.rm = TRUE), "High", "Low")
  dat <- data.frame(time = y_time, event = y_event, group = factor(group, levels = c("Low", "High")))
  fit <- survfit(Surv(time, event) ~ group, data = dat)
  plot(fit, col = palette_main[1:2], lwd = 2, xlab = "Days", ylab = "Overall survival",
       main = panel_title, mark.time = TRUE)
  grid(col = "gray88")
  legend("bottomleft", legend = c("Low", "High"), col = palette_main[1:2], lwd = 2, bty = "n", cex = 0.8)
  lr <- survdiff(Surv(time, event) ~ group, data = dat)
  p <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  mtext(sprintf("p = %.3g", p), side = 3, line = -1, adj = 1, cex = 0.8)
}

top_genes <- head(cox_res$gene[is.finite(cox_res$pvalue)], 4)
save_png("prognosis_top_gene_km.png", width = 2200, height = 1700, {
  par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
  for (gene in top_genes) km_gene_plot(gene)
})

# Exploratory weighted score from top univariate genes.
score_genes <- head(cox_res$gene[cox_res$FDR < 0.05 & is.finite(cox_res$beta)], 8)
if (length(score_genes) < 4) score_genes <- head(cox_res$gene[is.finite(cox_res$beta)], 8)
score_beta <- cox_res$beta[match(score_genes, cox_res$gene)]
z_expr <- scale(t(expr[score_genes, , drop = FALSE]))
risk_score <- as.vector(z_expr %*% score_beta)
risk_group <- ifelse(risk_score >= median(risk_score, na.rm = TRUE), "High", "Low")
risk_table <- data.frame(
  matrix_column = survival_sheet$matrix_column,
  sample_barcode = survival_sheet$sample_barcode,
  patient_barcode = survival_sheet$patient_barcode,
  os_time_days = y_time,
  os_event = y_event,
  risk_score = risk_score,
  risk_group = risk_group,
  stringsAsFactors = FALSE
)
write.table(risk_table, file.path(table_dir, "exploratory_univariate_risk_score.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(gene = score_genes, beta = score_beta),
            file.path(table_dir, "exploratory_risk_score_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

save_png("prognosis_exploratory_risk_score_km.png", {
  dat <- data.frame(time = y_time, event = y_event, group = factor(risk_group, levels = c("Low", "High")))
  fit <- survfit(Surv(time, event) ~ group, data = dat)
  plot(fit, col = palette_main[1:2], lwd = 2, xlab = "Days", ylab = "Overall survival probability",
       main = "Exploratory expression risk score", mark.time = TRUE)
  grid(col = "gray88")
  legend("bottomleft", legend = c("Low risk", "High risk"), col = palette_main[1:2], lwd = 2, bty = "n")
  lr <- survdiff(Surv(time, event) ~ group, data = dat)
  p <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  mtext(sprintf("Log-rank p = %.3g", p), side = 3, line = 0.2, adj = 1)
})

save_png("prognosis_risk_score_distribution.png", {
  ord <- order(risk_score)
  cols <- ifelse(y_event[ord] == 1, palette_main[2], palette_main[1])
  plot(seq_along(ord), risk_score[ord], pch = 16, col = adjustcolor(cols, alpha.f = 0.8),
       xlab = "Patients ordered by risk score", ylab = "Risk score",
       main = "Exploratory risk-score distribution")
  abline(h = median(risk_score, na.rm = TRUE), col = "gray40", lty = 2)
  legend("topleft", legend = c("Alive/censored", "Dead"), col = palette_main[1:2], pch = 16, bty = "n")
})

html <- c(
  "<!doctype html>",
  "<html><head><meta charset=\"utf-8\"><title>TCGA-BLCA prognosis screening</title>",
  "<style>body{font-family:Arial,sans-serif;margin:24px;color:#222}h1{font-size:24px}section{margin:24px 0 36px}img{max-width:100%;border:1px solid #ddd}code{background:#f3f3f3;padding:2px 4px}</style>",
  "</head><body>",
  "<h1>TCGA-BLCA prognosis screening</h1>",
  "<p>Exploratory univariate Cox screen. Confirm with penalized Cox and external validation before manuscript claims.</p>",
  "<section><h2>Univariate Cox volcano</h2><img src=\"figures/prognosis_univariate_cox_volcano.png\" alt=\"Univariate Cox volcano\"></section>",
  "<section><h2>Top gene KM curves</h2><img src=\"figures/prognosis_top_gene_km.png\" alt=\"Top gene KM curves\"></section>",
  "<section><h2>Exploratory risk score KM</h2><img src=\"figures/prognosis_exploratory_risk_score_km.png\" alt=\"Exploratory risk score KM\"></section>",
  "<section><h2>Risk score distribution</h2><img src=\"figures/prognosis_risk_score_distribution.png\" alt=\"Risk score distribution\"></section>",
  "</body></html>"
)
writeLines(html, file.path(project_root, "results", "prognosis_report.html"))

message("Top Cox gene: ", cox_res$gene[1], " p=", signif(cox_res$pvalue[1], 3), " FDR=", signif(cox_res$FDR[1], 3))
message("Risk genes: ", paste(score_genes, collapse = ", "))
message("Prognosis report complete: ", file.path(project_root, "results", "prognosis_report.html"))
