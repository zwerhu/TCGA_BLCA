#!/usr/bin/env Rscript
# Validate fixed TCGA model in GSE13507 and GSE31684.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE) else file.path(getwd(), "scripts", "validate_geo_fixed_model.R")
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

external_dir <- file.path(project_root, "data", "external", "processed")
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

library(survival)
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

roc_at_time <- function(time, event, score, horizon) {
  case <- time <= horizon & event == 1
  control <- time > horizon
  keep <- case | control
  y <- ifelse(case[keep], 1, 0)
  s <- score[keep]
  if (sum(y == 1) < 3 || sum(y == 0) < 3) return(NULL)
  ranks <- rank(s)
  auc <- (sum(ranks[y == 1]) - sum(y == 1) * (sum(y == 1) + 1) / 2) / (sum(y == 1) * sum(y == 0))
  cuts <- c(Inf, sort(unique(s), decreasing = TRUE), -Inf)
  sens <- spec <- numeric(length(cuts))
  for (i in seq_along(cuts)) {
    pred <- s >= cuts[i]
    sens[i] <- sum(pred & y == 1) / sum(y == 1)
    spec[i] <- sum(!pred & y == 0) / sum(y == 0)
  }
  data.frame(horizon_days = horizon, fpr = 1 - spec, tpr = sens, auc = as.numeric(auc))
}

stage_group <- function(cohort, stage) {
  s <- toupper(as.character(stage))
  if (cohort == "GSE13507") {
    ifelse(grepl("^TA|^T1", s), "Ta/T1", ifelse(grepl("^T2|^T3|^T4", s), "T2-T4", NA))
  } else {
    ifelse(grepl("PT1|PT2", s), "pT1/pT2", ifelse(grepl("PT3|PT4", s), "pT3/pT4", NA))
  }
}

cox_table <- function(fit, label) {
  s <- summary(fit)
  data.frame(
    cohort = label,
    variable = rownames(s$coefficients),
    beta = s$coefficients[, "coef"],
    HR = s$coefficients[, "exp(coef)"],
    lower95 = s$conf.int[, "lower .95"],
    upper95 = s$conf.int[, "upper .95"],
    pvalue = s$coefficients[, "Pr(>|z|)"],
    row.names = NULL
  )
}

read_cohort <- function(cohort) {
  dat <- read.delim(file.path(external_dir, paste0(cohort, "_tcga_fixed_risk_scores.tsv")), check.names = FALSE)
  dat$cohort_median_group <- ifelse(dat$tcga_fixed_risk_score >= median(dat$tcga_fixed_risk_score, na.rm = TRUE), "High", "Low")
  dat$cohort_median_group <- factor(dat$cohort_median_group, levels = c("Low", "High"))
  dat$stage_group <- factor(stage_group(cohort, dat$stage))
  dat$gender <- factor(dat$gender)
  dat
}

cohorts <- c("GSE13507", "GSE31684")
all_metrics <- list()
all_cox <- list()
all_roc <- list()

save_png("geo_validation_km.png", width = 1800, height = 850, {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  for (cohort in cohorts) {
    dat <- read_cohort(cohort)
    fit <- survfit(Surv(os_time_days, os_event) ~ cohort_median_group, data = dat)
    plot(fit, col = palette_main[1:2], lwd = 2, xlab = "Days", ylab = "Overall survival probability",
         main = cohort, mark.time = TRUE)
    grid(col = "gray88")
    legend("bottomleft", legend = c("Low risk", "High risk"), col = palette_main[1:2], lwd = 2, bty = "n")
    lr <- survdiff(Surv(os_time_days, os_event) ~ cohort_median_group, data = dat)
    p <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
    mtext(sprintf("Log-rank p = %.3g", p), side = 3, line = 0.2, adj = 1)
  }
})

save_png("geo_validation_risk_distribution.png", width = 1800, height = 850, {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  for (cohort in cohorts) {
    dat <- read_cohort(cohort)
    ord <- order(dat$tcga_fixed_risk_score)
    cols <- ifelse(dat$os_event[ord] == 1, palette_main[2], palette_main[1])
    plot(seq_along(ord), dat$tcga_fixed_risk_score[ord], pch = 16, col = adjustcolor(cols, alpha.f = 0.8),
         xlab = "Patients ordered by score", ylab = "TCGA fixed risk score", main = cohort)
    abline(h = median(dat$tcga_fixed_risk_score, na.rm = TRUE), lty = 2, col = "gray40")
    legend("topleft", legend = c("Alive/censored", "Dead"), col = palette_main[1:2], pch = 16, bty = "n")
  }
})

for (cohort in cohorts) {
  dat <- read_cohort(cohort)
  dat <- dat[!is.na(dat$os_time_days) & !is.na(dat$os_event), ]
  dat$negative_score <- -dat$tcga_fixed_risk_score
  cindex <- concordance(Surv(os_time_days, os_event) ~ negative_score, data = dat)$concordance
  lr <- survdiff(Surv(os_time_days, os_event) ~ cohort_median_group, data = dat)
  lr_p <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  all_metrics[[length(all_metrics) + 1]] <- data.frame(
    cohort = cohort,
    n = nrow(dat),
    events = sum(dat$os_event),
    c_index = cindex,
    km_logrank_p = lr_p,
    available_model_genes = unique(dat$n_model_genes_available)[1],
    missing_model_genes = unique(dat$missing_model_genes)[1]
  )
  all_cox[[length(all_cox) + 1]] <- cox_table(coxph(Surv(os_time_days, os_event) ~ tcga_fixed_risk_score, data = dat), paste0(cohort, "_univariate_score"))
  multi_dat <- dat[complete.cases(dat[, c("tcga_fixed_risk_score", "age", "gender", "stage_group")]), ]
  if (nrow(multi_dat) > 30 && length(unique(multi_dat$stage_group)) > 1) {
    fit <- coxph(Surv(os_time_days, os_event) ~ tcga_fixed_risk_score + age + gender + stage_group, data = multi_dat)
    all_cox[[length(all_cox) + 1]] <- cox_table(fit, paste0(cohort, "_multivariate"))
  }
  for (h in c(365, 1095, 1825)) {
    r <- roc_at_time(dat$os_time_days, dat$os_event, dat$tcga_fixed_risk_score, h)
    if (!is.null(r)) {
      r$cohort <- cohort
      all_roc[[length(all_roc) + 1]] <- r
    }
  }
}

metrics <- do.call(rbind, all_metrics)
cox_out <- do.call(rbind, all_cox)
roc_df <- do.call(rbind, all_roc)
auc_df <- unique(roc_df[, c("cohort", "horizon_days", "auc")])
write.table(metrics, file.path(table_dir, "geo_fixed_model_validation_metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cox_out, file.path(table_dir, "geo_fixed_model_cox_tables.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(roc_df, file.path(table_dir, "geo_fixed_model_roc_points.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(auc_df, file.path(table_dir, "geo_fixed_model_time_auc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

save_png("geo_validation_time_roc.png", width = 1800, height = 850, {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  for (cohort in cohorts) {
    plot(c(0, 1), c(0, 1), type = "n", xlab = "1 - specificity", ylab = "Sensitivity",
         main = paste(cohort, "time ROC"))
    abline(0, 1, col = "gray75", lty = 2)
    labs <- auc_df[auc_df$cohort == cohort, ]
    i <- 0
    for (h in c(365, 1095, 1825)) {
      i <- i + 1
      sub <- roc_df[roc_df$cohort == cohort & roc_df$horizon_days == h, ]
      if (nrow(sub) == 0) next
      ord <- order(sub$fpr, sub$tpr)
      lines(sub$fpr[ord], sub$tpr[ord], col = palette_main[i], lwd = 2)
    }
    legend("bottomright", legend = sprintf("%dy AUC %.3f", round(labs$horizon_days / 365), labs$auc),
           col = palette_main[seq_len(nrow(labs))], lwd = 2, bty = "n")
  }
})

html <- c(
  "<!doctype html>",
  "<html><head><meta charset=\"utf-8\"><title>GEO validation</title>",
  "<style>body{font-family:Arial,sans-serif;margin:24px;color:#222}h1{font-size:24px}section{margin:24px 0 36px}img{max-width:100%;border:1px solid #ddd}code{background:#f3f3f3;padding:2px 4px}</style>",
  "</head><body>",
  "<h1>GEO validation of fixed TCGA model</h1>",
  "<section><h2>Kaplan-Meier validation</h2><img src=\"figures/geo_validation_km.png\"></section>",
  "<section><h2>Time ROC</h2><img src=\"figures/geo_validation_time_roc.png\"></section>",
  "<section><h2>Risk-score distribution</h2><img src=\"figures/geo_validation_risk_distribution.png\"></section>",
  "</body></html>"
)
writeLines(html, file.path(project_root, "results", "geo_validation_report.html"))
message("GEO validation report complete: ", file.path(project_root, "results", "geo_validation_report.html"))
