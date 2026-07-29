#!/usr/bin/env Rscript
# Evaluate fixed TCGA LASSO-Cox model with KM, ROC, C-index, and multivariate Cox.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE) else file.path(getwd(), "scripts", "evaluate_tcga_fixed_model.R")
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

processed_dir <- file.path(project_root, "data", "processed")
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
model_dir <- file.path(project_root, "results", "model")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

library(survival)
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

stage_group <- function(x) {
  xu <- toupper(as.character(x))
  out <- ifelse(grepl("^STAGE III|^STAGE IV", xu), "Stage III-IV",
                ifelse(grepl("^STAGE I($|[ABC ]|\\s)|^STAGE II($|[ABC ]|\\s)", xu), "Stage I-II", NA))
  out
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
  cuts <- sort(unique(s), decreasing = TRUE)
  cuts <- c(Inf, cuts, -Inf)
  sens <- spec <- numeric(length(cuts))
  for (i in seq_along(cuts)) {
    pred <- s >= cuts[i]
    sens[i] <- sum(pred & y == 1) / sum(y == 1)
    spec[i] <- sum(!pred & y == 0) / sum(y == 0)
  }
  data.frame(horizon_days = horizon, fpr = 1 - spec, tpr = sens, auc = as.numeric(auc))
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

scores <- read.delim(file.path(table_dir, "tcga_lasso_risk_scores.tsv"), check.names = FALSE)
model <- read.delim(file.path(model_dir, "final_lasso_cox_model.tsv"), check.names = FALSE)
scores$stage_group <- stage_group(scores$ajcc_pathologic_stage)
scores$risk_group_factor <- factor(scores$lasso_risk_group, levels = c("Low", "High"))
scores$negative_lasso_risk_score <- -scores$lasso_risk_score
scores$gender <- factor(scores$gender)
scores$stage_group <- factor(scores$stage_group, levels = c("Stage I-II", "Stage III-IV"))

km_plot <- function(dat, title) {
  fit <- survfit(Surv(os_time_days, os_event) ~ risk_group_factor, data = dat)
  plot(fit, col = palette_main[1:2], lwd = 2, xlab = "Days", ylab = "Overall survival probability",
       main = title, mark.time = TRUE)
  grid(col = "gray88")
  legend("bottomleft", legend = c("Low risk", "High risk"), col = palette_main[1:2], lwd = 2, bty = "n")
  lr <- survdiff(Surv(os_time_days, os_event) ~ risk_group_factor, data = dat)
  p <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  mtext(sprintf("Log-rank p = %.3g", p), side = 3, line = 0.2, adj = 1)
}

save_png("tcga_fixed_model_km_train_test_all.png", width = 2400, height = 900, {
  par(mfrow = c(1, 3), mar = c(5, 5, 4, 1))
  km_plot(scores[scores$split == "train", ], "TCGA train")
  km_plot(scores[scores$split == "test", ], "TCGA test")
  km_plot(scores, "TCGA all")
})

save_png("tcga_fixed_model_coefficients.png", {
  vals <- model$coefficient
  names(vals) <- model$gene
  vals <- sort(vals)
  cols <- ifelse(vals > 0, palette_main[2], palette_main[1])
  par(mar = c(5, 9, 4, 1))
  pad <- diff(range(vals)) * 0.2
  bp <- barplot(vals, horiz = TRUE, las = 1, col = cols, border = NA,
                xlim = c(min(vals) - pad, max(vals) + pad),
                main = "Fixed LASSO-Cox coefficients", xlab = "Coefficient")
  abline(v = 0, col = "gray40")
  text(vals, bp, labels = sprintf("%.3f", vals), pos = ifelse(vals > 0, 4, 2), cex = 0.75)
})

save_png("tcga_fixed_model_risk_distribution.png", {
  ord <- order(scores$lasso_risk_score)
  cols <- ifelse(scores$os_event[ord] == 1, palette_main[2], palette_main[1])
  plot(seq_along(ord), scores$lasso_risk_score[ord], pch = 16, col = adjustcolor(cols, alpha.f = 0.8),
       xlab = "Patients ordered by fixed risk score", ylab = "Risk score",
       main = "Fixed TCGA risk-score distribution")
  abline(h = median(scores$lasso_risk_score[scores$split == "train"], na.rm = TRUE), lty = 2, col = "gray40")
  legend("topleft", legend = c("Alive/censored", "Dead"), col = palette_main[1:2], pch = 16, bty = "n")
})

roc_rows <- list()
for (cohort in c("train", "test", "all")) {
  dat <- if (cohort == "all") scores else scores[scores$split == cohort, ]
  for (h in c(365, 1095, 1825)) {
    r <- roc_at_time(dat$os_time_days, dat$os_event, dat$lasso_risk_score, h)
    if (!is.null(r)) {
      r$cohort <- cohort
      roc_rows[[length(roc_rows) + 1]] <- r
    }
  }
}
roc_df <- do.call(rbind, roc_rows)
write.table(roc_df, file.path(table_dir, "tcga_fixed_model_time_roc_points.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
auc_df <- unique(roc_df[, c("cohort", "horizon_days", "auc")])
write.table(auc_df, file.path(table_dir, "tcga_fixed_model_time_auc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

save_png("tcga_fixed_model_time_roc.png", width = 2400, height = 900, {
  par(mfrow = c(1, 3), mar = c(5, 5, 4, 1))
  for (cohort in c("train", "test", "all")) {
    plot(c(0, 1), c(0, 1), type = "n", xlab = "1 - specificity", ylab = "Sensitivity",
         main = paste("TCGA", cohort, "time ROC"))
    abline(0, 1, col = "gray75", lty = 2)
    i <- 0
    for (h in c(365, 1095, 1825)) {
      i <- i + 1
      sub <- roc_df[roc_df$cohort == cohort & roc_df$horizon_days == h, ]
      if (nrow(sub) == 0) next
      ord <- order(sub$fpr, sub$tpr)
      lines(sub$fpr[ord], sub$tpr[ord], col = palette_main[i], lwd = 2)
    }
    labs <- auc_df[auc_df$cohort == cohort, ]
    legend("bottomright", legend = sprintf("%dy AUC %.3f", round(labs$horizon_days / 365), labs$auc),
           col = palette_main[seq_len(nrow(labs))], lwd = 2, bty = "n")
  }
})

metric_rows <- list()
cox_rows <- list()
for (cohort in c("train", "test", "all")) {
  dat <- if (cohort == "all") scores else scores[scores$split == cohort, ]
  dat <- dat[!is.na(dat$os_time_days) & !is.na(dat$os_event), ]
  cindex <- concordance(Surv(os_time_days, os_event) ~ negative_lasso_risk_score, data = dat)$concordance
  metric_rows[[length(metric_rows) + 1]] <- data.frame(cohort = cohort, n = nrow(dat), events = sum(dat$os_event), c_index = cindex)
  cox_rows[[length(cox_rows) + 1]] <- cox_table(coxph(Surv(os_time_days, os_event) ~ lasso_risk_score, data = dat), paste0(cohort, "_univariate_score"))
  multi_dat <- dat[complete.cases(dat[, c("lasso_risk_score", "age_at_diagnosis_years", "gender", "stage_group")]), ]
  if (nrow(multi_dat) > 50 && length(unique(multi_dat$stage_group)) > 1) {
    fit <- coxph(Surv(os_time_days, os_event) ~ lasso_risk_score + age_at_diagnosis_years + gender + stage_group, data = multi_dat)
    cox_rows[[length(cox_rows) + 1]] <- cox_table(fit, paste0(cohort, "_multivariate"))
  }
}
metrics <- do.call(rbind, metric_rows)
write.table(metrics, file.path(table_dir, "tcga_fixed_model_cindex.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cox_out <- do.call(rbind, cox_rows)
write.table(cox_out, file.path(table_dir, "tcga_fixed_model_cox_tables.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

html <- c(
  "<!doctype html>",
  "<html><head><meta charset=\"utf-8\"><title>TCGA-BLCA fixed model</title>",
  "<style>body{font-family:Arial,sans-serif;margin:24px;color:#222}h1{font-size:24px}section{margin:24px 0 36px}img{max-width:100%;border:1px solid #ddd}code{background:#f3f3f3;padding:2px 4px}</style>",
  "</head><body>",
  "<h1>TCGA-BLCA fixed LASSO-Cox model</h1>",
  "<section><h2>Train/test/all survival</h2><img src=\"figures/tcga_fixed_model_km_train_test_all.png\"></section>",
  "<section><h2>Coefficients</h2><img src=\"figures/tcga_fixed_model_coefficients.png\"></section>",
  "<section><h2>Time ROC</h2><img src=\"figures/tcga_fixed_model_time_roc.png\"></section>",
  "<section><h2>Risk distribution</h2><img src=\"figures/tcga_fixed_model_risk_distribution.png\"></section>",
  "</body></html>"
)
writeLines(html, file.path(project_root, "results", "tcga_fixed_model_report.html"))
message("Fixed TCGA model report complete: ", file.path(project_root, "results", "tcga_fixed_model_report.html"))
