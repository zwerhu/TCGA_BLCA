#!/usr/bin/env Rscript
# Evaluate GEO-compatible reduced model in TCGA, GSE13507, and GSE31684.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE) else file.path(getwd(), "scripts", "evaluate_reduced_geo_compatible_model.R")
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

ext_dir <- file.path(project_root, "data", "external", "processed")
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

tcga_stage_group <- function(stage) {
  xu <- toupper(as.character(stage))
  ifelse(grepl("^STAGE III|^STAGE IV", xu), "Stage III-IV",
         ifelse(grepl("^STAGE I($|[ABC ]|\\s)|^STAGE II($|[ABC ]|\\s)", xu), "Stage I-II", NA))
}

geo_stage_group <- function(cohort, stage) {
  s <- toupper(as.character(stage))
  if (cohort == "GSE13507") {
    ifelse(grepl("^TA|^T1", s), "Ta/T1", ifelse(grepl("^T2|^T3|^T4", s), "T2-T4", NA))
  } else {
    ifelse(grepl("PT1|PT2", s), "pT1/pT2", ifelse(grepl("PT3|PT4", s), "pT3/pT4", NA))
  }
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

read_tcga <- function() {
  dat <- read.delim(file.path(table_dir, "tcga_reduced_geo_compatible_risk_scores.tsv"), check.names = FALSE)
  dat$cohort <- "TCGA"
  dat$sample_id <- dat$matrix_column
  dat$risk_score <- dat$reduced_risk_score
  dat$risk_group <- factor(dat$reduced_risk_group, levels = c("Low", "High"))
  dat$age <- dat$age_at_diagnosis_years
  dat$stage_group <- factor(tcga_stage_group(dat$ajcc_pathologic_stage))
  dat$gender <- factor(dat$gender)
  dat
}

read_geo <- function(cohort) {
  dat <- read.delim(file.path(ext_dir, paste0(cohort, "_reduced_geo_compatible_risk_scores.tsv")), check.names = FALSE)
  dat$risk_score <- dat$reduced_risk_score
  dat$risk_group <- factor(dat$reduced_risk_group, levels = c("Low", "High"))
  dat$stage_group <- factor(geo_stage_group(cohort, dat$stage))
  dat$gender <- factor(dat$gender)
  dat
}

km_plot <- function(dat, title) {
  fit <- survfit(Surv(os_time_days, os_event) ~ risk_group, data = dat)
  plot(fit, col = palette_main[1:2], lwd = 2, xlab = "Days", ylab = "Overall survival probability",
       main = title, mark.time = TRUE)
  grid(col = "gray88")
  legend("bottomleft", legend = c("Low risk", "High risk"), col = palette_main[1:2], lwd = 2, bty = "n")
  lr <- survdiff(Surv(os_time_days, os_event) ~ risk_group, data = dat)
  p <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  mtext(sprintf("Log-rank p = %.3g", p), side = 3, line = 0.2, adj = 1)
}

tcga <- read_tcga()
gse13507 <- read_geo("GSE13507")
gse31684 <- read_geo("GSE31684")
model <- read.delim(file.path(model_dir, "reduced_geo_compatible_model.tsv"), check.names = FALSE)

save_png("reduced_model_km_tcga_geo.png", width = 2400, height = 1700, {
  par(mfrow = c(2, 3), mar = c(5, 5, 4, 1))
  km_plot(tcga[tcga$split == "train", ], "TCGA train")
  km_plot(tcga[tcga$split == "test", ], "TCGA test")
  km_plot(tcga, "TCGA all")
  km_plot(gse13507, "GSE13507")
  km_plot(gse31684, "GSE31684")
  plot.new()
})

save_png("reduced_model_coefficients.png", {
  vals <- model$coefficient
  names(vals) <- model$gene
  vals <- sort(vals)
  cols <- ifelse(vals > 0, palette_main[2], palette_main[1])
  par(mar = c(5, 9, 4, 1))
  pad <- diff(range(vals)) * 0.25
  bp <- barplot(vals, horiz = TRUE, las = 1, col = cols, border = NA,
                xlim = c(min(vals) - pad, max(vals) + pad),
                main = "GEO-compatible reduced model coefficients",
                xlab = "Coefficient")
  abline(v = 0, col = "gray40")
  text(vals, bp, labels = sprintf("%.3f", vals), pos = ifelse(vals > 0, 4, 2), cex = 0.8)
})

cohorts <- list(
  "TCGA_train" = tcga[tcga$split == "train", ],
  "TCGA_test" = tcga[tcga$split == "test", ],
  "TCGA_all" = tcga,
  "GSE13507" = gse13507,
  "GSE31684" = gse31684
)

metric_rows <- list()
cox_rows <- list()
roc_rows <- list()
for (name in names(cohorts)) {
  dat <- cohorts[[name]]
  dat <- dat[!is.na(dat$os_time_days) & !is.na(dat$os_event), ]
  dat$negative_score <- -dat$risk_score
  cidx <- concordance(Surv(os_time_days, os_event) ~ negative_score, data = dat)$concordance
  lr <- survdiff(Surv(os_time_days, os_event) ~ risk_group, data = dat)
  lr_p <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  metric_rows[[length(metric_rows) + 1]] <- data.frame(cohort = name, n = nrow(dat), events = sum(dat$os_event), c_index = cidx, km_logrank_p = lr_p)
  cox_rows[[length(cox_rows) + 1]] <- cox_table(coxph(Surv(os_time_days, os_event) ~ risk_score, data = dat), paste0(name, "_univariate_score"))
  multi_dat <- dat[complete.cases(dat[, c("risk_score", "age", "gender", "stage_group")]), ]
  if (nrow(multi_dat) > 30 && length(unique(multi_dat$stage_group)) > 1) {
    fit <- coxph(Surv(os_time_days, os_event) ~ risk_score + age + gender + stage_group, data = multi_dat)
    cox_rows[[length(cox_rows) + 1]] <- cox_table(fit, paste0(name, "_multivariate"))
  }
  for (h in c(365, 1095, 1825)) {
    r <- roc_at_time(dat$os_time_days, dat$os_event, dat$risk_score, h)
    if (!is.null(r)) {
      r$cohort <- name
      roc_rows[[length(roc_rows) + 1]] <- r
    }
  }
}
metrics <- do.call(rbind, metric_rows)
cox_out <- do.call(rbind, cox_rows)
roc_df <- do.call(rbind, roc_rows)
auc_df <- unique(roc_df[, c("cohort", "horizon_days", "auc")])
write.table(metrics, file.path(table_dir, "reduced_geo_compatible_validation_metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cox_out, file.path(table_dir, "reduced_geo_compatible_cox_tables.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(roc_df, file.path(table_dir, "reduced_geo_compatible_roc_points.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(auc_df, file.path(table_dir, "reduced_geo_compatible_time_auc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

save_png("reduced_model_time_roc.png", width = 2400, height = 1700, {
  par(mfrow = c(2, 3), mar = c(5, 5, 4, 1))
  for (name in names(cohorts)) {
    plot(c(0, 1), c(0, 1), type = "n", xlab = "1 - specificity", ylab = "Sensitivity",
         main = paste(name, "time ROC"))
    abline(0, 1, col = "gray75", lty = 2)
    labs <- auc_df[auc_df$cohort == name, ]
    i <- 0
    for (h in c(365, 1095, 1825)) {
      i <- i + 1
      sub <- roc_df[roc_df$cohort == name & roc_df$horizon_days == h, ]
      if (nrow(sub) == 0) next
      ord <- order(sub$fpr, sub$tpr)
      lines(sub$fpr[ord], sub$tpr[ord], col = palette_main[i], lwd = 2)
    }
    if (nrow(labs) > 0) {
      legend("bottomright", legend = sprintf("%dy AUC %.3f", round(labs$horizon_days / 365), labs$auc),
             col = palette_main[seq_len(nrow(labs))], lwd = 2, bty = "n", cex = 0.8)
    }
  }
  plot.new()
})

save_png("reduced_model_risk_distribution_geo.png", width = 2400, height = 900, {
  par(mfrow = c(1, 3), mar = c(5, 5, 4, 1))
  for (name in c("TCGA_all", "GSE13507", "GSE31684")) {
    dat <- cohorts[[name]]
    ord <- order(dat$risk_score)
    cols <- ifelse(dat$os_event[ord] == 1, palette_main[2], palette_main[1])
    plot(seq_along(ord), dat$risk_score[ord], pch = 16, col = adjustcolor(cols, alpha.f = 0.8),
         xlab = "Patients ordered by score", ylab = "Reduced risk score", main = name)
    abline(h = median(dat$risk_score, na.rm = TRUE), lty = 2, col = "gray40")
    legend("topleft", legend = c("Alive/censored", "Dead"), col = palette_main[1:2], pch = 16, bty = "n", cex = 0.8)
  }
})

html <- c(
  "<!doctype html>",
  "<html><head><meta charset=\"utf-8\"><title>GEO-compatible reduced model</title>",
  "<style>body{font-family:Arial,sans-serif;margin:24px;color:#222}h1{font-size:24px}section{margin:24px 0 36px}img{max-width:100%;border:1px solid #ddd}code{background:#f3f3f3;padding:2px 4px}</style>",
  "</head><body>",
  "<h1>GEO-compatible reduced Cox model</h1>",
  "<p>Five genes shared by TCGA, GSE13507, and GSE31684. Trained in TCGA; GEO validation uses fixed coefficients on cohort-wise z-scored expression.</p>",
  "<section><h2>Kaplan-Meier validation</h2><img src=\"figures/reduced_model_km_tcga_geo.png\"></section>",
  "<section><h2>Coefficients</h2><img src=\"figures/reduced_model_coefficients.png\"></section>",
  "<section><h2>Time ROC</h2><img src=\"figures/reduced_model_time_roc.png\"></section>",
  "<section><h2>Risk distribution</h2><img src=\"figures/reduced_model_risk_distribution_geo.png\"></section>",
  "</body></html>"
)
writeLines(html, file.path(project_root, "results", "reduced_geo_compatible_model_report.html"))
message("Reduced GEO-compatible model report complete: ", file.path(project_root, "results", "reduced_geo_compatible_model_report.html"))
