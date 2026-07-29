#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survival)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- if (!is.na(script_arg)) sub("^--file=", "", script_arg) else "scripts/plot_immunotherapy_response_validation.R"
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
fig_dir <- file.path(root, "results", "figures")
table_dir <- file.path(root, "results", "tables")
processed_dir <- file.path(root, "data", "external", "immunotherapy", "processed")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

scores <- read.delim(
  file.path(processed_dir, "immunotherapy_validation_scores.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN", "nan")
)
model <- read.delim(
  file.path(root, "results", "model", "reduced_geo_compatible_model.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

scores$risk_group <- factor(scores$risk_group, levels = c("Low", "High"))
scores$response_group <- factor(scores$response_group, levels = c("Responder", "Non-responder"))
scores$sex <- factor(scores$sex)
scores$fgfr3_altered <- factor(scores$fgfr3_altered)
scores$immune_phenotype <- factor(scores$immune_phenotype)
scores$tcga_subtype <- factor(scores$tcga_subtype)
scores$cd8_tgfb_balance <- scores$sig_CD8_effector - scores$sig_TGFb_CAF

write.delim <- function(x, file) {
  write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE)
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return(formatC(p, format = "e", digits = 2))
  formatC(p, format = "f", digits = 3)
}

pal <- list(
  low = "#2B6CB0",
  high = "#C2410C",
  responder = "#2F855A",
  nonresponder = "#B7791F",
  gray = "#4A5568",
  light = "#E2E8F0"
)

cohort_order <- c("IMvigor210", "GSE176307")

safe_wilcox <- function(formula, data) {
  out <- tryCatch(wilcox.test(formula, data = data)$p.value, error = function(e) NA_real_)
  out
}

auc_rank <- function(label, score) {
  ok <- !is.na(label) & !is.na(score)
  label <- label[ok]
  score <- score[ok]
  n_pos <- sum(label == 1)
  n_neg <- sum(label == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[label == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

roc_points <- function(label, score) {
  ok <- !is.na(label) & !is.na(score)
  label <- label[ok]
  score <- score[ok]
  if (length(unique(label)) < 2) {
    return(data.frame(fpr = c(0, 1), tpr = c(0, 1), auc = NA_real_))
  }
  ord <- order(score, decreasing = TRUE)
  y <- label[ord]
  tpr <- c(0, cumsum(y == 1) / sum(y == 1))
  fpr <- c(0, cumsum(y == 0) / sum(y == 0))
  auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
  data.frame(fpr = fpr, tpr = tpr, auc = auc)
}

extract_cox <- function(fit, cohort, endpoint, model_type, n, events) {
  ss <- summary(fit)
  cf <- ss$coefficients
  ci <- ss$conf.int
  data.frame(
    cohort = cohort,
    endpoint = endpoint,
    model_type = model_type,
    n = n,
    events = events,
    term = rownames(cf),
    HR = ci[, "exp(coef)"],
    lower95 = ci[, "lower .95"],
    upper95 = ci[, "upper .95"],
    p_value = cf[, "Pr(>|z|)"],
    concordance = ss$concordance[1],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

choose_covariates <- function(d, cohort) {
  covars <- c("risk_score")
  if (sum(!is.na(d$tmb)) >= 20 && sd(d$tmb, na.rm = TRUE) > 0) covars <- c(covars, "tmb")
  if (sum(!is.na(d$ecog)) >= 20 && length(unique(na.omit(d$ecog))) > 1) covars <- c(covars, "ecog")
  if (length(unique(na.omit(d$sex))) > 1) covars <- c(covars, "sex")
  if (cohort == "IMvigor210" && length(unique(na.omit(d$immune_phenotype))) > 1) {
    covars <- c(covars, "immune_phenotype")
  }
  if (cohort == "GSE176307" && length(unique(na.omit(d$fgfr3_altered))) > 1) {
    covars <- c(covars, "fgfr3_altered")
  }
  covars
}

make_cox_tables <- function(scores) {
  rows <- list()
  idx <- 1
  for (cohort in cohort_order) {
    d0 <- scores[scores$cohort == cohort, ]
    for (ep in list(
      list(name = "OS", time = "os_time_months", event = "os_event"),
      list(name = "PFS", time = "pfs_time_months", event = "pfs_event")
    )) {
      d <- d0[!is.na(d0[[ep$time]]) & !is.na(d0[[ep$event]]) & !is.na(d0$risk_score), ]
      if (nrow(d) < 15 || sum(d[[ep$event]] == 1, na.rm = TRUE) < 5) next
      f1 <- as.formula(paste0("Surv(", ep$time, ", ", ep$event, ") ~ risk_score"))
      fit1 <- tryCatch(coxph(f1, data = d), error = function(e) NULL)
      if (!is.null(fit1)) {
        rows[[idx]] <- extract_cox(fit1, cohort, ep$name, "univariate", nrow(d), sum(d[[ep$event]] == 1, na.rm = TRUE))
        idx <- idx + 1
      }

      covars <- choose_covariates(d, cohort)
      f2 <- as.formula(paste0("Surv(", ep$time, ", ", ep$event, ") ~ ", paste(covars, collapse = " + ")))
      mf <- model.frame(f2, data = d, na.action = na.omit)
      if (nrow(mf) >= 25 && sum(model.response(mf)[, 2] == 1, na.rm = TRUE) >= 8) {
        fit2 <- tryCatch(coxph(f2, data = d), error = function(e) NULL)
        if (!is.null(fit2)) {
          rows[[idx]] <- extract_cox(fit2, cohort, ep$name, "multivariate", nrow(mf), sum(model.response(mf)[, 2] == 1, na.rm = TRUE))
          idx <- idx + 1
        }
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  write.delim(out, file.path(table_dir, "immunotherapy_cox_models.tsv"))
  out
}

extract_glm <- function(fit, cohort, model_type, n, nonresponders) {
  ss <- summary(fit)
  cf <- ss$coefficients
  data.frame(
    cohort = cohort,
    model_type = model_type,
    n = n,
    nonresponders = nonresponders,
    term = rownames(cf),
    OR = exp(cf[, "Estimate"]),
    lower95 = exp(cf[, "Estimate"] - 1.96 * cf[, "Std. Error"]),
    upper95 = exp(cf[, "Estimate"] + 1.96 * cf[, "Std. Error"]),
    p_value = cf[, "Pr(>|z|)"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

make_logistic_tables <- function(scores) {
  rows <- list()
  idx <- 1
  for (cohort in cohort_order) {
    d <- scores[scores$cohort == cohort & !is.na(scores$nonresponse_binary) & !is.na(scores$risk_score), ]
    if (nrow(d) < 20 || length(unique(d$nonresponse_binary)) < 2) next
    f1 <- nonresponse_binary ~ risk_score
    fit1 <- tryCatch(glm(f1, data = d, family = binomial()), warning = function(w) suppressWarnings(glm(f1, data = d, family = binomial())), error = function(e) NULL)
    if (!is.null(fit1)) {
      rows[[idx]] <- extract_glm(fit1, cohort, "univariate", nrow(d), sum(d$nonresponse_binary == 1))
      idx <- idx + 1
    }
    covars <- choose_covariates(d, cohort)
    f2 <- as.formula(paste("nonresponse_binary ~", paste(covars, collapse = " + ")))
    mf <- model.frame(f2, data = d, na.action = na.omit)
    if (nrow(mf) >= 25 && length(unique(model.response(mf))) > 1) {
      fit2 <- tryCatch(glm(f2, data = d, family = binomial()), warning = function(w) suppressWarnings(glm(f2, data = d, family = binomial())), error = function(e) NULL)
      if (!is.null(fit2)) {
        rows[[idx]] <- extract_glm(fit2, cohort, "multivariate", nrow(mf), sum(model.response(mf) == 1))
        idx <- idx + 1
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  write.delim(out, file.path(table_dir, "immunotherapy_logistic_response_models.tsv"))
  out
}

plot_response_box <- function() {
  png(file.path(fig_dir, "immunotherapy_response_boxplot.png"), width = 1700, height = 900, res = 170)
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1), oma = c(0, 0, 1, 0))
  for (cohort in cohort_order) {
    d <- scores[scores$cohort == cohort & !is.na(scores$response_group), ]
    boxplot(
      risk_score ~ response_group, data = d, outline = FALSE,
      col = c(grDevices::adjustcolor(pal$responder, 0.35), grDevices::adjustcolor(pal$nonresponder, 0.35)),
      border = c(pal$responder, pal$nonresponder), lwd = 1.5,
      ylab = "Fixed 5-gene risk score", xlab = "", main = cohort
    )
    stripchart(
      risk_score ~ response_group, data = d, vertical = TRUE, method = "jitter",
      add = TRUE, pch = 16, cex = 0.65,
      col = grDevices::adjustcolor("#1A202C", 0.45)
    )
    p <- safe_wilcox(risk_score ~ response_group, d)
    usr <- par("usr")
    text(1.5, usr[4] - 0.05 * diff(usr[3:4]), paste0("Wilcoxon p = ", fmt_p(p)), cex = 0.9)
  }
  mtext("Risk score is higher in ICB non-responders", outer = TRUE, cex = 1.1, font = 2)
  dev.off()
}

plot_response_roc <- function() {
  png(file.path(fig_dir, "immunotherapy_response_roc.png"), width = 1200, height = 1000, res = 170)
  par(mar = c(5, 5, 4, 2))
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - specificity", ylab = "Sensitivity", main = "Risk score predicts ICB non-response")
  abline(0, 1, col = "gray70", lty = 2)
  cols <- c(IMvigor210 = pal$high, GSE176307 = pal$low)
  labels <- c()
  for (cohort in cohort_order) {
    d <- scores[scores$cohort == cohort & !is.na(scores$nonresponse_binary), ]
    r <- roc_points(d$nonresponse_binary, d$risk_score)
    lines(r$fpr, r$tpr, lwd = 2.4, col = cols[[cohort]])
    labels <- c(labels, paste0(cohort, " AUC=", sprintf("%.3f", unique(r$auc)[1])))
  }
  legend("bottomright", legend = labels, col = cols[cohort_order], lwd = 2.4, bty = "n")
  dev.off()
}

plot_km_endpoint <- function(endpoint, time_col, event_col, file_name) {
  png(file.path(fig_dir, file_name), width = 1700, height = 900, res = 170)
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1), oma = c(0, 0, 1, 0))
  plotted <- FALSE
  for (cohort in cohort_order) {
    d <- scores[scores$cohort == cohort & !is.na(scores[[time_col]]) & !is.na(scores[[event_col]]) & !is.na(scores$risk_group), ]
    if (nrow(d) < 15 || sum(d[[event_col]] == 1, na.rm = TRUE) < 5 || length(unique(d$risk_group)) < 2) {
      plot.new()
      title(main = cohort)
      text(0.5, 0.5, "Endpoint unavailable or too few events")
      next
    }
    plotted <- TRUE
    surv_formula <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ risk_group"))
    fit <- survfit(surv_formula, data = d)
    plot(
      fit, col = c(pal$low, pal$high), lwd = 2.2, mark.time = TRUE,
      xlab = "Months", ylab = paste(endpoint, "probability"), main = cohort
    )
    grid(col = grDevices::adjustcolor("gray70", 0.35))
    legend("topright", legend = c("Low risk", "High risk"), col = c(pal$low, pal$high), lwd = 2.2, bty = "n")
    sd <- tryCatch(survdiff(surv_formula, data = d), error = function(e) NULL)
    p <- if (is.null(sd)) NA_real_ else pchisq(sd$chisq, length(sd$n) - 1, lower.tail = FALSE)
    hr_fit <- tryCatch(coxph(surv_formula, data = d), error = function(e) NULL)
    hr_label <- ""
    if (!is.null(hr_fit)) {
      ss <- summary(hr_fit)
      hr_label <- paste0("HR=", sprintf("%.2f", ss$conf.int[1, "exp(coef)"]))
    }
    legend("bottomleft", legend = c(paste0("Log-rank p=", fmt_p(p)), hr_label), bty = "n")
  }
  mtext(paste(endpoint, "by median risk group"), outer = TRUE, cex = 1.1, font = 2)
  dev.off()
  plotted
}

plot_response_stack <- function() {
  png(file.path(fig_dir, "immunotherapy_response_stacked_bar.png"), width = 1600, height = 900, res = 170)
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2), oma = c(0, 0, 1, 0))
  for (cohort in cohort_order) {
    d <- scores[scores$cohort == cohort & !is.na(scores$response_group), ]
    tab <- table(d$risk_group, d$response_group)
    prop <- prop.table(tab, 1)
    bp <- barplot(
      t(prop[, c("Responder", "Non-responder"), drop = FALSE]),
      beside = FALSE, ylim = c(0, 1), col = c(pal$responder, pal$nonresponder),
      ylab = "Fraction of response-evaluable samples", main = cohort
    )
    text(bp, prop[, "Responder"] / 2, paste0(round(100 * prop[, "Responder"]), "%"), col = "white", font = 2)
    text(bp, prop[, "Responder"] + prop[, "Non-responder"] / 2, paste0(round(100 * prop[, "Non-responder"]), "%"), col = "white", font = 2)
    legend("topright", legend = c("Responder", "Non-responder"), fill = c(pal$responder, pal$nonresponder), bty = "n")
  }
  mtext("Response composition by risk group", outer = TRUE, cex = 1.1, font = 2)
  dev.off()
}

plot_signature_corr <- function() {
  corr <- read.delim(file.path(table_dir, "immunotherapy_signature_correlations.tsv"), sep = "\t", stringsAsFactors = FALSE)
  sig_order <- c("CD8_effector", "Checkpoint", "TGFb_CAF", "EMT_stromal", "Endothelial", "Myeloid", "B_cell", "Basal_epithelial", "Luminal_epithelial", "Proliferation")
  mat <- matrix(NA_real_, nrow = length(sig_order), ncol = length(cohort_order), dimnames = list(sig_order, cohort_order))
  for (i in seq_len(nrow(corr))) {
    if (corr$signature[i] %in% sig_order && corr$cohort[i] %in% cohort_order) {
      mat[corr$signature[i], corr$cohort[i]] <- corr$spearman_r_with_risk[i]
    }
  }
  png(file.path(fig_dir, "immunotherapy_signature_correlation_heatmap.png"), width = 1100, height = 1200, res = 170)
  par(mar = c(6, 8, 4, 6))
  cols <- colorRampPalette(c("#2C7BB6", "#F7F7F7", "#D7191C"))(101)
  zlim <- c(-1, 1)
  image(
    x = seq_len(ncol(mat)), y = seq_len(nrow(mat)),
    z = t(mat[nrow(mat):1, , drop = FALSE]), col = cols, zlim = zlim,
    axes = FALSE, xlab = "", ylab = "", main = "Spearman correlation with risk score"
  )
  axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2)
  axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2)
  for (x in seq_len(ncol(mat))) {
    for (y in seq_len(nrow(mat))) {
      val <- mat[rev(seq_len(nrow(mat)))[y], x]
      if (!is.na(val)) text(x, y, sprintf("%.2f", val), cex = 0.85)
    }
  }
  legend_y <- seq(0.2, 0.8, length.out = 101)
  usr <- par("usr")
  x0 <- usr[2] + 0.35
  rect(x0, legend_y[-101] * nrow(mat), x0 + 0.12, legend_y[-1] * nrow(mat), col = cols[-101], border = NA, xpd = NA)
  text(x0 + 0.22, c(0.2, 0.5, 0.8) * nrow(mat), c("-1", "0", "1"), xpd = NA, adj = 0)
  dev.off()
}

plot_balance <- function() {
  png(file.path(fig_dir, "immunotherapy_effector_suppression_balance.png"), width = 1700, height = 900, res = 170)
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2), oma = c(0, 0, 1, 0))
  for (cohort in cohort_order) {
    d <- scores[scores$cohort == cohort & !is.na(scores$response_group), ]
    boxplot(
      cd8_tgfb_balance ~ response_group, data = d, outline = FALSE,
      col = c(grDevices::adjustcolor(pal$responder, 0.35), grDevices::adjustcolor(pal$nonresponder, 0.35)),
      border = c(pal$responder, pal$nonresponder), lwd = 1.5,
      ylab = "CD8 effector - TGFb/CAF score", xlab = "", main = cohort
    )
    stripchart(cd8_tgfb_balance ~ response_group, data = d, vertical = TRUE, method = "jitter", add = TRUE, pch = 16, cex = 0.65, col = grDevices::adjustcolor("#1A202C", 0.45))
    p <- safe_wilcox(cd8_tgfb_balance ~ response_group, d)
    usr <- par("usr")
    text(1.5, usr[4] - 0.05 * diff(usr[3:4]), paste0("Wilcoxon p = ", fmt_p(p)), cex = 0.9)
  }
  mtext("Effector/suppression balance by response", outer = TRUE, cex = 1.1, font = 2)
  dev.off()
}

plot_model_gene_heatmap <- function() {
  expr_path <- file.path(processed_dir, "immunotherapy_selected_log_expression.tsv.gz")
  expr <- read.delim(gzfile(expr_path), sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  genes <- model$gene[model$gene %in% rownames(expr)]
  expr <- as.matrix(expr[genes, , drop = FALSE])
  expr_z <- t(scale(t(expr)))
  expr_z[is.na(expr_z)] <- 0
  scores$full_id <- paste(scores$cohort, scores$analysis_id, sep = "::")
  ann <- scores[match(colnames(expr_z), scores$full_id), ]
  ord <- order(ann$cohort, ann$risk_group, ann$response_group, ann$risk_score, na.last = TRUE)
  expr_z <- expr_z[, ord, drop = FALSE]
  ann <- ann[ord, ]
  png(file.path(fig_dir, "immunotherapy_model_gene_heatmap.png"), width = 1800, height = 900, res = 170)
  layout(matrix(c(1, 2), nrow = 2), heights = c(1.5, 7))
  par(mar = c(1, 7, 3, 2))
  risk_cols <- ifelse(ann$risk_group == "High", pal$high, pal$low)
  resp_cols <- ifelse(ann$response_group == "Responder", pal$responder, ifelse(ann$response_group == "Non-responder", pal$nonresponder, "gray80"))
  plot.new()
  plot.window(xlim = c(0.5, ncol(expr_z) + 0.5), ylim = c(0.4, 2.7))
  rect(seq_len(ncol(expr_z)) - 0.5, 0.6, seq_len(ncol(expr_z)) + 0.5, 1.4, col = risk_cols, border = NA)
  rect(seq_len(ncol(expr_z)) - 0.5, 1.6, seq_len(ncol(expr_z)) + 0.5, 2.4, col = resp_cols, border = NA)
  text(0, c(1, 2), c("Risk", "Response"), adj = 1, xpd = NA)
  title("Sample ordering: cohort, risk group, response, risk score")
  par(mar = c(5, 7, 1, 2))
  cols <- colorRampPalette(c("#2C7BB6", "#F7F7F7", "#D7191C"))(101)
  image(seq_len(ncol(expr_z)), seq_len(nrow(expr_z)), t(expr_z[nrow(expr_z):1, , drop = FALSE]), axes = FALSE, col = cols, zlim = c(-2.5, 2.5), xlab = "Samples", ylab = "")
  axis(2, at = seq_len(nrow(expr_z)), labels = rev(rownames(expr_z)), las = 2)
  abline(v = which(diff(as.numeric(factor(ann$cohort))) != 0) + 0.5, lwd = 2)
  dev.off()
}

plot_sankey_style <- function() {
  d <- scores[!is.na(scores$response_group), ]
  d$biomarker <- ifelse(
    d$cohort == "IMvigor210",
    paste0("Phenotype: ", ifelse(is.na(d$immune_phenotype), "unknown", as.character(d$immune_phenotype))),
    paste0("FGFR3: ", ifelse(is.na(d$fgfr3_altered), "unknown", as.character(d$fgfr3_altered)))
  )
  axes <- c("cohort", "risk_group", "biomarker", "response_group")
  flow <- aggregate(rep(1, nrow(d)), d[, axes], sum)
  names(flow)[ncol(flow)] <- "n"
  node_info <- list()
  for (ax in axes) {
    counts <- sort(tapply(d[[ax]], d[[ax]], length), decreasing = TRUE)
    gap <- 0.025
    usable <- 1 - gap * (length(counts) - 1)
    heights <- usable * as.numeric(counts) / sum(counts)
    top <- 1
    info <- data.frame(node = names(counts), y0 = NA_real_, y1 = NA_real_, yc = NA_real_, stringsAsFactors = FALSE)
    for (i in seq_along(counts)) {
      info$y1[i] <- top
      info$y0[i] <- top - heights[i]
      info$yc[i] <- mean(c(info$y0[i], info$y1[i]))
      top <- info$y0[i] - gap
    }
    rownames(info) <- info$node
    node_info[[ax]] <- info
  }
  x_pos <- c(0.08, 0.33, 0.64, 0.90)
  png(file.path(fig_dir, "immunotherapy_sankey_response_flow.png"), width = 1800, height = 1100, res = 170)
  par(mar = c(2, 2, 4, 2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  max_n <- max(flow$n)
  for (i in seq_len(nrow(flow))) {
    yy <- numeric(length(axes))
    for (j in seq_along(axes)) {
      info <- node_info[[axes[j]]]
      yy[j] <- info[as.character(flow[i, axes[j]]), "yc"]
    }
    col <- ifelse(flow$response_group[i] == "Responder", pal$responder, pal$nonresponder)
    xspline(x_pos, yy, shape = -0.45, lwd = 2 + 18 * flow$n[i] / max_n, border = grDevices::adjustcolor(col, 0.28), col = NA)
  }
  for (j in seq_along(axes)) {
    info <- node_info[[axes[j]]]
    rect(x_pos[j] - 0.025, info$y0, x_pos[j] + 0.025, info$y1, col = grDevices::adjustcolor(pal$light, 0.9), border = "gray55")
    text(ifelse(j < 3, x_pos[j] + 0.035, x_pos[j] - 0.035), info$yc, info$node, adj = ifelse(j < 3, 0, 1), cex = 0.75)
    text(x_pos[j], 1.035, c("Cohort", "Risk", "Biomarker", "Response")[j], font = 2, cex = 0.95, xpd = NA)
  }
  title("ICB validation flow: cohort to risk, biomarker and response")
  dev.off()
}

plot_circos_style <- function() {
  d <- scores[!is.na(scores$response_group), ]
  d$balance_group <- ave(d$cd8_tgfb_balance, d$cohort, FUN = function(x) ifelse(x >= median(x, na.rm = TRUE), "High CD8/TGFb balance", "Low CD8/TGFb balance"))
  d$biomarker_group <- ifelse(
    d$cohort == "IMvigor210",
    paste0("Phenotype: ", ifelse(is.na(d$immune_phenotype), "unknown", as.character(d$immune_phenotype))),
    paste0("FGFR3: ", ifelse(is.na(d$fgfr3_altered), "unknown", as.character(d$fgfr3_altered)))
  )
  link_df <- rbind(
    data.frame(from = paste(d$risk_group, "risk"), to = as.character(d$response_group), stringsAsFactors = FALSE),
    data.frame(from = paste(d$risk_group, "risk"), to = d$balance_group, stringsAsFactors = FALSE),
    data.frame(from = paste(d$risk_group, "risk"), to = d$biomarker_group, stringsAsFactors = FALSE)
  )
  link_df <- link_df[!is.na(link_df$from) & !is.na(link_df$to), ]
  links <- aggregate(rep(1, nrow(link_df)), link_df, sum)
  names(links)[3] <- "n"
  links <- links[links$n >= 5, ]
  nodes <- unique(c(links$from, links$to))
  risk_nodes <- c("Low risk", "High risk")
  nodes <- c(risk_nodes[risk_nodes %in% nodes], setdiff(nodes, risk_nodes))
  theta <- seq(pi / 2, pi / 2 - 2 * pi + 2 * pi / length(nodes), length.out = length(nodes))
  pos <- data.frame(node = nodes, x = cos(theta), y = sin(theta), stringsAsFactors = FALSE)
  rownames(pos) <- pos$node
  png(file.path(fig_dir, "immunotherapy_circos_association.png"), width = 1300, height = 1300, res = 170)
  par(mar = c(1, 1, 4, 1))
  plot.new()
  plot.window(xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35), asp = 1)
  symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, fg = "gray85")
  max_n <- max(links$n)
  for (i in seq_len(nrow(links))) {
    p1 <- pos[links$from[i], ]
    p2 <- pos[links$to[i], ]
    col <- ifelse(grepl("^High", links$from[i]), pal$high, pal$low)
    xspline(c(p1$x * 0.92, 0, p2$x * 0.92), c(p1$y * 0.92, 0, p2$y * 0.92), shape = -0.65, lwd = 1 + 10 * links$n[i] / max_n, border = grDevices::adjustcolor(col, 0.25), col = NA)
  }
  node_counts <- table(c(link_df$from, link_df$to))
  for (node in nodes) {
    p <- pos[node, ]
    col <- ifelse(node == "High risk", pal$high, ifelse(node == "Low risk", pal$low, pal$gray))
    points(p$x, p$y, pch = 21, bg = grDevices::adjustcolor(col, 0.55), col = col, cex = 1.4 + 1.6 * node_counts[node] / max(node_counts))
    text(p$x * 1.17, p$y * 1.17, node, cex = 0.72, xpd = NA)
  }
  title("Circos-style associations among risk, response and immune context")
  dev.off()
}

make_key_tables <- function(cox_tbl, glm_tbl) {
  response_metrics <- read.delim(file.path(table_dir, "immunotherapy_response_metrics.tsv"), sep = "\t", stringsAsFactors = FALSE)
  survival_metrics <- read.delim(file.path(table_dir, "immunotherapy_survival_metrics.tsv"), sep = "\t", stringsAsFactors = FALSE)
  response_metrics$auc_risk_predicts_nonresponse <- round(response_metrics$auc_risk_predicts_nonresponse, 3)
  response_metrics$wilcoxon_p_nonresponder_vs_responder <- signif(response_metrics$wilcoxon_p_nonresponder_vs_responder, 3)
  survival_metrics$c_index_risk_higher_hazard <- round(survival_metrics$c_index_risk_higher_hazard, 3)
  write.delim(response_metrics, file.path(table_dir, "immunotherapy_response_metrics_formatted.tsv"))
  write.delim(survival_metrics, file.path(table_dir, "immunotherapy_survival_metrics_formatted.tsv"))
}

html_table <- function(df, n = 8) {
  if (nrow(df) == 0) return("<p>No rows.</p>")
  df <- head(df, n)
  df[] <- lapply(df, function(x) {
    if (is.numeric(x)) signif(x, 4) else as.character(x)
  })
  rows <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  paste0(
    "<table><thead><tr>", paste0("<th>", names(df), "</th>", collapse = ""), "</tr></thead><tbody>",
    paste(rows, collapse = "\n"), "</tbody></table>"
  )
}

write_report <- function(cox_tbl, glm_tbl) {
  response_metrics <- read.delim(file.path(table_dir, "immunotherapy_response_metrics.tsv"), sep = "\t", stringsAsFactors = FALSE)
  survival_metrics <- read.delim(file.path(table_dir, "immunotherapy_survival_metrics.tsv"), sep = "\t", stringsAsFactors = FALSE)
  figs <- c(
    "immunotherapy_response_boxplot.png",
    "immunotherapy_response_roc.png",
    "immunotherapy_response_stacked_bar.png",
    "immunotherapy_km_os.png",
    "immunotherapy_km_pfs.png",
    "immunotherapy_effector_suppression_balance.png",
    "immunotherapy_signature_correlation_heatmap.png",
    "immunotherapy_model_gene_heatmap.png",
    "immunotherapy_sankey_response_flow.png",
    "immunotherapy_circos_association.png"
  )
  fig_tags <- paste0(
    "<section><h3>", figs, "</h3><img src=\"figures/", figs, "\" alt=\"", figs, "\"></section>",
    collapse = "\n"
  )
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>ICB validation report</title>",
    "<style>body{font-family:Arial,sans-serif;margin:28px;color:#1f2933}h1,h2{margin-bottom:0.3em}img{max-width:100%;border:1px solid #ddd;margin:8px 0 24px 0}table{border-collapse:collapse;margin:10px 0 24px 0;font-size:13px}td,th{border:1px solid #ddd;padding:5px 7px;text-align:left}th{background:#f3f4f6}.note{color:#4b5563}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:18px}</style></head><body>",
    "<h1>IMvigor210 and GSE176307 immunotherapy validation</h1>",
    "<p class=\"note\">Fixed TCGA/GEO-compatible 5-gene formula: EMP1, AHNAK, TNFRSF14, CLEC2D, GSDMB. Higher score denotes higher TCGA Cox risk. Response ROC uses non-response (SD/PD) as the positive outcome.</p>",
    "<h2>Response metrics</h2>",
    html_table(response_metrics, n = 10),
    "<h2>Survival metrics</h2>",
    html_table(survival_metrics, n = 10),
    "<h2>Cox models</h2>",
    html_table(cox_tbl[cox_tbl$term == "risk_score", ], n = 10),
    "<h2>Logistic response models</h2>",
    html_table(glm_tbl[glm_tbl$term == "risk_score", ], n = 10),
    "<h2>Figures</h2><div class=\"grid\">",
    fig_tags,
    "</div></body></html>"
  )
  writeLines(html, file.path(root, "results", "immunotherapy_response_validation_report.html"))
}

cox_tbl <- make_cox_tables(scores)
glm_tbl <- make_logistic_tables(scores)
make_key_tables(cox_tbl, glm_tbl)
plot_response_box()
plot_response_roc()
plot_response_stack()
plot_km_endpoint("OS", "os_time_months", "os_event", "immunotherapy_km_os.png")
plot_km_endpoint("PFS", "pfs_time_months", "pfs_event", "immunotherapy_km_pfs.png")
plot_balance()
plot_signature_corr()
plot_model_gene_heatmap()
plot_sankey_style()
plot_circos_style()
write_report(cox_tbl, glm_tbl)

cat("wrote immunotherapy figures and report\n")
cat(file.path(root, "results", "immunotherapy_response_validation_report.html"), "\n")
