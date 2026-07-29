#!/usr/bin/env Rscript
# Plot TCGA fixed-risk pathway and immune-support analyses.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE) else file.path(getwd(), "scripts", "plot_risk_pathway_immune.R")
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
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

enrich <- read.delim(file.path(table_dir, "tcga_fixed_risk_pathway_enrichment.tsv"), check.names = FALSE)
immune <- read.delim(file.path(table_dir, "tcga_fixed_risk_immune_signature_stats.tsv"), check.names = FALSE)
screen <- read.delim(file.path(table_dir, "tcga_fixed_risk_high_vs_low_expression_screen.tsv"), check.names = FALSE)

save_png("tcga_fixed_risk_pathway_enrichment.png", width = 2100, height = 1600, {
  hall <- enrich[enrich$library == "Hallmark_2020", ]
  up <- head(hall[hall$direction == "High_risk_up", ], 15)
  down <- head(hall[hall$direction == "High_risk_down", ], 10)
  dat <- rbind(up, down)
  dat$score <- -log10(pmax(dat$FDR, .Machine$double.xmin))
  dat$label <- gsub("_", " ", dat$term)
  vals <- dat$score
  names(vals) <- dat$label
  vals <- rev(vals)
  cols <- ifelse(rev(dat$direction) == "High_risk_up", palette_main[2], palette_main[1])
  par(mar = c(5, 12, 4, 1))
  bp <- barplot(vals, horiz = TRUE, las = 1, col = cols, border = NA,
                main = "Hallmark enrichment in fixed-risk expression programs",
                xlab = "-log10(FDR)")
  legend("bottomright", legend = c("High-risk up", "High-risk down"), fill = c(palette_main[2], palette_main[1]), bty = "n")
})

save_png("tcga_fixed_risk_immune_signatures.png", width = 1800, height = 1200, {
  immune <- immune[order(immune$diff_high_minus_low), ]
  vals <- immune$diff_high_minus_low
  names(vals) <- immune$signature
  cols <- ifelse(vals > 0, palette_main[2], palette_main[1])
  par(mar = c(5, 10, 4, 1))
  pad <- diff(range(vals)) * 0.2
  bp <- barplot(vals, horiz = TRUE, las = 1, col = cols, border = NA,
                xlim = c(min(vals) - pad, max(vals) + pad),
                main = "Immune and stromal marker scores: high vs low risk",
                xlab = "Mean score difference, high - low")
  abline(v = 0, col = "gray40")
  text(vals, bp, labels = sprintf("FDR %.1e", immune$FDR), pos = ifelse(vals > 0, 4, 2), cex = 0.65)
})

save_png("tcga_fixed_risk_expression_volcano.png", {
  sig <- screen$FDR < 0.05 & abs(screen$diff_high_minus_low) > 0.5
  cols <- rep(adjustcolor("gray40", alpha.f = 0.35), nrow(screen))
  cols[sig & screen$diff_high_minus_low > 0] <- adjustcolor(palette_main[2], alpha.f = 0.75)
  cols[sig & screen$diff_high_minus_low < 0] <- adjustcolor(palette_main[1], alpha.f = 0.75)
  plot(screen$diff_high_minus_low, -log10(pmax(screen$FDR, .Machine$double.xmin)),
       pch = 16, cex = 0.45, col = cols,
       xlab = "High-risk - low-risk mean log2(TPM + 1)", ylab = "-log10(FDR)",
       main = "Risk-group expression program")
  abline(v = c(-0.5, 0.5), col = "gray60", lty = 2)
  abline(h = -log10(0.05), col = "gray60", lty = 2)
  top <- head(screen[order(screen$FDR), ], 12)
  text(top$diff_high_minus_low, -log10(pmax(top$FDR, .Machine$double.xmin)), labels = top$gene,
       cex = 0.65, pos = ifelse(top$diff_high_minus_low > 0, 4, 2))
})

html <- c(
  "<!doctype html>",
  "<html><head><meta charset=\"utf-8\"><title>Risk pathway and immune support</title>",
  "<style>body{font-family:Arial,sans-serif;margin:24px;color:#222}h1{font-size:24px}section{margin:24px 0 36px}img{max-width:100%;border:1px solid #ddd}</style>",
  "</head><body>",
  "<h1>TCGA fixed-risk pathway and immune support</h1>",
  "<section><h2>Hallmark enrichment</h2><img src=\"figures/tcga_fixed_risk_pathway_enrichment.png\"></section>",
  "<section><h2>Immune/stromal marker scores</h2><img src=\"figures/tcga_fixed_risk_immune_signatures.png\"></section>",
  "<section><h2>Risk-group expression program</h2><img src=\"figures/tcga_fixed_risk_expression_volcano.png\"></section>",
  "</body></html>"
)
writeLines(html, file.path(project_root, "results", "risk_pathway_immune_report.html"))
message("Risk pathway/immune report complete: ", file.path(project_root, "results", "risk_pathway_immune_report.html"))
