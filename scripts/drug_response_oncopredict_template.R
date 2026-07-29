#!/usr/bin/env Rscript
# Template for TCGA-BLCA drug-response prediction with oncoPredict.
#
# Required external training files, for example from the oncoPredict/GDSC data:
# - GDSC2_Expr.rds: genes x cell lines expression matrix
# - GDSC2_Res.rds: cell lines x drugs response matrix
#
# Example:
# Rscript scripts/drug_response_oncopredict_template.R data/external/GDSC2_Expr.rds data/external/GDSC2_Res.rds

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE) else file.path(getwd(), "scripts", "drug_response_oncopredict_template.R")
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

training_expr_path <- if (length(args) >= 1) args[1] else file.path(project_root, "data", "external", "GDSC2_Expr.rds")
training_resp_path <- if (length(args) >= 2) args[2] else file.path(project_root, "data", "external", "GDSC2_Res.rds")
test_expr_path <- file.path(project_root, "data", "processed", "blca_drug_prediction_expression_input.tsv.gz")
out_path <- file.path(project_root, "results", "tables", "oncoPredict_predicted_drug_response.tsv")

if (!requireNamespace("oncoPredict", quietly = TRUE)) {
  stop("Install oncoPredict first: install.packages('oncoPredict')", call. = FALSE)
}

if (!file.exists(training_expr_path) || !file.exists(training_resp_path)) {
  stop("Training RDS files not found. Provide GDSC/CTRP expression and response RDS paths.", call. = FALSE)
}

test_expr <- as.matrix(read.table(gzfile(test_expr_path), sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
training_expr <- readRDS(training_expr_path)
training_resp <- readRDS(training_resp_path)

common_genes <- intersect(rownames(training_expr), rownames(test_expr))
training_expr <- training_expr[common_genes, , drop = FALSE]
test_expr <- test_expr[common_genes, , drop = FALSE]

pred <- oncoPredict::calcPhenotype(
  trainingExprData = training_expr,
  trainingPtype = training_resp,
  testExprData = test_expr,
  batchCorrect = "standardize",
  powerTransformPhenotype = TRUE,
  removeLowVaryingGenes = 0.2,
  minNumSamples = 10,
  printOutput = TRUE,
  removeLowVaringGenesFrom = "rawData"
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.table(pred, out_path, sep = "\t", quote = FALSE, col.names = NA)
message("Wrote: ", out_path)
