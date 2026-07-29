$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent $PSScriptRoot
Set-Location $RepoRoot

New-Item -ItemType Directory -Force -Path "data", "results", "results\tables", "results\figures", "results\model", "logs" | Out-Null

if (Test-Path "models") {
  Copy-Item -Path "models\*.tsv", "models\*.json" -Destination "results\model" -Force -ErrorAction SilentlyContinue
}

Write-Host "[1/14] TCGA download and preprocessing"
python scripts\tcga_blca_pipeline.py all

Write-Host "[2/14] TCGA overview plots"
Rscript scripts\analysis_visualization_blca.R

Write-Host "[3/14] Differential expression"
Rscript scripts\differential_expression_screen.R

Write-Host "[4/14] Prognostic screening"
Rscript scripts\prognostic_screening.R

Write-Host "[5/14] TCGA fixed model fitting"
python scripts\fit_tcga_lasso_cox.py
Rscript scripts\evaluate_tcga_fixed_model.R

Write-Host "[6/14] GEO download and processing"
python scripts\download_geo_cohorts.py
python scripts\process_geo_validation.py
Rscript scripts\validate_geo_fixed_model.R

Write-Host "[7/14] GEO-compatible reduced model"
python scripts\fit_geo_compatible_reduced_model.py
Rscript scripts\evaluate_reduced_geo_compatible_model.R

Write-Host "[8/14] Downstream inputs"
python scripts\prepare_downstream_inputs.py

Write-Host "[9/14] Pathway and immune analyses"
python scripts\risk_pathway_immune_analysis.py
Rscript scripts\plot_risk_pathway_immune.R

Write-Host "[10/14] Multi-omics support"
python scripts\multiomics_support_analysis.py
Rscript scripts\plot_multiomics_support.R

Write-Host "[11/14] Drug sensitivity"
Rscript scripts\drug_sensitivity_gdsc2_ridge.R

Write-Host "[12/14] Single-cell GSE145137"
python scripts\download_single_cell_gse145137.py
python scripts\single_cell_gse145137_reduced_model_analysis.py
Rscript scripts\plot_single_cell_gse145137.R

Write-Host "[13/14] Candidate scRNA/spatial GEO sets"
python scripts\download_scRNA_candidate_geo_sets.py
python scripts\analyze_candidate_scRNA_geo_reduced_model.py
Rscript scripts\plot_candidate_scRNA_geo_comparison.R

Write-Host "[14/14] Immunotherapy validation and final HTML report"
python scripts\download_immunotherapy_cohorts.py
Rscript scripts\extract_imvigor210_core_data.R
python scripts\immunotherapy_response_validation.py
Rscript scripts\plot_immunotherapy_response_validation.R
Rscript scripts\generate_visual_report.R

Write-Host "Done. See results\figures, results\tables, results\model, and results\*.html."
