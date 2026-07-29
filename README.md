# TCGA-BLCA multi-omics risk-score analysis

This repository contains the reproducible analysis code for a bladder cancer
(BLCA) study integrating TCGA transcriptomics, survival modeling, GEO external
validation, pathway/immune analyses, multi-omics support, single-cell/spatial
validation, immunotherapy response validation, and drug sensitivity analysis.

The core manuscript model is a GEO-compatible five-gene risk score trained in
TCGA-BLCA and validated without refitting in GSE13507, GSE31684, IMvigor210,
and GSE176307 where data were available.

## Highlights

- TCGA-BLCA open-data download and preprocessing from the GDC API.
- Tumor-normal differential expression screening.
- Prognostic screening, train/test split, LASSO-Cox-style fitting, and
  multivariable Cox validation.
- Reduced GEO-compatible model using genes shared by TCGA, GSE13507, and
  GSE31684.
- External GEO validation with fixed TCGA-derived coefficients.
- Pathway enrichment, immune-signature scoring, mutation/CNV/RPPA support.
- Single-cell and spatial candidate dataset testing.
- Immunotherapy response validation in IMvigor210 and GSE176307.
- Drug sensitivity screening with a GDSC2 ridge predictor and an oncoPredict
  template.

## Fixed reduced model

The primary GEO-compatible model is stored in:

- `models/reduced_geo_compatible_model.tsv`
- `models/reduced_geo_compatible_model.json`

Risk score formula:

```text
risk_score =
  0.056583946477291314 * z(EMP1)
  + 0.011862678975374186 * z(AHNAK)
  - 0.01868590790890934 * z(TNFRSF14)
  - 0.03127689274218831 * z(CLEC2D)
  - 0.10035965453013068 * z(GSDMB)
```

For TCGA modeling, z-scores use the TCGA train-set scaling recorded in the
model JSON. For GEO and immunotherapy validation cohorts, each cohort is
z-scored gene-wise before multiplying by the fixed coefficients.

## Repository structure

```text
.
|-- README.md
|-- CODE_MANIFEST.tsv
|-- requirements.txt
|-- r_packages.txt
|-- .gitignore
|-- scripts/
|   |-- tcga_blca_pipeline.py
|   |-- differential_expression_screen.R
|   |-- prognostic_screening.R
|   |-- fit_tcga_lasso_cox.py
|   |-- fit_geo_compatible_reduced_model.py
|   |-- validate_geo_fixed_model.R
|   |-- evaluate_reduced_geo_compatible_model.R
|   |-- risk_pathway_immune_analysis.py
|   |-- multiomics_support_analysis.py
|   |-- immunotherapy_response_validation.py
|   |-- single_cell_gse145137_reduced_model_analysis.py
|   |-- analyze_candidate_scRNA_geo_reduced_model.py
|   `-- other plotting/download helper scripts
|-- models/
|   |-- final_lasso_cox_model.json
|   |-- final_lasso_cox_model.tsv
|   |-- reduced_geo_compatible_model.json
|   `-- reduced_geo_compatible_model.tsv
|-- workflow/
    `-- run_all_windows.ps1
```

Large raw data, processed matrices, result tables, figures, local package
caches, and logs are intentionally not committed. They are regenerated under
`data/`, `results/`, and `logs/`.

## Data sources

The workflow uses public/open-access resources:

- TCGA-BLCA RNA-seq STAR counts, clinical data, somatic mutation MAF, masked
  copy-number segments, and RPPA data from the Genomic Data Commons.
- GEO validation cohorts: GSE13507 and GSE31684.
- Candidate single-cell/spatial datasets: GSE145137, GSE222315, GSE171351,
  GSE269877, and GSE293189 where downloadable files were available.
- Immunotherapy cohorts: IMvigor210 and GSE176307 where expression, clinical,
  response, and survival endpoints were available.
- Drug sensitivity: GDSC2-style external training files for ridge/oncoPredict
  analyses.

Please follow each data provider's license, access, and citation requirements.
Controlled-access data are not included.

## Environment

Python 3.10+ is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

R 4.2+ is recommended. Install the packages listed in `r_packages.txt`.

```r
install.packages(c("survival", "ggplot2"))
```

`oncoPredict` is optional and is only needed for the publication-grade drug
response template in `scripts/drug_response_oncopredict_template.R`.

## Quick start

From the repository root:

```bash
python scripts/tcga_blca_pipeline.py all
Rscript scripts/analysis_visualization_blca.R
Rscript scripts/differential_expression_screen.R
Rscript scripts/prognostic_screening.R
python scripts/fit_tcga_lasso_cox.py
Rscript scripts/evaluate_tcga_fixed_model.R
python scripts/download_geo_cohorts.py
python scripts/process_geo_validation.py
python scripts/fit_geo_compatible_reduced_model.py
Rscript scripts/evaluate_reduced_geo_compatible_model.R
python scripts/risk_pathway_immune_analysis.py
Rscript scripts/plot_risk_pathway_immune.R
python scripts/multiomics_support_analysis.py
Rscript scripts/plot_multiomics_support.R
Rscript scripts/drug_sensitivity_gdsc2_ridge.R
```

Windows users can run the bundled workflow helper:

```powershell
.\workflow\run_all_windows.ps1
```

Some downstream analyses require external files that may be large or hosted
outside GDC. If those files are unavailable locally, run the relevant download
script first or skip that module.

## Workflow modules

1. `tcga_blca_pipeline.py`: query GDC, download core TCGA-BLCA files, and build
   processed expression/clinical/mutation/CNV/RPPA matrices.
2. `analysis_visualization_blca.R`: basic cohort and QC plots.
3. `differential_expression_screen.R`: tumor-normal Wilcoxon screening and
   volcano/heatmap outputs.
4. `prognostic_screening.R`: univariate Cox screening and exploratory survival
   plots.
5. `fit_tcga_lasso_cox.py`: train/test split and TCGA-only risk model fitting.
6. `evaluate_tcga_fixed_model.R`: TCGA fixed-model Kaplan-Meier, ROC, C-index,
   and Cox analyses.
7. `download_geo_cohorts.py`, `process_geo_validation.py`,
   `validate_geo_fixed_model.R`: GEO download, probe-to-symbol processing, and
   fixed-model validation.
8. `fit_geo_compatible_reduced_model.py` and
   `evaluate_reduced_geo_compatible_model.R`: reduced 4-8 gene model fitting and
   external validation.
9. `risk_pathway_immune_analysis.py` and `plot_risk_pathway_immune.R`: pathway
   and immune signature analysis.
10. `multiomics_support_analysis.py` and `plot_multiomics_support.R`: mutation,
    CNV, and RPPA support.
11. `download_single_cell_gse145137.py`,
    `single_cell_gse145137_reduced_model_analysis.py`, and
    `plot_single_cell_gse145137.R`: single-cell localization.
12. `download_scRNA_candidate_geo_sets.py`,
    `analyze_candidate_scRNA_geo_reduced_model.py`, and
    `plot_candidate_scRNA_geo_comparison.R`: candidate scRNA/spatial validation.
13. `download_immunotherapy_cohorts.py`, `extract_imvigor210_core_data.R`,
    `immunotherapy_response_validation.py`, and
    `plot_immunotherapy_response_validation.R`: immunotherapy response,
    survival, Sankey/circos-style association, and immune context analysis.
14. `drug_sensitivity_gdsc2_ridge.R` and
    `drug_response_oncopredict_template.R`: drug sensitivity analysis.

See `CODE_MANIFEST.tsv` for a script-by-script index.

## Outputs

Generated outputs are written to:

```text
results/tables/
results/figures/
results/model/
results/*.html
```

The manuscript draft and figure/table plan are provided under `manuscript/` and
`docs/` for editing, but they are not a substitute for final manual biological
review.

## Reproducibility notes

- The reduced model uses seed `20260745` and a 70/30 TCGA train/test split.
- GEO validation applies fixed TCGA coefficients without refitting.
- Microarray/RNA-seq platform differences are handled by cohort-wise gene
  z-scoring for validation cohorts.
- The GDSC2 ridge predictor is a lightweight screening implementation; rerun
  `drug_response_oncopredict_template.R` with official oncoPredict training
  files for publication-grade drug response claims.
- Single-cell candidate datasets vary in format and annotation quality. Treat
  those outputs as localization support rather than independent clinical
  validation.

## Data and GitHub policy

Do not commit:

- `data/`
- `results/`
- `logs/`
- local virtual environments
- downloaded GDC/GEO/GDSC files
- local package caches such as `py_libs/`

The `.gitignore` in this repository excludes those files by default.

## Citation

If using this workflow, cite the relevant data sources and methods, including
TCGA/GDC, GEO accessions, IMvigor210, GSE176307, GDSC/oncoPredict where used,
and standard survival/differential-expression methods.

