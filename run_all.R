source("R/common.R")

set_reproducible_rng(2026L)
ensure_directories()

analysis_scripts <- c(
  "scripts/01_posthoc_verification.R",
  "scripts/02_experiment1_covariance.R",
  "scripts/03_experiment2_bootstrap_pca.R",
  "scripts/04_sensitivity_analyses.R",
  "scripts/04b_dimension_to_sample_size_sensitivity.R",
  "scripts/05_experiment3_prediction_conformal.R",
  "scripts/06_experiment4_kmeans.R",
  "scripts/07_ridge_alpha_robustness.R",
  "scripts/08_all_gene_expression.R"
)

for (script_path in analysis_scripts) {
  message("Running ", script_path)
  source(script_path, local = new.env(parent = globalenv()))
}

message("All analyses completed.")