# ==============================================================================
# Real-Data Gene-Expression Identity Illustration
#
# Online Appendix 2, Table 3:
# Global covariance-scale calibration using the Bioconductor ALL package.
#
# The analysis:
#   - selects observations whose ALL package metadata variable BT begins with B;
#   - uses the first 50 expression probes in the package feature order;
#   - verifies the covariance-trace, PCA-eigenvalue, and percentile-bootstrap
#     interval-width scaling identities under surrogate preprocessing.
#
# This is an empirical identity illustration only. It does not estimate
# population covariance bias, predictive performance, or interval coverage.
#
# Outputs:
#   results/all_real_data_covariance_pca.rds
#   results/all_real_data_covariance_pca_results.csv
#   results/all_real_data_configuration.csv
#   results/all_dataset_citation.txt
#   results/all_real_data_session_info.txt
#   tables/app2_table3_all_identity.tex
#
# Run from the repository root:
#   source("scripts/08_all_gene_expression.R")
# ==============================================================================

source("R/common.R")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (package_name in c("ALL", "Biobase")) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    BiocManager::install(
      package_name,
      ask = FALSE,
      update = FALSE
    )
  }
}

suppressPackageStartupMessages({
  library(ALL)
  library(Biobase)
})

set_reproducible_rng(20260811L)
ensure_directories()

B <- 2000L
panel_size <- 50L
alpha <- 0.05
master_seed <- 20260811L

data(ALL)

bt <- as.character(Biobase::pData(ALL)$BT)

if (anyNA(bt)) {
  stop(
    "The ALL metadata field 'BT' contains missing values.",
    call. = FALSE
  )
}

b_cell_index <- grepl("^B", bt)

if (!any(b_cell_index)) {
  stop(
    "No observations beginning with 'B' were found in ALL metadata field 'BT'.",
    call. = FALSE
  )
}

all_b <- ALL[, b_cell_index]
expression_matrix <- Biobase::exprs(all_b)

if (nrow(expression_matrix) < panel_size) {
  stop(
    "The requested probe panel is larger than the number of available probes.",
    call. = FALSE
  )
}

probe_panel <- Biobase::featureNames(all_b)[seq_len(panel_size)]

expression_matrix <- expression_matrix[
  probe_panel,
  ,
  drop = FALSE
]

X <- t(expression_matrix)
n <- nrow(X)
p <- ncol(X)

if (n <= 1L) {
  stop(
    "At least two selected observations are required.",
    call. = FALSE
  )
}

if (p != panel_size) {
  stop(
    "The selected expression panel does not have the requested number of probes.",
    call. = FALSE
  )
}

Y <- surrogate_population(X)

S_X <- empirical_covariance_n(X)
S_Y <- empirical_covariance_n(Y)

eigenvalues_X <- eigen(
  S_X,
  symmetric = TRUE,
  only.values = TRUE
)$values

eigenvalues_Y <- eigen(
  S_Y,
  symmetric = TRUE,
  only.values = TRUE
)$values

if (length(eigenvalues_X) < 2L || length(eigenvalues_Y) < 2L) {
  stop(
    "At least two covariance eigenvalues are required for this analysis.",
    call. = FALSE
  )
}

theoretical_ratio <- n / (n - 1L)

trace_ratio <- sum(diag(S_Y)) / sum(diag(S_X))

eigen_ratios <- eigenvalues_Y[seq_len(2L)] /
  eigenvalues_X[seq_len(2L)]

bootstrap_eigenvalues_X <- matrix(
  NA_real_,
  nrow = B,
  ncol = 2L
)

bootstrap_eigenvalues_Y <- matrix(
  NA_real_,
  nrow = B,
  ncol = 2L
)

for (b in seq_len(B)) {
  index <- sample.int(n, size = n, replace = TRUE)

  S_X_boot <- empirical_covariance_n(
    X[index, , drop = FALSE]
  )

  S_Y_boot <- empirical_covariance_n(
    Y[index, , drop = FALSE]
  )

  bootstrap_eigenvalues_X[b, ] <- eigen(
    S_X_boot,
    symmetric = TRUE,
    only.values = TRUE
  )$values[seq_len(2L)]

  bootstrap_eigenvalues_Y[b, ] <- eigen(
    S_Y_boot,
    symmetric = TRUE,
    only.values = TRUE
  )$values[seq_len(2L)]
}

interval_width <- function(values, alpha = 0.05) {
  endpoints <- stats::quantile(
    values,
    probs = c(alpha / 2, 1 - alpha / 2),
    names = FALSE,
    type = 7
  )

  unname(diff(endpoints))
}

width_X <- apply(
  bootstrap_eigenvalues_X,
  2L,
  interval_width,
  alpha = alpha
)

width_Y <- apply(
  bootstrap_eigenvalues_Y,
  2L,
  interval_width,
  alpha = alpha
)

width_ratios <- width_Y / width_X

real_data_results <- data.frame(
  quantity = c(
    "Covariance trace",
    "First PCA eigenvalue",
    "Second PCA eigenvalue",
    "First percentile-bootstrap interval width",
    "Second percentile-bootstrap interval width"
  ),
  ratio_y_x = c(
    trace_ratio,
    eigen_ratios[1],
    eigen_ratios[2],
    width_ratios[1],
    width_ratios[2]
  ),
  theoretical_ratio = rep(theoretical_ratio, 5L),
  row.names = NULL
)

tolerance <- 1e-10

if (any(abs(
  real_data_results$ratio_y_x -
    real_data_results$theoretical_ratio
) > tolerance)) {
  stop(
    "At least one real-data covariance-scale identity check exceeded tolerance.",
    call. = FALSE
  )
}

print(real_data_results, row.names = FALSE)

latex_num <- function(x, digits = 6L) {
  formatC(
    as.numeric(x),
    format = "f",
    digits = digits
  )
}

latex_rows <- vapply(
  seq_len(nrow(real_data_results)),
  function(i) {
    paste0(
      real_data_results$quantity[i], " & ",
      latex_num(real_data_results$ratio_y_x[i]), " & ",
      latex_num(real_data_results$theoretical_ratio[i]),
      " \\\\"
    )
  },
  character(1)
)

latex_table <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\small",
  paste0(
    "\\caption{Real-data illustration of global covariance-scale calibration ",
    "using $n=", n, "$ observations designated as B-cell acute lymphoblastic ",
    "leukemia in the Bioconductor \\texttt{ALL} package metadata and ",
    "$p=", p, "$ pre-specified expression probes. The probe panel consists of ",
    "the first ", panel_size, " probes in package feature order. All ratios ",
    "equal $n/(n-1)$ up to floating-point precision.}"
  ),
  "\\label{tab:real_data_covariance_pca}",
  "\\begin{tabular}{lcc}",
  "\\hline\\hline",
  "\\textbf{Quantity} & \\textbf{Observed ratio: $Y/X$} & ",
  "\\textbf{Theoretical ratio: $n/(n-1)$} \\\\",
  "\\hline",
  latex_rows,
  "\\hline\\hline",
  "\\end{tabular}",
  "\\end{table}"
)

analysis_configuration <- data.frame(
  parameter = c(
    "master_seed",
    "bootstrap_replicates_B",
    "nominal_alpha",
    "metadata_field",
    "metadata_selection_rule",
    "selected_observations_n",
    "probe_panel_size_p",
    "probe_selection_rule",
    "theoretical_ratio_n_over_n_minus_1"
  ),
  value = c(
    as.character(master_seed),
    as.character(B),
    as.character(alpha),
    "BT",
    "grepl('^B', as.character(pData(ALL)$BT))",
    as.character(n),
    as.character(p),
    paste0(
      "First ", panel_size,
      " probes in Biobase::featureNames(ALL) package order"
    ),
    format(theoretical_ratio, digits = 16)
  ),
  row.names = NULL
)

saveRDS(
  list(
    real_data_results = real_data_results,
    selected_sample_names = Biobase::sampleNames(all_b),
    selected_probe_panel = probe_panel,
    bootstrap_eigenvalues_X = bootstrap_eigenvalues_X,
    bootstrap_eigenvalues_Y = bootstrap_eigenvalues_Y,
    configuration = analysis_configuration
  ),
  file = "results/all_real_data_covariance_pca.rds"
)

write.csv(
  real_data_results,
  file = "results/all_real_data_covariance_pca_results.csv",
  row.names = FALSE
)

write.csv(
  analysis_configuration,
  file = "results/all_real_data_configuration.csv",
  row.names = FALSE
)

writeLines(
  latex_table,
  con = "tables/app2_table3_all_identity.tex"
)

writeLines(
  capture.output(citation("ALL")),
  con = "results/all_dataset_citation.txt"
)

writeLines(
  capture.output(sessionInfo()),
  con = "results/all_real_data_session_info.txt"
)

message("Real-data gene-expression identity analysis completed.")
message("Saved results/all_real_data_covariance_pca.rds")
message("Saved results/all_real_data_covariance_pca_results.csv")
message("Saved results/all_real_data_configuration.csv")
message("Saved results/all_dataset_citation.txt")
message("Saved results/all_real_data_session_info.txt")
message("Saved tables/app2_table3_all_identity.tex")