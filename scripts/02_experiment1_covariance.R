# Experiment 1: High-Dimensional Covariance Estimation
#
# Generates the six-row Experiment 1 table for the manuscript using paired
# Monte Carlo simulations. Default manuscript configuration:
# n = 30, p = 100, rho = 0.7, M = 1000.
#
# Required project file: R/common.R
# Required package: MASS
# Outputs:
#   results/experiment1_covariance.rds
#   results/experiment1_covariance_summary.csv
#   results/experiment1_covariance_paired_summary.csv
#   tables/table4_experiment1.tex

source("R/common.R")

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required. Install it before running this script.",
       call. = FALSE)
}

set_reproducible_rng(2026L)
ensure_directories()

make_ar1_covariance <- function(p, rho = 0.7) {
  if (length(p) != 1L || !is.finite(p) || p < 1L || p != as.integer(p)) {
    stop("p must be a positive integer.", call. = FALSE)
  }

  if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) {
    stop("rho must be a finite scalar strictly between -1 and 1.",
         call. = FALSE)
  }

  rho ^ abs(outer(seq_len(as.integer(p)), seq_len(as.integer(p)), "-"))
}

safe_inverse <- function(Sigma_hat, inversion_ridge = 1e-4) {
  if (length(inversion_ridge) != 1L || !is.finite(inversion_ridge) ||
      inversion_ridge < 0) {
    stop("inversion_ridge must be a nonnegative finite scalar.",
         call. = FALSE)
  }

  p <- ncol(Sigma_hat)

  tryCatch(
    solve(Sigma_hat + inversion_ridge * diag(p)),
    error = function(e) {
      stop(
        "Unable to invert the shrinkage covariance matrix. Increase " ,
        "inversion_ridge or inspect the covariance estimate.",
        call. = FALSE
      )
    }
  )
}

shrinkage_covariance <- function(X, estimator = c("LW", "OAS"),
                                 inversion_ridge = 1e-4) {
  estimator <- match.arg(estimator)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  if (n <= 1L || p < 1L) {
    stop("X must have at least two rows and one column.", call. = FALSE)
  }

  X_c <- sweep(X, 2L, colMeans(X), FUN = "-")
  S_emp <- crossprod(X_c) / n
  mu <- sum(diag(S_emp)) / p
  target <- mu * diag(p)
  d2 <- sum((S_emp - target)^2)

  if (d2 <= .Machine$double.eps) {
    shrinkage <- 1
  } else if (estimator == "LW") {
    beta_hat <- sum(vapply(
      seq_len(n),
      function(i) {
        x_i <- X_c[i, ]
        sum((tcrossprod(x_i) - S_emp)^2)
      },
      numeric(1)
    )) / n^2

    shrinkage <- min(1, max(0, beta_hat / d2))
  } else {
    tr_s <- sum(diag(S_emp))
    tr_s2 <- sum(S_emp^2)
    numerator <- (1 - 2 / p) * tr_s2 + tr_s^2
    denominator <- (n + 1 - 2 / p) * (tr_s2 - tr_s^2 / p)

    shrinkage <- if (denominator <= .Machine$double.eps) {
      1
    } else {
      min(1, max(0, numerator / denominator))
    }
  }

  Sigma_hat <- (1 - shrinkage) * S_emp + shrinkage * target
  Omega_hat <- safe_inverse(Sigma_hat, inversion_ridge = inversion_ridge)

  list(
    Sigma = Sigma_hat,
    Precision = Omega_hat,
    shrinkage = shrinkage,
    empirical_covariance = S_emp
  )
}

covariance_metrics <- function(Sigma_hat, Sigma_true, Omega_hat, Omega_true) {
  c(
    frob_cov = norm(Sigma_hat - Sigma_true, type = "F"),
    frob_prec = norm(Omega_hat - Omega_true, type = "F"),
    trace_bias = sum(diag(Sigma_hat)) - sum(diag(Sigma_true))
  )
}

one_replication <- function(seed, n, p, Sigma_true, Omega_true,
                            inversion_ridge = 1e-4) {
  set.seed(as.integer(seed))

  X <- MASS::mvrnorm(
    n = n,
    mu = rep(0, p),
    Sigma = Sigma_true
  )

  Y <- surrogate_population(X)

  fits <- list(
    LWX = shrinkage_covariance(X, "LW", inversion_ridge),
    LWY = shrinkage_covariance(Y, "LW", inversion_ridge),
    OASX = shrinkage_covariance(X, "OAS", inversion_ridge),
    OASY = shrinkage_covariance(Y, "OAS", inversion_ridge)
  )

  metrics_lw_x <- covariance_metrics(
    fits$LWX$Sigma, Sigma_true, fits$LWX$Precision, Omega_true
  )
  metrics_lw_y <- covariance_metrics(
    fits$LWY$Sigma, Sigma_true, fits$LWY$Precision, Omega_true
  )
  metrics_oas_x <- covariance_metrics(
    fits$OASX$Sigma, Sigma_true, fits$OASX$Precision, Omega_true
  )
  metrics_oas_y <- covariance_metrics(
    fits$OASY$Sigma, Sigma_true, fits$OASY$Precision, Omega_true
  )

  c(
    LW_X_frob_cov = unname(metrics_lw_x["frob_cov"]),
    LW_Y_frob_cov = unname(metrics_lw_y["frob_cov"]),
    LW_X_frob_prec = unname(metrics_lw_x["frob_prec"]),
    LW_Y_frob_prec = unname(metrics_lw_y["frob_prec"]),
    LW_X_trace_bias = unname(metrics_lw_x["trace_bias"]),
    LW_Y_trace_bias = unname(metrics_lw_y["trace_bias"]),
    OAS_X_frob_cov = unname(metrics_oas_x["frob_cov"]),
    OAS_Y_frob_cov = unname(metrics_oas_y["frob_cov"]),
    OAS_X_frob_prec = unname(metrics_oas_x["frob_prec"]),
    OAS_Y_frob_prec = unname(metrics_oas_y["frob_prec"]),
    OAS_X_trace_bias = unname(metrics_oas_x["trace_bias"]),
    OAS_Y_trace_bias = unname(metrics_oas_y["trace_bias"])
  )
}

run_experiment1 <- function(n = 30L, p = 100L, rho = 0.7, M = 1000L,
                            inversion_ridge = 1e-4,
                            master_seed = 2026L) {
  n <- as.integer(n)
  p <- as.integer(p)
  M <- as.integer(M)

  if (n <= 1L || p <= 0L || M <= 1L) {
    stop("n must exceed 1, p must be positive, and M must exceed 1.",
         call. = FALSE)
  }

  Sigma_true <- make_ar1_covariance(p, rho)
  Omega_true <- solve(Sigma_true)
  replication_seeds <- as.integer(master_seed + seq_len(M) - 1L)

  message(sprintf(
    "Running Experiment 1: n = %d, p = %d, rho = %.1f, M = %d",
    n, p, rho, M
  ))

  raw <- do.call(
    rbind,
    lapply(
      replication_seeds,
      one_replication,
      n = n,
      p = p,
      Sigma_true = Sigma_true,
      Omega_true = Omega_true,
      inversion_ridge = inversion_ridge
    )
  )

  paired_raw <- cbind(
    LW_frob_cov_Y_minus_X = raw[, "LW_Y_frob_cov"] - raw[, "LW_X_frob_cov"],
    LW_frob_prec_Y_minus_X = raw[, "LW_Y_frob_prec"] - raw[, "LW_X_frob_prec"],
    LW_trace_bias_Y_minus_X = raw[, "LW_Y_trace_bias"] - raw[, "LW_X_trace_bias"],
    OAS_frob_cov_Y_minus_X = raw[, "OAS_Y_frob_cov"] - raw[, "OAS_X_frob_cov"],
    OAS_frob_prec_Y_minus_X = raw[, "OAS_Y_frob_prec"] - raw[, "OAS_X_frob_prec"],
    OAS_trace_bias_Y_minus_X = raw[, "OAS_Y_trace_bias"] - raw[, "OAS_X_trace_bias"]
  )

  list(
    n = n,
    p = p,
    rho = rho,
    M = M,
    inversion_ridge = inversion_ridge,
    master_seed = master_seed,
    raw = raw,
    paired_raw = paired_raw,
    summaries = summarize_columns(raw),
    paired_summaries = summarize_columns(paired_raw)
  )
}

make_experiment1_table <- function(result, label = "tab:exp1_results") {
  s <- result$summaries
  d <- result$paired_summaries

  row_text <- c(
    paste0(
      "Ledoit--Wolf & Frob. covariance error & ",
      format_mean_mcse(s, "LW_X_frob_cov"), " & ",
      format_mean_mcse(s, "LW_Y_frob_cov"), " & ",
      format_mean_mcse(d, "LW_frob_cov_Y_minus_X", 3L, 3L), " \\\\"
    ),
    paste0(
      "Ledoit--Wolf & Frob. precision error & ",
      format_mean_mcse(s, "LW_X_frob_prec"), " & ",
      format_mean_mcse(s, "LW_Y_frob_prec"), " & ",
      format_mean_mcse(d, "LW_frob_prec_Y_minus_X", 3L, 3L), " \\\\"
    ),
    paste0(
      "Ledoit--Wolf & Trace bias & ",
      format_mean_mcse(s, "LW_X_trace_bias"), " & ",
      format_mean_mcse(s, "LW_Y_trace_bias"), " & ",
      format_mean_mcse(d, "LW_trace_bias_Y_minus_X", 3L, 3L), " \\\\"
    ),
    "\\hline",
    paste0(
      "OAS & Frob. covariance error & ",
      format_mean_mcse(s, "OAS_X_frob_cov"), " & ",
      format_mean_mcse(s, "OAS_Y_frob_cov"), " & ",
      format_mean_mcse(d, "OAS_frob_cov_Y_minus_X", 3L, 3L), " \\\\"
    ),
    paste0(
      "OAS & Frob. precision error & ",
      format_mean_mcse(s, "OAS_X_frob_prec"), " & ",
      format_mean_mcse(s, "OAS_Y_frob_prec"), " & ",
      format_mean_mcse(d, "OAS_frob_prec_Y_minus_X", 3L, 3L), " \\\\"
    ),
    paste0(
      "OAS & Trace bias & ",
      format_mean_mcse(s, "OAS_X_trace_bias"), " & ",
      format_mean_mcse(s, "OAS_Y_trace_bias"), " & ",
      format_mean_mcse(d, "OAS_trace_bias_Y_minus_X", 3L, 3L), " \\\\"
    )
  )

  sprintf(
    paste0(
      "\\begin{table}[htbp]\n",
      "\\centering\n",
      "\\small\n",
      "\\setlength{\\tabcolsep}{4pt}\n",
      "\\caption{Experiment 1 results for $(n,p)=(%d,%d)$ and $\\rho=%.1f$, ",
      "based on $M=%d$ paired Monte Carlo replications. Entries for standard ",
      "and surrogate data are Monte Carlo means (MCSEs). The final column is ",
      "the paired difference $Y-X$ (MCSE); negative values favor the surrogate ",
      "for Frobenius-error metrics. For the precision-error metric, ",
      "$\\widehat{\\Omega}=(\\widehat{\\Sigma}+10^{-4}\\mathbf{I}_p)^{-1}$ ",
      "was used for both members of each paired replication.}\n",
      "\\label{%s}\n",
      "\\begin{tabular}{llccc}\n",
      "\\hline\\hline\n",
      "\\textbf{Estimator} & \\textbf{Metric} & ",
      "\\textbf{Standard ($\\mathbf{X}$)} & ",
      "\\textbf{Surrogate ($\\mathbf{Y}$)} & ",
      "\\textbf{Difference ($Y-X$)} \\\\\n",
      "\\hline\n",
      "%s\n",
      "\\hline\\hline\n",
      "\\end{tabular}\n",
      "\\end{table}\n"
    ),
    result$n,
    result$p,
    result$rho,
    result$M,
    label,
    paste(row_text, collapse = "\n")
  )
}

result_exp1 <- run_experiment1(
  n = 30L,
  p = 100L,
  rho = 0.7,
  M = 1000L,
  inversion_ridge = 1e-4,
  master_seed = 2026L
)

saveRDS(
  result_exp1,
  "results/experiment1_covariance.rds"
)

write_summary_csv(
  result_exp1$summaries,
  "results/experiment1_covariance_summary.csv"
)

write_summary_csv(
  result_exp1$paired_summaries,
  "results/experiment1_covariance_paired_summary.csv"
)

latex_exp1 <- make_experiment1_table(
  result_exp1,
  label = "tab:exp1_results"
)

writeLines(
  latex_exp1,
  "tables/table4_experiment1.tex"
)
cat(latex_exp1, "\n")
message("Experiment 1 completed.")
message("Saved results/experiment1_covariance.rds")
message("Saved results/experiment1_covariance_summary.csv")
message("Saved results/experiment1_covariance_paired_summary.csv")
message("Saved tables/table4_experiment1.tex")