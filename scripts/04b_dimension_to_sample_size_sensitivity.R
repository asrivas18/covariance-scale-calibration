# ==============================================================================
# Dimension-to-Sample-Size Sensitivity Analysis
#
# Generates manuscript Tables 9 and 10:
#   Panel A: fixed n = 100, varying p = 25, 50, 100, 200, 500
#   Panel B: fixed p = 100, varying n = 25, 50, 100, 200, 500
#
# Model:
#   X_i ~ N_p(0, Sigma), with Sigma_jk = rho^{|j-k|}, rho = 0.7.
#
# Estimator:
#   Sigma_hat_alpha = (1-alpha) S + alpha {tr(S)/p} I_p,
#   alpha = 0.10,
# where S is the denominator-n empirical covariance. The associated precision
# estimator is Omega_hat_alpha = Sigma_hat_alpha^{-1}.
#
# Outputs:
#   results/dimension_to_sample_size_sensitivity.rds
#   results/dimension_sensitivity_fixed_n.csv
#   results/dimension_sensitivity_fixed_p.csv
#   tables/table9_dimension_sensitivity_fixed_n.tex
#   tables/table10_dimension_sensitivity_fixed_p.tex
#
# Run from the repository root:
#   source("scripts/04b_dimension_to_sample_size_sensitivity.R")
# ==============================================================================

source("R/common.R")

set_reproducible_rng(2026L)
ensure_directories()

make_ar1_covariance <- function(p, rho = 0.7) {
  p <- as.integer(p)

  if (length(p) != 1L || !is.finite(p) || p < 1L) {
    stop("p must be one positive finite integer.", call. = FALSE)
  }

  if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) {
    stop(
      "rho must be a finite scalar strictly between -1 and 1.",
      call. = FALSE
    )
  }

  rho ^ abs(outer(seq_len(p), seq_len(p), "-"))
}

ridge_shrinkage_covariance <- function(S, alpha = 0.10) {
  S <- as.matrix(S)
  p <- ncol(S)

  if (nrow(S) != p) {
    stop("S must be a square matrix.", call. = FALSE)
  }

  if (length(alpha) != 1L || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("alpha must lie in [0, 1].", call. = FALSE)
  }

  trace_target <- sum(diag(S)) / p

  Sigma_hat <- (1 - alpha) * S +
    alpha * trace_target * diag(p)

  Omega_hat <- tryCatch(
    solve(Sigma_hat),
    error = function(e) {
      stop(
        "The ridge-shrinkage covariance estimate could not be inverted.",
        call. = FALSE
      )
    }
  )

  list(
    Sigma = Sigma_hat,
    Precision = Omega_hat
  )
}

relative_frobenius_error <- function(A, B) {
  denominator <- norm(B, type = "F")

  if (!is.finite(denominator) || denominator <= 0) {
    stop("Reference matrix has nonpositive Frobenius norm.", call. = FALSE)
  }

  norm(A - B, type = "F") / denominator
}

one_replication <- function(seed, n, p, rho = 0.7, alpha = 0.10) {
  set.seed(as.integer(seed))

  Sigma_true <- make_ar1_covariance(p = p, rho = rho)
  Omega_true <- solve(Sigma_true)
  chol_sigma <- chol(Sigma_true)

  X <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% chol_sigma
  Y <- surrogate_population(X)

  S_X <- empirical_covariance_n(X)
  S_Y <- empirical_covariance_n(Y)

  fit_X <- ridge_shrinkage_covariance(S_X, alpha = alpha)
  fit_Y <- ridge_shrinkage_covariance(S_Y, alpha = alpha)

  c(
    trace_bias_X =
      sum(diag(S_X)) - sum(diag(Sigma_true)),

    trace_bias_Y =
      sum(diag(S_Y)) - sum(diag(Sigma_true)),

    trace_ratio_Y_X =
      sum(diag(S_Y)) / sum(diag(S_X)),

    covariance_error_X =
      relative_frobenius_error(fit_X$Sigma, Sigma_true),

    covariance_error_Y =
      relative_frobenius_error(fit_Y$Sigma, Sigma_true),

    precision_error_X =
      relative_frobenius_error(fit_X$Precision, Omega_true),

    precision_error_Y =
      relative_frobenius_error(fit_Y$Precision, Omega_true)
  )
}

run_configuration <- function(n, p, rho = 0.7, alpha = 0.10,
                              M = 1000L, master_seed = 2026L) {
  n <- as.integer(n)
  p <- as.integer(p)
  M <- as.integer(M)

  if (n <= 1L || p <= 0L || M <= 1L) {
    stop(
      "n must exceed 1, p must be positive, and M must exceed 1.",
      call. = FALSE
    )
  }

  replication_seeds <- as.integer(
    master_seed + seq_len(M) - 1L
  )

  raw <- do.call(
    rbind,
    lapply(
      replication_seeds,
      one_replication,
      n = n,
      p = p,
      rho = rho,
      alpha = alpha
    )
  )

  paired_raw <- cbind(
    trace_bias_Y_minus_X =
      raw[, "trace_bias_Y"] - raw[, "trace_bias_X"],

    covariance_error_Y_minus_X =
      raw[, "covariance_error_Y"] - raw[, "covariance_error_X"],

    precision_error_Y_minus_X =
      raw[, "precision_error_Y"] - raw[, "precision_error_X"]
  )

  list(
    raw = raw,
    paired_raw = paired_raw,
    summaries = summarize_columns(raw),
    paired_summaries = summarize_columns(paired_raw)
  )
}

make_panel_summary <- function(configurations, panel_name) {
  metric_specs <- list(
    list(
      metric = "Trace bias",
      x_name = "trace_bias_X",
      y_name = "trace_bias_Y",
      difference_name = "trace_bias_Y_minus_X"
    ),
    list(
      metric = "Covariance error",
      x_name = "covariance_error_X",
      y_name = "covariance_error_Y",
      difference_name = "covariance_error_Y_minus_X"
    ),
    list(
      metric = "Precision error",
      x_name = "precision_error_X",
      y_name = "precision_error_Y",
      difference_name = "precision_error_Y_minus_X"
    )
  )

  rows <- lapply(configurations, function(config) {
    s <- config$result$summaries
    d <- config$result$paired_summaries

    do.call(
      rbind,
      lapply(metric_specs, function(spec) {
        data.frame(
          panel = panel_name,
          n = as.integer(config$n),
          p = as.integer(config$p),
          p_over_n = config$p / config$n,
          rho = config$rho,
          alpha = config$alpha,
          M = config$M,
          metric = spec$metric,
          x_mean = s[spec$x_name, "mean"],
          x_mcse = s[spec$x_name, "mcse"],
          y_mean = s[spec$y_name, "mean"],
          y_mcse = s[spec$y_name, "mcse"],
          difference_y_minus_x_mean =
            d[spec$difference_name, "mean"],
          difference_y_minus_x_mcse =
            d[spec$difference_name, "mcse"],
          row.names = NULL,
          check.names = FALSE
        )
      })
    )
  })

  do.call(rbind, rows)
}

format_entry <- function(mean_value, mcse_value,
                         digits_mean = 4L,
                         digits_mcse = 4L) {
  sprintf(
    paste0("%.", digits_mean, "f (%.", digits_mcse, "f)"),
    as.numeric(mean_value),
    as.numeric(mcse_value)
  )
}

make_latex_panel <- function(summary_df,
                             panel = c("fixed_n", "fixed_p"),
                             label,
                             caption) {
  panel <- match.arg(panel)

  required_columns <- c(
    "n",
    "p",
    "p_over_n",
    "metric",
    "x_mean",
    "x_mcse",
    "y_mean",
    "y_mcse",
    "difference_y_minus_x_mean",
    "difference_y_minus_x_mcse"
  )

  missing_columns <- setdiff(required_columns, names(summary_df))

  if (length(missing_columns) > 0L) {
    stop(
      "Summary data frame is missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  metric_order <- c(
    "Trace bias",
    "Covariance error",
    "Precision error"
  )

  if (panel == "fixed_n") {
    fixed_n <- unique(summary_df$n)

    if (length(fixed_n) != 1L) {
      stop("Panel A must have one fixed n value.", call. = FALSE)
    }

    p_values <- sort(unique(summary_df$p))
    configuration_keys <- paste(fixed_n, p_values, sep = "_")
  } else {
    fixed_p <- unique(summary_df$p)

    if (length(fixed_p) != 1L) {
      stop("Panel B must have one fixed p value.", call. = FALSE)
    }

    n_values <- sort(unique(summary_df$n))
    configuration_keys <- paste(n_values, fixed_p, sep = "_")
  }

  summary_df$config_key <- paste(
    summary_df$n,
    summary_df$p,
    sep = "_"
  )

  blocks <- lapply(configuration_keys, function(key) {
    block <- summary_df[
      summary_df$config_key == key,
      ,
      drop = FALSE
    ]

    block <- block[
      match(metric_order, block$metric),
      ,
      drop = FALSE
    ]

    if (nrow(block) != length(metric_order) ||
        anyNA(block$metric)) {
      stop(
        "Missing or duplicated metric rows for configuration ",
        key,
        call. = FALSE
      )
    }

    block
  })

  block_text <- vapply(blocks, function(block) {
    rows <- vapply(seq_len(nrow(block)), function(i) {
      digits <- if (block$metric[i] == "Trace bias") 3L else 4L

      configuration_fields <- if (i == 1L) {
        sprintf(
          "%d & %d & %.2f",
          as.integer(block$n[i]),
          as.integer(block$p[i]),
          as.numeric(block$p_over_n[i])
        )
      } else {
        " &  & "
      }

      sprintf(
        "%s & %s & %s & %s & %s \\\\",
        configuration_fields,
        block$metric[i],
        format_entry(
          block$x_mean[i],
          block$x_mcse[i],
          digits_mean = digits,
          digits_mcse = digits
        ),
        format_entry(
          block$y_mean[i],
          block$y_mcse[i],
          digits_mean = digits,
          digits_mcse = digits
        ),
        format_entry(
          block$difference_y_minus_x_mean[i],
          block$difference_y_minus_x_mcse[i],
          digits_mean = digits,
          digits_mcse = digits
        )
      )
    }, character(1))

    paste(rows, collapse = "\n")
  }, character(1))

  paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\begin{tabular}{rrrlrrr}\n",
    "\\hline\\hline\n",
    "\\textbf{$n$} & \\textbf{$p$} & \\textbf{$p/n$} & ",
    "\\textbf{Metric} & \\textbf{$X$} & \\textbf{$Y$} & ",
    "\\textbf{$Y-X$} \\\\\n",
    "\\hline\n",
    paste(block_text, collapse = "\n\\hline\n"),
    "\n\\hline\\hline\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
}

run_panel_a <- function(p_grid = c(25L, 50L, 100L, 200L, 500L),
                        n = 100L,
                        rho = 0.7,
                        alpha = 0.10,
                        M = 1000L,
                        master_seed = 2026L) {
  lapply(seq_along(p_grid), function(i) {
    p <- as.integer(p_grid[i])

    message(sprintf(
      "Panel A: n = %d, p = %d, p/n = %.2f, M = %d",
      n, p, p / n, M
    ))

    list(
      n = as.integer(n),
      p = p,
      rho = rho,
      alpha = alpha,
      M = as.integer(M),
      result = run_configuration(
        n = n,
        p = p,
        rho = rho,
        alpha = alpha,
        M = M,
        master_seed = master_seed + 100000L * i
      )
    )
  })
}

run_panel_b <- function(n_grid = c(25L, 50L, 100L, 200L, 500L),
                        p = 100L,
                        rho = 0.7,
                        alpha = 0.10,
                        M = 1000L,
                        master_seed = 4026L) {
  lapply(seq_along(n_grid), function(i) {
    n <- as.integer(n_grid[i])

    message(sprintf(
      "Panel B: n = %d, p = %d, p/n = %.2f, M = %d",
      n, p, p / n, M
    ))

    list(
      n = n,
      p = as.integer(p),
      rho = rho,
      alpha = alpha,
      M = as.integer(M),
      result = run_configuration(
        n = n,
        p = p,
        rho = rho,
        alpha = alpha,
        M = M,
        master_seed = master_seed + 100000L * i
      )
    )
  })
}

M_final <- 1000L
rho_final <- 0.7
alpha_final <- 0.10

panel_a_runs <- run_panel_a(
  p_grid = c(25L, 50L, 100L, 200L, 500L),
  n = 100L,
  rho = rho_final,
  alpha = alpha_final,
  M = M_final,
  master_seed = 2026L
)

panel_b_runs <- run_panel_b(
  n_grid = c(25L, 50L, 100L, 200L, 500L),
  p = 100L,
  rho = rho_final,
  alpha = alpha_final,
  M = M_final,
  master_seed = 4026L
)

panel_a_summary <- make_panel_summary(
  panel_a_runs,
  panel_name = "fixed_n_varying_p"
)

panel_b_summary <- make_panel_summary(
  panel_b_runs,
  panel_name = "fixed_p_varying_n"
)

saveRDS(
  list(
    panel_a_runs = panel_a_runs,
    panel_b_runs = panel_b_runs,
    panel_a_summary = panel_a_summary,
    panel_b_summary = panel_b_summary
  ),
  file = "results/dimension_to_sample_size_sensitivity.rds"
)

write.csv(
  panel_a_summary,
  file = "results/dimension_sensitivity_fixed_n.csv",
  row.names = FALSE
)

write.csv(
  panel_b_summary,
  file = "results/dimension_sensitivity_fixed_p.csv",
  row.names = FALSE
)

table9_caption <- paste0(
  "Panel A: dimension-to-sample-size sensitivity with fixed $n=100$ and ",
  "varying $p$. Trace bias is reported on the covariance-trace scale; ",
  "covariance and precision errors are relative Frobenius-norm errors. ",
  "For every row, the deterministic trace ratio is ",
  "$\\operatorname{tr}(\\mathbf{S}_Y)/\\operatorname{tr}(\\mathbf{S}_X)",
  "=100/99=1.01010$. Precision-error comparisons are specific to the fixed ",
  "ridge-shrinkage parameter $\\alpha=0.10$ and should not be interpreted ",
  "as a universal precision-estimation result."
)

table10_caption <- paste0(
  "Panel B: dimension-to-sample-size sensitivity with fixed $p=100$ and ",
  "varying $n$. Trace bias is reported on the covariance-trace scale; ",
  "covariance and precision errors are relative Frobenius-norm errors. ",
  "The trace ratio is deterministic conditional on each simulated sample and ",
  "equals $n/(n-1)$. Precision-error comparisons are specific to the fixed ",
  "ridge-shrinkage parameter $\\alpha=0.10$ and should not be interpreted ",
  "as a universal precision-estimation result."
)

latex_table9 <- make_latex_panel(
  panel_a_summary,
  panel = "fixed_n",
  label = "tab:pn_sensitivity_fixed_n",
  caption = table9_caption
)

latex_table10 <- make_latex_panel(
  panel_b_summary,
  panel = "fixed_p",
  label = "tab:pn_sensitivity_fixed_p",
  caption = table10_caption
)

writeLines(
  latex_table9,
  con = "tables/table9_dimension_sensitivity_fixed_n.tex"
)

writeLines(
  latex_table10,
  con = "tables/table10_dimension_sensitivity_fixed_p.tex"
)

cat("\n--- Table 9: fixed n, varying p ---\n\n")
cat(latex_table9, "\n")

cat("\n--- Table 10: fixed p, varying n ---\n\n")
cat(latex_table10, "\n")

message("Dimension-to-sample-size sensitivity analysis completed.")
message("Saved results/dimension_to_sample_size_sensitivity.rds")
message("Saved results/dimension_sensitivity_fixed_n.csv")
message("Saved results/dimension_sensitivity_fixed_p.csv")
message("Saved tables/table9_dimension_sensitivity_fixed_n.tex")
message("Saved tables/table10_dimension_sensitivity_fixed_p.tex")