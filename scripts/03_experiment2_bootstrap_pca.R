# ==============================================================================
# Experiment 2: Bootstrap PCA and Eigenvalue Intervals
#
# Ordinary percentile bootstrap on X and surrogate Y.
# BCa interval benchmark on X.
#
# Main-manuscript outputs:
#   Tables 5 and 6
#
# Online Appendix 2 output:
#   Table 1
#
# Default manuscript configuration:
#   n = 25, p = 5, B = 2000, M = 1000
#   population eigenvalues = (10, 5, 2, 1, 0.5)
#
# Outputs:
#   results/experiment2_bootstrap_pca.rds
#   results/experiment2_bootstrap_pca_summary.csv
#   results/experiment2_bootstrap_pca_paired_summary.csv
#   results/experiment2_bootstrap_pca_configuration.csv
#   results/experiment2_bootstrap_pca_session_info.txt
#   tables/table5_experiment2_eigenvalues.tex
#   tables/table6_experiment2_coverage.tex
#   tables/app2_table1_pca_widths.tex
#
# Run from repository root:
#   source("scripts/03_experiment2_bootstrap_pca.R")
# ==============================================================================

source("R/common.R")

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required.", call. = FALSE)
}

set_reproducible_rng(1234L)
ensure_directories()

master_seed <- 1234L

generate_pca_data <- function(
  n,
  p = 5L,
  true_eigenvalues = c(10, 5, 2, 1, 0.5)
) {
  n <- as.integer(n)
  p <- as.integer(p)

  if (n <= 1L) {
    stop("n must exceed one.", call. = FALSE)
  }

  if (length(true_eigenvalues) != p) {
    stop("length(true_eigenvalues) must equal p.", call. = FALSE)
  }

  if (any(!is.finite(true_eigenvalues)) || any(true_eigenvalues <= 0)) {
    stop(
      "All population eigenvalues must be finite and positive.",
      call. = FALSE
    )
  }

  Z <- matrix(rnorm(p * p), nrow = p, ncol = p)
  qr_Z <- qr(Z)
  Q <- qr.Q(qr_Z)

  diagonal_R <- diag(qr.R(qr_Z))
  signs <- ifelse(diagonal_R >= 0, 1, -1)
  Q <- sweep(Q, 2L, signs, FUN = "*")

  Sigma_true <- Q %*%
    diag(true_eigenvalues) %*%
    t(Q)

  X <- MASS::mvrnorm(
    n = n,
    mu = rep(0, p),
    Sigma = Sigma_true
  )

  list(
    X = X,
    Sigma = Sigma_true,
    true_eigenvalues = sort(true_eigenvalues, decreasing = TRUE)
  )
}

pca_eigenvalues <- function(X, k = 2L) {
  values <- eigen(
    empirical_covariance_n(X),
    symmetric = TRUE,
    only.values = TRUE
  )$values

  if (k > length(values)) {
    stop("k exceeds the number of PCA eigenvalues.", call. = FALSE)
  }

  values[seq_len(k)]
}

percentile_interval_local <- function(values, conf = 0.95) {
  alpha <- 1 - conf

  as.numeric(
    stats::quantile(
      values,
      probs = c(alpha / 2, 1 - alpha / 2),
      names = FALSE,
      type = 7
    )
  )
}

bca_interval <- function(theta_hat, bootstrap_statistics, jackknife_statistics,
                         conf = 0.95) {
  B <- length(bootstrap_statistics)
  alpha <- 1 - conf

  proportion_less <- (
    sum(bootstrap_statistics < theta_hat) + 0.5
  ) / (B + 1)

  z0 <- stats::qnorm(proportion_less)

  jackknife_mean <- mean(jackknife_statistics)
  influence <- jackknife_mean - jackknife_statistics

  denominator <- 6 * sum(influence^2)^(3 / 2)

  acceleration <- if (denominator <= .Machine$double.eps) {
    0
  } else {
    sum(influence^3) / denominator
  }

  z_alpha <- stats::qnorm(c(alpha / 2, 1 - alpha / 2))

  adjusted_probabilities <- stats::pnorm(
    z0 + (z0 + z_alpha) /
      (1 - acceleration * (z0 + z_alpha))
  )

  adjusted_probabilities <- pmin(
    pmax(adjusted_probabilities, 1 / (B + 1)),
    B / (B + 1)
  )

  as.numeric(
    stats::quantile(
      bootstrap_statistics,
      probs = adjusted_probabilities,
      names = FALSE,
      type = 7
    )
  )
}

one_exp2_replication <- function(
  seed,
  n = 25L,
  p = 5L,
  true_eigenvalues = c(10, 5, 2, 1, 0.5),
  B = 2000L,
  k = 2L,
  conf = 0.95
) {
  set.seed(as.integer(seed))

  dat <- generate_pca_data(
    n = n,
    p = p,
    true_eigenvalues = true_eigenvalues
  )

  X <- dat$X
  Y <- surrogate_population(X)

  target_eigenvalues <- dat$true_eigenvalues[seq_len(k)]

  theta_X <- pca_eigenvalues(X, k = k)
  theta_Y <- pca_eigenvalues(Y, k = k)

  bootstrap_X <- matrix(NA_real_, nrow = B, ncol = k)
  bootstrap_Y <- matrix(NA_real_, nrow = B, ncol = k)

  bootstrap_indices <- matrix(
    NA_integer_,
    nrow = n,
    ncol = B
  )

  for (b in seq_len(B)) {
    index <- sample.int(n, size = n, replace = TRUE)
    bootstrap_indices[, b] <- index

    bootstrap_X[b, ] <- pca_eigenvalues(
      X[index, , drop = FALSE],
      k = k
    )

    bootstrap_Y[b, ] <- pca_eigenvalues(
      Y[index, , drop = FALSE],
      k = k
    )
  }

  jackknife_X <- t(vapply(
    seq_len(n),
    function(i) {
      pca_eigenvalues(
        X[-i, , drop = FALSE],
        k = k
      )
    },
    numeric(k)
  ))

  percentile_X <- lapply(
    seq_len(k),
    function(j) percentile_interval_local(bootstrap_X[, j], conf = conf)
  )

  percentile_Y <- lapply(
    seq_len(k),
    function(j) percentile_interval_local(bootstrap_Y[, j], conf = conf)
  )

  bca_X <- lapply(
    seq_len(k),
    function(j) {
      bca_interval(
        theta_hat = theta_X[j],
        bootstrap_statistics = bootstrap_X[, j],
        jackknife_statistics = jackknife_X[, j],
        conf = conf
      )
    }
  )

  output <- numeric(0)

  for (j in seq_len(k)) {
    bootstrap_mean_X <- mean(bootstrap_X[, j])
    bootstrap_mean_Y <- mean(bootstrap_Y[, j])

    output <- c(
      output,

      setNames(theta_X[j], paste0("X_sample_l", j)),
      setNames(theta_Y[j], paste0("Y_sample_l", j)),

      setNames(
        theta_X[j] - target_eigenvalues[j],
        paste0("X_sample_bias_l", j)
      ),
      setNames(
        theta_Y[j] - target_eigenvalues[j],
        paste0("Y_sample_bias_l", j)
      ),

      setNames(
        bootstrap_mean_X,
        paste0("X_boot_mean_l", j)
      ),
      setNames(
        bootstrap_mean_Y,
        paste0("Y_boot_mean_l", j)
      ),

      setNames(
        bootstrap_mean_X - theta_X[j],
        paste0("X_cond_boot_bias_l", j)
      ),
      setNames(
        bootstrap_mean_Y - theta_Y[j],
        paste0("Y_cond_boot_bias_l", j)
      ),

      setNames(
        bootstrap_mean_X - target_eigenvalues[j],
        paste0("X_boot_pop_error_l", j)
      ),
      setNames(
        bootstrap_mean_Y - target_eigenvalues[j],
        paste0("Y_boot_pop_error_l", j)
      ),

      setNames(
        inside_interval(target_eigenvalues[j], percentile_X[[j]]),
        paste0("X_pct_cov_l", j)
      ),
      setNames(
        inside_interval(target_eigenvalues[j], percentile_Y[[j]]),
        paste0("Y_pct_cov_l", j)
      ),
      setNames(
        inside_interval(target_eigenvalues[j], bca_X[[j]]),
        paste0("X_bca_cov_l", j)
      ),

      setNames(
        diff(percentile_X[[j]]),
        paste0("X_pct_width_l", j)
      ),
      setNames(
        diff(percentile_Y[[j]]),
        paste0("Y_pct_width_l", j)
      ),
      setNames(
        diff(bca_X[[j]]),
        paste0("X_bca_width_l", j)
      )
    )
  }

  output
}

run_experiment2 <- function(
  n = 25L,
  p = 5L,
  true_eigenvalues = c(10, 5, 2, 1, 0.5),
  B = 2000L,
  M = 1000L,
  k = 2L,
  conf = 0.95,
  master_seed = 1234L
) {
  n <- as.integer(n)
  p <- as.integer(p)
  B <- as.integer(B)
  M <- as.integer(M)
  k <- as.integer(k)

  if (M <= 1L) {
    stop("M must exceed one to calculate MCSEs.", call. = FALSE)
  }

  replication_seeds <- as.integer(
    master_seed + seq_len(M) - 1L
  )

  message(sprintf(
    "Running Experiment 2: n = %d, p = %d, M = %d, B = %d",
    n, p, M, B
  ))

  raw <- do.call(
    rbind,
    lapply(
      replication_seeds,
      one_exp2_replication,
      n = n,
      p = p,
      true_eigenvalues = true_eigenvalues,
      B = B,
      k = k,
      conf = conf
    )
  )

  paired_raw <- cbind(
    Y_minus_X_pct_cov_l1 =
      raw[, "Y_pct_cov_l1"] - raw[, "X_pct_cov_l1"],

    BCa_minus_X_pct_cov_l1 =
      raw[, "X_bca_cov_l1"] - raw[, "X_pct_cov_l1"],

    Y_minus_X_pct_cov_l2 =
      raw[, "Y_pct_cov_l2"] - raw[, "X_pct_cov_l2"],

    BCa_minus_X_pct_cov_l2 =
      raw[, "X_bca_cov_l2"] - raw[, "X_pct_cov_l2"],

    Y_minus_X_pct_width_l1 =
      raw[, "Y_pct_width_l1"] - raw[, "X_pct_width_l1"],

    BCa_minus_X_pct_width_l1 =
      raw[, "X_bca_width_l1"] - raw[, "X_pct_width_l1"],

    Y_minus_X_pct_width_l2 =
      raw[, "Y_pct_width_l2"] - raw[, "X_pct_width_l2"],

    BCa_minus_X_pct_width_l2 =
      raw[, "X_bca_width_l2"] - raw[, "X_pct_width_l2"]
  )

  list(
    n = n,
    p = p,
    true_eigenvalues = true_eigenvalues,
    B = B,
    M = M,
    k = k,
    conf = conf,
    master_seed = master_seed,
    raw = raw,
    paired_raw = paired_raw,
    summaries = summarize_columns(raw),
    paired_summaries = summarize_columns(paired_raw)
  )
}

make_exp2_eigenvalue_table <- function(
  result,
  label = "tab:exp2_eigenvalues"
) {
  s <- result$summaries

  rows <- vapply(
    seq_len(result$k),
    function(j) {
      sprintf(
        "$\\lambda_%d$ & %.1f & %s & %s & %s & %s \\\\",
        j,
        result$true_eigenvalues[j],
        format_mean_mcse(s, paste0("X_sample_l", j)),
        format_mean_mcse(s, paste0("Y_sample_l", j)),
        format_mean_mcse(s, paste0("X_boot_mean_l", j)),
        format_mean_mcse(s, paste0("Y_boot_mean_l", j))
      )
    },
    character(1)
  )

  sprintf(
    paste0(
      "\\begin{table}[htbp]\n",
      "\\centering\n",
      "\\small\n",
      "\\caption{Experiment 2 eigenvalue summaries for $n=%d$, $p=%d$, ",
      "$B=%d$, and $M=%d$. Entries are Monte Carlo means (MCSEs).}\n",
      "\\label{%s}\n",
      "\\begin{tabular}{lccccc}\n",
      "\\hline\\hline\n",
      "\\textbf{Eigenvalue} & \\textbf{Population} & ",
      "\\textbf{Sample: $X$} & \\textbf{Sample: $Y$} & ",
      "\\textbf{Bootstrap Mean: $X$} & ",
      "\\textbf{Bootstrap Mean: $Y$} \\\\\n",
      "\\hline\n",
      "%s\n",
      "\\hline\\hline\n",
      "\\end{tabular}\n",
      "\\end{table}\n"
    ),
    result$n,
    result$p,
    result$B,
    result$M,
    label,
    paste(rows, collapse = "\n")
  )
}

make_exp2_coverage_table <- function(
  result,
  label = "tab:exp2_coverage"
) {
  s <- result$summaries
  d <- result$paired_summaries

  rows <- vapply(
    seq_len(result$k),
    function(j) {
      sprintf(
        "$\\lambda_%d$ & %s & %s & %s & %s & %s \\\\",
        j,
        format_pct_mcse(s, paste0("X_pct_cov_l", j)),
        format_pct_mcse(s, paste0("Y_pct_cov_l", j)),
        format_pct_mcse(s, paste0("X_bca_cov_l", j)),
        format_pct_mcse(d, paste0("Y_minus_X_pct_cov_l", j)),
        format_pct_mcse(d, paste0("BCa_minus_X_pct_cov_l", j))
      )
    },
    character(1)
  )

  sprintf(
    paste0(
      "\\begin{table}[htbp]\n",
      "\\centering\n",
      "\\small\n",
      "\\caption{Experiment 2 empirical 95\\%% interval coverage for ",
      "$n=%d$, $p=%d$, $B=%d$, and $M=%d$. Entries are coverage (MCSE). ",
      "The final two columns are paired coverage differences, in percentage ",
      "points, relative to the ordinary percentile interval based on $X$.}\n",
      "\\label{%s}\n",
      "\\begin{tabular}{lccccc}\n",
      "\\hline\\hline\n",
      "\\textbf{Eigenvalue} & \\textbf{Percentile: $X$} & ",
      "\\textbf{Percentile: $Y$} & \\textbf{BCa: $X$} & ",
      "\\textbf{$Y-X$} & \\textbf{BCa--Percentile: $X$} \\\\\n",
      "\\hline\n",
      "%s\n",
      "\\hline\\hline\n",
      "\\end{tabular}\n",
      "\\end{table}\n"
    ),
    result$n,
    result$p,
    result$B,
    result$M,
    label,
    paste(rows, collapse = "\n")
  )
}

make_exp2_width_table <- function(
  result,
  label = "tab:exp2_widths"
) {
  s <- result$summaries
  d <- result$paired_summaries

  rows <- vapply(
    seq_len(result$k),
    function(j) {
      sprintf(
        "$\\lambda_%d$ & %s & %s & %s & %s & %s \\\\",
        j,
        format_mean_mcse(s, paste0("X_pct_width_l", j)),
        format_mean_mcse(s, paste0("Y_pct_width_l", j)),
        format_mean_mcse(s, paste0("X_bca_width_l", j)),
        format_mean_mcse(
          d,
          paste0("Y_minus_X_pct_width_l", j),
          digits_mean = 3L,
          digits_mcse = 3L
        ),
        format_mean_mcse(
          d,
          paste0("BCa_minus_X_pct_width_l", j),
          digits_mean = 3L,
          digits_mcse = 3L
        )
      )
    },
    character(1)
  )

  sprintf(
    paste0(
      "\\begin{table}[htbp]\n",
      "\\centering\n",
      "\\small\n",
      "\\caption{Experiment 2 mean 95\\%% interval widths for $n=%d$, ",
      "$p=%d$, $B=%d$, and $M=%d$. Entries are Monte Carlo means (MCSEs). ",
      "The final two columns are paired width differences relative to the ",
      "ordinary percentile interval based on $X$; positive values indicate ",
      "wider intervals.}\n",
      "\\label{%s}\n",
      "\\begin{tabular}{lccccc}\n",
      "\\hline\\hline\n",
      "\\textbf{Eigenvalue} & \\textbf{Percentile: $X$} & ",
      "\\textbf{Percentile: $Y$} & \\textbf{BCa: $X$} & ",
      "\\textbf{$Y-X$} & \\textbf{BCa--Percentile: $X$} \\\\\n",
      "\\hline\n",
      "%s\n",
      "\\hline\\hline\n",
      "\\end{tabular}\n",
      "\\end{table}\n"
    ),
    result$n,
    result$p,
    result$B,
    result$M,
    label,
    paste(rows, collapse = "\n")
  )
}

result_exp2 <- run_experiment2(
  n = 25L,
  p = 5L,
  true_eigenvalues = c(10, 5, 2, 1, 0.5),
  B = 2000L,
  M = 1000L,
  k = 2L,
  conf = 0.95,
  master_seed = master_seed
)

experiment2_configuration <- data.frame(
  parameter = c(
    "master_seed",
    "n",
    "p",
    "B",
    "M",
    "k",
    "confidence_level",
    "population_eigenvalues",
    "covariance_denominator",
    "bootstrap_matching"
  ),
  value = c(
    as.character(result_exp2$master_seed),
    as.character(result_exp2$n),
    as.character(result_exp2$p),
    as.character(result_exp2$B),
    as.character(result_exp2$M),
    as.character(result_exp2$k),
    as.character(result_exp2$conf),
    paste(result_exp2$true_eigenvalues, collapse = ", "),
    "n",
    "Matched bootstrap index vectors for X and Y"
  ),
  row.names = NULL
)

latex_eigenvalues <- make_exp2_eigenvalue_table(
  result_exp2,
  label = "tab:exp2_eigenvalues"
)

latex_coverage <- make_exp2_coverage_table(
  result_exp2,
  label = "tab:exp2_coverage"
)

latex_widths <- make_exp2_width_table(
  result_exp2,
  label = "tab:exp2_widths"
)

saveRDS(
  result_exp2,
  file = "results/experiment2_bootstrap_pca.rds"
)

write_summary_csv(
  result_exp2$summaries,
  "results/experiment2_bootstrap_pca_summary.csv"
)

write_summary_csv(
  result_exp2$paired_summaries,
  "results/experiment2_bootstrap_pca_paired_summary.csv"
)

write.csv(
  experiment2_configuration,
  "results/experiment2_bootstrap_pca_configuration.csv",
  row.names = FALSE
)

writeLines(
  latex_eigenvalues,
  "tables/table5_experiment2_eigenvalues.tex"
)

writeLines(
  latex_coverage,
  "tables/table6_experiment2_coverage.tex"
)

writeLines(
  latex_widths,
  "tables/app2_table1_pca_widths.tex"
)

writeLines(
  capture.output(sessionInfo()),
  "results/experiment2_bootstrap_pca_session_info.txt"
)

cat("\n--- Table 5: Experiment 2 Eigenvalue Summaries ---\n\n")
cat(latex_eigenvalues, "\n")

cat("\n--- Table 6: Experiment 2 Interval Coverage ---\n\n")
cat(latex_coverage, "\n")

cat("\n--- Online Appendix 2, Table 1: Interval Widths ---\n\n")
cat(latex_widths, "\n")

message("Experiment 2 completed.")
message("Saved results/experiment2_bootstrap_pca.rds")
message("Saved results/experiment2_bootstrap_pca_summary.csv")
message("Saved results/experiment2_bootstrap_pca_paired_summary.csv")
message("Saved results/experiment2_bootstrap_pca_configuration.csv")
message("Saved results/experiment2_bootstrap_pca_session_info.txt")
message("Saved tables/table5_experiment2_eigenvalues.tex")
message("Saved tables/table6_experiment2_coverage.tex")
message("Saved tables/app2_table1_pca_widths.tex")