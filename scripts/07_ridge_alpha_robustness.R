# ==============================================================================
# Ridge-Shrinkage Robustness Analysis
#
# Online Appendix 2, Table 2:
# Robustness of relative precision-estimation error to the ridge-shrinkage
# parameter alpha.
#
# This script is not the source of main-manuscript Tables 9--10. It evaluates:
#   (n, p) in {(100, 25), (100, 100), (100, 500)}
#   alpha in {0.05, 0.10, 0.20}
#
# Model:
#   X_i ~ N_p(0, Sigma),  Sigma_jk = rho^{|j-k|}, rho = 0.7.
#
# Estimator:
#   Sigma_hat_alpha =
#     (1 - alpha) S_X + alpha {tr(S_X)/p} I_p,
# where S_X is the denominator-n empirical covariance.
#
# Under surrogate calibration:
#   S_Y = {n/(n-1)} S_X
# and therefore:
#   Omega_hat_alpha,Y = {(n-1)/n} Omega_hat_alpha,X.
#
# Outputs:
#   results/ridge_alpha_robustness.rds
#   results/ridge_alpha_robustness_results_M1000.csv
#   results/ridge_alpha_robustness_configuration.csv
#   tables/app2_table2_ridge_robustness.tex
#
# Run from the repository root:
#   source("scripts/07_ridge_alpha_robustness.R")
# ==============================================================================

source("R/common.R")

if (!requireNamespace("parallel", quietly = TRUE)) {
  stop("Package 'parallel' is required.", call. = FALSE)
}

set_reproducible_rng(20260811L)
ensure_directories()

rho <- 0.7
M <- 1000L
master_seed <- 20260811L

alpha_grid <- c(0.05, 0.10, 0.20)

design <- data.frame(
  n = c(100L, 100L, 100L),
  p = c(25L, 100L, 500L)
)

ar1_covariance <- function(p, rho = 0.7) {
  p <- as.integer(p)

  if (length(p) != 1L || !is.finite(p) || p < 1L) {
    stop("p must be a positive integer.", call. = FALSE)
  }

  if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) {
    stop(
      "rho must be a finite scalar strictly between -1 and 1.",
      call. = FALSE
    )
  }

  rho ^ abs(outer(seq_len(p), seq_len(p), "-"))
}

summarize_alpha <- function(x, y, n, p, alpha) {
  difference <- y - x

  data.frame(
    n = as.integer(n),
    p = as.integer(p),
    ratio_p_n = p / n,
    alpha = alpha,
    precision_x = mean(x),
    precision_x_mcse = stats::sd(x) / sqrt(length(x)),
    precision_y = mean(y),
    precision_y_mcse = stats::sd(y) / sqrt(length(y)),
    difference_y_minus_x = mean(difference),
    difference_mcse = stats::sd(difference) / sqrt(length(difference)),
    row.names = NULL
  )
}

one_configuration <- function(seed, n, p, alpha_grid, rho) {
  set.seed(as.integer(seed))

  Sigma_true <- ar1_covariance(p, rho)
  chol_sigma <- chol(Sigma_true)
  Omega_true <- chol2inv(chol_sigma)

  omega_norm_sq <- sum(Omega_true * Omega_true)
  scale_factor <- n / (n - 1L)

  n_alpha <- length(alpha_grid)

  precision_x <- numeric(n_alpha)
  precision_y <- numeric(n_alpha)

  X <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% chol_sigma
  X_c <- sweep(X, 2L, colMeans(X), FUN = "-")
  S_x <- crossprod(X_c) / n

  eig_sx <- eigen(S_x, symmetric = TRUE)
  eigenvalues_x <- pmax(eig_sx$values, 0)
  eigenvectors <- eig_sx$vectors

  mean_eigenvalue_x <- sum(eigenvalues_x) / p
  omega_in_eigenbasis_diag <- colSums(
    eigenvectors * (Omega_true %*% eigenvectors)
  )

  for (a in seq_along(alpha_grid)) {
    alpha <- alpha_grid[a]

    ridge_eigenvalues_x <-
      (1 - alpha) * eigenvalues_x +
      alpha * mean_eigenvalue_x

    if (any(ridge_eigenvalues_x <= 0)) {
      stop(
        "Nonpositive ridge-shrinkage eigenvalue encountered.",
        call. = FALSE
      )
    }

    inverse_eigenvalues_x <- 1 / ridge_eigenvalues_x
    inverse_eigenvalues_y <- inverse_eigenvalues_x / scale_factor

    omega_hat_x_sq <- sum(inverse_eigenvalues_x^2)
    cross_term_x <- sum(
      inverse_eigenvalues_x * omega_in_eigenbasis_diag
    )

    omega_hat_y_sq <- sum(inverse_eigenvalues_y^2)
    cross_term_y <- sum(
      inverse_eigenvalues_y * omega_in_eigenbasis_diag
    )

    error_x_sq <- omega_hat_x_sq +
      omega_norm_sq -
      2 * cross_term_x

    error_y_sq <- omega_hat_y_sq +
      omega_norm_sq -
      2 * cross_term_y

    precision_x[a] <- sqrt(max(error_x_sq, 0) / omega_norm_sq)
    precision_y[a] <- sqrt(max(error_y_sq, 0) / omega_norm_sq)
  }

  list(
    precision_x = precision_x,
    precision_y = precision_y
  )
}

run_configuration_all_alpha <- function(n, p, alpha_grid, rho,
                                        M = 1000L,
                                        master_seed = 20260811L) {
  replication_seeds <- as.integer(
    master_seed + seq_len(M) - 1L
  )

  n_alpha <- length(alpha_grid)

  precision_x <- matrix(
    NA_real_,
    nrow = M,
    ncol = n_alpha
  )

  precision_y <- matrix(
    NA_real_,
    nrow = M,
    ncol = n_alpha
  )

  for (m in seq_len(M)) {
    one_result <- one_configuration(
      seed = replication_seeds[m],
      n = n,
      p = p,
      alpha_grid = alpha_grid,
      rho = rho
    )

    precision_x[m, ] <- one_result$precision_x
    precision_y[m, ] <- one_result$precision_y
  }

  summary_df <- do.call(
    rbind,
    lapply(seq_along(alpha_grid), function(a) {
      summarize_alpha(
        x = precision_x[, a],
        y = precision_y[, a],
        n = n,
        p = p,
        alpha = alpha_grid[a]
      )
    })
  )

  list(
    summary = summary_df,
    precision_x = precision_x,
    precision_y = precision_y,
    replication_seeds = replication_seeds
  )
}

message(
  sprintf(
    "Running %d configurations over %d alpha values with M = %d serial replications.",
    nrow(design),
    length(alpha_grid),
    M
  )
)

configuration_results <- vector("list", nrow(design))

for (i in seq_len(nrow(design))) {
  current_n <- design$n[i]
  current_p <- design$p[i]

  message(
    sprintf(
      "Ridge robustness: n = %d, p = %d, M = %d",
      current_n,
      current_p,
      M
    )
  )

  configuration_results[[i]] <- run_configuration_all_alpha(
    n = current_n,
    p = current_p,
    alpha_grid = alpha_grid,
    rho = rho,
    M = M,
    master_seed = master_seed + 100000L * i
  )
}

ridge_robustness <- do.call(
  rbind,
  lapply(configuration_results, `[[`, "summary")
)

ridge_robustness <- ridge_robustness[
  order(ridge_robustness$p, ridge_robustness$alpha),
  ,
  drop = FALSE
]

print(ridge_robustness, row.names = FALSE)

latex_num <- function(x, digits = 4L) {
  formatC(as.numeric(x), format = "f", digits = digits)
}

latex_mcse <- function(x, digits = 4L, threshold = 0.00005) {
  x <- abs(as.numeric(x))

  if (x < threshold) {
    return(paste0("$<", formatC(threshold, format = "f", digits = 5L), "$"))
  }

  latex_num(x, digits = digits)
}

latex_entry <- function(mean_value, mcse_value, digits = 4L) {
  paste0(
    latex_num(mean_value, digits = digits),
    " (",
    latex_mcse(mcse_value, digits = digits),
    ")"
  )
}

latex_rows <- vapply(
  seq_len(nrow(ridge_robustness)),
  function(i) {
    row <- ridge_robustness[i, ]

    paste0(
      as.integer(row$n), " & ",
      as.integer(row$p), " & ",
      latex_num(row$ratio_p_n, digits = 2L), " & ",
      latex_num(row$alpha, digits = 2L), " & ",
      latex_entry(
        row$precision_x,
        row$precision_x_mcse,
        digits = 4L
      ), " & ",
      latex_entry(
        row$precision_y,
        row$precision_y_mcse,
        digits = 4L
      ), " & ",
      latex_entry(
        row$difference_y_minus_x,
        row$difference_mcse,
        digits = 4L
      ), " \\\\"
    )
  },
  character(1)
)

latex_table <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\small",
  paste0(
    "\\caption{Robustness of relative precision-estimation error to the ",
    "ridge-shrinkage parameter $\\alpha$. Entries are Monte Carlo means ",
    "(MCSEs) over $M=", M, "$ paired replications. Negative $Y-X$ values ",
    "favor the surrogate under relative precision Frobenius loss.}"
  ),
  "\\label{tab:ridge_alpha_robustness}",
  "\\begin{tabular}{rrrrccc}",
  "\\hline\\hline",
  "\\textbf{$n$} & \\textbf{$p$} & \\textbf{$p/n$} & ",
  "\\textbf{$\\alpha$} & \\textbf{$X$} & \\textbf{$Y$} & ",
  "\\textbf{$Y-X$} \\\\",
  "\\hline",
  latex_rows,
  "\\hline\\hline",
  "\\end{tabular}",
  "\\end{table}"
)

saveRDS(
  list(
    configuration_results = configuration_results,
    ridge_robustness = ridge_robustness,
    rho = rho,
    M = M,
    alpha_grid = alpha_grid,
    design = design,
    master_seed = master_seed
  ),
  file = "results/ridge_alpha_robustness.rds"
)

write.csv(
  ridge_robustness,
  file = "results/ridge_alpha_robustness_results_M1000.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    parameter = c(
      "master_seed",
      "rho",
      "M",
      "alpha_grid",
      "design"
    ),
    value = c(
      as.character(master_seed),
      as.character(rho),
      as.character(M),
      paste(alpha_grid, collapse = ", "),
      paste(
        paste0("(n=", design$n, ", p=", design$p, ")"),
        collapse = "; "
      )
    )
  ),
  file = "results/ridge_alpha_robustness_configuration.csv",
  row.names = FALSE
)

writeLines(
  latex_table,
  con = "tables/app2_table2_ridge_robustness.tex"
)

writeLines(
  capture.output(sessionInfo()),
  con = "results/ridge_alpha_robustness_session_info.txt"
)

message("Ridge-shrinkage robustness analysis completed.")
message("Saved results/ridge_alpha_robustness.rds")
message("Saved results/ridge_alpha_robustness_results_M1000.csv")
message("Saved results/ridge_alpha_robustness_configuration.csv")
message("Saved results/ridge_alpha_robustness_session_info.txt")
message("Saved tables/app2_table2_ridge_robustness.tex")