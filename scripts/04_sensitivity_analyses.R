## ==============================================================================
## Sample-Size Sensitivity: Covariance Estimation and Bootstrap PCA
## High-Performance, Cross-Platform R Script Optimized for RStudio / macOS & Windows
## ==============================================================================

library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)

n_grid <- c(10L, 20L, 30L, 50L, 100L, 250L)

## 1. Core transformation and utilities
get_surrogate_population <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n <= 1L) stop("Sample size n must be greater than 1.")
  x_bar <- colMeans(X)
  sqrt(n / (n - 1)) * sweep(X, 2L, x_bar, "-") +
    matrix(x_bar, nrow = n, ncol = ncol(X), byrow = TRUE)
}

cov_ml <- function(X) {
  X_c <- sweep(as.matrix(X), 2L, colMeans(X), "-")
  crossprod(X_c) / nrow(X)
}

make_ar1_cov <- function(p, rho) {
  rho ^ abs(outer(seq_len(p), seq_len(p), "-"))
}

mc_summary <- function(x) {
  n_eff <- sum(is.finite(x))
  c(mean = mean(x, na.rm = TRUE), mcse = sd(x, na.rm = TRUE) / sqrt(n_eff))
}

summarize_matrix <- function(raw, metadata) {
  summary <- t(apply(raw, 2L, mc_summary))
  data.frame(
    metadata,
    metric = rownames(summary),
    mean = summary[, "mean"],
    mcse = summary[, "mcse"],
    row.names = NULL,
    check.names = FALSE
  )
}

## 2. Experiment 1: Optimized Shrinkage & Trace Identity
shrinkage_covariance <- function(X, estimator = c("LW", "OAS"), ridge = 1e-4) {
  estimator <- match.arg(estimator)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  X_c <- sweep(X, 2L, colMeans(X), "-")
  S_emp <- crossprod(X_c) / n
  mu <- sum(diag(S_emp)) / p
  target <- mu * diag(p)
  d2 <- sum((S_emp - target)^2)
  
  if (d2 <= .Machine$double.eps) {
    shrinkage <- 1
  } else if (estimator == "LW") {
    row_norm4 <- sum(rowSums(X_c^2)^2)
    quadratic <- sum(rowSums((X_c %*% S_emp) * X_c))
    b_bar2 <- (row_norm4 - 2 * quadratic + n * sum(S_emp^2)) / n^2
    shrinkage <- min(1, max(0, b_bar2 / d2))
  } else {
    tr_s <- sum(diag(S_emp))
    tr_s2 <- sum(S_emp^2)
    numerator <- (1 - 2 / p) * tr_s2 + tr_s^2
    denominator <- (n + 1 - 2 / p) * (tr_s2 - tr_s^2 / p)
    shrinkage <- if (denominator <= .Machine$double.eps) 1 else {
      min(1, max(0, numerator / denominator))
    }
  }
  
  Sigma_hat <- (1 - shrinkage) * S_emp + shrinkage * target
  Omega_hat <- chol2inv(chol(Sigma_hat + ridge * diag(p)))
  list(Sigma = Sigma_hat, Precision = Omega_hat)
}

exp1_single <- function(n, p, chol_sigma, Omega_true) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% chol_sigma
  Y <- get_surrogate_population(X)
  
  S_X <- cov_ml(X)
  S_Y <- cov_ml(Y)
  target_trace <- p
  
  lw_X <- shrinkage_covariance(X, "LW")
  lw_Y <- shrinkage_covariance(Y, "LW")
  oas_X <- shrinkage_covariance(X, "OAS")
  oas_Y <- shrinkage_covariance(Y, "OAS")
  
  true_cov <- crossprod(chol_sigma)
  
  c(
    raw_trace_X = sum(diag(S_X)),
    raw_trace_Y = sum(diag(S_Y)),
    raw_trace_bias_X = sum(diag(S_X)) - target_trace,
    raw_trace_bias_Y = sum(diag(S_Y)) - target_trace,
    raw_trace_ratio_Y_X = sum(diag(S_Y)) / sum(diag(S_X)),
    LW_frob_cov_X = norm(lw_X$Sigma - true_cov, "F"),
    LW_frob_cov_Y = norm(lw_Y$Sigma - true_cov, "F"),
    LW_frob_prec_X = norm(lw_X$Precision - Omega_true, "F"),
    LW_frob_prec_Y = norm(lw_Y$Precision - Omega_true, "F"),
    LW_trace_bias_X = sum(diag(lw_X$Sigma)) - target_trace,
    LW_trace_bias_Y = sum(diag(lw_Y$Sigma)) - target_trace,
    OAS_frob_cov_X = norm(oas_X$Sigma - true_cov, "F"),
    OAS_frob_cov_Y = norm(oas_Y$Sigma - true_cov, "F"),
    OAS_frob_prec_X = norm(oas_X$Precision - Omega_true, "F"),
    OAS_frob_prec_Y = norm(oas_Y$Precision - Omega_true, "F"),
    OAS_trace_bias_X = sum(diag(oas_X$Sigma)) - target_trace,
    OAS_trace_bias_Y = sum(diag(oas_Y$Sigma)) - target_trace
  )
}

run_exp1_sample_size <- function(ns = n_grid, p = 100L, rho = 0.7,
                                 M = 1000L, ncores = max(1L, detectCores() - 1L)) {
  Sigma_true <- make_ar1_cov(p, rho)
  chol_sigma <- chol(Sigma_true)
  Omega_true <- chol2inv(chol(Sigma_true))
  output <- vector("list", length(ns))
  
  windows_cluster <- NULL
  if (.Platform$OS.type == "windows" && ncores > 1L) {
    windows_cluster <- makeCluster(ncores)
    on.exit(stopCluster(windows_cluster), add = TRUE)
    clusterSetRNGStream(windows_cluster, iseed = 2026)
    clusterExport(windows_cluster,
                  c("get_surrogate_population", "cov_ml", "shrinkage_covariance", "exp1_single"),
                  envir = environment())
  }
  
  for (ii in seq_along(ns)) {
    n <- ns[ii]
    cat(sprintf("Experiment 1 sensitivity: n = %d, p = %d, M = %d using %d core(s)...\n", 
                n, p, M, ncores))
    one_rep <- function(i) exp1_single(n, p, chol_sigma, Omega_true)
    
    if (!is.null(windows_cluster)) {
      reps <- parLapply(windows_cluster, seq_len(M), one_rep)
    } else if (.Platform$OS.type != "windows" && ncores > 1L) {
      reps <- mclapply(seq_len(M), one_rep, mc.cores = ncores, mc.set.seed = TRUE)
    } else {
      reps <- lapply(seq_len(M), one_rep)
    }
    
    raw <- do.call(rbind, reps)
    paired <- cbind(
      LW_frob_cov_Y_minus_X = raw[, "LW_frob_cov_Y"] - raw[, "LW_frob_cov_X"],
      LW_frob_prec_Y_minus_X = raw[, "LW_frob_prec_Y"] - raw[, "LW_frob_prec_X"],
      LW_trace_bias_Y_minus_X = raw[, "LW_trace_bias_Y"] - raw[, "LW_trace_bias_X"],
      OAS_frob_cov_Y_minus_X = raw[, "OAS_frob_cov_Y"] - raw[, "OAS_frob_cov_X"],
      OAS_frob_prec_Y_minus_X = raw[, "OAS_frob_prec_Y"] - raw[, "OAS_frob_prec_X"],
      OAS_trace_bias_Y_minus_X = raw[, "OAS_trace_bias_Y"] - raw[, "OAS_trace_bias_X"]
    )
    metadata <- data.frame(
      experiment = "covariance_sensitivity", n = n, p = p, rho = rho, M = M,
      theoretical_factor = n / (n - 1)
    )
    output[[ii]] <- rbind(
      transform(summarize_matrix(raw, metadata), panel = "raw_and_estimators"),
      transform(summarize_matrix(paired, metadata), panel = "paired_Y_minus_X")
    )
  }
  do.call(rbind, output)
}

## 3. Experiment 2: PCA scaling and percentile-bootstrap coverage sensitivity
random_covariance <- function(true_evals) {
  p <- length(true_evals)
  Q <- qr.Q(qr(matrix(rnorm(p * p), nrow = p, ncol = p)))
  list(Sigma = Q %*% diag(true_evals) %*% t(Q), evals = sort(true_evals, decreasing = TRUE))
}

ordered_eigenvalues_ml <- function(X, k = 2L) {
  eigen(cov_ml(X), symmetric = TRUE, only.values = TRUE)$values[seq_len(k)]
}

percentile_interval <- function(x, alpha = 0.05) {
  as.numeric(quantile(x, probs = c(alpha / 2, 1 - alpha / 2),
                      names = FALSE, type = 7))
}

inside_interval <- function(value, interval) {
  as.numeric(value >= interval[1] && value <= interval[2])
}

exp2_single <- function(n, chol_sigma, target_evals, B = 2000L,
                        k = 2L, alpha = 0.05) {
  p <- length(target_evals)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% chol_sigma
  Y <- get_surrogate_population(X)
  target <- target_evals[seq_len(k)]
  
  eig_X <- ordered_eigenvalues_ml(X, k)
  eig_Y <- ordered_eigenvalues_ml(Y, k)
  
  boot_X <- matrix(NA_real_, B, k)
  boot_Y <- matrix(NA_real_, B, k)
  
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    boot_X[b, ] <- ordered_eigenvalues_ml(X[idx, , drop = FALSE], k)
    boot_Y[b, ] <- ordered_eigenvalues_ml(Y[idx, , drop = FALSE], k)
  }
  
  output <- numeric(0)
  for (j in seq_len(k)) {
    ci_X <- percentile_interval(boot_X[, j], alpha)
    ci_Y <- percentile_interval(boot_Y[, j], alpha)
    cov_X <- inside_interval(target[j], ci_X)
    cov_Y <- inside_interval(target[j], ci_Y)
    
    output <- c(
      output,
      setNames(eig_X[j], paste0("sample_eig_X_l", j)),
      setNames(eig_Y[j], paste0("sample_eig_Y_l", j)),
      setNames(eig_Y[j] / eig_X[j], paste0("sample_ratio_Y_X_l", j)),
      setNames(mean(boot_X[, j]), paste0("boot_mean_X_l", j)),
      setNames(mean(boot_Y[, j]), paste0("boot_mean_Y_l", j)),
      setNames(mean(boot_Y[, j]) / mean(boot_X[, j]), paste0("boot_ratio_Y_X_l", j)),
      setNames(cov_X, paste0("pct_cov_X_l", j)),
      setNames(cov_Y, paste0("pct_cov_Y_l", j)),
      setNames(cov_Y - cov_X, paste0("pct_cov_diff_Y_minus_X_l", j)),
      setNames(diff(ci_X), paste0("pct_width_X_l", j)),
      setNames(diff(ci_Y), paste0("pct_width_Y_l", j)),
      setNames(diff(ci_Y) / diff(ci_X), paste0("pct_width_ratio_Y_X_l", j))
    )
  }
  output
}

run_exp2_sample_size <- function(ns = n_grid,
                                 true_evals = c(10, 5, 2, 1, 0.5),
                                 B = 2000L, M = 1000L, ncores = max(1L, detectCores() - 1L),
                                 k = 2L, alpha = 0.05) {
  spec <- random_covariance(true_evals)
  chol_sigma <- chol(spec$Sigma)
  target_evals <- spec$evals
  output <- vector("list", length(ns))
  
  windows_cluster <- NULL
  if (.Platform$OS.type == "windows" && ncores > 1L) {
    windows_cluster <- makeCluster(ncores)
    on.exit(stopCluster(windows_cluster), add = TRUE)
    clusterSetRNGStream(windows_cluster, iseed = 3026)
    clusterExport(windows_cluster,
                  c("get_surrogate_population", "cov_ml", "ordered_eigenvalues_ml",
                    "percentile_interval", "inside_interval", "exp2_single"),
                  envir = environment())
  }
  
  for (ii in seq_along(ns)) {
    n <- ns[ii]
    cat(sprintf("Experiment 2 sensitivity: n = %d, p = %d, M = %d, B = %d using %d core(s)...\n",
                n, length(target_evals), M, B, ncores))
    one_rep <- function(i) exp2_single(n, chol_sigma, target_evals, B, k, alpha)
    
    if (!is.null(windows_cluster)) {
      reps <- parLapply(windows_cluster, seq_len(M), one_rep)
    } else if (.Platform$OS.type != "windows" && ncores > 1L) {
      reps <- mclapply(seq_len(M), one_rep, mc.cores = ncores, mc.set.seed = TRUE)
    } else {
      reps <- lapply(seq_len(M), one_rep)
    }
    
    raw <- do.call(rbind, reps)
    metadata <- data.frame(
      experiment = "pca_sensitivity", n = n, p = length(target_evals), B = B, M = M,
      theoretical_factor = n / (n - 1)
    )
    output[[ii]] <- summarize_matrix(raw, metadata)
  }
  do.call(rbind, output)
}

## 4. LaTeX table helpers
format_mean_mcse <- function(mean, mcse, digits_mean = 4L, digits_mcse = 4L) {
  sprintf(paste0("%.", digits_mean, "f (%.", digits_mcse, "f)"), mean, mcse)
}

get_stat <- function(df, n, metric, panel = NULL) {
  row <- df[df$n == n & df$metric == metric, , drop = FALSE]
  if (!is.null(panel)) row <- row[row$panel == panel, , drop = FALSE]
  if (nrow(row) != 1L) stop("Missing or duplicate statistic: ", metric, " at n = ", n)
  row
}

make_exp1_trace_latex <- function(df, label = "tab:sample_size_trace") {
  ns <- sort(unique(df$n))
  rows <- vapply(ns, function(n) {
    factor <- n / (n - 1)
    x_bias <- get_stat(df, n, "raw_trace_bias_X", "raw_and_estimators")
    y_bias <- get_stat(df, n, "raw_trace_bias_Y", "raw_and_estimators")
    ratio <- get_stat(df, n, "raw_trace_ratio_Y_X", "raw_and_estimators")
    sprintf("%d & %.5f & %s & %s & %s \\\\", n, factor,
            format_mean_mcse(x_bias$mean, x_bias$mcse, 3L, 3L),
            format_mean_mcse(y_bias$mean, y_bias$mcse, 3L, 3L),
            format_mean_mcse(ratio$mean, ratio$mcse, 5L, 5L))
  }, character(1))
  
  sprintf("\\begin{table}[htbp]\\centering\\small\\caption{Sample-size verification of the exact empirical covariance-scale identity in the AR(1) design. Entries are Monte Carlo means (MCSEs).}\\label{%s}\\begin{tabular}{rcccc}\\hline\\hline$n$ & $n/(n-1)$ & Trace bias: $X$ & Trace bias: $Y$ & $\\operatorname{tr}(S_Y)/\\operatorname{tr}(S_X)$ \\\\\\hline%s\\hline\\hline\\end{tabular}\\end{table}", label, paste(rows, collapse = "\n"))
}

make_exp2_pca_latex <- function(df, label = "tab:sample_size_pca") {
  ns <- sort(unique(df$n))
  rows <- vapply(ns, function(n) {
    factor <- n / (n - 1)
    ratio_1 <- get_stat(df, n, "sample_ratio_Y_X_l1")
    ratio_2 <- get_stat(df, n, "sample_ratio_Y_X_l2")
    cov_x <- get_stat(df, n, "pct_cov_X_l2")
    cov_y <- get_stat(df, n, "pct_cov_Y_l2")
    cov_diff <- get_stat(df, n, "pct_cov_diff_Y_minus_X_l2")
    sprintf("%d & %.5f & %s & %s & %.1f\\%% (%.1f\\%%) & %.1f\\%% (%.1f\\%%) & %.1f\\%% (%.1f\\%%) \\\\",
            n, factor,
            format_mean_mcse(ratio_1$mean, ratio_1$mcse, 5L, 5L),
            format_mean_mcse(ratio_2$mean, ratio_2$mcse, 5L, 5L),
            100 * cov_x$mean, 100 * cov_x$mcse,
            100 * cov_y$mean, 100 * cov_y$mcse,
            100 * cov_diff$mean, 100 * cov_diff$mcse)
  }, character(1))
  
  sprintf("\\begin{table}[htbp]\\centering\\scriptsize\\caption{Sample-size sensitivity of bootstrap PCA. The surrogate-to-standard eigenvalue ratios verify the exact factor $n/(n-1)$. Coverage is shown for the second population eigenvalue, $\\lambda_2=5$; the final column is the paired coverage difference.}\\label{%s}\\begin{tabular}{rcccccc}\\hline\\hline$n$ & $n/(n-1)$ & $\\lambda_1$: $Y/X$ & $\\lambda_2$: $Y/X$ & Coverage $X$ & Coverage $Y$ & $Y-X$ \\\\\\hline%s\\hline\\hline\\end{tabular}\\end{table}", label, paste(rows, collapse = "\n"))
}

## ==============================================================================
## Execution in RStudio (Multi-core parallel enabled)
## ==============================================================================
M_final <- 1000L
B_final <- 2000L
ncores  <- max(1L, detectCores() - 1L)

exp1_sensitivity <- run_exp1_sample_size(
  ns = n_grid, p = 100L, rho = 0.7, M = M_final, ncores = ncores
)

exp2_sensitivity <- run_exp2_sample_size(
  ns = n_grid, true_evals = c(10, 5, 2, 1, 0.5),
  B = B_final, M = M_final, ncores = ncores, k = 2L, alpha = 0.05
)

saveRDS(
  list(
    covariance_sample_size = exp1_sensitivity,
    pca_sample_size = exp2_sensitivity
  ),
  file = "results/sensitivity_analyses.rds"
)

write.csv(
  exp1_sensitivity,
  "results/sample_size_sensitivity_exp1.csv",
  row.names = FALSE
)

write.csv(
  exp2_sensitivity,
  "results/sample_size_sensitivity_exp2.csv",
  row.names = FALSE
)

latex_exp1_trace <- make_exp1_trace_latex(
  exp1_sensitivity,
  label = "tab:sample_size_trace"
)

latex_exp2_pca <- make_exp2_pca_latex(
  exp2_sensitivity,
  label = "tab:sample_size_pca"
)

writeLines(
  latex_exp1_trace,
  "tables/table7_sample_size_trace.tex"
)

writeLines(
  latex_exp2_pca,
  "tables/table8_sample_size_pca.tex"
)

message("Saved results/sensitivity_analyses.rds")
message("Saved results/sample_size_sensitivity_exp1.csv")
message("Saved results/sample_size_sensitivity_exp2.csv")
message("Saved tables/table7_sample_size_trace.tex")
message("Saved tables/table8_sample_size_pca.tex")


cat("\n--- Experiment 1 trace-scale table ---\n\n")
cat(latex_exp1_trace, "\n\n")
cat("\n--- Experiment 2 PCA sensitivity table ---\n\n")
cat(latex_exp2_pca, "\n")