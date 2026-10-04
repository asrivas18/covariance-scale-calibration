options(stringsAsFactors = FALSE)

set_reproducible_rng <- function(seed = 2026L) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(seed))
  invisible(seed)
}

ensure_directories <- function() {
  dir.create("results", showWarnings = FALSE, recursive = TRUE)
  dir.create("tables", showWarnings = FALSE, recursive = TRUE)
}

surrogate_population <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)

  if (n <= 1L) {
    stop("Sample size n must be greater than 1.", call. = FALSE)
  }

  x_bar <- colMeans(X)
  s <- sqrt(n / (n - 1L))

  s * sweep(X, 2L, x_bar, FUN = "-") +
    matrix(x_bar, nrow = n, ncol = ncol(X), byrow = TRUE)
}

transform_with_training_map <- function(X_new, x_bar_train, n_train) {
  X_new <- as.matrix(X_new)
  x_bar_train <- as.numeric(x_bar_train)

  if (ncol(X_new) != length(x_bar_train)) {
    stop("X_new and x_bar_train have incompatible dimensions.", call. = FALSE)
  }

  if (n_train <= 1L) {
    stop("n_train must be greater than 1.", call. = FALSE)
  }

  s <- sqrt(n_train / (n_train - 1L))

  s * sweep(X_new, 2L, x_bar_train, FUN = "-") +
    matrix(
      x_bar_train,
      nrow = nrow(X_new),
      ncol = ncol(X_new),
      byrow = TRUE
    )
}

empirical_covariance_n <- function(X) {
  X <- as.matrix(X)

  if (nrow(X) <= 1L) {
    stop("At least two observations are required.", call. = FALSE)
  }

  X_c <- sweep(X, 2L, colMeans(X), FUN = "-")
  crossprod(X_c) / nrow(X)
}

mc_summary <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    return(c(mean = NA_real_, mcse = NA_real_))
  }

  c(
    mean = mean(x),
    mcse = stats::sd(x) / sqrt(length(x))
  )
}

summarize_columns <- function(X) {
  out <- t(apply(as.matrix(X), 2L, mc_summary))
  colnames(out) <- c("mean", "mcse")
  out
}

format_mean_mcse <- function(summary_matrix, name,
                             digits_mean = 2L,
                             digits_mcse = 2L) {
  value <- summary_matrix[name, "mean"]
  se <- summary_matrix[name, "mcse"]

  if (!is.finite(value) || !is.finite(se)) {
    return("NA")
  }

  sprintf(
    paste0("%.", digits_mean, "f (%.", digits_mcse, "f)"),
    value,
    se
  )
}

format_pct_mcse <- function(summary_matrix, name,
                            digits_pct = 1L,
                            digits_mcse = 1L) {
  value <- summary_matrix[name, "mean"]
  se <- summary_matrix[name, "mcse"]

  if (!is.finite(value) || !is.finite(se)) {
    return("NA")
  }

  sprintf(
    paste0("%.", digits_pct, "f\\\\%% (%.", digits_mcse, "f\\\\%%)"),
    100 * value,
    100 * se
  )
}

write_summary_csv <- function(summary_matrix, path) {
  utils::write.csv(
    data.frame(
      statistic = rownames(summary_matrix),
      mean = summary_matrix[, "mean"],
      mcse = summary_matrix[, "mcse"],
      row.names = NULL
    ),
    file = path,
    row.names = FALSE
  )
}

make_ar1_covariance <- function(p, rho = 0.7) {
  rho ^ abs(outer(seq_len(p), seq_len(p), "-"))
}

ordered_pca_eigenvalues <- function(X, k = 2L) {
  values <- eigen(
    empirical_covariance_n(X),
    symmetric = TRUE,
    only.values = TRUE
  )$values

  values[seq_len(k)]
}

percentile_interval <- function(values, alpha = 0.05) {
  as.numeric(
    stats::quantile(
      values,
      probs = c(alpha / 2, 1 - alpha / 2),
      names = FALSE,
      type = 7
    )
  )
}

inside_interval <- function(value, interval) {
  interval[1] <= value && value <= interval[2]
}