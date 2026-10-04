source("R/common.R")

set_reproducible_rng(2026L)
ensure_directories()

fit_ridge <- function(X, y, lambda = 1) {
  X <- as.matrix(X)
  y <- as.numeric(y)

  if (!is.numeric(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda <= 0) {
    stop("lambda must be a positive finite scalar.", call. = FALSE)
  }

  if (nrow(X) != length(y)) {
    stop("X and y have incompatible dimensions.", call. = FALSE)
  }

  x_bar <- colMeans(X)
  y_bar <- mean(y)
  X_c <- sweep(X, 2L, x_bar, FUN = "-")

  ridge_matrix <- crossprod(X_c) + lambda * diag(ncol(X))

  R <- tryCatch(
    chol(ridge_matrix),
    error = function(e) {
      stop(
        "Ridge normal-equation matrix is not numerically positive definite.",
        call. = FALSE
      )
    }
  )

  beta_hat <- backsolve(
    R,
    forwardsolve(t(R), crossprod(X_c, y - y_bar))
  )

  list(beta = beta_hat, x_bar = x_bar, y_bar = y_bar)
}

predict_ridge <- function(fit, X_new) {
  X_new <- as.matrix(X_new)

  as.vector(
    fit$y_bar +
      sweep(X_new, 2L, fit$x_bar, FUN = "-") %*% fit$beta
  )
}

predict_knn <- function(X_train, y_train, X_test, k = 3L) {
  X_train <- as.matrix(X_train)
  X_test <- as.matrix(X_test)
  y_train <- as.numeric(y_train)

  if (nrow(X_train) != length(y_train)) {
    stop("X_train and y_train have incompatible dimensions.", call. = FALSE)
  }

  if (ncol(X_train) != ncol(X_test)) {
    stop("Training and test matrices have incompatible dimensions.", call. = FALSE)
  }

  k <- min(as.integer(k), nrow(X_train))

  if (k < 1L) {
    stop("k must be at least one.", call. = FALSE)
  }

  d2 <- outer(rowSums(X_test^2), rowSums(X_train^2), "+") -
    2 * tcrossprod(X_test, X_train)

  vapply(
    seq_len(nrow(X_test)),
    function(i) {
      neighbors <- order(d2[i, ])[seq_len(k)]
      mean(y_train[neighbors])
    },
    numeric(1)
  )
}

generate_exp3_data <- function(n = 30L, p = 10L,
                               noise_sd = 0.5, n_test = 100L) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)

  beta <- runif(p, min = 1, max = 3)
  beta <- beta / sqrt(sum(beta^2))

  y <- as.vector(X %*% beta + rnorm(n, sd = noise_sd))

  X_test <- matrix(rnorm(n_test * p), nrow = n_test, ncol = p)
  y_test <- as.vector(X_test %*% beta + rnorm(n_test, sd = noise_sd))

  list(X = X, y = y, X_test = X_test, y_test = y_test)
}

bootstrap_interval <- function(prediction_draws, alpha = 0.05) {
  interval <- t(
    apply(
      prediction_draws,
      2L,
      stats::quantile,
      probs = c(alpha / 2, 1 - alpha / 2),
      names = FALSE,
      type = 7
    )
  )

  colnames(interval) <- c("lower", "upper")
  interval
}

widen_interval <- function(interval, center, factor) {
  half_width <- 0.5 * (interval[, "upper"] - interval[, "lower"])

  cbind(
    lower = center - factor * half_width,
    upper = center + factor * half_width
  )
}

mean_interval_width <- function(interval) {
  mean(interval[, "upper"] - interval[, "lower"])
}

mean_coverage <- function(interval, y) {
  y <- as.numeric(y)

  if (nrow(interval) != length(y)) {
    stop("Interval rows and outcome length do not match.", call. = FALSE)
  }

  mean(interval[, "lower"] <= y & y <= interval[, "upper"])
}

conformal_quantile <- function(scores, alpha = 0.05) {
  scores <- abs(as.numeric(scores))
  m <- length(scores)

  if (m < 1L) {
    stop("At least one calibration score is required.", call. = FALSE)
  }

  k_order <- min(m, ceiling((m + 1) * (1 - alpha)))
  sort(scores, partial = k_order)[k_order]
}

split_conformal <- function(X, y, X_test,
                            learner = c("ridge", "knn"),
                            lambda = 1,
                            k = 3L,
                            alpha = 0.05,
                            train_fraction = 0.7) {
  learner <- match.arg(learner)

  n <- nrow(X)
  n_fit <- floor(train_fraction * n)
  n_fit <- max(2L, min(n - 1L, n_fit))

  fit_index <- sample.int(n, size = n_fit, replace = FALSE)
  calibration_index <- setdiff(seq_len(n), fit_index)

  X_fit <- X[fit_index, , drop = FALSE]
  y_fit <- y[fit_index]
  X_cal <- X[calibration_index, , drop = FALSE]
  y_cal <- y[calibration_index]

  if (learner == "ridge") {
    model <- fit_ridge(X_fit, y_fit, lambda = lambda)
    prediction_cal <- predict_ridge(model, X_cal)
    prediction_test <- predict_ridge(model, X_test)
  } else {
    prediction_cal <- predict_knn(X_fit, y_fit, X_cal, k = k)
    prediction_test <- predict_knn(X_fit, y_fit, X_test, k = k)
  }

  q_hat <- conformal_quantile(y_cal - prediction_cal, alpha = alpha)

  cbind(
    lower = prediction_test - q_hat,
    upper = prediction_test + q_hat,
    prediction = prediction_test
  )
}

one_exp3_replication <- function(seed,
                                 n = 30L,
                                 p = 10L,
                                 B = 500L,
                                 lambda = 1,
                                 k = 3L,
                                 noise_sd = 0.5,
                                 n_test = 100L,
                                 alpha = 0.05,
                                 conformal_train_fraction = 0.7) {
  set.seed(seed)

  dat <- generate_exp3_data(
    n = n,
    p = p,
    noise_sd = noise_sd,
    n_test = n_test
  )

  X <- dat$X
  y <- dat$y
  X_test <- dat$X_test
  y_test <- dat$y_test

  Y <- surrogate_population(X)
  X_test_Y <- transform_with_training_map(
    X_new = X_test,
    x_bar_train = colMeans(X),
    n_train = n
  )

  ridge_draws_X <- matrix(NA_real_, nrow = B, ncol = n_test)
  ridge_draws_Y <- matrix(NA_real_, nrow = B, ncol = n_test)
  knn_draws_X <- matrix(NA_real_, nrow = B, ncol = n_test)
  knn_draws_Y <- matrix(NA_real_, nrow = B, ncol = n_test)

  for (b in seq_len(B)) {
    index <- sample.int(n, size = n, replace = TRUE)

    X_b <- X[index, , drop = FALSE]
    Y_b <- Y[index, , drop = FALSE]
    y_b <- y[index]

    ridge_X <- fit_ridge(X_b, y_b, lambda = lambda)
    ridge_Y <- fit_ridge(Y_b, y_b, lambda = lambda)

    ridge_draws_X[b, ] <- predict_ridge(ridge_X, X_test)
    ridge_draws_Y[b, ] <- predict_ridge(ridge_Y, X_test_Y)

    knn_draws_X[b, ] <- predict_knn(X_b, y_b, X_test, k = k)
    knn_draws_Y[b, ] <- predict_knn(Y_b, y_b, X_test_Y, k = k)
  }

  ci_ridge_X <- bootstrap_interval(ridge_draws_X, alpha = alpha)
  ci_ridge_Y <- bootstrap_interval(ridge_draws_Y, alpha = alpha)
  ci_knn_X <- bootstrap_interval(knn_draws_X, alpha = alpha)
  ci_knn_Y <- bootstrap_interval(knn_draws_Y, alpha = alpha)

  scale_factor <- sqrt(n / (n - 1L))

  ci_ridge_post <- widen_interval(
    ci_ridge_X,
    center = colMeans(ridge_draws_X),
    factor = scale_factor
  )

  ci_knn_post <- widen_interval(
    ci_knn_X,
    center = colMeans(knn_draws_X),
    factor = scale_factor
  )

  ci_ridge_conformal <- split_conformal(
    X, y, X_test,
    learner = "ridge",
    lambda = lambda,
    k = k,
    alpha = alpha,
    train_fraction = conformal_train_fraction
  )

  ci_knn_conformal <- split_conformal(
    X, y, X_test,
    learner = "knn",
    lambda = lambda,
    k = k,
    alpha = alpha,
    train_fraction = conformal_train_fraction
  )

  ridge_mean_X <- colMeans(ridge_draws_X)
  ridge_mean_Y <- colMeans(ridge_draws_Y)
  knn_mean_X <- colMeans(knn_draws_X)
  knn_mean_Y <- colMeans(knn_draws_Y)

  c(
    ridge_X_mse = mean((ridge_mean_X - y_test)^2),
    ridge_Y_mse = mean((ridge_mean_Y - y_test)^2),
    ridge_conformal_mse = mean((ci_ridge_conformal[, "prediction"] - y_test)^2),

    ridge_X_width = mean_interval_width(ci_ridge_X),
    ridge_Y_width = mean_interval_width(ci_ridge_Y),
    ridge_X_post_width = mean_interval_width(ci_ridge_post),
    ridge_conformal_width = mean_interval_width(ci_ridge_conformal),

    ridge_X_coverage = mean_coverage(ci_ridge_X, y_test),
    ridge_Y_coverage = mean_coverage(ci_ridge_Y, y_test),
    ridge_X_post_coverage = mean_coverage(ci_ridge_post, y_test),
    ridge_conformal_coverage = mean_coverage(ci_ridge_conformal, y_test),

    knn_X_mse = mean((knn_mean_X - y_test)^2),
    knn_Y_mse = mean((knn_mean_Y - y_test)^2),
    knn_conformal_mse = mean((ci_knn_conformal[, "prediction"] - y_test)^2),

    knn_X_width = mean_interval_width(ci_knn_X),
    knn_Y_width = mean_interval_width(ci_knn_Y),
    knn_X_post_width = mean_interval_width(ci_knn_post),
    knn_conformal_width = mean_interval_width(ci_knn_conformal),

    knn_X_coverage = mean_coverage(ci_knn_X, y_test),
    knn_Y_coverage = mean_coverage(ci_knn_Y, y_test),
    knn_X_post_coverage = mean_coverage(ci_knn_post, y_test),
    knn_conformal_coverage = mean_coverage(ci_knn_conformal, y_test)
  )
}

run_exp3 <- function(M = 500L, seed = 2026L,
                     n = 30L, p = 10L, B = 500L,
                     lambda = 1, k = 3L,
                     noise_sd = 0.5, n_test = 100L,
                     alpha = 0.05,
                     conformal_train_fraction = 0.7) {
  replication_seeds <- seed + seq_len(M) - 1L

  raw <- do.call(
    rbind,
    lapply(
      replication_seeds,
      one_exp3_replication,
      n = n,
      p = p,
      B = B,
      lambda = lambda,
      k = k,
      noise_sd = noise_sd,
      n_test = n_test,
      alpha = alpha,
      conformal_train_fraction = conformal_train_fraction
    )
  )

  list(
    raw = raw,
    summaries = summarize_columns(raw),
    n = n,
    p = p,
    B = B,
    M = M,
    n_test = n_test,
    lambda = lambda,
    k = k,
    noise_sd = noise_sd,
    alpha = alpha,
    conformal_train_fraction = conformal_train_fraction
  )
}

make_exp3_main_table <- function(res, label = "tab:exp3results") {
  s <- res$summaries

  row_block <- function(prefix, display_name) {
    paste0(
      display_name, " & Test MSE & ",
      format_mean_mcse(s, paste0(prefix, "_X_mse")), " & ",
      format_mean_mcse(s, paste0(prefix, "_Y_mse")), " & -- \\\\\n",
      " & 95\\% interval width & ",
      format_mean_mcse(s, paste0(prefix, "_X_width")), " & ",
      format_mean_mcse(s, paste0(prefix, "_Y_width")), " & ",
      format_mean_mcse(s, paste0(prefix, "_X_post_width")), " \\\\\n",
      " & 95\\% empirical coverage & ",
      format_pct_mcse(s, paste0(prefix, "_X_coverage")), " & ",
      format_pct_mcse(s, paste0(prefix, "_Y_coverage")), " & ",
      format_pct_mcse(s, paste0(prefix, "_X_post_coverage")), " \\\\"
    )
  }

  body <- paste(
    row_block("ridge", "Bagged Ridge"),
    "\\hline",
    row_block("knn", "Bagged $K$-NN"),
    sep = "\n"
  )

  sprintf(
    paste0(
      "\\begin{table}[htbp]\n",
      "\\centering\n",
      "\\small\n",
      "\\caption{Experiment 3 predictive performance for ",
      "$n=%d$, $p=%d$, $B=%d$, $M=%d$ paired Monte Carlo replications, ",
      "and $n_{\\mathrm{test}}=%d$. Entries are Monte Carlo means with ",
      "MCSEs in parentheses. Bootstrap intervals are quantiles of the ",
      "distributions of bootstrap fitted predictions.}\n",
      "\\label{%s}\n",
      "\\begin{tabular}{llccc}\n",
      "\\hline\\hline\n",
      "\\textbf{Algorithm} & \\textbf{Metric} & ",
      "\\textbf{Bootstrap: $X$} & \\textbf{Bootstrap: $Y$} & ",
      "\\textbf{Post-hoc: $X$} \\\\\n",
      "\\hline\n",
      "%s\n",
      "\\hline\\hline\n",
      "\\end{tabular}\n",
      "\\end{table}\n"
    ),
    res$n, res$p, res$B, res$M, res$n_test, label, body
  )
}

make_exp3_conformal_table <- function(
  res,
  label = "tab:split_conformal_reference"
) {
  s <- res$summaries

  sprintf(
    paste0(
      "\\begin{table}[htbp]\n",
      "\\centering\n",
      "\\small\n",
      "\\caption{Split-conformal reference results for Experiment 3. ",
      "Entries are Monte Carlo means (MCSEs) over $M=%d$ replications. ",
      "These values are not directly compared statistically with the ",
      "bootstrap prediction-quantile results because the procedures use ",
      "different training--calibration constructions and target different ",
      "inferential objects.}\n",
      "\\label{%s}\n",
      "\\begin{tabular}{lccc}\n",
      "\\hline\\hline\n",
      "\\textbf{Algorithm} & \\textbf{Test MSE} & ",
      "\\textbf{95\\%% interval width} & ",
      "\\textbf{95\\%% empirical coverage} \\\\\n",
      "\\hline\n",
      "Ridge & %s & %s & %s \\\\\n",
      "$K$-NN & %s & %s & %s \\\\\n",
      "\\hline\\hline\n",
      "\\end{tabular}\n",
      "\\end{table}\n"
    ),
    res$M,
    label,
    format_mean_mcse(s, "ridge_conformal_mse"),
    format_mean_mcse(s, "ridge_conformal_width"),
    format_pct_mcse(s, "ridge_conformal_coverage"),
    format_mean_mcse(s, "knn_conformal_mse"),
    format_mean_mcse(s, "knn_conformal_width"),
    format_pct_mcse(s, "knn_conformal_coverage")
  )
}

result <- run_exp3(
  M = 500L,
  n = 30L,
  p = 10L,
  B = 500L,
  lambda = 1,
  k = 3L,
  noise_sd = 0.5,
  n_test = 100L,
  alpha = 0.05,
  conformal_train_fraction = 0.7
)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)

saveRDS(
  result,
  "results/experiment3_prediction_conformal.rds"
)

write_summary_csv(
  result$summaries,
  "results/experiment3_prediction_conformal_summary.csv"
)

main_exp3_table <- make_exp3_main_table(
  result,
  label = "tab:exp3results"
)

conformal_table <- make_exp3_conformal_table(
  result,
  label = "tab:split_conformal_reference"
)

writeLines(
  main_exp3_table,
  "tables/table11_experiment3.tex"
)

writeLines(
  conformal_table,
  "tables/app2_table4_split_conformal.tex"
)

message("Saved results/experiment3_prediction_conformal.rds")
message("Saved results/experiment3_prediction_conformal_summary.csv")
message("Saved tables/table11_experiment3.tex")
message("Saved tables/app2_table4_split_conformal.tex")