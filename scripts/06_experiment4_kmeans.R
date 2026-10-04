## =============================================================================
## Experiment 4: Bootstrap K-Means Clustering and Centroid Variability
## Matched initialization for exact affine-invariance verification
## =============================================================================

library(parallel)

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)

## 1. Surrogate population transformation
get_surrogate_population <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n <= 1L) stop("Sample size n must be greater than 1.")
  x_bar <- colMeans(X)
  sqrt(n / (n - 1)) * sweep(X, 2L, x_bar, "-") +
    matrix(x_bar, nrow = n, ncol = ncol(X), byrow = TRUE)
}

## 2. Three isotropic Gaussian clusters in R^2
## The present verification experiment is intentionally fixed at K = 3 and p = 2.
gen_exp4_data <- function(n_per_cluster = 25L, K = 3L, p = 2L, delta = 2) {
  if (K != 3L || p != 2L) stop("This experiment is specified for K = 3 and p = 2.")
  centers <- rbind(
    c(-delta, 0),
    c( delta, 0),
    c(0, delta * sqrt(3))
  )
  X <- do.call(rbind, lapply(seq_len(K), function(g) {
    matrix(rnorm(n_per_cluster * p), nrow = n_per_cluster, ncol = p) +
      matrix(centers[g, ], nrow = n_per_cluster, ncol = p, byrow = TRUE)
  }))
  list(X = X, true_labels = rep(seq_len(K), each = n_per_cluster), centers = centers)
}

## 3. Adjusted Rand index
calc_ari <- function(labels1, labels2) {
  tab <- table(labels1, labels2)
  n <- sum(tab)
  if (n < 2L) return(1)
  sum_ij <- sum(choose(tab, 2))
  sum_i <- sum(choose(rowSums(tab), 2))
  sum_j <- sum(choose(colSums(tab), 2))
  expected <- sum_i * sum_j / choose(n, 2)
  maximum <- 0.5 * (sum_i + sum_j)
  if (abs(maximum - expected) <= .Machine$double.eps) return(1)
  (sum_ij - expected) / (maximum - expected)
}

## 4. K-means with matched affine initializations
## Supplying corresponding initial centers ensures that X and Y follow the same
## Lloyd iterations up to the surrogate affine transformation.
matched_kmeans <- function(X_b, Y_b, K, initial_rows, iter.max = 100L) {
  starts_X <- X_b[initial_rows, , drop = FALSE]
  starts_Y <- Y_b[initial_rows, , drop = FALSE]

  list(
    X = stats::kmeans(
      X_b,
      centers = starts_X,
      iter.max = iter.max,
      nstart = 1L,
      algorithm = "Lloyd"
    ),
    Y = stats::kmeans(
      Y_b,
      centers = starts_Y,
      iter.max = iter.max,
      nstart = 1L,
      algorithm = "Lloyd"
    )
  )
}
## 5. One paired Monte Carlo replication
exp4_single <- function(n_per_cluster = 25L, K = 3L, p = 2L,
                        B = 500L, delta = 2) {
  dat <- gen_exp4_data(n_per_cluster, K, p, delta)
  X <- dat$X
  true_labels <- dat$true_labels
  n <- nrow(X)
  Y <- get_surrogate_population(X)

  ari_X <- numeric(B)
  ari_Y <- numeric(B)
  centers_X <- matrix(NA_real_, nrow = B, ncol = K * p)
  centers_Y <- matrix(NA_real_, nrow = B, ncol = K * p)
  assignment_agreement <- numeric(B)

  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    X_b <- X[idx, , drop = FALSE]
    Y_b <- Y[idx, , drop = FALSE]

    ## Bootstrap duplicates may occur, but continuous Gaussian observations make
    ## at least K distinct sampled original observations overwhelmingly likely.
    distinct_positions <- match(unique(idx), idx)
    if (length(distinct_positions) < K) next
    initial_rows <- sample(distinct_positions, K, replace = FALSE)

    fit <- matched_kmeans(X_b, Y_b, K, initial_rows)
    km_X <- fit$X
    km_Y <- fit$Y

    ari_X[b] <- calc_ari(true_labels[idx], km_X$cluster)
    ari_Y[b] <- calc_ari(true_labels[idx], km_Y$cluster)

if (!identical(km_X$cluster, km_Y$cluster)) {
  stop(
    "Matched K-means assignments disagree. Check the ordered initial centers, ",
    "tie-breaking, numerical tolerances, stopping rules, and implementation.",
    call. = FALSE
  )
}

assignment_agreement[b] <- 1

## The same ordered bootstrap-row positions define the initial centers in the
## X and Y runs. Under matched Lloyd iterations and deterministic tie-breaking,
## this initialization order defines the common cluster-label correspondence;
## no post-hoc label matching is applied.
centers_X[b, ] <- as.vector(km_X$centers)
centers_Y[b, ] <- as.vector(km_Y$centers)

  }

  keep <- is.finite(ari_X) & is.finite(ari_Y) &
    apply(is.finite(centers_X), 1L, all) & apply(is.finite(centers_Y), 1L, all)
  if (sum(keep) < 2L) stop("Too few valid bootstrap clustering replicates.")

  centers_X <- centers_X[keep, , drop = FALSE]
  centers_Y <- centers_Y[keep, , drop = FALSE]
  ari_X <- ari_X[keep]
  ari_Y <- ari_Y[keep]
  assignment_agreement <- assignment_agreement[keep]

  total_centroid_var_X <- sum(apply(centers_X, 2L, var))
  total_centroid_var_Y <- sum(apply(centers_Y, 2L, var))

  c(
    valid_bootstraps = sum(keep),
    ari_X_mean = mean(ari_X),
    ari_Y_mean = mean(ari_Y),
    ari_X_sd = sd(ari_X),
    ari_Y_sd = sd(ari_Y),
    assignment_agreement = mean(assignment_agreement),
    centroid_var_X = total_centroid_var_X,
    centroid_var_Y = total_centroid_var_Y,
    variance_ratio = total_centroid_var_Y / total_centroid_var_X
  )
}

mc_summary <- function(x) {
  n_eff <- sum(is.finite(x))
  c(mean = mean(x, na.rm = TRUE), mcse = sd(x, na.rm = TRUE) / sqrt(n_eff))
}

## 6. Cross-platform Monte Carlo runner
run_exp4_grid <- function(n_per_cluster = 25L, K = 3L, p = 2L,
                          B = 500L, M = 1000L, ncores = 1L, delta = 2) {
  n_total <- n_per_cluster * K
  cat(sprintf("Running Experiment 4: n = %d, K = %d, p = %d, M = %d, B = %d\n",
              n_total, K, p, M, B))
  one_rep <- function(i) exp4_single(n_per_cluster, K, p, B, delta)

  if (.Platform$OS.type == "windows" && ncores > 1L) {
    cl <- makeCluster(ncores)
    on.exit(stopCluster(cl), add = TRUE)
    clusterSetRNGStream(cl, iseed = 2026)
    clusterExport(cl, c("get_surrogate_population", "gen_exp4_data", "calc_ari",
                        "matched_kmeans", "exp4_single"), envir = environment())
    res_list <- parLapply(cl, seq_len(M), one_rep)
  } else if (.Platform$OS.type != "windows" && ncores > 1L) {
    res_list <- mclapply(seq_len(M), one_rep, mc.cores = ncores, mc.set.seed = TRUE)
  } else {
    res_list <- lapply(seq_len(M), one_rep)
  }

  raw <- do.call(rbind, res_list)
  list(
    n = n_total, K = K, p = p, B = B, M = M, delta = delta,
    summaries = t(apply(raw, 2L, mc_summary)),
    raw = raw
  )
}

format_mean_mcse <- function(s, name, digits_mean = 3L, digits_mcse = 3L) {
  sprintf(paste0("%.", digits_mean, "f (%.", digits_mcse, "f)"),
          s[name, "mean"], s[name, "mcse"])
}

## 7. LaTeX table generator
## The target variance ratio is n/(n-1), because centroid deviations under the
## surrogate are multiplied by sqrt(n/(n-1)).
generate_exp4_latex_table <- function(res_config, label = "tab:exp4results") {
  s <- res_config$summaries
  target_ratio <- res_config$n / (res_config$n - 1)

  latex <- sprintf(
"\\begin{table}[htbp]
\\centering
\\small
\\caption{Experiment 4: matched-initialization bootstrap $K$-means verification ($n=%d$, $K=%d$, $p=%d$, $B=%d$, $M=%d$). Entries are Monte Carlo means (MCSEs).}
\\label{%s}
\\begin{tabular}{lccc}
\\hline\\hline
\\textbf{Metric} & \\textbf{Theoretical target} & \\textbf{Standard ($\\mathbf{X}$)} & \\textbf{Surrogate ($\\mathbf{Y}$)} \\\\
\\hline
Mean bootstrap ARI & -- & %s & %s \\\\
Bootstrap ARI SD & -- & %s & %s \\\\
Cluster-assignment agreement & 1.000 & \\multicolumn{2}{c}{%s} \\\\
Total centroid variance & -- & %s & %s \\\\
Centroid-variance ratio $Y/X$ & %.4f & \\multicolumn{2}{c}{%s} \\\\
\\hline\\hline
\\end{tabular}
\\end{table}",
    res_config$n, res_config$K, res_config$p, res_config$B, res_config$M, label,
    format_mean_mcse(s, "ari_X_mean"),
    format_mean_mcse(s, "ari_Y_mean"),
    format_mean_mcse(s, "ari_X_sd"),
    format_mean_mcse(s, "ari_Y_sd"),
    format_mean_mcse(s, "assignment_agreement", 4L, 4L),
    format_mean_mcse(s, "centroid_var_X"),
    format_mean_mcse(s, "centroid_var_Y"),
    target_ratio,
    format_mean_mcse(s, "variance_ratio", 4L, 4L)
  )
  cat(latex, "\n")
  invisible(latex)
}

## =============================================================================
## Execution
## =============================================================================
## Use M = 200 for a quick diagnostic; use M = 1000 for final reporting.
res_exp4 <- run_exp4_grid(
  n_per_cluster = 25L,
  K = 3L,
  p = 2L,
  B = 500L,
  M = 1000L,
  ncores = 1L,
  delta = 2
)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)

saveRDS(
  res_exp4,
  "results/experiment4_kmeans.rds"
)

write.csv(
  data.frame(
    statistic = rownames(res_exp4$summaries),
    mean = res_exp4$summaries[, "mean"],
    mcse = res_exp4$summaries[, "mcse"],
    row.names = NULL
  ),
  "results/experiment4_kmeans_summary.csv",
  row.names = FALSE
)

latex_code_exp4 <- generate_exp4_latex_table(
  res_exp4,
  label = "tab:exp4results"
)

writeLines(
  latex_code_exp4,
  "tables/table12_experiment4.tex"
)

message("Saved results/experiment4_kmeans.rds")
message("Saved results/experiment4_kmeans_summary.csv")
message("Saved tables/table12_experiment4.tex")