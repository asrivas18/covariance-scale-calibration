# ==============================================================================
# Table 3: Numerical Verification of Analytical Post-Hoc Rescaling Relations
#
# This script is separate from the primary empirical experiments.
#
# It uses M_verify = 100 independent Monte Carlo replications to verify:
#   1. Denominator-n empirical covariance scaling;
#   2. PCA eigenvalue scaling;
#   3. Matched bootstrap PCA draw scaling;
#   4. Matched percentile-bootstrap endpoint scaling;
#   5. K-nearest-neighbor ordering invariance;
#   6. Matched-initialization K-means assignment invariance;
#   7. K-means centroid-deviation scaling;
#   8. K-means centroid-variance scaling.
#
# Primary manuscript experiments use:
#   - M = 1000 for Experiments 1, 2, 4, and sensitivity analyses;
#   - M = 500 for Experiment 3.
#
# Outputs:
#   results/table3_posthoc_verification.rds
#   results/table3_posthoc_verification.csv
#   results/table3_posthoc_configuration.csv
#   results/table3_posthoc_session_info.txt
#   tables/table3_posthoc_benchmarks.tex
#
# Run from repository root:
#   source("scripts/01_posthoc_verification.R")
# ==============================================================================

source("R/common.R")

set_reproducible_rng(2026L)
ensure_directories()

M_verify <- 100L
verification_seed <- 2026L
verification_tolerance <- 1e-10

surrogate_scale <- function(n) {
  n <- as.integer(n)

  if (length(n) != 1L || !is.finite(n) || n <= 1L) {
    stop(
      "n must be one finite integer greater than 1.",
      call. = FALSE
    )
  }

  sqrt(n / (n - 1L))
}

transform_with_training_map_local <- function(
  X_new,
  x_bar_train,
  n_train
) {
  X_new <- as.matrix(X_new)
  x_bar_train <- as.numeric(x_bar_train)

  if (ncol(X_new) != length(x_bar_train)) {
    stop(
      "X_new and x_bar_train have incompatible dimensions.",
      call. = FALSE
    )
  }

  s <- surrogate_scale(n_train)

  s * sweep(X_new, 2L, x_bar_train, FUN = "-") +
    matrix(
      x_bar_train,
      nrow = nrow(X_new),
      ncol = ncol(X_new),
      byrow = TRUE
    )
}

nearest_neighbor_order <- function(X_train, x_query) {
  X_train <- as.matrix(X_train)
  x_query <- as.numeric(x_query)

  if (ncol(X_train) != length(x_query)) {
    stop("Training and query dimensions do not match.", call. = FALSE)
  }

  d2 <- rowSums(
    (
      X_train -
        matrix(
          x_query,
          nrow = nrow(X_train),
          ncol = ncol(X_train),
          byrow = TRUE
        )
    )^2
  )

  order(d2)
}

kmeans_fixed <- function(X, centers, max_iter = 100L) {
  X <- as.matrix(X)
  centers <- as.matrix(centers)

  K <- nrow(centers)

  if (ncol(X) != ncol(centers)) {
    stop(
      "X and centers have incompatible dimensions.",
      call. = FALSE
    )
  }

  cluster <- rep.int(NA_integer_, nrow(X))

  for (iteration in seq_len(max_iter)) {
    squared_distances <- vapply(
      seq_len(K),
      function(k) {
        rowSums(
          (
            X -
              matrix(
                centers[k, ],
                nrow = nrow(X),
                ncol = ncol(X),
                byrow = TRUE
              )
          )^2
        )
      },
      numeric(nrow(X))
    )

    new_cluster <- max.col(
      -squared_distances,
      ties.method = "first"
    )

    new_centers <- centers

    for (k in seq_len(K)) {
      cluster_members <- X[new_cluster == k, , drop = FALSE]

      if (nrow(cluster_members) > 0L) {
        new_centers[k, ] <- colMeans(cluster_members)
      }
    }

    if (
      identical(new_cluster, cluster) &&
      isTRUE(all.equal(new_centers, centers, tolerance = 1e-12))
    ) {
      cluster <- new_cluster
      centers <- new_centers
      break
    }

    cluster <- new_cluster
    centers <- new_centers
  }

  list(
    cluster = cluster,
    centers = centers
  )
}

total_centroid_variance <- function(centers) {
  if (length(dim(centers)) != 3L) {
    stop(
      "centers must be a three-dimensional array: bootstrap x cluster x dimension.",
      call. = FALSE
    )
  }

  total <- 0

  for (k in seq_len(dim(centers)[2])) {
    center_matrix <- centers[, k, , drop = FALSE]
    center_matrix <- matrix(
      center_matrix,
      nrow = dim(centers)[1],
      ncol = dim(centers)[3]
    )

    total <- total + sum(diag(stats::cov(center_matrix)))
  }

  total
}

one_table3_replication <- function(
  seed,
  n_cov = 30L,
  p_cov = 10L,
  n_pca = 25L,
  p_pca = 5L,
  B_pca = 200L,
  n_knn = 30L,
  p_knn = 10L,
  n_kmeans = 75L,
  p_kmeans = 2L,
  K = 3L,
  B_kmeans = 100L,
  alpha = 0.05
) {
  set.seed(as.integer(seed))

  # ---------------------------------------------------------------------------
  # 1. Exact denominator-n covariance scaling.
  # ---------------------------------------------------------------------------
  X_cov <- matrix(rnorm(n_cov * p_cov), nrow = n_cov, ncol = p_cov)
  Y_cov <- surrogate_population(X_cov)

  covariance_error <- max(abs(
    empirical_covariance_n(Y_cov) -
      (n_cov / (n_cov - 1L)) * empirical_covariance_n(X_cov)
  ))

  # ---------------------------------------------------------------------------
  # 2. Exact PCA eigenvalue scaling and matched bootstrap scaling.
  # ---------------------------------------------------------------------------
  X_pca <- matrix(rnorm(n_pca * p_pca), nrow = n_pca, ncol = p_pca)
  Y_pca <- surrogate_population(X_pca)
  pca_factor <- n_pca / (n_pca - 1L)

  eigenvalues_X <- eigen(
    empirical_covariance_n(X_pca),
    symmetric = TRUE,
    only.values = TRUE
  )$values

  eigenvalues_Y <- eigen(
    empirical_covariance_n(Y_pca),
    symmetric = TRUE,
    only.values = TRUE
  )$values

  pca_eigenvalue_error <- max(abs(
    eigenvalues_Y - pca_factor * eigenvalues_X
  ))

  bootstrap_draw_error <- 0
  bootstrap_indices <- matrix(
    NA_integer_,
    nrow = n_pca,
    ncol = B_pca
  )

  for (b in seq_len(B_pca)) {
    index <- sample.int(n_pca, size = n_pca, replace = TRUE)
    bootstrap_indices[, b] <- index

    eigenvalues_X_b <- eigen(
      empirical_covariance_n(X_pca[index, , drop = FALSE]),
      symmetric = TRUE,
      only.values = TRUE
    )$values

    eigenvalues_Y_b <- eigen(
      empirical_covariance_n(Y_pca[index, , drop = FALSE]),
      symmetric = TRUE,
      only.values = TRUE
    )$values

    bootstrap_draw_error <- max(
      bootstrap_draw_error,
      max(abs(eigenvalues_Y_b - pca_factor * eigenvalues_X_b))
    )
  }

  bootstrap_lambda_X <- vapply(
    seq_len(B_pca),
    function(b) {
      index <- bootstrap_indices[, b]

      eigen(
        empirical_covariance_n(X_pca[index, , drop = FALSE]),
        symmetric = TRUE,
        only.values = TRUE
      )$values[1]
    },
    numeric(1)
  )

  bootstrap_lambda_Y <- vapply(
    seq_len(B_pca),
    function(b) {
      index <- bootstrap_indices[, b]

      eigen(
        empirical_covariance_n(Y_pca[index, , drop = FALSE]),
        symmetric = TRUE,
        only.values = TRUE
      )$values[1]
    },
    numeric(1)
  )

  pca_percentile_endpoint_error <- max(abs(
    percentile_interval(bootstrap_lambda_Y, alpha = alpha) -
      pca_factor * percentile_interval(bootstrap_lambda_X, alpha = alpha)
  ))

  # ---------------------------------------------------------------------------
  # 3. K-NN ordering under a common training-derived transformation.
  # ---------------------------------------------------------------------------
  X_knn <- matrix(rnorm(n_knn * p_knn), nrow = n_knn, ncol = p_knn)
  X_query <- matrix(rnorm(20L * p_knn), nrow = 20L, ncol = p_knn)

  Y_knn <- surrogate_population(X_knn)

  Y_query <- transform_with_training_map_local(
    X_new = X_query,
    x_bar_train = colMeans(X_knn),
    n_train = n_knn
  )

  knn_neighbor_order_error <- 0

  for (q in seq_len(nrow(X_query))) {
    order_X <- nearest_neighbor_order(X_knn, X_query[q, ])
    order_Y <- nearest_neighbor_order(Y_knn, Y_query[q, ])

    knn_neighbor_order_error <- max(
      knn_neighbor_order_error,
      as.numeric(!identical(order_X, order_Y))
    )
  }

  # ---------------------------------------------------------------------------
  # 4. Matched-initialization K-means identities.
  # ---------------------------------------------------------------------------
  centers_true <- rbind(
    c(-2, 0),
    c(2, 0),
    c(0, 2.5)
  )

  cluster_sizes <- rep.int(n_kmeans %/% K, K)
  cluster_sizes[seq_len(n_kmeans %% K)] <-
    cluster_sizes[seq_len(n_kmeans %% K)] + 1L

  X_kmeans <- do.call(
    rbind,
    lapply(seq_len(K), function(k) {
      matrix(
        rnorm(cluster_sizes[k] * p_kmeans, sd = 0.5),
        ncol = p_kmeans
      ) +
        matrix(
          centers_true[k, ],
          nrow = cluster_sizes[k],
          ncol = p_kmeans,
          byrow = TRUE
        )
    })
  )

  Y_kmeans <- surrogate_population(X_kmeans)
  s_kmeans <- surrogate_scale(n_kmeans)

  initial_X <- centers_true +
    matrix(
      rnorm(K * p_kmeans, sd = 0.1),
      nrow = K,
      ncol = p_kmeans
    )

  initial_Y <- transform_with_training_map_local(
    X_new = initial_X,
    x_bar_train = colMeans(X_kmeans),
    n_train = n_kmeans
  )

  assignment_error <- 0
  centroid_deviation_error <- 0

  centers_X <- array(
    NA_real_,
    dim = c(B_kmeans, K, p_kmeans)
  )

  centers_Y <- array(
    NA_real_,
    dim = c(B_kmeans, K, p_kmeans)
  )

  for (b in seq_len(B_kmeans)) {
    index <- sample.int(n_kmeans, size = n_kmeans, replace = TRUE)

    fit_X <- kmeans_fixed(
      X_kmeans[index, , drop = FALSE],
      initial_X
    )

    fit_Y <- kmeans_fixed(
      Y_kmeans[index, , drop = FALSE],
      initial_Y
    )

    assignment_error <- max(
      assignment_error,
      as.numeric(!identical(fit_X$cluster, fit_Y$cluster))
    )

    x_bar_b <- colMeans(X_kmeans[index, , drop = FALSE])
    y_bar_b <- colMeans(Y_kmeans[index, , drop = FALSE])

    deviation_error <- max(abs(
      (
        fit_Y$centers -
          matrix(
            y_bar_b,
            nrow = K,
            ncol = p_kmeans,
            byrow = TRUE
          )
      ) -
        s_kmeans *
          (
            fit_X$centers -
              matrix(
                x_bar_b,
                nrow = K,
                ncol = p_kmeans,
                byrow = TRUE
              )
          )
    ))

    centroid_deviation_error <- max(
      centroid_deviation_error,
      deviation_error
    )

    centers_X[b, , ] <- fit_X$centers
    centers_Y[b, , ] <- fit_Y$centers
  }

  variance_X <- total_centroid_variance(centers_X)
  variance_Y <- total_centroid_variance(centers_Y)

  centroid_variance_error <- abs(
    variance_Y - s_kmeans^2 * variance_X
  )

  c(
    covariance = covariance_error,
    pca_eigenvalue = pca_eigenvalue_error,
    pca_bootstrap_draw = bootstrap_draw_error,
    pca_percentile_endpoint = pca_percentile_endpoint_error,
    knn_neighbor_order = knn_neighbor_order_error,
    kmeans_assignment = assignment_error,
    kmeans_centroid_deviation = centroid_deviation_error,
    kmeans_centroid_variance = centroid_variance_error
  )
}

run_table3_verification <- function(
  M = M_verify,
  master_seed = verification_seed,
  tolerance = verification_tolerance
) {
  replication_seeds <- as.integer(
    master_seed + seq_len(M) - 1L
  )

  raw_discrepancies <- t(vapply(
    replication_seeds,
    one_table3_replication,
    numeric(8)
  ))

  maximum_discrepancies <- apply(
    raw_discrepancies,
    2L,
    max
  )

  list(
    M = M,
    master_seed = master_seed,
    tolerance = tolerance,
    max_abs_discrepancies = maximum_discrepancies,
    passed = maximum_discrepancies <= tolerance,
    raw_discrepancies = raw_discrepancies
  )
}

format_scientific_latex <- function(x) {
  x <- abs(as.numeric(x))

  if (x == 0) {
    return("$0$")
  }

  exponent <- floor(log10(x))
  mantissa <- x / (10^exponent)

  sprintf(
    "$%.2f\\times 10^{%d}$",
    mantissa,
    exponent
  )
}

make_table3_latex <- function(
  checks,
  label = "tab:posthoc_benchmarks"
) {
  discrepancies <- checks$max_abs_discrepancies

  sprintf(
    paste0(
      "\\begin{table}[htbp]\n",
      "\\centering\n",
      "\\small\n",
      "\\setlength{\\tabcolsep}{4pt}\n",
      "\\caption{Numerical verification of exact analytical post-hoc ",
      "rescaling relations over $M=%d$ independent Monte Carlo replications. ",
      "Entries report maximum absolute discrepancies. The checks use matched ",
      "stochastic inputs and the stated matching conditions.}\n",
      "\\label{%s}\n",
      "\\begin{tabular}{p{3.7cm}p{4.6cm}c}\n",
      "\\hline\\hline\n",
      "\\textbf{Quantity} & \\textbf{Analytical relation} & ",
      "\\textbf{Maximum discrepancy} \\\\\n",
      "\\hline\n",
      "Denominator-$n$ covariance & ",
      "$\\mathbf{S}_{Y,\\mathrm{emp}}=\\frac{n}{n-1}",
      "\\mathbf{S}_{X,\\mathrm{emp}}$ & %s \\\\\n",
      "Matched PCA eigenvalues & ",
      "$\\lambda_j(\\mathbf{S}_Y)=\\frac{n}{n-1}",
      "\\lambda_j(\\mathbf{S}_X)$ & %s \\\\\n",
      "Matched PCA bootstrap draws & ",
      "$\\lambda_{j,Y}^{*(b)}=\\frac{n}{n-1}",
      "\\lambda_{j,X}^{*(b)}$ & %s \\\\\n",
      "Matched PCA percentile endpoints & ",
      "$q_{\\alpha,Y}=\\frac{n}{n-1}q_{\\alpha,X}$ & %s \\\\\n",
      "Matched $K$-NN neighbor ordering & ",
      "identical distance ordering under matched scaling & %s \\\\\n",
      "Matched $K$-means assignments & ",
      "identical assignments under stated matching & %s \\\\\n",
      "Matched $K$-means centroid deviations & ",
      "$\\widehat{\\boldsymbol{\\mu}}_Y-\\bar{\\mathbf{Y}}=",
      "\\sqrt{\\frac{n}{n-1}}",
      "(\\widehat{\\boldsymbol{\\mu}}_X-\\bar{\\mathbf{X}})$ & %s \\\\\n",
      "Matched $K$-means centroid variance & ",
      "$V_Y=\\frac{n}{n-1}V_X$ & %s \\\\\n",
      "\\hline\\hline\n",
      "\\end{tabular}\n",
      "\\end{table}\n"
    ),
    checks$M,
    label,
    format_scientific_latex(discrepancies["covariance"]),
    format_scientific_latex(discrepancies["pca_eigenvalue"]),
    format_scientific_latex(discrepancies["pca_bootstrap_draw"]),
    format_scientific_latex(discrepancies["pca_percentile_endpoint"]),
    format_scientific_latex(discrepancies["knn_neighbor_order"]),
    format_scientific_latex(discrepancies["kmeans_assignment"]),
    format_scientific_latex(discrepancies["kmeans_centroid_deviation"]),
    format_scientific_latex(discrepancies["kmeans_centroid_variance"])
  )
}

checks_table3 <- run_table3_verification()

if (!all(checks_table3$passed)) {
  stop(
    "At least one Table 3 verification relation exceeded the stated tolerance.",
    call. = FALSE
  )
}

table3_latex <- make_table3_latex(
  checks_table3,
  label = "tab:posthoc_benchmarks"
)

table3_summary <- data.frame(
  quantity = names(checks_table3$max_abs_discrepancies),
  max_abs_discrepancy = unname(checks_table3$max_abs_discrepancies),
  passed = unname(checks_table3$passed),
  row.names = NULL
)

table3_configuration <- data.frame(
  parameter = c(
    "master_seed",
    "M_verify",
    "verification_tolerance",
    "n_cov",
    "p_cov",
    "n_pca",
    "p_pca",
    "B_pca",
    "n_knn",
    "p_knn",
    "n_kmeans",
    "p_kmeans",
    "K",
    "B_kmeans",
    "nominal_alpha"
  ),
  value = c(
    as.character(verification_seed),
    as.character(M_verify),
    format(verification_tolerance, scientific = TRUE),
    "30",
    "10",
    "25",
    "5",
    "200",
    "30",
    "10",
    "75",
    "2",
    "3",
    "100",
    "0.05"
  ),
  row.names = NULL
)

saveRDS(
  checks_table3,
  file = "results/table3_posthoc_verification.rds"
)

write.csv(
  table3_summary,
  file = "results/table3_posthoc_verification.csv",
  row.names = FALSE
)

write.csv(
  table3_configuration,
  file = "results/table3_posthoc_configuration.csv",
  row.names = FALSE
)

writeLines(
  table3_latex,
  con = "tables/table3_posthoc_benchmarks.tex"
)

writeLines(
  capture.output(sessionInfo()),
  con = "results/table3_posthoc_session_info.txt"
)

print(checks_table3$max_abs_discrepancies)
print(checks_table3$passed)

cat("\n--- Table 3: Post-hoc verification ---\n\n")
cat(table3_latex, "\n")

message("Table 3 post-hoc verification completed.")
message("Saved results/table3_posthoc_verification.rds")
message("Saved results/table3_posthoc_verification.csv")
message("Saved results/table3_posthoc_configuration.csv")
message("Saved results/table3_posthoc_session_info.txt")
message("Saved tables/table3_posthoc_benchmarks.tex")