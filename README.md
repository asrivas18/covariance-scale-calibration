# Covariance-Scale Calibration of Empirical Populations

Code and reproducibility materials for:

> Apurv Srivastav and Sudesh K. Srivastav.  
> *Covariance-Scale Calibration of Empirical Populations: Theory and Applications to Resampling and Machine Learning.*

## Overview

This repository reproduces the numerical analyses reported in the accompanying
main manuscript and its two online appendices.

The central deterministic empirical-support transformation is

$\mathbf{y}_i =
\bar{\mathbf{x}}
+
\sqrt{\frac{n}{n-1}}
(\mathbf{x}_i-\bar{\mathbf{x}}),
\qquad
i=1,\ldots,n.$

It preserves the sample mean and makes the denominator-\(n\) covariance of the
transformed empirical support equal the denominator-\((n-1)\) sample covariance
of the original observations.

The repository includes scripts for:

- Exact analytical post-hoc identity checks
- Experiment 1: high-dimensional covariance and precision estimation
- Experiment 2: bootstrap PCA and eigenvalue intervals
- Sample-size sensitivity analyses
- Dimension-to-sample-size sensitivity analyses
- Experiment 3: bagged Ridge, bagged K-NN, and split-conformal reference analyses
- Experiment 4: bootstrap K-means and centroid variability
- Ridge-shrinkage robustness
- A gene-expression identity illustration using the Bioconductor `ALL` package

## Repository layout

```text
R/          Shared R functions
scripts/    Authoritative analysis scripts
results/    Generated numerical outputs and session information
tables/     Generated LaTeX tables
manuscript/ Final main manuscript and online appendix PDFs
archive/    Historical scripts not used to generate final reported results
```

The authoritative pipeline is:

```text
run_all.R
R/common.R
scripts/
```

Files in `archive/legacy/`, if present, are retained only for historical
reference and are not used to generate the final manuscript tables.

## Requirements

The analyses were run in R. Required packages include:

- `MASS`
- `parallel`
- `BiocManager`
- `Biobase`
- `ALL`

Install the Bioconductor dependencies with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  c("ALL", "Biobase"),
  ask = FALSE,
  update = FALSE
)
```

The frozen repository release includes:

```text
results/session_info.txt
```

and analysis-specific session-information files for selected scripts. These
record the R version, package versions, operating system, and attached package
environment used to generate the archived outputs.

## Reproduce all analyses

From the repository root, run:

```r
source("run_all.R")
```

The pipeline creates or overwrites files in:

```text
results/
tables/
```

The complete run can be computationally intensive. In particular:

- Experiment 2 uses \(M=1000\) Monte Carlo replications and \(B=2000\)
  bootstrap resamples per replication.
- Sample-size sensitivity repeats the bootstrap PCA calculation over
  multiple values of \(n\).
- Dimension-to-sample-size sensitivity runs 1,000 paired replications for
  each configuration.
- Experiment 4 uses \(M=1000\) Monte Carlo replications and \(B=500\)
  bootstrap clusterings per replication.

For a quick diagnostic run, edit the relevant script locally to reduce
\(M\) and \(B\). Do not treat quick-run results as manuscript results.

## Script-to-manuscript map

| Manuscript item | Script | Generated output |
|---|---|---|
| Table 3: Post-hoc identity verification | `scripts/01_posthoc_verification.R` | `tables/table3_posthoc_benchmarks.tex` |
| Table 4: Experiment 1 covariance estimation | `scripts/02_experiment1_covariance.R` | `tables/table4_experiment1.tex` |
| Table 5: Experiment 2 eigenvalues | `scripts/03_experiment2_bootstrap_pca.R` | `tables/table5_experiment2_eigenvalues.tex` |
| Table 6: Experiment 2 interval coverage | `scripts/03_experiment2_bootstrap_pca.R` | `tables/table6_experiment2_coverage.tex` |
| Table 7: Covariance sample-size sensitivity | `scripts/04_sensitivity_analyses.R` | `tables/table7_sample_size_trace.tex` |
| Table 8: Bootstrap PCA sample-size sensitivity | `scripts/04_sensitivity_analyses.R` | `tables/table8_sample_size_pca.tex` |
| Table 9: Fixed-\(n\), varying-\(p\) sensitivity | `scripts/04b_dimension_to_sample_size_sensitivity.R` | `tables/table9_dimension_sensitivity_fixed_n.tex` |
| Table 10: Fixed-\(p\), varying-\(n\) sensitivity | `scripts/04b_dimension_to_sample_size_sensitivity.R` | `tables/table10_dimension_sensitivity_fixed_p.tex` |
| Table 11: Experiment 3 prediction analysis | `scripts/05_experiment3_prediction_conformal.R` | `tables/table11_experiment3.tex` |
| Table 12: Experiment 4 bootstrap K-means | `scripts/06_experiment4_kmeans.R` | `tables/table12_experiment4.tex` |
| Online Appendix 2, Table 1: PCA widths | `scripts/03_experiment2_bootstrap_pca.R` | `tables/app2_table1_pca_widths.tex` |
| Online Appendix 2, Table 2: Ridge robustness | `scripts/07_ridge_alpha_robustness.R` | `tables/app2_table2_ridge_robustness.tex` |
| Online Appendix 2, Table 3: ALL data identity check | `scripts/08_all_gene_expression.R` | `tables/app2_table3_all_identity.tex` |
| Online Appendix 2, Table 4: Split conformal | `scripts/05_experiment3_prediction_conformal.R` | `tables/app2_table4_split_conformal.tex` |

## Reproducibility conventions

- Each named Monte Carlo analysis uses a fixed random-number generator and
  its own fixed replication-level random-number stream.
- Ordinary and surrogate procedures are evaluated on matched simulated data and,
  where applicable, matched bootstrap-resampling index vectors.
- Consequently, overlapping parameter configurations in separately seeded
  analyses need not yield identical Monte Carlo summaries across scripts.
- The surrogate test covariates in the prediction experiment are transformed
  only with training-sample parameters.
- In Experiment 4, corresponding ordered initialization centers are supplied
  to the ordinary and surrogate K-means fits. Cluster labels are compared under
  that initialization-induced correspondence; no post-hoc label matching is
  performed.
- The gene-expression illustration uses observations designated as B-cell acute
  lymphoblastic leukemia in the `ALL` package metadata and the first 50 probes
  in package feature order.

## Generated outputs

The repository stores both machine-readable results and manuscript-ready
LaTeX tables.

Examples include:

```text
results/experiment1_covariance.rds
results/experiment1_covariance_summary.csv
results/experiment2_bootstrap_pca.rds
results/experiment3_prediction_conformal.rds
results/experiment4_kmeans.rds
results/ridge_alpha_robustness_results_M1000.csv
results/all_real_data_covariance_pca_results.csv

tables/table3_posthoc_benchmarks.tex
tables/table4_experiment1.tex
tables/table5_experiment2_eigenvalues.tex
tables/table6_experiment2_coverage.tex
tables/table7_sample_size_trace.tex
tables/table8_sample_size_pca.tex
tables/table9_dimension_sensitivity_fixed_n.tex
tables/table10_dimension_sensitivity_fixed_p.tex
tables/table11_experiment3.tex
tables/table12_experiment4.tex
tables/app2_table1_pca_widths.tex
tables/app2_table2_ridge_robustness.tex
tables/app2_table3_all_identity.tex
tables/app2_table4_split_conformal.tex
```

## Manuscript files

The `manuscript/` directory contains the submitted paper and online appendices:

- `ms_JMLR.pdf`: Main manuscript
- `ms_appendix-1.pdf`: Online Appendix 1, technical proofs
- `ms_appendix-2.pdf`: Online Appendix 2, additional empirical and computational results

## License

The R code and repository documentation are released under the MIT License.
The manuscript and online appendices are governed by their respective
publication and copyright terms. Third-party software and data, including the
Bioconductor `ALL` package, remain subject to their original licenses.

## Citation

If you use this code, please cite the accompanying manuscript. Release-specific
citation metadata are provided in `CITATION.cff`.