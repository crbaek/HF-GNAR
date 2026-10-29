# Homogeneity Fusion for Sparse Grouped Network VAR Models

**Empirical analysis code** accompanying the manuscript by **Panagiotis Andreou,
Changryong Baek, and Vladas Pipiras** (2026).

Repository: <https://github.com/crbaek/GNAR>

## Purpose

The paper studies grouped network vector autoregression when the sparse network,
coefficient groups, and model parameters are unknown. Its computational procedure
alternates between network selection and joint estimation of the groups and
coefficients using mixed integer programming (MIP).

This directory implements the paper's empirical application to daily realized
volatilities of 31 global equity indices. It fits two coefficient groups under
nine specifications: three initial network methods (SIS, BS, and AL) crossed with
three sparsity settings (RV10, RV15, and RV20). It also compares the estimated
coefficient groups with communities detected from the fitted coefficient-support
graphs. The estimation code is based on the completed `all9_v4` analysis.

The supplied input contains 1,569 observations dated January 4, 2016, through
December 31, 2021. It is the prepared analysis input used by the estimation
routine. Original data collection and preprocessing are described in the paper.
This distribution covers the empirical application.

## Files

| File | Purpose |
|---|---|
| `README.md` | Paper context, installation, execution, and output definitions |
| `run_analysis.R` | Settings, input checks, nine fits, checkpointing, and exports |
| `library.R` | Estimation, solver model construction, and reporting functions |
| `data/volatility.rds` | Prepared empirical input and market labels |
| `packages/sparseVAR_0.4.4.tar.gz` | Source package used for adaptive-Lasso initialization |
| `.gitignore` | Excludes generated results and local R session files |

All paths are resolved relative to `run_analysis.R`. The estimation flow is
written at the beginning of `library.R` in `fit_gnar()`, followed by its helpers.
Sourcing `library.R` only defines functions.

## Requirements and installation

The completed empirical run used R 4.3.3, Matrix 1.6-5, cluster 2.1.6,
Gurobi 11.0.3 (R package 11.0-3), and sparseVAR 0.4.4 on Windows. A Gurobi
installation and a license supporting the model sizes in this application are
required. After installing the dependencies below, install the Gurobi R package
from its installation's `R/` directory
using `install.packages("/path/to/gurobi-package-file", repos = NULL)`.
See the [official Gurobi R installation guide](https://docs.gurobi.com/projects/optimizer/en/current/reference/r/setup.html).

From this directory, install any missing dependencies and the bundled package:

```r
dependencies <- c("Matrix", "cluster", "slam", "glmnet", "Rcpp", "RcppArmadillo")
missing <- dependencies[!vapply(dependencies, requireNamespace,
                                logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
install.packages("packages/sparseVAR_0.4.4.tar.gz",
                 repos = NULL, type = "source")
```

The source package requires a C++17 compiler. On Windows, install the Rtools
version matching your R version before installing it
([Rtools downloads](https://cran.r-project.org/bin/windows/Rtools/)). Use dependency versions
compatible with your R installation. The runner checks that sparseVAR is exactly
version 0.4.4. The bundled package retains its own GPL-3 license.

## Run the analysis

From the repository root, in R or RStudio:

```r
source("run_analysis.R", encoding = "UTF-8")
```

An absolute path to `run_analysis.R` also works from another working directory.

Alternatively, from this directory in a terminal:

```sh
Rscript run_analysis.R
```

Settings are collected at the top of `run_analysis.R`:

| Setting | Default |
|---|---|
| Specification order | RV10: SIS, BS, AL; then RV15; then RV20 |
| Off-diagonal edge budgets | 655, 446, and 259 |
| Number of groups | 2 |
| Gurobi time limit | 200 seconds per optimization |
| Gurobi threads | 60, as in the completed empirical run |
| Maximum refinement attempts | 2 |
| Acceptance rule | At least a 0.1% reduction in uncentered input-data RSS |
| Group/network MIP relative-gap target | 0.01 |
| Master and community seeds | 20261011 |

The budgets are computed by thresholding absolute marginal VAR(1) coefficients
at 0.10, 0.15, and 0.20. SIS and BS use the specified support budget; AL selects
its initial support by five-fold cross-validation. BS uses binary support
indicators. Each fit then follows the common initialization and refinement
routine. The network-update objective retains the full quadratic cross-products.
BS uses Gurobi's default relative-gap target; the 0.01 target applies to the
subsequent group and network MIPs.

In `fit_gnar()`, initialization is followed by grouped least squares and a
fixed-network group MIP. Each refinement attempt then runs a network MIP and
another fixed-network group MIP. The RSS rule determines whether the proposal
is retained, and the first rejected proposal ends refinement.

The specifications run sequentially. `threads` controls Gurobi's threads within
each solve. The runner records settings, package versions, solver status, and
iteration history. With time-limited MIPs, feasible solutions can vary with the
solver version, hardware, and thread count, so numerical equality across
environments is not guaranteed.

## Results and restart

Outputs are written to `results/`:

| File | Contents |
|---|---|
| `results.rds` | Run metadata and the initial, LS, group-MIP, and final retained fits |
| `coefficients.csv` | Final group sizes, GNAR coefficients, RSS, and fit diagnostics |
| `market_groups.csv` | Market-level group and community memberships |
| `block_effects.csv` | Connectivity-weighted within- and between-group effects |
| `AL_RV20_crossclassification.csv` | Group/community cross-classification for AL at RV20 |
| `iteration_flow.csv` | Attempted updates, RSS changes, acceptance, solver status, and gaps |
| `groupstructure.pdf`, `communitydetection.pdf` | Membership figures with specification labels below |
| `progress.log`, `sessionInfo.txt` | Execution log and R environment |

Group 1 is displayed as the group containing `us.SPX`, with coefficients permuted
consistently. Community labels are aligned within each specification to maximize
agreement with its coefficient groups. These display conventions do not change
the fitted models. `coefficients.csv` reports RSS on centered observations;
`iteration_flow.csv` reports the uncentered-input RSS used for acceptance.

The runner saves `results.rds` after each completed specification. Restart the
same script to skip completed fits and continue. A checkpoint must match the
code, data, settings, R version/platform, and recorded package versions. To start another
run, set `output_subdir` to a new directory name. The distribution starts with
no fitted results. Checkpoints from a different release should be kept in
their original directory; they can be used for the export-only workflow below.

To regenerate the tables and figures from your completed checkpoint without
refitting, run the following from this directory:

```r
source("library.R", encoding = "UTF-8")
# Or use an absolute path to a completed all9_v4 checkpoint.
checkpoint_file <- "results/results.rds"
saved <- readRDS(checkpoint_file)
write_results(saved$results, saved$labels, "results",
              community_seed = saved$run_info$settings$community_seed)
```

## Citation

Andreou, P., Baek, C., and Pipiras, V. (2026). *Homogeneity Fusion for Sparse
Grouped Network VAR Models*. Manuscript.
