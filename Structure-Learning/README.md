# Structure Learning — PARNI for DAG Model Selection

R implementation of the PARNI sampler adapted for Bayesian structure learning of Directed Acyclic Graphs (DAGs). Includes comparison with order MCMC and partition MCMC via the BiDAG package.

## Overview

This module extends the PARNI scheme from variable selection to DAG structure learning. Instead of selecting covariates, the sampler explores the space of adjacency matrices (DAGs) representing causal relationships between variables, using a node-wise linear regression likelihood with g-prior.

## Requirements

- R (>= 3.6)
- `pcalg` — for skeleton estimation via PC algorithm
- `BiDAG` — for partition MCMC and order MCMC comparisons
- `rlist`, `tidyverse` — for data handling

```r
install.packages(c("pcalg", "BiDAG", "rlist", "tidyverse"))
```

## Quick Start

```r
setwd("Structure-Learning/")
source("0_examples.R")  # sources all algorithm files and runs example
```

Or step by step:

```r
# Source algorithms
setwd("algorithms/")
source("PARNI.R")
source("ADR.R")
source("Compute_LA_DAG.R")
source("sample_ind_DAG.R")
source("logit_e.R")
source("is.DAG_adjmat.R")
source("log_llh_DAG.R")
source("log_llh_DAG_update_table.R")
source("other_functions.R")
source("DAG_heatmap.R")
source("marPIPs_DAG_H.R")
source("log_llh_bge_table.R")

# Prepare data
X <- scale(your_data)
X <- t(X)  # variables in rows, observations in columns

# Skeleton estimation (restricts search space)
library(pcalg)
suffStat <- list(C = cor(t(X)), n = ncol(X))
skel <- skeleton(suffStat, indepTest = gaussCItest,
                 labels = as.character(1:nrow(X)), alpha = 0.05)
skel.W <- as(skel@graph, "matrix")

# Hyperparameters
hyper_par <- list(
  n = ncol(X), p = nrow(X),
  g = 10, h = 1 / nrow(X),
  X = X, XtX = X %*% t(X)
)

# Run PARNI for DAGs
alg_par <- list(
  N = 50000, Nb = 10000,
  n_chain = 1, verbose = TRUE,
  kappa = 0.01,
  omega_adap = "rm", omega_par = c(-0.7, 20),
  eps = 1 / (hyper_par$p * (hyper_par$p - 1)),
  use_logit_e = TRUE,
  store_chains = FALSE,
  bal_fun = function(x) { pmin(1, x) },
  H = skel.W  # skeleton constraint
)

results <- PARNI(alg_par, hyper_par)

# Visualise results
DAG_heatmap(results$estm_PIPs, text_bound = 0.7)
```

## File Structure

### Core Algorithms (`algorithms/`)

| File | Description |
|------|-------------|
| `PARNI.R` | PARNI sampler for DAG space |
| `ADR.R` | Add-Delete-Reverse baseline sampler |
| `Compute_LA_DAG.R` | Approximate likelihood computation for DAGs |
| `sample_ind_DAG.R` | Independent edge sampling utilities |
| `log_llh_DAG.R` | Node-wise log-likelihood under g-prior |
| `log_llh_DAG_update_table.R` | Incremental likelihood update table |
| `log_llh_bge_table.R` | BGe score table computation |
| `update_LA_DAG_sawp.R` | Likelihood update after edge swap |
| `is.DAG_adjmat.R` | DAG acyclicity check on adjacency matrix |
| `marPIPs_DAG_H.R` | Marginal PIP estimation for DAG edges |
| `marPIPs_bge_H.R` | Marginal PIPs under BGe score |
| `logit_e.R` | Epsilon-bounded logit transform |
| `other_functions.R` | Miscellaneous helper functions |
| `DAG_heatmap.R` | Heatmap visualisation for adjacency/PIP matrices |

### Experiments

| Folder / File | Description |
|---------------|-------------|
| `0_examples.R` | Main example script: PARNI + partition MCMC comparison |
| `Section 4.2/` | Experiments from paper Section 4.2 (Ecoli70, arth150, magic-niab, protein) |
| `gsim100/` | Gaussian simulation with p=100 nodes — PARNI, ADR, order, partition MCMC |

### Section 4.2 Experiments

| File | Description |
|------|-------------|
| `Ecoli70.R` | E. coli gene regulatory network (p=70) |
| `arth150.R` | Arabidopsis gene expression network (p=150) |
| `magic-niab.R` | MAGIC wheat dataset |
| `Protein_Data_BCDnets.R` | Protein signalling network |
| `*_pen.R` | Penalised variants of the above experiments |

### gsim100 Experiments

| File | Description |
|------|-------------|
| `PARNI.R` / `ADR.R` | Single-run PARNI and ADR |
| `iterative_PARNI.R` / `iterative_ADR.R` | Repeated runs for performance comparison |
| `PC_PARNI.R` / `PC_ADR.R` | Runs with PC-algorithm skeleton constraint |
| `order.R` / `partition.R` | Order MCMC and partition MCMC via BiDAG |
| `Best_skeleton.R` | Oracle skeleton comparison |
| `traces.R` | Trace plot generation |

## Key Differences from Variable Selection

| Aspect | Variable Selection (`R/`) | Structure Learning |
|--------|---------------------------|--------------------|
| Search space | Binary vector gamma in {0,1}^p | Adjacency matrix W in {0,1}^(p x p), constrained to DAGs |
| Likelihood | Linear regression marginal likelihood | Node-wise linear regression with g-prior |
| Constraint | None | Acyclicity (DAG) |
| Skeleton | Not used | PC-algorithm skeleton restricts search space |
| Baseline | ADS (Add-Delete-Swap) | ADR (Add-Delete-Reverse) |
| Output | PIP vector (length p) | PIP matrix (p x p) |

## Outputs

The sampler returns a list containing:
- `estm_PIPs` — edge inclusion probability matrix (p x p)
- `ad_PIPs` — adaptive PIP estimates used during sampling
- `omega` — adapted neighbourhood size parameter
- `acc_rate` — acceptance rate
- `ESJD` — expected squared jumping distance
- `CPU_time` — total computation time
- `log_post_trace` — log-posterior trace
- `k_sizes` — neighbourhood sizes over iterations
