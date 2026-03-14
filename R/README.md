# R Implementation — Linear Regression BVS

R implementation of the PARNI, ASI, and ADS samplers for Bayesian variable selection in linear regression with Zellner's g-prior or independent normal prior.

## Requirements

- Base R (>= 3.6)
- `tidyverse` (for data handling in examples)
- `rlist` (for loading real datasets)

## Quick Start

```r
setwd("R/")
source("main functions.R")
source("simulation/linear regression random sample generator.R")
source("make_hyper_par.R")

# Generate simulated data
dataset <- lrrsg(n = 500, p = 5000, rho = 0.6, SNR = 2)

# Prepare hyperparameters
hyper_par <- make_hyper_par(
  y = dataset$y, X = dataset$X,
  g = 9, h = 10 / ncol(dataset$X),
  Z = NULL,          # NULL removes intercept via projection
  prior_type = 1     # 1 = independent prior, 2 = g-prior
)

# Run PARNI
alg_par <- list(
  N = 3000, Nb = 1000,
  n_chain = 25, n_temp = 1,
  full_adap = FALSE, verbose = TRUE,
  kappa = 0.001, eps = 0.1 / hyper_par$p,
  omega_adap = "kw", omega_par = c(-1, -0.5),
  omega_init = 0.5,
  store_chains = FALSE,
  bal_fun = function(x) { min(x, 1) },
  use_rb = TRUE
)

results <- PARNI(alg_par = alg_par, hyper_par = hyper_par)
```

## File Structure

### Core

| File | Description |
|------|-------------|
| `main functions.R` | Entry point — sources all sampler files |
| `make_hyper_par.R` | Constructs `hyper_par` list (data preprocessing, prior setup) |

### Samplers (`sampler/`)

| File | Description |
|------|-------------|
| `PARNI.R` | PARNI sampler with parallel tempering, RM/KW omega adaptation |
| `ASI.R` | Adaptively Scaled Individual sampler with zeta adaptation |
| `ADS.R` | Add-Delete-Swap baseline sampler |
| `other_supportive_functions.R` | Bayes factors, log-likelihood, and helper functions |

### Simulation (`simulation/`)

| File | Description |
|------|-------------|
| `linear regression random sample generator.R` | `lrrsg()` — simulated data with correlated covariates |
| `simulation_example.R` | Full example: PARNI, ADS, and ASI on simulated data |
| `simulation_ADS_runs.R` | Batch ADS runs for comparison |
| `load_simulation.R` | Load and process simulation results |
| `dataset_generation.R` | Dataset generation scripts |
| `dataset/` | Saved simulation datasets |

### Real Data (`real dataset/`)

| File | Description |
|------|-------------|
| `real_dataset_example.R` | Example running PARNI/ADS/ASI on real datasets |
| `real_datasets.csv` | Index of available real datasets |
| `dataset/` | Stored real dataset files |

## Key Parameters

### `make_hyper_par(y, X, g, h, Z, prior_type)`

| Parameter | Description |
|-----------|-------------|
| `y` | Response vector |
| `X` | Covariate matrix (n x p) |
| `g` | Prior variance scaling |
| `h` | Prior inclusion probability (scalar for fixed-h, length-2 vector for beta-binomial) |
| `Z` | Fixed effects matrix (`NULL` to project out intercept) |
| `prior_type` | `1` = independent normal prior, `2` = Zellner's g-prior |

### Algorithm Parameters (`alg_par`)

| Parameter | Description |
|-----------|-------------|
| `N` | Total MCMC iterations |
| `Nb` | Burn-in iterations |
| `n_chain` | Number of parallel chains |
| `n_temp` | Number of temperature levels (1 = no tempering) |
| `omega_adap` | PARNI adaptation: `"kw"` (Kiefer-Wolfowitz) or `"rm"` (Robbins-Monro) |
| `omega_par` | Adaptation step-size parameters |
| `bal_fun` | Balancing function: `function(x) min(x, 1)` (Hastings) or `function(x) x/(1+x)` (Barker) |
| `kappa` | PIP floor/ceiling (e.g., 0.001) |
| `use_rb` | Use Rao-Blackwellised PIP estimates |
| `store_chains` | Whether to store full chain history |

## Outputs

The sampler returns a list containing:
- `estm_PIPs` — estimated posterior inclusion probabilities
- `acc_rate` — acceptance rates
- `CPU_time` — total computation time
- `log_post_trace` — log-posterior trace (for convergence diagnostics)
