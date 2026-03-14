# GLM_PARNI — PARNI Sampler for Generalised Linear Models

R implementation of the PARNI and ADS samplers for Bayesian variable selection in GLMs, supporting logistic regression, Cox proportional hazards, and Weibull survival models.

## Supported Models

| Model | Likelihood Approximation | Key Files |
|-------|--------------------------|-----------|
| **Logistic** | ALA, Laplace, CPM, Data Augmentation (Polya-Gamma) | `logistic_*.R` |
| **Cox PH** | ALA, Laplace, CPM | `cox_*.R` |
| **Weibull** | ALA, Laplace, CPM | `weibull_*.R` |

**Likelihood approximation methods:**
- **ALA** — Approximate Laplace Approximation (fast, used within MCMC)
- **Laplace** — Full Laplace approximation
- **CPM** — Conditional Predictive Marginal
- **DA** — Data Augmentation via Polya-Gamma (logistic only)

## File Structure

### Core Samplers

| File | Description |
|------|-------------|
| `main.R` | Entry point — sources all required files |
| `GLM_PARNI.R` | PARNI sampler with parallel tempering |
| `GLM_ADS.R` | Add-Delete-Swap baseline sampler |
| `PARNI_proposal.R` | PARNI proposal mechanism (sequential locally-informed flips) |

### Logistic Regression

| File | Description |
|------|-------------|
| `logistic_ALA.R` | ALA about zero |
| `logistic_ALA2.R` | ALA about current mode |
| `logistic_ALA_rb.R` | ALA with Rao-Blackwellised Bayes factors |
| `logistic_DA_rb.R` | Data augmentation (Polya-Gamma) with Rao-Blackwellised BFs |
| `logistic_IRLS.R` | Iteratively Reweighted Least Squares for MLE |
| `logistic_laplace.R` | Full Laplace approximation |
| `logistic_cpm.R` | Conditional Predictive Marginal |
| `logistic_log_post_beta.R` | Log posterior of regression coefficients |
| `logistic_update_matrices.R` | Design matrix updates for active model |
| `logistic_update_omega.R` | Polya-Gamma auxiliary variable update |

### Cox Proportional Hazards

| File | Description |
|------|-------------|
| `cox_ALA.R` | ALA about zero |
| `cox_ALA2.R` | ALA about current mode |
| `cox_ALA_rb.R` | ALA with Rao-Blackwellised Bayes factors |
| `cox_NR_coxph.R` | Newton-Raphson via `coxph` for MLE |
| `cox_compute_ALA_W.R` | Compute ALA weight matrix |
| `cox_functions.R` | Cox model helper functions |
| `cox_laplace.R` | Full Laplace approximation |
| `cox_cpm.R` | Conditional Predictive Marginal |
| `cox_log_post_beta.R` | Log posterior of regression coefficients |

### Weibull Survival

| File | Description |
|------|-------------|
| `weibull_ALA.R` | ALA about zero |
| `weibull_ALA2.R` | ALA about current mode |
| `weibull_ALA_rb.R` | ALA with Rao-Blackwellised Bayes factors |
| `weibull_NR.R` | Newton-Raphson for MLE |
| `weibull_laplace.R` | Full Laplace approximation |
| `weibull_cpm.R` | Conditional Predictive Marginal |
| `weibull_log_post_beta.R` | Log posterior of regression coefficients |
| `weibull_update_k.R` | Shape parameter update |

### Shared Utilities

| File | Description |
|------|-------------|
| `logit_e.R` | Epsilon-bounded logit transform |
| `logit.R` | Standard logit/inverse-logit |
| `bounded_x.R` | Bounded parameter transforms |
| `mvnorm_trans.R` | Multivariate normal sampling via Cholesky |
| `sample_ind.R` | Independent sampling utilities |
| `make_fixed_variances.R` | Prior variance construction for fixed effects |
| `rpolyagamma.R` | Polya-Gamma random variate generator |

### Random g-Prior

| File | Description |
|------|-------------|
| `hc_update_g.R` | MH update for g under half-Cauchy prior |
| `rhcauchy.R` | Half-Cauchy random variate generator |
| `log_half_cauchy_wj.R` | Log density of half-Cauchy (Wand-Jones parameterisation) |

### Mixing Comparison (`mixing_comparison/`)

| File | Description |
|------|-------------|
| `informed_mixing.R` | Simulation study comparing informed vs uninformed proposals |
| `plot.R` | Plotting utilities for comparison results |

## Usage

```r
# Source all functions
source("main.R")

# Set up hyperparameters
hyper_par <- list(
  n = nrow(X), p = ncol(X), p_z = 1,
  X = X, y = y, Z = Z,
  g = 100, h = ncol(X) / 2,  # or h = c(alpha, beta) for beta-binomial
  model = "logistic",         # "logistic", "cox", or "weibull"
  t = NULL, d = NULL          # survival time and event indicator (Cox/Weibull)
)

# Set up algorithm parameters
alg_par <- list(
  N = 5000, Nb = 2000,
  n_chain = 10,
  verbose = TRUE,
  method = "ALA",             # "ALA", "LA", "CPM", or "DA"
  bal_fun = "hastings",       # "hastings" or "barker"
  kappa = 0.001,
  store_chains = TRUE
)

# Run PARNI sampler
result <- GLM_PARNI(alg_par, hyper_par)

# Or run ADS baseline
result <- GLM_ADS(alg_par, hyper_par)
```

## Dependencies

- Base R (>= 3.6)
- `survival` (for Cox model via `coxph`)
- `rlist` (for mixing comparison scripts)
