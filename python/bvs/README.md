# BVS Package — Bayesian Variable Selection via Adaptive MCMC

## Overview

This package implements MCMC algorithms for Bayesian variable selection in linear regression, with a modular design that makes it straightforward to add new likelihood models and proposal mechanisms.

## Module Reference

| Module | Role |
|--------|------|
| `sampler.py` | **BVS_MCMC** orchestrator — assembles components and runs the main loop |
| `gamma_par.py` | **GammaPar** dataclass — stores a model indicator and computed statistics |
| `model_prior.py` | **ModelPrior** — fixed-h or beta-binomial prior on model inclusion |
| `likelihoods/base.py` | **BaseLikelihood** ABC — interface for log-likelihood and Bayes factors |
| `likelihoods/linear.py` | **LinearGPrior**, **LinearIndPrior** — linear regression implementations |
| `proposals/base.py` | **BaseProposal** ABC — interface for MCMC proposals |
| `proposals/adaptation.py` | **PIPAdaptation** mixin — shared PIP adaptation logic |
| `proposals/ads.py` | **ADSProposal** — Add-Delete-Swap baseline |
| `proposals/asi.py` | **ASIProposal** — Adaptively Scaled Individual proposal |
| `proposals/parni.py` | **PARNIProposal** — PARNI proposal with RM or KW adaptation |
| `parallel_tempering.py` | **ParallelTempering** / **NullTempering** — temperature-based exploration |
| `g_prior.py` | **RandomGSampler** — MH sampler for random g with half-Cauchy prior |
| `logit_utils.py` | Epsilon-bounded logit/inverse-logit transforms |
| `plotting.py` | `plot_temperatures()`, `plot_pips()` — visualisation utilities |

## How Components Connect

```
BVS_MCMC (sampler.py)
├── likelihood  (likelihoods/)  — evaluates log p(y | gamma)
├── model_prior (model_prior.py) — evaluates log p(gamma)
├── proposal    (proposals/)     — generates candidate models
│   └── uses PIPAdaptation mixin for adaptive PIPs
├── tempering   (parallel_tempering.py) — manages temperature ladder
└── g_sampler   (g_prior.py)     — updates g if random
```

The main sampling loop in `sampler.py` calls these components in sequence:
1. **Swap** between temperature levels (parallel tempering)
2. **Update g** if using random g-prior
3. **Propose** a new model
4. **Accept/reject** via Metropolis-Hastings
5. **Adapt** proposal parameters (PIPs, omega/zeta)
6. **Adapt temperatures** at the end of each sweep

## Adding a New Likelihood Model

To add support for a new regression model (e.g., logistic, Cox):

1. Create a new file `likelihoods/logistic.py`.
2. Subclass `BaseLikelihood`:
   ```python
   from .base import BaseLikelihood

   class LogisticLikelihood(BaseLikelihood):
       def __init__(self, X, y, ...):
           self.n, self.p = X.shape
           self.random_g = False
           # Pre-compute sufficient statistics...

       def log_llh(self, gamma_par):
           """Compute log marginal likelihood, store in gamma_par.log_llh."""
           ...

       def compute_bf(self, gamma_par):
           """Compute Bayes factors, store in gamma_par.BF."""
           ...
   ```
3. Register it in `likelihoods/__init__.py`.
4. Add a branch in `sampler._build_likelihood()` for `model="logistic"`.

## Adding a New Proposal

1. Create a new file `proposals/my_proposal.py`.
2. Subclass `BaseProposal`:
   ```python
   from .base import BaseProposal
   from ..gamma_par import GammaPar, update_gamma_par

   class MyProposal(BaseProposal):
       def __init__(self, p, likelihood, model_prior, ...):
           ...

       def propose(self, gamma_par, temp=1.0):
           """Return (proposed_gamma_par, log_prop_odd, change_size)."""
           ...

       # If adaptive:
       @property
       def is_adaptive(self):
           return True

       def update_params(self, gamma_par, acc_rate, change,
                         t=0, i=0, epoch=0, temp=1.0):
           ...
   ```
3. If it needs PIP adaptation, use the `PIPAdaptation` mixin:
   ```python
   from .adaptation import PIPAdaptation

   class MyProposal(PIPAdaptation, BaseProposal):
       def __init__(self, ...):
           self.init_pip_adaptation(n_temp, n_chain, p, h_exp, kappa,
                                    N_adapt_PIPs, N_rb, likelihood, model_prior)
           ...
   ```
4. Register it in `proposals/__init__.py`.
5. Add a branch in `sampler._build_proposal()`.

## Quick Start

```python
import sys
sys.path.insert(0, '..')  # if running from python/ directory
from bvs import BVS_MCMC
from utils import lrrsg

# Generate data
beta, y, X = lrrsg(n=1000, p=20)

# Fit
fit = BVS_MCMC(X, y, g=100, prior_type="ind")
fit.set_alg_par(sampler="PARNI", N_iter=2000, N_burnin=1000,
                n_chain=5, n_temp=3, verbose=True)
fit.sample_now()

print(f"PIPs: {fit.estm_PIPs[0]}")
print(f"Time: {fit.time_total:.1f}s")
```
