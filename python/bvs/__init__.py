"""BVS: Bayesian Variable Selection via adaptive MCMC.

This package implements the PARNI, ASI, and ADS samplers for
high-dimensional Bayesian variable selection with optional
parallel tempering and random g-prior.

Quick start
-----------
>>> from bvs import BVS_MCMC
>>> fit = BVS_MCMC(X, y, g=100, prior_type="ind")
>>> fit.set_alg_par(sampler="PARNI", N_iter=2000, N_burnin=1000)
>>> fit.sample_now()
>>> print(fit.estm_PIPs[0])
"""

from .sampler import BVS_MCMC
from .gamma_par import GammaPar
from .plotting import plot_temperatures, plot_pips

__all__ = ["BVS_MCMC", "GammaPar", "plot_temperatures", "plot_pips"]
