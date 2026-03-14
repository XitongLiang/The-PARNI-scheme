"""Likelihood models for Bayesian variable selection.

Each likelihood class implements log_llh() and compute_bf() for a specific
model type and prior structure. To add a new model, subclass BaseLikelihood.
"""

from .linear import LinearGPrior, LinearIndPrior

__all__ = ["LinearGPrior", "LinearIndPrior"]
