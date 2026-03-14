"""Model indicator dataclass for Bayesian variable selection.

Stores the binary inclusion vector gamma and derived quantities
(model size, included/excluded indices) along with computed
statistics (log-likelihood, log model prior, Bayes factors).
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class GammaPar:
    """State object for a single variable selection model.

    Parameters
    ----------
    gamma : np.ndarray
        Binary inclusion vector of length p. gamma[j] = 1 means
        variable j is included in the model.

    Attributes (computed on init)
    ----------------------------
    p_gam : int
        Model size (number of included variables).
    includes : np.ndarray
        Indices of included variables.
    excludes : np.ndarray
        Indices of excluded variables.

    Attributes (set by likelihood / prior / Bayes factor methods)
    -------------------------------------------------------------
    log_llh : float or None
        Log marginal likelihood of this model.
    log_mp : float or None
        Log model prior probability.
    inv_V_gam : np.ndarray or None
        Inverse of the posterior precision sub-matrix for included variables.
    BF : np.ndarray or None
        Bayes factor vector of length p for single-variable changes.
    h_til : float or np.ndarray or None
        Tilted model prior parameter used in Rao-Blackwellised PIP estimation.
    """

    gamma: np.ndarray
    p_gam: int = field(init=False)
    includes: np.ndarray = field(init=False)
    excludes: np.ndarray = field(init=False)
    log_llh: Optional[float] = None
    log_mp: Optional[float] = None
    inv_V_gam: Optional[np.ndarray] = None
    BF: Optional[np.ndarray] = None
    h_til: object = None  # float or np.ndarray

    def __post_init__(self):
        self.p_gam = int(self.gamma.sum())
        self.includes = np.where(self.gamma)[0]
        self.excludes = np.where(self.gamma == 0)[0]


def update_gamma_par(gamma_par, likelihood, model_prior):
    """Compute both log-likelihood and log model prior for a model.

    Parameters
    ----------
    gamma_par : GammaPar
        Model indicator to update. Modified in place.
    likelihood : BaseLikelihood
        Likelihood object providing log_llh().
    model_prior : ModelPrior
        Model prior object providing log_m_prior().

    Side Effects
    ------------
    Sets gamma_par.log_llh and gamma_par.log_mp.
    """
    likelihood.log_llh(gamma_par)
    model_prior.log_m_prior(gamma_par)
