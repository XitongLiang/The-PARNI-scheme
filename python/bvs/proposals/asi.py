"""Adaptively Scaled Individual (ASI) proposal.

The ASI proposal samples a random subset of variables to flip simultaneously,
with per-variable probabilities derived from adaptive PIPs. The scaling
parameter zeta is adapted via Robbins-Monro to achieve a target acceptance rate.

Reference: Griffin et al. (2021), Biometrika.
"""

import math
import numpy as np

from .base import BaseProposal
from .adaptation import PIPAdaptation
from ..gamma_par import GammaPar, update_gamma_par
from ..logit_utils import logit_e, inv_logit_e


class ASIProposal(PIPAdaptation, BaseProposal):
    """Adaptively Scaled Individual proposal with PIP adaptation.

    Parameters
    ----------
    p : int
        Number of candidate variables.
    n_temp : int
        Number of temperature levels.
    n_chain : int
        Number of parallel chains.
    N_total : int
        Total iterations (burn-in + sampling).
    likelihood : BaseLikelihood
        Likelihood for evaluating proposals and computing Bayes factors.
    model_prior : ModelPrior
        Model prior.
    h_exp : float
        Expected prior inclusion probability.
    N_adapt_PIPs : int
        Iterations for PIP adaptation.
    N_rb : int
        Iterations for Rao-Blackwellised PIP estimation.
    adapt_phi : float
        Step-size exponent for zeta adaptation.
    target_acc : float
        Target acceptance rate for zeta.
    zeta_init : float
        Initial scaling parameter zeta in (0, 1).
    kappa : float
        PIP floor/ceiling.
    """

    def __init__(self, p, n_temp, n_chain, N_total,
                 likelihood, model_prior, h_exp,
                 N_adapt_PIPs, N_rb,
                 adapt_phi=-0.7, target_acc=0.234, zeta_init=0.5, kappa=0.001):
        self._p = p
        self._adapt_phi = adapt_phi
        self._target_acc = target_acc
        self._logit_eps = 0.1 / p

        self.init_pip_adaptation(n_temp, n_chain, p, h_exp, kappa,
                                 N_adapt_PIPs, N_rb, likelihood, model_prior)

        # ASI-specific: A_til, D_til are unscaled; A, D = zeta * A_til, zeta * D_til
        ratio = h_exp / (1 - h_exp) if h_exp < 1 else 1.0
        self._A_til = np.ones((n_temp, p)) * min(1, ratio)
        inv_ratio = (1 - h_exp) / h_exp if h_exp > 0 else 1.0
        self._D_til = np.ones((n_temp, p)) * min(1, inv_ratio)
        self.A = self._A_til.copy()
        self.D = self._A_til.copy()

        self._temp_acc_rates = np.zeros((n_temp, n_chain))

        logitzeta_init = logit_e(zeta_init, self._logit_eps)
        self._logitzeta = np.ones(n_temp) * logitzeta_init
        self.zetas = np.zeros((n_temp, N_total + 1))
        self.zetas[:, 0] = zeta_init

    @property
    def is_adaptive(self):
        return True

    def _ad_sample(self, new_sample, AD_prob, sample=None):
        """Sample which variables to flip, or compute log-probability of a given set.

        Parameters
        ----------
        new_sample : bool
            If True, draw a new sample. If False, compute probability of `sample`.
        AD_prob : np.ndarray, shape (p,)
            Per-variable flip probabilities.
        sample : np.ndarray or None
            Indices of variables to evaluate (required if new_sample=False).

        Returns
        -------
        sample : np.ndarray
            Indices of variables selected.
        log_prob : float
            Log probability of the sample.
        """
        if new_sample:
            sample = np.where(np.random.random(self._p) < AD_prob)[0]
        log_prob = np.log(AD_prob[sample]).sum()
        return sample, log_prob

    def propose(self, gamma_par, temp=1.0):
        """Propose a new model by flipping a random subset of variables.

        Each variable j is independently selected for flipping with
        probability A[t,j] (if excluded) or D[t,j] (if included),
        where A and D are scaled by zeta.

        Parameters
        ----------
        gamma_par : GammaPar
            Current model state.
        temp : float
            Temperature (unused directly; adaptation uses it).

        Returns
        -------
        prop_gamma_par : GammaPar
            Proposed model.
        log_prop_odd : float
            Log reverse/forward proposal ratio.
        change : int
            Number of variables flipped.
        """
        # Use temperature index stored during update_params context
        t = self._current_t if hasattr(self, '_current_t') else 0

        AD_prob = (1 - gamma_par.gamma) * self.A[t] + gamma_par.gamma * self.D[t]

        which_change, log_prob_prop = self._ad_sample(True, AD_prob)

        gamma_prop = gamma_par.gamma.copy()
        gamma_prop[which_change] = 1 - gamma_prop[which_change]

        # Reverse probability
        AD_prob[which_change] = ((1 - gamma_prop[which_change]) * self.A[t, which_change]
                                 + gamma_prop[which_change] * self.D[t, which_change])
        _, log_prob_rev = self._ad_sample(False, AD_prob, which_change)

        prop_gamma_par = GammaPar(gamma_prop)
        update_gamma_par(prop_gamma_par, self._likelihood, self._model_prior)

        return prop_gamma_par, log_prob_rev - log_prob_prop, which_change.size

    def update_params(self, gamma_par, acc_rate, change,
                      t=0, i=0, epoch=0, temp=1.0):
        """Adapt PIPs and zeta after an MH step.

        Parameters
        ----------
        gamma_par : GammaPar
            Accepted model state.
        acc_rate : float
            MH acceptance probability.
        change : int
            Number of flipped variables.
        t : int
            Temperature index.
        i : int
            Chain index.
        epoch : int
            Current iteration.
        temp : float
            Temperature value.

        Side Effects
        ------------
        Updates adapt_PIPs, A, D, and zeta.
        """
        self._current_t = t
        self._temp_acc_rates[t, i] = acc_rate

        self._update_pips_for_chain(gamma_par, t, i, epoch, temp)

        # End-of-sweep aggregation
        if i == self._n_chain - 1:
            self._aggregate_pips(t, epoch)

            if epoch <= self._N_adapt_PIPs:
                self._A_til[t] = np.minimum(1, self.adapt_PIPs[t] / (1 - self.adapt_PIPs[t]))
                self._D_til[t] = np.minimum(1, (1 - self.adapt_PIPs[t]) / self.adapt_PIPs[t])

            self._logitzeta[t] += epoch ** self._adapt_phi * (
                self._temp_acc_rates[t].mean() - self._target_acc
            )

            zeta = inv_logit_e(self._logitzeta[t], self._logit_eps)
            self.zetas[t, epoch] = zeta

            self.A[t] = zeta * self._A_til[t]
            self.D[t] = zeta * self._D_til[t]
