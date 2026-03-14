"""Shared PIP (Posterior Inclusion Probability) adaptation logic.

Both ASI and PARNI proposals adapt their neighbourhood probabilities
using estimates of the marginal posterior inclusion probabilities.
During the Rao-Blackwellisation phase, PIPs are estimated via Bayes
factors; afterwards, the ergodic average of gamma is used.

This mixin provides the common initialisation and update logic.
"""

import numpy as np


class PIPAdaptation:
    """Mixin providing adaptive PIP estimation for proposals.

    Maintains running estimates of posterior inclusion probabilities
    and derives the add/delete probability vectors A and D.

    Attributes
    ----------
    adapt_PIPs : np.ndarray, shape (n_temp, p)
        Current PIP estimates (clipped by kappa).
    A : np.ndarray, shape (n_temp, p)
        Add probabilities: min(1, PIP / (1-PIP)).
    D : np.ndarray, shape (n_temp, p)
        Delete probabilities: min(1, (1-PIP) / PIP).
    """

    def init_pip_adaptation(self, n_temp, n_chain, p, h_exp, kappa,
                            N_adapt_PIPs, N_rb, likelihood, model_prior):
        """Initialise PIP adaptation state.

        Parameters
        ----------
        n_temp : int
            Number of temperature levels.
        n_chain : int
            Number of parallel chains.
        p : int
            Number of candidate variables.
        h_exp : float
            Expected prior inclusion probability.
        kappa : float
            Floor/ceiling for PIPs to prevent degeneracy (e.g., 0.001).
        N_adapt_PIPs : int
            Number of iterations to adapt PIPs.
        N_rb : int
            Number of iterations to use Rao-Blackwellised estimates
            (afterwards switches to ergodic average).
        likelihood : BaseLikelihood
            Likelihood object for Bayes factor computation.
        model_prior : ModelPrior
            Model prior for h_til computation.
        """
        self._n_temp = n_temp
        self._n_chain = n_chain
        self._p = p
        self._kappa = kappa
        self._N_adapt_PIPs = N_adapt_PIPs
        self._N_rb = N_rb
        self._likelihood = likelihood
        self._model_prior = model_prior

        self.adapt_PIPs = np.ones((n_temp, p)) * h_exp
        self._sum_adapt_PIPs = np.zeros((n_temp, p))

        # A = min(1, PIP/(1-PIP)), D = min(1, (1-PIP)/PIP)
        ratio = h_exp / (1 - h_exp) if h_exp < 1 else 1.0
        self.A = np.ones((n_temp, p)) * min(1, ratio)
        inv_ratio = (1 - h_exp) / h_exp if h_exp > 0 else 1.0
        self.D = np.ones((n_temp, p)) * min(1, inv_ratio)

        self._temp_PIPs = np.zeros((n_temp, n_chain, p))

    def _update_pips_for_chain(self, gamma_par, t, i, epoch, temp):
        """Update the per-chain PIP estimate for one (t, i) pair.

        During Rao-Blackwellisation phase (epoch <= N_rb), computes
        Bayes-factor-based PIP estimates. Afterwards, uses the raw
        gamma vector as the PIP estimate.

        Parameters
        ----------
        gamma_par : GammaPar
            Current model state.
        t : int
            Temperature index.
        i : int
            Chain index.
        epoch : int
            Current iteration.
        temp : float
            Temperature value at level t.
        """
        if epoch > self._N_adapt_PIPs:
            return

        if epoch <= self._N_rb:
            if gamma_par.BF is None:
                self._model_prior.h_til(gamma_par)
                self._likelihood.compute_bf(gamma_par)

                if t == 0:
                    self._temp_PIPs[t, i] = ((gamma_par.h_til * gamma_par.BF)
                                             / (1 - gamma_par.h_til + gamma_par.h_til * gamma_par.BF))
                else:
                    BF_temp = gamma_par.BF ** temp
                    self._temp_PIPs[t, i] = ((gamma_par.h_til * BF_temp)
                                             / (1 - gamma_par.h_til + gamma_par.h_til * BF_temp))
        else:
            self._temp_PIPs[t, i] = gamma_par.gamma

    def _aggregate_pips(self, t, epoch):
        """Aggregate per-chain PIPs into the adaptive PIP estimate.

        Called once per sweep (after all chains processed) to update
        adapt_PIPs, A, and D for temperature t.

        Parameters
        ----------
        t : int
            Temperature index.
        epoch : int
            Current iteration.
        """
        if epoch > self._N_adapt_PIPs:
            return

        self._sum_adapt_PIPs[t] += self._temp_PIPs[t].sum(axis=0)
        self.adapt_PIPs[t] = (self._kappa
                              + (1 - 2 * self._kappa)
                              * self._sum_adapt_PIPs[t] / (epoch * self._n_chain))

        self.A[t] = np.minimum(1, self.adapt_PIPs[t] / (1 - self.adapt_PIPs[t]))
        self.D[t] = np.minimum(1, (1 - self.adapt_PIPs[t]) / self.adapt_PIPs[t])
