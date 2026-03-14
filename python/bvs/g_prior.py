"""Random g-prior sampler via Metropolis-Hastings.

When g is treated as a random variable with a half-Cauchy prior,
this module provides the MH update step that proposes new g values
on the log scale and adapts the proposal standard deviation.
"""

import math
import numpy as np
from scipy import stats

from .gamma_par import GammaPar


class RandomGSampler:
    """Metropolis-Hastings sampler for the g parameter.

    Uses a log-normal random walk proposal with adaptive step size
    targeting a specified acceptance rate.

    Parameters
    ----------
    n_temp : int
        Number of temperature levels.
    n_chain : int
        Number of parallel chains.
    N_total : int
        Total iterations (burn-in + sampling).
    g_init : float or None
        Initial g value. If None, draws from half-Cauchy.
    adapt_phi : float
        Step-size exponent for Robbins-Monro adaptation.
    target_acc : float
        Target acceptance rate for the g update.

    Attributes
    ----------
    g_matrix : np.ndarray, shape (n_temp, n_chain)
        Current g values for each (temperature, chain).
    gs : np.ndarray, shape (n_temp, n_chain, N_total + 1)
        Full trace of g values.
    acc_times : np.ndarray, shape (n_temp,)
        Accepted g updates per temperature (after burn-in).
    """

    def __init__(self, n_temp, n_chain, N_total,
                 g_init=None, adapt_phi=-0.7, target_acc=0.234):
        self.n_temp = n_temp
        self.n_chain = n_chain
        self.adapt_phi = adapt_phi
        self.target_acc = target_acc

        if isinstance(g_init, (int, float)):
            self.g_matrix = g_init * np.ones((n_temp, n_chain))
        else:
            self.g_matrix = stats.halfcauchy.rvs(size=(n_temp, n_chain))

        self.gs = np.zeros((n_temp, n_chain, N_total + 1))
        self.gs[:, :, 0] = self.g_matrix
        self.logsd = np.ones(n_temp)
        self.acc_times = np.ones(n_temp)

    @staticmethod
    def log_half_cauchy(g):
        """Log density of the half-Cauchy prior (with Jacobian for log-transform).

        Parameters
        ----------
        g : float
            Positive g value.

        Returns
        -------
        float
            log p(g) + log |dg/d(log g)| = log(g) - log(1 + g^2).
        """
        return math.log(g) - math.log(1 + g ** 2)

    def update(self, gamma_par, likelihood, t, i, epoch, temp, N_burnin):
        """Perform one MH step for g at temperature t, chain i.

        Proposes a new g via log-normal random walk, evaluates the
        acceptance ratio using the tempered likelihood and half-Cauchy prior,
        and adapts the proposal standard deviation.

        Parameters
        ----------
        gamma_par : GammaPar
            Current model state. May be replaced if g update changes likelihood.
        likelihood : BaseLikelihood
            Likelihood object whose .g attribute will be modified on acceptance.
        t : int
            Temperature index.
        i : int
            Chain index.
        epoch : int
            Current iteration number.
        temp : float
            Temperature value at level t.
        N_burnin : int
            Number of burn-in iterations.

        Returns
        -------
        GammaPar
            Updated gamma_par (may be a new object if g was accepted).

        Side Effects
        ------------
        Modifies likelihood.g on acceptance. Updates self.g_matrix, self.gs,
        self.logsd, and self.acc_times.
        """
        old_g = self.g_matrix[t, i]
        likelihood.g = old_g

        # Propose new g on log scale
        new_g = old_g * math.exp(math.exp(self.logsd[t]) * np.random.normal())

        if new_g == 0 or np.isinf(new_g):
            g_acc_rate = 0
        else:
            likelihood.g = new_g
            new_gamma_par = GammaPar(gamma_par.gamma.copy())
            new_gamma_par.log_mp = gamma_par.log_mp
            likelihood.log_llh(new_gamma_par)

            g_acc_rate = min(1, math.exp(
                temp * (new_gamma_par.log_llh - gamma_par.log_llh)
                + self.log_half_cauchy(new_g) - self.log_half_cauchy(old_g)
            ))

        if np.random.random() < g_acc_rate:
            gamma_par = new_gamma_par
            self.g_matrix[t, i] = new_g
            if epoch > N_burnin:
                self.acc_times[t] += 1
        else:
            likelihood.g = old_g
            self.g_matrix[t, i] = old_g

        self.gs[t, i, epoch] = self.g_matrix[t, i]
        self.logsd[t] += epoch ** self.adapt_phi * (g_acc_rate - self.target_acc)

        return gamma_par
