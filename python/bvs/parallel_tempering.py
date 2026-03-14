"""Parallel tempering for multi-modal posterior exploration.

Maintains a ladder of temperatures and adaptively tunes the temperature
schedule to achieve a target swap acceptance rate.
"""

import math
import numpy as np


class ParallelTempering:
    """Adaptive parallel tempering with temperature schedule learning.

    Parameters
    ----------
    n_temp : int
        Number of temperature levels (must be >= 2).
    n_chain : int
        Number of parallel chains per temperature.
    N_total : int
        Total number of MCMC iterations (burn-in + sampling).
    PT_taus : np.ndarray or None
        Initial log-log temperature spacings. Shape (n_temp - 1,).
        Defaults to -4 * ones.
    adapt_phi : float
        Step-size exponent for Robbins-Monro adaptation of taus.
    target_acc : float
        Target swap acceptance rate.

    Attributes
    ----------
    temps : np.ndarray
        Current temperature schedule, shape (n_temp,). temps[0] = 1.0.
    temperatures : np.ndarray
        Full trace of temperatures, shape (n_temp, N_total + 1).
    acc_times : int
        Number of accepted swaps (after burn-in).
    """

    def __init__(self, n_temp, n_chain, N_total,
                 PT_taus=None, adapt_phi=-0.7, target_acc=0.234):
        self.n_temp = n_temp
        self.n_chain = n_chain
        self.target_acc = target_acc
        self.adapt_phi = adapt_phi

        if PT_taus is not None:
            self.taus = PT_taus.copy()
        else:
            self.taus = -4 * np.ones(n_temp - 1)

        self._update_temps()
        self.temperatures = np.zeros((n_temp, N_total + 1))
        self.temperatures[:, 0] = self.temps

        self.acc_times = 0
        self.H = np.zeros((n_chain, n_temp - 1))

    def _update_temps(self):
        """Recompute temperature schedule from current taus."""
        self.temps = np.append(1, np.exp(-np.exp(self.taus))).cumprod()

    def swap_chains(self, agg_gamma_par, i, epoch, N_burnin):
        """Attempt a swap between adjacent temperature levels for chain i.

        Randomly selects an adjacent pair and performs a Metropolis-Hastings
        accept/reject step based on the likelihood difference.

        Parameters
        ----------
        agg_gamma_par : dict
            Dictionary mapping (t, i) -> GammaPar for all temperatures and chains.
        i : int
            Chain index.
        epoch : int
            Current iteration number.
        N_burnin : int
            Number of burn-in iterations.

        Side Effects
        ------------
        May swap agg_gamma_par entries between adjacent temperatures.
        Increments self.acc_times if swap is accepted and epoch > N_burnin.
        """
        swap_idx = np.random.randint(self.n_temp - 1)
        swap_acc_rate = math.exp(
            (self.temps[swap_idx] - self.temps[swap_idx + 1])
            * (agg_gamma_par[swap_idx + 1, i].log_llh - agg_gamma_par[swap_idx, i].log_llh)
        )

        if np.random.random() < swap_acc_rate:
            temp_gp = agg_gamma_par[swap_idx + 1, i]
            agg_gamma_par[swap_idx + 1, i] = agg_gamma_par[swap_idx, i]
            agg_gamma_par[swap_idx, i] = temp_gp

            if epoch > N_burnin:
                self.acc_times += 1

    def update_H(self, i, t, log_llh):
        """Accumulate log-likelihood differences for temperature adaptation.

        Called once per (chain, temperature) during the main loop to build
        the swap statistic H used in adapt_temperatures().

        Parameters
        ----------
        i : int
            Chain index.
        t : int
            Temperature index.
        log_llh : float
            Log-likelihood of the current model at temperature t.
        """
        if t == 0:
            self.H[i, t] = -log_llh
        elif t == self.n_temp - 1:
            self.H[i, t - 1] += log_llh
        else:
            self.H[i, t] = -log_llh
            self.H[i, t - 1] += log_llh

    def adapt_temperatures(self, epoch):
        """Adapt the temperature schedule using Robbins-Monro.

        Parameters
        ----------
        epoch : int
            Current iteration number (used for step-size decay).

        Side Effects
        ------------
        Updates self.taus, self.temps, and records in self.temperatures.
        """
        self.H *= self.temps[:(self.n_temp - 1)] - self.temps[1:]
        self.taus += epoch ** self.adapt_phi * (
            np.minimum(1, np.exp(self.H)).mean(axis=0) - self.target_acc
        )
        self._update_temps()
        self.temperatures[:, epoch] = self.temps


class NullTempering:
    """No-op tempering for single-temperature runs.

    Provides the same interface as ParallelTempering but does nothing.

    Attributes
    ----------
    n_temp : int
        Always 1.
    temps : np.ndarray
        [1.0].
    acc_times : int
        Always 0.
    """

    def __init__(self):
        self.n_temp = 1
        self.temps = np.ones(1)
        self.acc_times = 0
        self.temperatures = None

    def swap_chains(self, agg_gamma_par, i, epoch, N_burnin):
        """No-op: no swaps with a single temperature."""
        pass

    def update_H(self, i, t, log_llh):
        """No-op."""
        pass

    def adapt_temperatures(self, epoch):
        """No-op."""
        pass
