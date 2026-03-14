"""Point-wise Adaptive Random Neighbourhood Informed (PARNI) proposal.

PARNI selects a random neighbourhood of variables (via adaptive PIPs),
then sequentially decides whether to flip each selected variable using
a locally informed balancing function that incorporates the posterior
odds ratio. The proposal scaling parameter omega controls the
change-vs-keep trade-off and is adapted via either Robbins-Monro (RM)
or Kiefer-Wolfowitz (KW) stochastic approximation.

Reference: Liang, Livingstone, Griffin (2022), Statistics and Computing.
"""

import math
import numpy as np

from .base import BaseProposal
from .adaptation import PIPAdaptation
from ..gamma_par import GammaPar, update_gamma_par
from ..logit_utils import logit_e, inv_logit_e, inv_logit_e_vec


def _hastings(x):
    """Hastings balancing function: min(1, x).

    Parameters
    ----------
    x : float
        Non-negative input.

    Returns
    -------
    float
    """
    return min(1, x)


def _barker(x):
    """Barker balancing function: x / (1 + x).

    Parameters
    ----------
    x : float
        Non-negative input.

    Returns
    -------
    float
    """
    return x / (1 + x)


class PARNIProposal(PIPAdaptation, BaseProposal):
    """PARNI proposal with adaptive omega and PIPs.

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
        Likelihood for evaluating proposals and Bayes factors.
    model_prior : ModelPrior
        Model prior.
    h_exp : float
        Expected prior inclusion probability.
    N_adapt_PIPs : int
        Iterations for PIP adaptation.
    N_rb : int
        Iterations for Rao-Blackwellised estimation.
    omega_init : float
        Initial omega value in (0, 1).
    bal_fun : str or callable
        Balancing function. "hastings", "barker", "sqrt", or a callable.
    omega_adapt : str
        Adaptation scheme: "RM" (Robbins-Monro) or "KW" (Kiefer-Wolfowitz).
    adapt_phi : float or list of two floats
        For RM: single exponent. For KW: [phi_a, phi_c].
    target_acc : float
        Target acceptance rate (RM only).
    kappa : float
        PIP floor/ceiling.
    """

    def __init__(self, p, n_temp, n_chain, N_total,
                 likelihood, model_prior, h_exp,
                 N_adapt_PIPs, N_rb,
                 omega_init=0.5, bal_fun="hastings", omega_adapt="RM",
                 adapt_phi=-0.7, target_acc=0.65, kappa=0.001):
        self._p = p
        self._logit_eps = 0.1 / p
        self._omega_adapt = omega_adapt

        # Balancing function
        if bal_fun == "hastings":
            self._bal_fun = _hastings
        elif bal_fun == "barker":
            self._bal_fun = _barker
        elif bal_fun == "sqrt":
            self._bal_fun = math.sqrt
        else:
            self._bal_fun = bal_fun

        self.init_pip_adaptation(n_temp, n_chain, p, h_exp, kappa,
                                 N_adapt_PIPs, N_rb, likelihood, model_prior)

        if omega_adapt == "RM":
            self._init_rm(n_temp, n_chain, N_total, omega_init, adapt_phi, target_acc)
        elif omega_adapt == "KW":
            if n_chain == 1:
                raise ValueError("KW adaptation requires n_chain > 1.")
            self._init_kw(n_temp, n_chain, N_total, omega_init, adapt_phi)

    @property
    def is_adaptive(self):
        return True

    # ----------------------------------------------------------------
    # Robbins-Monro initialisation
    # ----------------------------------------------------------------

    def _init_rm(self, n_temp, n_chain, N_total, omega_init, adapt_phi, target_acc):
        """Initialise Robbins-Monro omega adaptation.

        Parameters
        ----------
        n_temp, n_chain, N_total : int
        omega_init : float
            Initial omega.
        adapt_phi : float
            Step-size exponent.
        target_acc : float
            Target acceptance rate.
        """
        self._adapt_phi = adapt_phi
        self._target_acc = target_acc
        self._temp_acc_rates = np.zeros((n_temp, n_chain))

        logitomega_init = logit_e(omega_init, self._logit_eps)
        self._logitomega = np.ones(n_temp) * logitomega_init

        self.omega = np.ones((n_temp, n_chain)) * omega_init
        self.omegas = np.zeros((n_temp, N_total + 1))
        self.omegas[:, 0] = omega_init

    # ----------------------------------------------------------------
    # Kiefer-Wolfowitz initialisation
    # ----------------------------------------------------------------

    def _init_kw(self, n_temp, n_chain, N_total, omega_init, adapt_phi):
        """Initialise Kiefer-Wolfowitz omega adaptation.

        Uses two groups of chains with perturbed omega values to
        estimate the gradient of ASJD with respect to omega.

        Parameters
        ----------
        n_temp, n_chain, N_total : int
        omega_init : float
            Initial omega.
        adapt_phi : list of two floats
            [phi_a, phi_c] for KW step sizes.
        """
        self._adapt_phi_a = adapt_phi[0]
        self._adapt_phi_c = adapt_phi[1]

        logitomega_init = logit_e(omega_init, self._logit_eps)
        self._logitomega = np.ones(n_temp) * logitomega_init
        omega_vec = inv_logit_e_vec(logitomega_init + np.array([0, -1, 1]), self._logit_eps)

        self._n_pos = math.floor(n_chain / 2)
        pos_idx = np.zeros(n_chain, dtype=bool)
        pos_idx[0:self._n_pos] = True
        self._pos_idx = np.random.permutation(pos_idx)

        self.omega = np.ones((n_temp, n_chain))
        self.omega[:, self._pos_idx] = omega_vec[2]
        self.omega[:, ~self._pos_idx] = omega_vec[1]

        self.omegas = np.zeros((n_temp, N_total + 1))
        self.omegas[:, 0] = omega_vec[0]

        self._ASJD_temp = np.zeros((n_temp, n_chain))

    # ----------------------------------------------------------------
    # Neighbourhood sampling
    # ----------------------------------------------------------------

    def _ad_sample(self, new_sample, AD_prob, sample=None):
        """Sample or evaluate a variable subset.

        Parameters
        ----------
        new_sample : bool
            If True, draw new sample. If False, compute log-prob of `sample`.
        AD_prob : np.ndarray, shape (p,)
            Per-variable selection probabilities.
        sample : np.ndarray or None
            Variable indices (required if new_sample=False).

        Returns
        -------
        sample : np.ndarray
            Selected variable indices.
        log_prob : float
            Log probability of the selection.
        """
        if new_sample:
            sample = np.where(np.random.random(self._p) < AD_prob)[0]
        log_prob = np.log(AD_prob[sample]).sum()
        return sample, log_prob

    # ----------------------------------------------------------------
    # Propose
    # ----------------------------------------------------------------

    def propose(self, gamma_par, temp=1.0):
        """Propose a new model via sequential locally-informed flips.

        For each variable in the randomly selected neighbourhood:
        1. Evaluate the posterior odds of flipping variable j.
        2. Weight by the adaptive PIP odds and balancing function.
        3. Flip with probability proportional to omega * bal(odds).
        4. Accumulate the normalising constants for the MH ratio.

        Parameters
        ----------
        gamma_par : GammaPar
            Current model state.
        temp : float
            Temperature for tempered posterior evaluation.

        Returns
        -------
        prop_gamma_par : GammaPar
            Proposed model.
        log_prop_odd : float
            Log ratio for MH correction (includes balancing constants
            and posterior difference, cancelling in the standard ratio).
        change : int
            Number of variables actually flipped.
        """
        t = self._current_t if hasattr(self, '_current_t') else 0
        i = self._current_i if hasattr(self, '_current_i') else 0
        omega = self.omega[t, i]

        AD_prob = (1 - gamma_par.gamma) * self.A[t] + gamma_par.gamma * self.D[t]
        k, _ = self._ad_sample(True, AD_prob)
        k = np.random.permutation(k)

        JD = 0
        prod_bal_const = 0
        rev_prod_bal_const = 0

        temp_gamma_par = gamma_par
        log_post_curr = temp * temp_gamma_par.log_llh + temp_gamma_par.log_mp
        log_post_temp = log_post_curr

        for kj in k:
            temp_gamma_prop = temp_gamma_par.gamma.copy()
            temp_kj = temp_gamma_prop[kj]
            temp_gamma_prop[kj] = 1 - temp_kj

            temp_prop_gp = GammaPar(temp_gamma_prop)
            update_gamma_par(temp_prop_gp, self._likelihood, self._model_prior)

            log_post_temp_prop = temp * temp_prop_gp.log_llh + temp_prop_gp.log_mp
            post_prop_temp = math.exp(log_post_temp_prop - log_post_temp)

            mar_eff = self.adapt_PIPs[t, kj]
            odd_kj = (mar_eff / (1 - mar_eff)) ** (2 * temp_kj - 1)

            change_prob = omega * self._bal_fun(post_prop_temp * odd_kj)
            keep_prob = (1 - omega) * self._bal_fun(1)

            bal_const = change_prob + keep_prob
            change_prob /= bal_const
            keep_prob /= bal_const

            if np.random.random() < change_prob:
                rev_change_prob = omega * self._bal_fun(1 / post_prop_temp / odd_kj)
                rev_keep_prob = (1 - omega) * self._bal_fun(1)

                rev_bal_const = rev_change_prob + rev_keep_prob

                temp_gamma_par = temp_prop_gp
                log_post_temp = log_post_temp_prop
                JD += 1

            else:
                rev_bal_const = bal_const

            prod_bal_const += math.log(bal_const)
            rev_prod_bal_const += math.log(rev_bal_const)

        return (temp_gamma_par,
                prod_bal_const - rev_prod_bal_const + log_post_curr - log_post_temp,
                JD)

    # ----------------------------------------------------------------
    # Parameter updates
    # ----------------------------------------------------------------

    def update_params(self, gamma_par, acc_rate, change,
                      t=0, i=0, epoch=0, temp=1.0):
        """Adapt omega and PIPs after an MH step.

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
        Updates adapt_PIPs, A, D, omega, and omegas trace.
        """
        self._current_t = t
        self._current_i = i

        self._update_pips_for_chain(gamma_par, t, i, epoch, temp)

        if self._omega_adapt == "RM":
            self._update_rm(gamma_par, acc_rate, change, t, i, epoch)
        elif self._omega_adapt == "KW":
            self._update_kw(gamma_par, acc_rate, change, t, i, epoch)

    def _update_rm(self, gamma_par, acc_rate, change, t, i, epoch):
        """Robbins-Monro omega adaptation.

        Adjusts logit(omega) towards achieving the target acceptance rate.

        Parameters
        ----------
        gamma_par : GammaPar
        acc_rate : float
        change : int
        t, i, epoch : int
        """
        self._temp_acc_rates[t, i] = acc_rate

        if i == self._n_chain - 1:
            self._aggregate_pips(t, epoch)

            self._logitomega[t] += epoch ** self._adapt_phi * (
                self._temp_acc_rates[t].mean() - self._target_acc
            )

            omega = inv_logit_e(self._logitomega[t], self._logit_eps)
            self.omegas[t, epoch] = omega
            self.omega[t, :] = omega

    def _update_kw(self, gamma_par, acc_rate, change, t, i, epoch):
        """Kiefer-Wolfowitz omega adaptation.

        Compares ASJD between positive and negative perturbation groups
        to estimate the gradient and adjust omega.

        Parameters
        ----------
        gamma_par : GammaPar
        acc_rate : float
        change : int
        t, i, epoch : int
        """
        self._ASJD_temp[t, i] = change * acc_rate

        if i == self._n_chain - 1:
            self._aggregate_pips(t, epoch)

            self._logitomega[t] += (
                epoch ** self._adapt_phi_a
                * (self._ASJD_temp[t, self._pos_idx].mean()
                   - self._ASJD_temp[t, ~self._pos_idx].mean())
                / (2 * epoch ** self._adapt_phi_c)
            )
            omega_vec = inv_logit_e_vec(
                self._logitomega[t] + np.array([0,
                                                -(epoch + 1) ** self._adapt_phi_c,
                                                (epoch + 1) ** self._adapt_phi_c]),
                self._logit_eps
            )
            self.omegas[t, epoch] = omega_vec[0]

            if t == 0:
                self._new_pos_idx = np.random.permutation(self._pos_idx)

            self.omega[t, self._new_pos_idx] = omega_vec[2]
            self.omega[t, ~self._new_pos_idx] = omega_vec[1]

            if t == self._n_temp - 1:
                self._pos_idx = self._new_pos_idx
