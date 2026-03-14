"""BVS_MCMC: Main orchestrator for Bayesian Variable Selection via MCMC.

This module assembles the modular components (likelihood, model prior,
proposal, parallel tempering, random g) into a single user-facing class.
The sampling loop is a clean orchestration of these components.

Usage
-----
>>> from bvs import BVS_MCMC
>>> fit = BVS_MCMC(X, y, g=100, prior_type="ind")
>>> fit.set_alg_par(sampler="PARNI", N_iter=2000, N_burnin=1000, n_chain=5)
>>> fit.sample_now()
>>> print(fit.estm_PIPs[0])
"""

import math
import time
import numpy as np
from tqdm import tqdm

from .gamma_par import GammaPar, update_gamma_par
from .model_prior import ModelPrior
from .likelihoods.linear import LinearGPrior, LinearIndPrior
from .proposals.ads import ADSProposal
from .proposals.asi import ASIProposal
from .proposals.parni import PARNIProposal
from .parallel_tempering import ParallelTempering, NullTempering
from .g_prior import RandomGSampler


def _build_likelihood(model, prior_type, X, y, g, **kwargs):
    """Construct the appropriate likelihood object.

    Parameters
    ----------
    model : str
        Model type ("linear"). Future: "logistic", "cox", "weibull".
    prior_type : str
        Prior structure ("g" or "ind").
    X : np.ndarray
        Design matrix.
    y : np.ndarray
        Response vector.
    g : float or str
        Prior variance parameter.
    **kwargs
        Passed to the likelihood constructor (Z, scale, ddof).

    Returns
    -------
    BaseLikelihood
        Configured likelihood object.
    """
    Z = kwargs.get('Z', None)
    scale = kwargs.get('scale', True)
    ddof = kwargs.get('ddof', 0)

    if model == "linear":
        if prior_type == "g":
            return LinearGPrior(X, y, g, Z=Z, scale=scale, ddof=ddof)
        elif prior_type == "ind":
            return LinearIndPrior(X, y, g, Z=Z, scale=scale, ddof=ddof)
    raise ValueError(f"Unknown model={model!r} or prior_type={prior_type!r}")


def _build_proposal(sampler, p, n_temp, n_chain, N_total,
                    likelihood, model_prior, h_exp, N_adapt_PIPs, N_rb, **kwargs):
    """Construct the appropriate proposal object.

    Parameters
    ----------
    sampler : str
        Proposal type ("PARNI", "ASI", "ADS").
    p : int
        Number of variables.
    n_temp, n_chain, N_total : int
        Dimensions.
    likelihood : BaseLikelihood
    model_prior : ModelPrior
    h_exp : float
        Expected inclusion probability.
    N_adapt_PIPs, N_rb : int
        Adaptation horizons.
    **kwargs
        Proposal-specific parameters.

    Returns
    -------
    BaseProposal
    """
    if sampler == "PARNI":
        return PARNIProposal(
            p, n_temp, n_chain, N_total, likelihood, model_prior, h_exp,
            N_adapt_PIPs, N_rb,
            omega_init=kwargs.get('PARNI_omega_init', 0.5),
            bal_fun=kwargs.get('PARNI_bal_fun', 'hastings'),
            omega_adapt=kwargs.get('PARNI_omega_adapt', 'RM'),
            adapt_phi=kwargs.get('PARNI_adapt_phi', -0.7),
            target_acc=kwargs.get('PARNI_target_acc', 0.65),
            kappa=kwargs.get('kappa', 0.001),
        )
    elif sampler == "ASI":
        return ASIProposal(
            p, n_temp, n_chain, N_total, likelihood, model_prior, h_exp,
            N_adapt_PIPs, N_rb,
            adapt_phi=kwargs.get('ASI_adapt_phi', -0.7),
            target_acc=kwargs.get('ASI_target_acc', 0.234),
            zeta_init=kwargs.get('ASI_zeta_init', 0.5),
            kappa=kwargs.get('kappa', 0.001),
        )
    elif sampler == "ADS":
        return ADSProposal(p, likelihood, model_prior)
    raise ValueError(f"Unknown sampler={sampler!r}")


class BVS_MCMC:
    """Bayesian Variable Selection via MCMC.

    Assembles a likelihood model, model prior, MCMC proposal, and optional
    parallel tempering and random-g sampling into a complete sampler.

    Parameters
    ----------
    X : np.ndarray, shape (n, p)
        Design matrix.
    y : np.ndarray, shape (n,) or (n, 1)
        Response vector.
    g : float or str
        Prior variance parameter, or "random" for half-Cauchy random g.
    model : str
        Model type. Currently supports "linear".
    prior_type : str
        Prior structure for coefficients: "g" (Zellner) or "ind" (independent).
    h_exp_size : float
        Expected model size for the model prior.
    h_type : int
        Model prior type (1 = fixed, 2 = beta-binomial).
    h : float, list, or None
        Explicit model prior parameter. Overrides h_type/h_exp_size if given.
    Z : np.ndarray or None
        Nuisance covariates to project out.
    scale : bool
        Whether to standardise X columns.
    ddof : int
        Degrees of freedom correction for standardisation.

    Attributes (available after sample_now)
    ----------------------------------------
    estm_PIPs : np.ndarray, shape (n_temp, p)
        Estimated posterior inclusion probabilities.
    acc_times : np.ndarray, shape (n_temp,)
        Acceptance rates per temperature.
    ASJD : np.ndarray, shape (n_temp,)
        Average squared jumping distance per temperature.
    log_posts : np.ndarray, shape (n_temp, n_chain, N_total + 1)
        Log-posterior trace.
    model_sizes : np.ndarray, shape (n_temp, n_chain, N_total + 1)
        Model size trace.
    time_total, time_burnin, time_sample : float
        CPU time in seconds.
    """

    def __init__(self, X, y, g,
                 model="linear", prior_type="ind",
                 h_exp_size=5, h_type=2, h=None,
                 Z=None, scale=True, ddof=0):
        self.likelihood = _build_likelihood(
            model, prior_type, X, y, g, Z=Z, scale=scale, ddof=ddof
        )
        self.model_prior = ModelPrior(
            p=X.shape[1], h=h, h_type=h_type, h_exp_size=h_exp_size
        )
        self.p = X.shape[1]
        self.n = X.shape[0]

    def set_alg_par(self, sampler="PARNI", N_iter=500, N_burnin=500, n_chain=1,
                    N_adapt_PIPs=None, N_rb=None, use_rb=True, kappa=0.001,
                    # PARNI parameters
                    PARNI_omega_init=0.5, PARNI_bal_fun="hastings",
                    PARNI_omega_adapt="RM", PARNI_adapt_phi=-0.7, PARNI_target_acc=0.65,
                    # ASI parameters
                    ASI_adapt_phi=-0.7, ASI_target_acc=0.234, ASI_zeta_init=0.5,
                    # g parameters
                    g_adapt_phi=-0.7, g_target_acc=0.234, g_init=None,
                    # Parallel tempering
                    n_temp=1, PT_taus=None, PT_adapt_phi=-0.7, PT_target_acc=0.234,
                    verbose=False, store_chains=False, gamma_init=None, f=None):
        """Configure algorithm parameters and initialise all components.

        Parameters
        ----------
        sampler : str
            Proposal type: "PARNI", "ASI", or "ADS".
        N_iter : int
            Number of post-burn-in iterations.
        N_burnin : int
            Number of burn-in iterations.
        n_chain : int
            Number of parallel chains per temperature.
        N_adapt_PIPs : int or None
            PIP adaptation horizon. Defaults to N_burnin.
        N_rb : int or None
            Rao-Blackwellisation horizon. Defaults to N_burnin.
        use_rb : bool
            Whether to use Rao-Blackwellised PIP estimates.
        kappa : float
            PIP floor/ceiling.
        PARNI_omega_init : float
            Initial omega for PARNI.
        PARNI_bal_fun : str or callable
            Balancing function for PARNI.
        PARNI_omega_adapt : str
            "RM" or "KW" for PARNI omega adaptation.
        PARNI_adapt_phi : float or list
            Adaptation step-size exponent(s).
        PARNI_target_acc : float
            Target acceptance rate for PARNI (RM only).
        ASI_adapt_phi : float
            Step-size exponent for ASI.
        ASI_target_acc : float
            Target acceptance rate for ASI.
        ASI_zeta_init : float
            Initial zeta for ASI.
        g_adapt_phi : float
            Step-size exponent for g adaptation.
        g_target_acc : float
            Target acceptance rate for g.
        g_init : float or None
            Initial g value (for random g).
        n_temp : int
            Number of temperature levels.
        PT_taus : np.ndarray or None
            Initial temperature spacings.
        PT_adapt_phi : float
            Step-size exponent for temperature adaptation.
        PT_target_acc : float
            Target swap acceptance rate.
        verbose : bool
            Show progress bar.
        store_chains : bool
            Store full chain history.
        gamma_init : np.ndarray or None
            Fixed initial model. If None, random initialisation.
        f : callable or None
            Function to evaluate on each model (e.g., for posterior expectations).
        """
        self.N_iter = N_iter
        self.N_burnin = N_burnin
        self.N_total = N_iter + N_burnin
        self.n_chain = n_chain
        self.n_temp = n_temp
        self.verbose = verbose
        self.store_chains = store_chains
        self._gamma_init = gamma_init

        N_adapt_PIPs = N_adapt_PIPs or N_burnin
        N_rb = N_rb or N_burnin

        # Build proposal
        self.proposal = _build_proposal(
            sampler, self.p, n_temp, n_chain, self.N_total,
            self.likelihood, self.model_prior, self.model_prior.h_exp,
            N_adapt_PIPs, N_rb,
            PARNI_omega_init=PARNI_omega_init, PARNI_bal_fun=PARNI_bal_fun,
            PARNI_omega_adapt=PARNI_omega_adapt, PARNI_adapt_phi=PARNI_adapt_phi,
            PARNI_target_acc=PARNI_target_acc,
            ASI_adapt_phi=ASI_adapt_phi, ASI_target_acc=ASI_target_acc,
            ASI_zeta_init=ASI_zeta_init, kappa=kappa,
        )

        # Build tempering
        if n_temp > 1:
            self.tempering = ParallelTempering(
                n_temp, n_chain, self.N_total,
                PT_taus=PT_taus, adapt_phi=PT_adapt_phi, target_acc=PT_target_acc,
            )
        else:
            self.tempering = NullTempering()

        # Build g sampler
        if self.likelihood.random_g:
            self.g_sampler = RandomGSampler(
                n_temp, n_chain, self.N_total,
                g_init=g_init, adapt_phi=g_adapt_phi, target_acc=g_target_acc,
            )
        else:
            self.g_sampler = None

        # Function evaluation
        if f is not None:
            self._eval_f = True
            self._f = f
            self._f_sum = np.zeros(f(np.zeros(self.p)).shape)
            self.fs = np.zeros((n_chain, self.N_total + 1) + self._f_sum.shape)
        else:
            self._eval_f = False

        # Storage arrays
        if store_chains:
            self.chains = np.zeros((n_chain, self.N_total + 1, self.p))

        self.log_posts = np.zeros((n_temp, n_chain, self.N_total + 1))
        self.model_sizes = np.zeros((n_temp, n_chain, self.N_total + 1))
        self.acc_times = np.zeros(n_temp)
        self.ASJD = np.zeros(n_temp)
        self.estm_PIPs = np.zeros((n_temp, self.p))

        # Initialise starting models
        self.agg_gamma_par = {}
        h_exp = self.model_prior.h_exp

        for t_idx in range(n_temp):
            for i_idx in range(n_chain):
                if self.g_sampler is not None:
                    self.likelihood.g = self.g_sampler.g_matrix[t_idx, i_idx]

                if gamma_init is None:
                    gamma = np.random.random(self.p) < h_exp
                else:
                    gamma = gamma_init.copy()

                gp = GammaPar(gamma)
                update_gamma_par(gp, self.likelihood, self.model_prior)
                self.agg_gamma_par[t_idx, i_idx] = gp

                temp = self.tempering.temps[t_idx]
                self.log_posts[t_idx, i_idx, 0] = gp.log_mp + temp * gp.log_llh
                self.model_sizes[t_idx, i_idx, 0] = gp.p_gam

                if t_idx == 0 and store_chains:
                    self.chains[i_idx, 0] = gp.gamma
                if t_idx == 0 and self._eval_f:
                    self.fs[i_idx, 0] = self._f(gp.gamma)

    def sample_now(self):
        """Run the MCMC sampler.

        Executes N_burnin + N_iter iterations of the configured algorithm.
        Results are stored in instance attributes (estm_PIPs, acc_times,
        ASJD, log_posts, model_sizes, chains, time_total, etc.).
        """
        time_1 = time.time()

        if self.verbose:
            pbar = tqdm(total=self.N_total, position=0, leave=True)

        for epoch in range(1, self.N_total + 1):
            if self.verbose:
                pbar.update(1)

            for i in range(self.n_chain):

                # Parallel tempering swap
                if self.n_temp > 1:
                    self.tempering.swap_chains(self.agg_gamma_par, i, epoch, self.N_burnin)

                for t in range(self.n_temp):
                    temp = self.tempering.temps[t]
                    curr = self.agg_gamma_par[t, i]

                    # Random g update
                    if self.g_sampler is not None:
                        curr = self.g_sampler.update(
                            curr, self.likelihood, t, i, epoch, temp, self.N_burnin
                        )
                        self.agg_gamma_par[t, i] = curr

                    # Propose
                    prop, log_prop_odd, change = self.proposal.propose(curr, temp=temp)

                    # Accept / reject
                    curr_log_post = curr.log_mp + temp * curr.log_llh
                    prop_log_post = prop.log_mp + temp * prop.log_llh
                    acc_rate = min(1, math.exp(prop_log_post - curr_log_post + log_prop_odd))

                    if np.random.random() < acc_rate:
                        curr = prop
                        self.agg_gamma_par[t, i] = prop
                        accepted = True
                    else:
                        accepted = False

                    # Adapt proposal
                    if self.proposal.is_adaptive:
                        self.proposal.update_params(
                            curr, acc_rate, change,
                            t=t, i=i, epoch=epoch, temp=temp,
                        )

                    # Record
                    self.log_posts[t, i, epoch] = curr.log_mp + temp * curr.log_llh
                    self.model_sizes[t, i, epoch] = curr.p_gam

                    if epoch > self.N_burnin:
                        self.acc_times[t] += accepted
                        self.ASJD[t] += change * acc_rate
                        self.estm_PIPs[t] += curr.gamma

                    if t == 0:
                        if self.store_chains:
                            self.chains[i, epoch] = curr.gamma
                        if self._eval_f:
                            f_curr = self._f(curr.gamma)
                            self.fs[i, epoch] = f_curr
                            self._f_sum += f_curr

                    # Temperature adaptation bookkeeping
                    if self.n_temp > 1:
                        self.tempering.update_H(i, t, curr.log_llh)

            # Adapt temperatures
            if self.n_temp > 1:
                self.tempering.adapt_temperatures(epoch)

            if epoch == self.N_burnin:
                time_2 = time.time()

        time_3 = time.time()

        if self.verbose:
            pbar.close()

        # Normalise
        N_const = self.n_chain * self.N_iter
        self.ASJD /= N_const
        self.acc_times /= N_const
        self.estm_PIPs /= N_const

        if self.n_temp > 1:
            self.tempering.acc_times /= N_const

        if self._eval_f:
            self.f_sum = self._f_sum / N_const

        if self.g_sampler is not None:
            self.g_sampler.acc_times /= N_const

        # Timing
        self.time_total = time_3 - time_1
        self.time_burnin = time_2 - time_1
        self.time_sample = time_3 - time_2
