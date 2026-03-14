"""Add-Delete-Swap (ADS) proposal for Bayesian variable selection.

A non-adaptive baseline proposal that randomly selects one of three
move types with equal probability:
- Add: include a randomly chosen excluded variable.
- Delete: exclude a randomly chosen included variable.
- Swap: simultaneously add one and delete one variable.

At boundary models (empty or full), only add or delete is available.
"""

import math
import numpy as np

from .base import BaseProposal
from ..gamma_par import GammaPar, update_gamma_par


class ADSProposal(BaseProposal):
    """Add-Delete-Swap proposal (non-adaptive).

    Parameters
    ----------
    p : int
        Number of candidate variables.
    likelihood : BaseLikelihood
        Likelihood object for evaluating proposed models.
    model_prior : ModelPrior
        Model prior for computing log model prior of proposals.
    """

    def __init__(self, p, likelihood, model_prior):
        self._p = p
        self._likelihood = likelihood
        self._model_prior = model_prior

    def propose(self, gamma_par, temp=1.0):
        """Propose a new model via add, delete, or swap.

        Parameters
        ----------
        gamma_par : GammaPar
            Current model state.
        temp : float
            Temperature (unused by ADS but required by interface).

        Returns
        -------
        prop_gamma_par : GammaPar
            Proposed model with log_llh and log_mp computed.
        log_prop_odd : float
            Log ratio of reverse-to-forward proposal probabilities.
        change : int
            Number of variables flipped (1 for add/delete, 2 for swap).
        """
        p = self._p

        if 0 < gamma_par.p_gam < p:

            if np.random.random() < 1 / 3:
                # swap
                del_index = gamma_par.includes[np.random.randint(gamma_par.p_gam)]
                add_index = gamma_par.excludes[np.random.randint(p - gamma_par.p_gam)]

                gamma_prop = gamma_par.gamma.copy()
                gamma_prop[del_index] -= 1
                gamma_prop[add_index] += 1

                log_prob_prop = 0
                log_prob_rev = 0
                change = 2

            elif np.random.random() < 1 / 2:
                # add
                add_index = gamma_par.excludes[np.random.randint(p - gamma_par.p_gam)]

                gamma_prop = gamma_par.gamma.copy()
                gamma_prop[add_index] += 1

                log_prob_prop = -math.log(p - gamma_par.p_gam)
                log_prob_rev = -math.log(gamma_par.p_gam + 1)
                change = 1

            else:
                # delete
                del_index = gamma_par.includes[np.random.randint(gamma_par.p_gam)]

                gamma_prop = gamma_par.gamma.copy()
                gamma_prop[del_index] -= 1

                log_prob_prop = -math.log(gamma_par.p_gam)
                log_prob_rev = -math.log(p - gamma_par.p_gam + 1)
                change = 1

        elif gamma_par.p_gam == 0:
            # add (only option)
            add_index = gamma_par.excludes[np.random.randint(p)]

            gamma_prop = gamma_par.gamma.copy()
            gamma_prop[add_index] += 1

            log_prob_prop = -math.log(p)
            log_prob_rev = 0
            change = 1

        elif gamma_par.p_gam == p:
            # delete (only option)
            del_index = gamma_par.includes[np.random.randint(p)]

            gamma_prop = gamma_par.gamma.copy()
            gamma_prop[del_index] -= 1

            log_prob_prop = -math.log(p)
            log_prob_rev = 0
            change = 1

        prop_gamma_par = GammaPar(gamma_prop)
        update_gamma_par(prop_gamma_par, self._likelihood, self._model_prior)

        return prop_gamma_par, log_prob_rev - log_prob_prop, change
