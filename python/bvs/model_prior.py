"""Model prior for Bayesian variable selection.

Supports two parameterisation schemes:
1. Fixed h: each variable is included independently with probability h.
2. Beta-binomial: the inclusion probability follows a Beta(alpha, beta) prior,
   leading to a beta-binomial distribution on model size.
"""

import math
from scipy import special


class ModelPrior:
    """Model prior over the inclusion vector gamma.

    Parameters
    ----------
    p : int
        Total number of candidate variables.
    h : float, list of two floats, or None
        Prior inclusion parameter.
        - float: fixed inclusion probability.
        - [alpha, beta]: beta-binomial hyperparameters.
        - None: determined by h_type and h_exp_size.
    h_type : int
        Prior type when h is None. 1 = fixed, 2 = beta-binomial.
    h_exp_size : float
        Expected model size used to derive h when h is None.

    Attributes
    ----------
    h : float or list
        Stored prior parameter(s).
    h_exp : float
        Expected inclusion probability per variable.
    """

    def __init__(self, p, h=None, h_type=2, h_exp_size=5):
        self.p = p

        if h is not None:
            self.h = h
            if isinstance(h, float):
                self._type = 1
                self.h_exp = h
            elif isinstance(h, list) and len(h) == 2:
                self._type = 2
                self.h_exp = h[1] / (h[1] + h[2])
        else:
            if h_type == 1:
                self.h_exp = h_exp_size / p
                self.h = self.h_exp
                self._type = 1
            elif h_type == 2:
                self.h_exp = h_exp_size / p
                self.h = [1, (1 - self.h_exp) / self.h_exp]
                self._type = 2

    def h_odd(self, I, gamma_par=None):
        """Compute the prior odds for adding or removing a variable.

        Parameters
        ----------
        I : bool
            True if the variable is currently included (computing deletion odds),
            False if excluded (computing addition odds).
        gamma_par : GammaPar or None
            Current model state. Required for beta-binomial (type 2).

        Returns
        -------
        float
            Prior odds ratio for the inclusion/exclusion change.
        """
        if self._type == 1:
            return self.h / (1 - self.h) if I else (1 - self.h) / self.h
        else:
            return ((gamma_par.p_gam + self.h[0]) / (self.p - gamma_par.p_gam - 1 + self.h[1]) if I
                    else (self.p - gamma_par.p_gam + self.h[1]) / (gamma_par.p_gam - 1 + self.h[0]))

    def log_m_prior(self, gamma_par):
        """Compute the log model prior and store in gamma_par.log_mp.

        Parameters
        ----------
        gamma_par : GammaPar
            Model indicator. Must have p_gam set.

        Side Effects
        ------------
        Sets gamma_par.log_mp (float).
        """
        if self._type == 1:
            if self.h >= 1.0:
                gamma_par.log_mp = 0.0  # uniform prior when h=1
            else:
                gamma_par.log_mp = gamma_par.p_gam * (math.log(self.h) - math.log(1 - self.h))
        else:
            gamma_par.log_mp = special.betaln(gamma_par.p_gam + self.h[0],
                                              self.p - gamma_par.p_gam + self.h[1])

    def h_til(self, gamma_par):
        """Compute the tilted prior parameter for Rao-Blackwellised PIP estimation.

        Parameters
        ----------
        gamma_par : GammaPar
            Model indicator. Must have gamma and p_gam set.

        Side Effects
        ------------
        Sets gamma_par.h_til (float or np.ndarray).
        """
        if self._type == 1:
            gamma_par.h_til = self.h
        else:
            gamma_par.h_til = ((gamma_par.p_gam - gamma_par.gamma + self.h[0])
                               / (self.p + self.h[0] + self.h[1] - 1))
