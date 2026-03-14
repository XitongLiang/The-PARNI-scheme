"""Epsilon-bounded logit and inverse-logit transformations.

These utilities map values in (eps, 1-eps) to the real line and back,
preventing numerical issues at the boundaries. Used for adapting
proposal parameters (omega, zeta) that must remain in (0, 1).
"""

import math
import numpy as np


def logit_e(y, eps):
    """Epsilon-bounded logit transform.

    Maps y in (eps, 1-eps) to the real line.

    Parameters
    ----------
    y : float
        Value to transform, must satisfy eps < y < 1 - eps.
    eps : float
        Boundary buffer (typically 0.1/p).

    Returns
    -------
    float
        logit_e(y) = log(y - eps) - log(1 - y - eps).
    """
    return math.log(y - eps) - math.log(1 - y - eps)


def inv_logit_e(x, eps):
    """Epsilon-bounded inverse logit (scalar).

    Maps a real number x back to (eps, 1-eps).

    Parameters
    ----------
    x : float
        Value on the real line.
    eps : float
        Boundary buffer.

    Returns
    -------
    float
        Value in (eps, 1-eps).
    """
    if x < 0:
        y = (eps + (1 - eps) * math.exp(x)) / (1 + math.exp(x))
    else:
        y = (eps * math.exp(-x) + 1 - eps) / (1 + math.exp(-x))
    return y


def inv_logit_e_vec(x, eps):
    """Epsilon-bounded inverse logit (vectorised).

    Parameters
    ----------
    x : np.ndarray
        Array of values on the real line.
    eps : float
        Boundary buffer.

    Returns
    -------
    np.ndarray
        Array of values in (eps, 1-eps), same shape as x.
    """
    y = np.empty_like(x)
    mask_neg = x <= 0
    mask_pos = ~mask_neg
    y[mask_neg] = (eps + (1 - eps) * np.exp(x[mask_neg])) / (1 + np.exp(x[mask_neg]))
    y[mask_pos] = (eps * np.exp(-x[mask_pos]) + 1 - eps) / (1 + np.exp(-x[mask_pos]))
    return y
