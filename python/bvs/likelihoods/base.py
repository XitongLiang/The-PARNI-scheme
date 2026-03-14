"""Abstract base class for likelihood models."""

from abc import ABC, abstractmethod


class BaseLikelihood(ABC):
    """Interface for likelihood computation in Bayesian variable selection.

    Subclasses must implement log_llh() and compute_bf(). These methods
    operate on a GammaPar object, reading the model indicator (gamma) and
    writing computed statistics (log_llh, inv_V_gam, BF) back to it.

    Attributes
    ----------
    n : int
        Number of observations.
    p : int
        Number of candidate variables.
    g : float
        Prior variance scaling parameter.
    random_g : bool
        Whether g is treated as a random variable.
    """

    @abstractmethod
    def log_llh(self, gamma_par):
        """Compute the log marginal likelihood for the given model.

        Parameters
        ----------
        gamma_par : GammaPar
            Model indicator. Reads gamma, p_gam, includes.

        Side Effects
        ------------
        Sets gamma_par.log_llh (float).
        May set gamma_par.inv_V_gam (np.ndarray) as a by-product.
        """

    @abstractmethod
    def compute_bf(self, gamma_par):
        """Compute Bayes factors for all single-variable changes.

        BF[j] is the Bayes factor for adding variable j (if excluded)
        or removing variable j (if included), relative to the current model.

        Parameters
        ----------
        gamma_par : GammaPar
            Model indicator. Reads gamma, p_gam, includes, inv_V_gam.

        Side Effects
        ------------
        Sets gamma_par.BF (np.ndarray of length p).
        """
