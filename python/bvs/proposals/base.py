"""Abstract base class for MCMC proposals."""

from abc import ABC, abstractmethod


class BaseProposal(ABC):
    """Interface for variable selection proposals.

    A proposal generates a candidate model from the current model and
    computes the log proposal odds ratio needed for the MH acceptance step.

    Subclasses must implement propose(). Adaptive proposals should also
    override update_params() and is_adaptive.
    """

    @abstractmethod
    def propose(self, gamma_par, temp=1.0):
        """Generate a proposal from the current model.

        Parameters
        ----------
        gamma_par : GammaPar
            Current model indicator with log_llh and log_mp set.
        temp : float
            Current temperature (1.0 for the target distribution).

        Returns
        -------
        prop_gamma_par : GammaPar
            Proposed model indicator with log_llh and log_mp set.
        log_prop_odd : float
            Log ratio of reverse-to-forward proposal probabilities,
            i.e. log q(current | proposed) - log q(proposed | current).
        change : int
            Number of variables flipped (Hamming distance).
        """

    def update_params(self, gamma_par, acc_rate, change,
                      t=0, i=0, epoch=0, temp=1.0):
        """Adapt proposal parameters after an MH step.

        Default implementation is a no-op. Override in adaptive proposals.

        Parameters
        ----------
        gamma_par : GammaPar
            Accepted model state (current after accept/reject).
        acc_rate : float
            MH acceptance probability for this step.
        change : int
            Number of variables flipped in the proposal.
        t : int
            Temperature index.
        i : int
            Chain index.
        epoch : int
            Current iteration number.
        temp : float
            Temperature value.
        """
        pass

    @property
    def is_adaptive(self):
        """Whether this proposal adapts its parameters.

        Returns
        -------
        bool
        """
        return False
