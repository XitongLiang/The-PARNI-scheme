"""MCMC proposal mechanisms for Bayesian variable selection.

Available proposals:
- ADSProposal: Add-Delete-Swap (non-adaptive baseline).
- ASIProposal: Adaptively Scaled Individual proposal.
- PARNIProposal: Point-wise Adaptive Random Neighbourhood Informed proposal.

To add a new proposal, subclass BaseProposal and implement propose().
"""

from .ads import ADSProposal
from .asi import ASIProposal
from .parni import PARNIProposal

__all__ = ["ADSProposal", "ASIProposal", "PARNIProposal"]
