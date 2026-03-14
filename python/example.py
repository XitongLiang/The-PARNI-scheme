"""Example: PARNI sampler on simulated linear regression data."""

import sys
import os

# Add parent directory to path so bvs package is importable
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from utils import lrrsg
from bvs import BVS_MCMC


# Generate simulated data
beta, y, X = lrrsg(n=10000, p=10, beta0=np.array([2, -3, 2, 2, -3]), rho=0)

# Set up BVS_MCMC
fit = BVS_MCMC(X=X, y=y, Z=None, model="linear", ddof=1,
               g=100, prior_type="ind",
               h_exp_size=2, h_type=1, h=None, scale=True)

# Configure algorithm parameters
fit.set_alg_par(sampler="PARNI", N_iter=2000, N_burnin=1000,
                PARNI_bal_fun="hastings",
                PARNI_omega_adapt="KW", PARNI_adapt_phi=[-1, -0.5],
                ASI_zeta_init=0.5,
                n_chain=25, n_temp=5, verbose=True)

# Run sampler
fit.sample_now()

print(f"Total time:   {fit.time_total:.2f}s")
print(f"Burn-in time: {fit.time_burnin:.2f}s")
print(f"Sample time:  {fit.time_sample:.2f}s")
print(f"Estimated PIPs: {fit.estm_PIPs[0]}")
