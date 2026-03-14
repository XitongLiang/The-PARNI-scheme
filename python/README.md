# Python Implementation

Bayesian Variable Selection via adaptive MCMC (PARNI, ASI, ADS samplers).

## Requirements

- Python 3.8+
- NumPy
- Matplotlib (optional, for plotting)

Install dependencies:

```bash
pip install numpy matplotlib
```

## Quick Start

```bash
cd python
python example.py
```

Or in your own script:

```python
import numpy as np
from utils import lrrsg
from bvs import BVS_MCMC

# Generate simulated data
beta, y, X = lrrsg(n=1000, p=20, beta0=np.array([2, -3, 2, 2, -3]))

# Set up and run
fit = BVS_MCMC(X, y, g=100, prior_type="ind")
fit.set_alg_par(sampler="PARNI", N_iter=2000, N_burnin=1000,
                n_chain=5, n_temp=3, verbose=True)
fit.sample_now()

print(f"PIPs: {fit.estm_PIPs[0]}")
print(f"Time: {fit.time_total:.1f}s")
```

## File Structure

| File / Folder | Description |
|----------------|-------------|
| `example.py` | End-to-end example with simulated data |
| `utils.py` | Data generation utility (`lrrsg`) |
| `bvs/` | Main package — see [`bvs/README.md`](bvs/README.md) for architecture details |

## Key Parameters

### `BVS_MCMC(X, y, ...)`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `g` | `100` | Prior variance scaling (integer or `"random"` for half-Cauchy) |
| `prior_type` | `"ind"` | Likelihood: `"g"` (Zellner's g-prior) or `"ind"` (independent normal) |
| `h_type` | `1` | Model prior: `1` = fixed-h, `2` = beta-binomial |
| `h_exp_size` | `2` | Expected model size for computing inclusion probability |
| `scale` | `True` | Standardise covariates before fitting |

### `fit.set_alg_par(...)`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sampler` | `"PARNI"` | Proposal: `"PARNI"`, `"ASI"`, or `"ADS"` |
| `N_iter` | `2000` | Total MCMC iterations |
| `N_burnin` | `1000` | Burn-in iterations |
| `n_chain` | `25` | Number of parallel chains |
| `n_temp` | `5` | Number of temperature levels (1 = no tempering) |
| `PARNI_omega_adapt` | `"KW"` | Adaptation: `"KW"` (Kiefer-Wolfowitz) or `"RM"` (Robbins-Monro) |
| `PARNI_bal_fun` | `"hastings"` | Balancing function: `"hastings"` or `"barker"` |

## Outputs

After `fit.sample_now()`:

- `fit.estm_PIPs[0]` — estimated posterior inclusion probabilities (array of length p)
- `fit.time_total`, `fit.time_burnin`, `fit.time_sample` — timing in seconds
- `fit.gamma_history` — full model indicator trace

## Plotting

```python
from bvs import plot_temperatures, plot_pips

plot_temperatures(fit)  # temperature ladder over iterations
plot_pips(fit)           # PIP bar chart
```
