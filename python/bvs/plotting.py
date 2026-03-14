"""Visualisation utilities for BVS_MCMC results.

Free functions that take a fitted BVS_MCMC object and produce
diagnostic plots.
"""

import matplotlib.pyplot as plt


def plot_temperatures(sampler):
    """Plot the temperature schedule evolution across iterations.

    Parameters
    ----------
    sampler : BVS_MCMC
        Fitted sampler. Must have tempering.temperatures available
        (requires n_temp > 1).
    """
    if sampler.tempering.temperatures is None:
        print("No temperature trace available (n_temp=1).")
        return

    for t in range(sampler.n_temp):
        plt.plot(sampler.tempering.temperatures[t])
    plt.ylim(0, 1.1)
    plt.xlabel("Iteration")
    plt.ylabel("Temperature")
    plt.title("Temperature schedule")
    plt.show()


def plot_pips(sampler):
    """Bar plot of estimated posterior inclusion probabilities.

    Produces one subplot per temperature level.

    Parameters
    ----------
    sampler : BVS_MCMC
        Fitted sampler with estm_PIPs attribute.
    """
    n_temp = sampler.n_temp
    p = sampler.p

    if n_temp == 1:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.bar(range(p), sampler.estm_PIPs[0])
        ax.set_ylim(0, 1)
        ax.set_xlabel("Variable index")
        ax.set_ylabel("PIP")
        ax.set_title("Posterior Inclusion Probabilities")
    else:
        fig, axs = plt.subplots(n_temp, figsize=(7, 4 * n_temp))
        for t in range(n_temp):
            axs[t].bar(range(p), sampler.estm_PIPs[t])
            axs[t].set_ylim(0, 1)
            axs[t].set_xlabel("Variable index")
            axs[t].set_ylabel("PIP")
            axs[t].set_title(f"Temperature level {t}")
    fig.tight_layout()
    plt.show()
