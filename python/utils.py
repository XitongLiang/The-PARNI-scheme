import numpy as np


def lrrsg(n, p, a=0, beta0=np.array([2, -3, 2, 2, -3, 3, -2, 3, -2, 3]),
          rho=0, SNR=1, sigma2=1, seed=None):
    """Linear regression random sample generator."""

    if seed is not None:
        np.random.seed(seed)

    beta = np.zeros((p, 1))
    b_n = beta0.size

    beta0 = beta0.reshape(b_n, 1)
    beta[range(b_n), ] = beta0

    beta = SNR * np.sqrt(sigma2 * np.log(p) / n) * beta

    b = np.sqrt(1 - rho**2)

    X = np.random.normal(size=(n, p))

    for j in range(1, p):
        X[:, j] = rho * X[:, j - 1] + b * X[:, j]

    # responses
    y = a + X.dot(beta) + np.random.normal(size=(n, 1))

    return beta, y, X
