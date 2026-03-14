"""Linear regression likelihood models.

Provides two prior structures for the regression coefficients:
- g-prior (Zellner): shared shrinkage via scalar g.
- Independent normal prior: diagonal prior precision 1/g per coefficient.
"""

import math
import numpy as np
from .base import BaseLikelihood


class _LinearBase(BaseLikelihood):
    """Shared data preprocessing for linear regression likelihoods.

    Parameters
    ----------
    X : np.ndarray, shape (n, p)
        Design matrix.
    y : np.ndarray, shape (n,) or (n, 1)
        Response vector.
    g : float or str
        Prior variance parameter. Pass "random" for half-Cauchy random g.
    Z : np.ndarray or None, shape (n, q)
        Nuisance covariates to project out (e.g., intercept columns).
    scale : bool
        Whether to standardise columns of X to unit variance.
    ddof : int
        Delta degrees of freedom for std computation (0 or 1).
    """

    def __init__(self, X, y, g, Z=None, scale=True, ddof=0):
        self.n, self.p = X.shape

        if isinstance(g, (int, float)):
            self.g = g
            self.random_g = False
        elif g == "random":
            self.g = 1.0  # initial value, updated by RandomGSampler
            self.random_g = True

        # Centre and optionally project out nuisance covariates
        if Z is None:
            X = X - X.mean(axis=0)
            y = y - y.mean()
            if scale:
                X = X / X.std(axis=0, ddof=ddof)
        else:
            Z = np.hstack((np.ones((self.n, 1)), Z))
            P_Z = Z.dot(np.linalg.inv(Z.T.dot(Z))).dot(Z.T)
            proj = np.eye(self.n) - P_Z
            y = proj.dot(y)
            X = proj.dot(X)
            if scale:
                X = X / X.std(axis=0, ddof=ddof)

        self.X = X
        self.y = y
        self.yty = np.sum(y ** 2)
        self.ytX = y.T.dot(X)
        self.diag_XtX = np.sum(X ** 2, axis=0)


class LinearGPrior(_LinearBase):
    """Linear regression with Zellner's g-prior.

    The coefficient prior is beta | gamma ~ N(0, g * (X_gamma' X_gamma)^{-1} sigma^2).
    The marginal likelihood has a closed-form expression.
    """

    def log_llh(self, gamma_par):
        """Compute log marginal likelihood under the g-prior.

        Parameters
        ----------
        gamma_par : GammaPar
            Model indicator with gamma, p_gam, includes set.

        Side Effects
        ------------
        Sets gamma_par.log_llh (float).
        Sets gamma_par.inv_V_gam (np.ndarray) when p_gam > 0.
        """
        if gamma_par.p_gam == 0:
            gamma_par.log_llh = -self.n * math.log(self.yty) / 2
        else:
            X_gam = self.X[:, gamma_par.includes]
            V_gam = X_gam.T.dot(X_gam)
            gamma_par.inv_V_gam = np.linalg.inv(V_gam)

            ytX_gam = self.ytX[0, gamma_par.includes]
            ytXgFXgty = ytX_gam.dot(gamma_par.inv_V_gam).dot(ytX_gam.T)

            A = self.yty - self.g / (1 + self.g) * ytXgFXgty
            gamma_par.log_llh = (-gamma_par.p_gam / 2 * math.log(1 + self.g)
                                 - self.n * math.log(A) / 2)

    def compute_bf(self, gamma_par):
        """Compute Bayes factors for all single-variable changes under g-prior.

        BF[j] measures the evidence for toggling variable j relative to the
        current model gamma.

        Parameters
        ----------
        gamma_par : GammaPar
            Model indicator with gamma, p_gam, includes, inv_V_gam set.

        Side Effects
        ------------
        Sets gamma_par.BF (np.ndarray of length p).
        """
        g_ratio = self.g / (self.g + 1)
        inv_sqrtg1 = 1 / math.sqrt(self.g + 1)
        n_power = self.n / 2

        if gamma_par.p_gam == 0:
            BF = (self.yty / (self.yty - g_ratio * self.ytX[0] ** 2 / self.diag_XtX)) ** n_power * inv_sqrtg1
        else:
            X_gam = self.X[:, gamma_par.includes]
            XgtX = X_gam.T.dot(self.X)

            ytXg = self.ytX[0, gamma_par.includes]
            ytXgFXgty = ytXg.dot(gamma_par.inv_V_gam).dot(ytXg.T)

            A = self.yty - ytXgFXgty * g_ratio
            d_vec = 1 / (self.diag_XtX - np.einsum('ij,ji->i', XgtX.T.dot(gamma_par.inv_V_gam), XgtX))

            ytXFXtxj_vec = ytXg.dot(gamma_par.inv_V_gam).dot(XgtX)
            tilda_A_vec = A - d_vec * (ytXFXtxj_vec - self.ytX) ** 2 * g_ratio

            BF = (A / tilda_A_vec.reshape(self.p)) ** n_power * inv_sqrtg1

            if gamma_par.p_gam == 1:
                BF[gamma_par.includes] = (self.yty / A) ** n_power * inv_sqrtg1
            else:
                zj = -1
                for j in gamma_par.includes:
                    zj += 1
                    A_ratio = 1 + (ytXg.dot(gamma_par.inv_V_gam[:, zj])) ** 2 * g_ratio / (A * gamma_par.inv_V_gam[zj, zj])
                    BF[j] = A_ratio ** n_power * inv_sqrtg1

        gamma_par.BF = BF


class LinearIndPrior(_LinearBase):
    """Linear regression with independent normal prior.

    The coefficient prior is beta_j | gamma ~ N(0, g * sigma^2) independently.
    """

    def log_llh(self, gamma_par):
        """Compute log marginal likelihood under the independent prior.

        Parameters
        ----------
        gamma_par : GammaPar
            Model indicator with gamma, p_gam, includes set.

        Side Effects
        ------------
        Sets gamma_par.log_llh (float).
        Sets gamma_par.inv_V_gam (np.ndarray) when p_gam > 0.
        """
        if gamma_par.p_gam == 0:
            gamma_par.log_llh = -self.n * math.log(self.yty) / 2
        else:
            X_gam = self.X[:, gamma_par.includes]
            V_gam = X_gam.T.dot(X_gam)
            V_gam[np.diag_indices(gamma_par.p_gam)] += 1 / self.g
            gamma_par.inv_V_gam = np.linalg.inv(V_gam)

            ytX_gam = self.ytX[0, gamma_par.includes]
            ytXgFXgty = ytX_gam.dot(gamma_par.inv_V_gam).dot(ytX_gam.T)

            A = self.yty - ytXgFXgty
            sqrt_det_Vg = math.log(np.linalg.det(V_gam)) / 2

            gamma_par.log_llh = (-gamma_par.p_gam / 2 * math.log(self.g)
                                 - sqrt_det_Vg
                                 - self.n * math.log(A) / 2)

    def compute_bf(self, gamma_par):
        """Compute Bayes factors for all single-variable changes under independent prior.

        Parameters
        ----------
        gamma_par : GammaPar
            Model indicator with gamma, p_gam, includes, inv_V_gam set.

        Side Effects
        ------------
        Sets gamma_par.BF (np.ndarray of length p).
        """
        inv_g = 1 / self.g
        inv_sqrt_g = math.sqrt(inv_g)
        n_power = self.n / 2

        diag_V = self.diag_XtX + 1 / self.g

        if gamma_par.p_gam == 0:
            BF = ((self.yty / (self.yty - self.ytX[0] ** 2 / diag_V)) ** n_power
                  * np.sqrt(1 / diag_V) * inv_sqrt_g)
        else:
            X_gam = self.X[:, gamma_par.includes]
            XgtX = X_gam.T.dot(self.X)

            ytXg = self.ytX[0, gamma_par.includes]
            ytXgFXgty = ytXg.dot(gamma_par.inv_V_gam).dot(ytXg.T)

            A = self.yty - ytXgFXgty
            d_vec = 1 / (diag_V - np.einsum('ij,ji->i', XgtX.T.dot(gamma_par.inv_V_gam), XgtX))
            d_vec[gamma_par.includes] = 0

            ytXFXtxj_vec = ytXg.dot(gamma_par.inv_V_gam).dot(XgtX)
            tilda_A_vec = A - d_vec * (ytXFXtxj_vec - self.ytX) ** 2

            BF = np.sqrt(d_vec) * (A / tilda_A_vec.reshape(self.p)) ** n_power * inv_sqrt_g

            if gamma_par.p_gam == 1:
                BF[gamma_par.includes] = (math.sqrt(gamma_par.inv_V_gam)
                                          * (self.yty / A) ** n_power * inv_sqrt_g)
            else:
                zj = -1
                for j in gamma_par.includes:
                    zj += 1
                    A_ratio = (A + (ytXg.dot(gamma_par.inv_V_gam[:, zj])) ** 2 / gamma_par.inv_V_gam[zj, zj]) / A
                    BF[j] = (math.sqrt(gamma_par.inv_V_gam[zj, zj])
                             * A_ratio ** n_power * inv_sqrt_g)

        gamma_par.BF = BF
