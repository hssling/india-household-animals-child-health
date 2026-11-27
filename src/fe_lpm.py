"""Weighted linear model with absorbed group fixed effects and a cluster-robust covariance.

The fixed effects (here, one per primary sampling unit) are removed by weighted within-group demeaning
of the outcome and every regressor; weighted least squares on the demeaned data gives the same slopes as
including one dummy per group (Frisch–Waugh–Lovell). The covariance is the cluster sandwich

    V = (X'WX)^-1 [ sum_c (X_c' W_c e_c)(X_c' W_c e_c)' ] (X'WX)^-1 * G/(G-1) * (n-1)/(n-k)

where k counts only the non-absorbed regressors. That is the usual convention when the fixed effects are
nested within the clusters, as they are when both are the primary sampling unit.

tests/test_fe_lpm.py checks the slopes and standard errors against statsmodels WLS with explicit group
dummies, so this module is verified rather than assumed.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats


@dataclass
class FEFit:
    params: pd.Series
    se: pd.Series
    p: pd.Series
    lo: pd.Series
    hi: pd.Series
    n: int
    n_groups: int
    n_clusters: int

    def frame(self) -> pd.DataFrame:
        return pd.DataFrame({"coef": self.params, "se": self.se, "lo": self.lo, "hi": self.hi, "p": self.p})


def _demean(a: np.ndarray, groups: np.ndarray, w: np.ndarray) -> np.ndarray:
    df = pd.DataFrame(a)
    wsum = pd.Series(w).groupby(groups).transform("sum").to_numpy()
    means = (df.mul(w, axis=0)).groupby(groups).transform("sum").to_numpy() / wsum[:, None]
    return a - means


def fit_fe_lpm(y: np.ndarray, X: pd.DataFrame, groups: np.ndarray, clusters: np.ndarray,
               weights: np.ndarray | None = None, fe_nested_in_clusters: bool = True) -> FEFit:
    y = np.asarray(y, dtype=float)
    Xa = np.asarray(X, dtype=float)
    w = np.ones(len(y)) if weights is None else np.asarray(weights, dtype=float)
    groups = np.asarray(groups)
    clusters = np.asarray(clusters)

    yd = _demean(y[:, None], groups, w).ravel()
    Xd = _demean(Xa, groups, w)
    # Regressors that are constant within every group carry no information once demeaned.
    keep = np.abs(Xd).sum(axis=0) > 1e-10
    names = [c for c, k in zip(X.columns, keep) if k]
    Xd = Xd[:, keep]

    XtWX = Xd.T @ (Xd * w[:, None])
    bread = np.linalg.pinv(XtWX)
    beta = bread @ (Xd.T @ (w * yd))
    resid = yd - Xd @ beta

    scores = Xd * (w * resid)[:, None]
    cs = pd.DataFrame(scores).groupby(clusters, sort=False).sum().to_numpy()
    G = cs.shape[0]
    n, k = Xd.shape
    n_groups = len(pd.unique(groups))
    k_eff = k if fe_nested_in_clusters else k + n_groups
    corr = (G / (G - 1)) * ((n - 1) / (n - k_eff))
    V = bread @ (cs.T @ cs) @ bread * corr

    params = pd.Series(beta, index=names)
    se = pd.Series(np.sqrt(np.diag(V)), index=names)
    tcrit = stats.t.ppf(0.975, G - 1)
    p = pd.Series(2 * stats.t.sf(np.abs(params / se), G - 1), index=names)
    return FEFit(params, se, p, params - tcrit * se, params + tcrit * se, n, n_groups, G)
