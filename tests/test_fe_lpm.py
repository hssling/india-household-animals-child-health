"""fe_lpm must reproduce statsmodels WLS with explicit group dummies and cluster-robust covariance."""
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from fe_lpm import fit_fe_lpm  # noqa: E402


def test_matches_statsmodels_with_dummies():
    rng = np.random.default_rng(7)
    G, m = 60, 12
    g = np.repeat(np.arange(G), m)
    a = rng.normal(size=G)[g]
    x1 = rng.binomial(1, 0.4, size=G * m) + 0.3 * a
    x2 = rng.normal(size=G * m)
    w = np.repeat(rng.uniform(0.5, 2.0, size=G), m)          # weights constant within cluster, as in NFHS
    y = 0.2 * x1 - 0.1 * x2 + a + rng.normal(size=G * m)
    X = pd.DataFrame({"x1": x1, "x2": x2})

    ours = fit_fe_lpm(y, X, groups=g, clusters=g, weights=w)

    D = pd.get_dummies(pd.Series(g), prefix="g", drop_first=False).astype(float)
    Xs = pd.concat([X, D], axis=1)
    ref = sm.WLS(y, Xs, weights=w).fit(cov_type="cluster", cov_kwds={"groups": g, "use_correction": True})

    assert np.allclose(ours.params.values, ref.params[["x1", "x2"]].values, atol=1e-10)
    # statsmodels counts the dummies in k; with FE nested in clusters the convention drops them, so compare
    # after putting both on the same small-sample factor.
    n, k_sm = Xs.shape
    k_ours = 2
    ratio = np.sqrt(((n - 1) / (n - k_ours)) / ((n - 1) / (n - k_sm)))
    assert np.allclose(ours.se.values, ref.bse[["x1", "x2"]].values * ratio, rtol=1e-8)


if __name__ == "__main__":
    test_matches_statsmodels_with_dummies()
    print("fe_lpm matches statsmodels: PASS")
