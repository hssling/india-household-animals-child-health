"""Sensitivity of the State pooling to the meta-analysis estimator (post-hoc, added 28 November 2026).

DerSimonian-Laird can understate between-State variance and the width of the pooled interval. This re-pools the
State-specific estimates from 04_review_analyses.py with restricted maximum likelihood (REML) for tau-squared and
the Hartung-Knapp-Sidik-Jonkman (HKSJ) interval (t distribution, k - 1 degrees of freedom).

Usage:  python -I src/05_meta_sensitivity.py
"""
import os

import numpy as np
import pandas as pd
from scipy import optimize, stats

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "outputs")


def reml_tau2(b, v):
    def negll(t2):
        w = 1 / (v + t2)
        mu = np.sum(w * b) / np.sum(w)
        return 0.5 * (np.sum(np.log(v + t2)) + np.log(np.sum(w)) + np.sum(w * (b - mu) ** 2))
    res = optimize.minimize_scalar(negll, bounds=(0, max(10 * np.var(b), 1e-6)), method="bounded")
    return max(0.0, res.x)


def main():
    st = pd.read_csv(os.path.join(OUT, "review_state_estimates_5.csv"))
    rows = []
    for e, g in st.groupby("exposure"):
        b, v = g["pp"].to_numpy(), g["se"].to_numpy() ** 2
        k = len(b)
        t2 = reml_tau2(b, v)
        w = 1 / (v + t2)
        mu = np.sum(w * b) / np.sum(w)
        q_hk = np.sum(w * (b - mu) ** 2) / (k - 1)
        se_hk = np.sqrt(q_hk / np.sum(w))
        tcrit = stats.t.ppf(0.975, k - 1)
        rows.append({"kind": "post-hoc: REML with Hartung-Knapp-Sidik-Jonkman interval", "outcome": "stunted", "exposure": e,
                     "states": k, "pooled_reml_pp": mu, "lo": mu - tcrit * se_hk, "hi": mu + tcrit * se_hk,
                     "tau2_reml": t2, "p": 2 * stats.t.sf(abs(mu / se_hk), k - 1)})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "review_state_pooled_reml_hksj_5.csv"), index=False)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
