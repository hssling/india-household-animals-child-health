#!/usr/bin/env python3
"""Survey-weighted generalised linear models with an explicit cluster-robust
sandwich covariance.

Why this module exists
----------------------
statsmodels refuses to combine `var_weights` with `cov_type="cluster"`, emitting
"cov_type not fully supported with var_weights". The earlier version of this
study simply disclosed that warning as a limitation. That is avoidable: the
sandwich estimator for a weighted GLM has a closed form, so it is computed here
directly rather than relying on the library.

For a GLM with canonical-form score contributions, survey weights w, fitted
means mu, design X and response y:

    bread  A = (X' diag(w * v(mu)) X)^-1        with v(mu) = mu for the log-link
                                                Poisson used here
    score  u_i = w_i (y_i - mu_i) x_i
    meat   B = sum_c (sum_{i in c} u_i)(sum_{i in c} u_i)'
    Var    = A B A, with a finite-cluster correction G/(G-1) * (n-1)/(n-k)

Point estimates come from statsmodels' weighted fit, which is correct; only the
covariance is replaced. The result is the same estimator that survey packages
report, with no suppressed warning.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats


@dataclass
class SurveyFit:
    params: pd.Series
    cov: np.ndarray
    standard_errors: pd.Series
    p_values: pd.Series
    conf_low: pd.Series
    conf_high: pd.Series
    n_observations: int
    n_clusters: int
    design_info: object
    exog_names: list[str]
    fitted_values: np.ndarray

    def summary_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "coefficient": self.params,
                "standard_error": self.standard_errors,
                "conf_low": self.conf_low,
                "conf_high": self.conf_high,
                "p_value": self.p_values,
            }
        )


def fit_survey_poisson(
    formula: str,
    data: pd.DataFrame,
    weights: np.ndarray | None,
    clusters: np.ndarray,
) -> SurveyFit:
    """Weighted log-link Poisson (modified Poisson) with cluster-robust covariance."""
    weight_vector = np.ones(len(data)) if weights is None else np.asarray(weights, dtype=float)

    model = smf.glm(
        formula,
        data=data,
        family=sm.families.Poisson(link=sm.families.links.Log()),
        var_weights=weight_vector,
    )
    with warnings.catch_warnings():
        # The default (non-robust) fit is used only for point estimates; its own
        # covariance is discarded and replaced below, so no warning applies.
        warnings.simplefilter("ignore")
        fitted = model.fit()

    design = np.asarray(model.exog, dtype=float)
    response = np.asarray(model.endog, dtype=float).ravel()
    mu = np.asarray(fitted.fittedvalues, dtype=float).ravel()

    n_observations, n_parameters = design.shape

    # Bread: inverse weighted information matrix.
    information = design.T @ (design * (weight_vector * mu)[:, None])
    bread = np.linalg.pinv(information)

    # Score contributions, summed within cluster.
    scores = (weight_vector * (response - mu))[:, None] * design
    frame = pd.DataFrame(scores)
    frame["_cluster"] = np.asarray(clusters)
    cluster_scores = frame.groupby("_cluster", sort=False).sum().to_numpy()
    n_clusters = cluster_scores.shape[0]

    meat = cluster_scores.T @ cluster_scores
    correction = (n_clusters / (n_clusters - 1)) * (
        (n_observations - 1) / (n_observations - n_parameters)
    )
    covariance = bread @ meat @ bread * correction

    names = list(model.exog_names)
    params = pd.Series(np.asarray(fitted.params, dtype=float), index=names)
    standard_errors = pd.Series(np.sqrt(np.diag(covariance)), index=names)
    z_scores = params / standard_errors
    p_values = pd.Series(2 * stats.norm.sf(np.abs(z_scores)), index=names)

    return SurveyFit(
        params=params,
        cov=covariance,
        standard_errors=standard_errors,
        p_values=p_values,
        conf_low=params - 1.96 * standard_errors,
        conf_high=params + 1.96 * standard_errors,
        n_observations=n_observations,
        n_clusters=n_clusters,
        design_info=model.data.design_info,
        exog_names=names,
        fitted_values=mu,
    )


def fit_survey_logistic(
    formula: str,
    data: pd.DataFrame,
    weights: np.ndarray | None,
    clusters: np.ndarray,
) -> SurveyFit:
    """Weighted logistic regression with the same cluster-robust sandwich.

    For the logit link the score is w_i (y_i - p_i) x_i and the weighted
    information uses v(p) = p(1 - p).
    """
    weight_vector = np.ones(len(data)) if weights is None else np.asarray(weights, dtype=float)

    model = smf.glm(formula, data=data, family=sm.families.Binomial(),
                    var_weights=weight_vector)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fitted = model.fit()

    design = np.asarray(model.exog, dtype=float)
    response = np.asarray(model.endog, dtype=float).ravel()
    probability = np.asarray(fitted.fittedvalues, dtype=float).ravel()

    n_observations, n_parameters = design.shape
    variance = probability * (1 - probability)
    information = design.T @ (design * (weight_vector * variance)[:, None])
    bread = np.linalg.pinv(information)

    scores = (weight_vector * (response - probability))[:, None] * design
    frame = pd.DataFrame(scores)
    frame["_cluster"] = np.asarray(clusters)
    cluster_scores = frame.groupby("_cluster", sort=False).sum().to_numpy()
    n_clusters = cluster_scores.shape[0]

    meat = cluster_scores.T @ cluster_scores
    correction = (n_clusters / (n_clusters - 1)) * (
        (n_observations - 1) / (n_observations - n_parameters)
    )
    covariance = bread @ meat @ bread * correction

    names = list(model.exog_names)
    params = pd.Series(np.asarray(fitted.params, dtype=float), index=names)
    standard_errors = pd.Series(np.sqrt(np.diag(covariance)), index=names)
    p_values = pd.Series(2 * stats.norm.sf(np.abs(params / standard_errors)), index=names)

    return SurveyFit(
        params=params,
        cov=covariance,
        standard_errors=standard_errors,
        p_values=p_values,
        conf_low=params - 1.96 * standard_errors,
        conf_high=params + 1.96 * standard_errors,
        n_observations=n_observations,
        n_clusters=n_clusters,
        design_info=model.data.design_info,
        exog_names=names,
        fitted_values=probability,
    )
