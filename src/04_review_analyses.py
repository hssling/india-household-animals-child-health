"""Analyses added at internal review (27 November 2026), after the primary results were seen.

Two kinds, labelled in every output file:
  PRE-SPECIFIED, previously omitted (plan sections 5 and 8.6): antibiotic among children with non-bloody
    diarrhoea and among children with fever or cough; conditional logit with cluster fixed effects.
  POST-HOC robustness (docs/DEVIATIONS.md): extended adjustment (maternal height and BMI, land area, sex of
    household head); continuous height-for-age and haemoglobin; occupation-adjusted subsample; Oster
    coefficient-stability bounds; State-level heterogeneity of the stunting associations.

Usage:  python -I src/04_review_analyses.py
"""
import importlib.util
import os
import sys
import warnings

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(REPO, "outputs")
sys.path.insert(0, HERE)
from fe_lpm import fit_fe_lpm  # noqa: E402

spec = importlib.util.spec_from_file_location("analysis", os.path.join(HERE, "02_analysis.py"))
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)


def load(r):
    d = pd.read_parquet(os.path.join(REPO, "data_local", f"nfhs{r}_children.parquet"))
    d = d[(d["age_m"] >= 0) & (d["age_m"] <= 59)].copy()
    d["district"] = d["district"].fillna(-1)
    c = A.covariates(d)
    c["assetq"] = A.asset_index(d)
    # extended covariates (post-hoc)
    c["mheight"] = pd.cut(d["mheight_cm"], [0, 144.99, 149.99, 154.99, 999], labels=["<145", "145-149", "150-154", "155+"])
    c["mbmi"] = pd.cut(d["mbmi"], [0, 18.49, 24.99, 99], labels=["<18.5", "18.5-24.9", "25+"])
    ha = (d["hv245"] / 10).where(d["hv245"] < 950)
    ha = ha.where(c["land"] == 1, 0.0)
    c["land_ha"] = pd.cut(ha, [-0.01, 0, 0.99, 1.99, 999], labels=["none", "<1", "1-2", "2+"])
    c["female_head"] = (d["v151"] == 2).astype(float).where(d["v151"].isin([1, 2]))
    occ = {0: "not_working", 1: "professional", 2: "professional", 3: "clerical_sales", 4: "clerical_sales",
           5: "services", 6: "agricultural", 7: "manual", 8: "manual", 9: "other"}
    c["occ_husband"] = d["v705"].map(occ)
    c["occ_mother"] = d["v717"].map(occ)
    d["haz_cont"] = d["haz"]
    d["hb_cont"] = d["hb"].where(d["age_m"] >= 6)
    return d, c


EXT = A.ADJ + ["mheight", "mbmi", "land_ha", "female_head"]


def species_of(d):
    return [s for s in A.SPECIES if d[s].notna().any()]


def tag(df, kind):
    df.insert(0, "kind", kind)
    return df


def conditional_logit(d, c, outcomes, sets):
    from statsmodels.discrete.conditional_models import ConditionalLogit
    rows = []
    for y in outcomes:
        for expo in sets:
            frame = pd.concat([d[[y, "psu"]], d[expo], c[A.ADJ]], axis=1).dropna()
            g = frame.groupby("psu")[y].transform(lambda s: s.nunique() > 1)
            frame = frame[g]
            X = A.design(frame.drop(columns=[y, "psu"]), expo + A.ADJ)
            Xd = X - X.groupby(frame["psu"].to_numpy()).transform("mean")
            X = X.loc[:, (Xd.abs().sum() > 1e-9).to_numpy()]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m = ConditionalLogit(frame[y].to_numpy(), X.to_numpy(), groups=frame["psu"].to_numpy())
                res = m.fit(disp=False)
            names = list(X.columns)
            for e in expo:
                if e in names:
                    i = names.index(e)
                    b, se = res.params[i], res.bse[i]
                    rows.append({"outcome": y, "exposure": e, "or": np.exp(b), "lo": np.exp(b - 1.96 * se),
                                 "hi": np.exp(b + 1.96 * se), "p": res.pvalues[i], "n": len(frame),
                                 "clusters_informative": frame["psu"].nunique()})
    return pd.DataFrame(rows)


def oster(d, c, y, species):
    """Oster (2019) bounds within PSU: short model = species + PSU fixed effects; long = primary model."""
    out = []
    frame = pd.concat([d[[y, "psu", "weight"]], d[species], c[A.ADJ]], axis=1).dropna()
    yv, g, w = frame[y].to_numpy(), frame["psu"].to_numpy(), frame["weight"].to_numpy()
    short = fit_fe_lpm(yv, frame[species].astype(float), g, g, w)
    long = fit_fe_lpm(yv, A.design(frame.drop(columns=[y, "psu", "weight"]), species + A.ADJ), g, g, w)
    rmax = min(1.0, 1.3 * long.r2_within)
    for e in species:
        bs, bl = short.params[e], long.params[e]
        rs, rl = short.r2_within, long.r2_within
        delta = bl * (rl - rs) / ((bs - bl) * (rmax - rl)) if (bs - bl) != 0 else np.inf
        bstar = bl - (bs - bl) * (rmax - rl) / (rl - rs)
        out.append({"outcome": y, "exposure": e, "beta_short_pp": 100 * bs, "beta_long_pp": 100 * bl,
                    "r2_within_short": rs, "r2_within_long": rl, "r2_max": rmax,
                    "delta_for_zero": delta, "beta_star_delta1_pp": 100 * bstar, "n": long.n})
    return pd.DataFrame(out)


def state_heterogeneity(d, c, y, exposures):
    rows = []
    for st, idx in d.groupby("state").groups.items():
        sub = d.loc[idx]
        frame = pd.concat([sub[[y, "psu", "weight"]], sub[exposures], c.loc[idx, A.ADJ]], axis=1).dropna()
        if frame["psu"].nunique() < 50:
            continue
        try:
            X = A.design(frame.drop(columns=[y, "psu", "weight"]), exposures + A.ADJ)
            f = fit_fe_lpm(frame[y].to_numpy(), X, frame["psu"].to_numpy(), frame["psu"].to_numpy(), frame["weight"].to_numpy())
        except Exception:
            continue
        for e in exposures:
            if e in f.params.index and np.isfinite(f.se[e]) and f.se[e] > 0:
                rows.append({"state": int(st), "exposure": e, "pp": 100 * f.params[e], "se": 100 * f.se[e], "n": f.n})
    st = pd.DataFrame(rows)
    pooled = []
    for e, g in st.groupby("exposure"):
        b, v = g["pp"].to_numpy(), g["se"].to_numpy() ** 2
        wi = 1 / v
        fixed = np.sum(wi * b) / np.sum(wi)
        q = np.sum(wi * (b - fixed) ** 2)
        k = len(b)
        tau2 = max(0.0, (q - (k - 1)) / (np.sum(wi) - np.sum(wi ** 2) / np.sum(wi)))
        wr = 1 / (v + tau2)
        re = np.sum(wr * b) / np.sum(wr)
        se_re = np.sqrt(1 / np.sum(wr))
        i2 = max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0
        pooled.append({"outcome": y, "exposure": e, "states": k, "pooled_re_pp": re, "lo": re - 1.96 * se_re,
                       "hi": re + 1.96 * se_re, "tau2": tau2, "i2_pct": 100 * i2,
                       "states_same_sign": int(np.sum(np.sign(b) == np.sign(re))),
                       "pi_lo": re - 1.96 * np.sqrt(tau2 + se_re ** 2), "pi_hi": re + 1.96 * np.sqrt(tau2 + se_re ** 2)})
    return st, pd.DataFrame(pooled)


def main():
    res = {}
    for r in (5, 4):
        d, c = load(r)
        sp = species_of(d)
        sets = [["any_livestock"], sp]
        res[f"secondary_{r}"] = tag(A.fe_models(d, c, ["abx_nonbloody_diarrhoea", "abx_fever_cough"], sets, A.ADJ,
                                                f"nfhs{r}"), "pre-specified secondary outcome")
        res[f"extended_{r}"] = tag(A.fe_models(d, c, A.PRIMARY_OUTCOMES, sets, EXT, f"nfhs{r}_extended"),
                                   "post-hoc: + maternal height and BMI, land area, sex of household head")
        res[f"continuous_{r}"] = tag(A.fe_models(d, c, ["haz_cont", "hb_cont"], sets, A.ADJ, f"nfhs{r}_continuous"),
                                     "post-hoc: continuous outcome (coefficient x100: HAZ in 0.01 SD, Hb in 0.01 g/dl)")
        res[f"oster_{r}"] = tag(oster(d, c, "stunted", sp).assign(round=r), "post-hoc: Oster coefficient stability")
        if r == 5:
            res["clogit_5"] = tag(conditional_logit(d, c, A.PRIMARY_OUTCOMES, sets), "pre-specified sensitivity 6: conditional logit")
            occ_ok = c["occ_husband"].notna() & c["occ_mother"].notna()
            base = A.fe_models(d, c, ["stunted"], [sp], A.ADJ, "occupation_subsample_base", sample=occ_ok)
            adj = A.fe_models(d, c, ["stunted"], [sp], A.ADJ + ["occ_husband", "occ_mother"], "occupation_subsample_adjusted", sample=occ_ok)
            res["occupation_5"] = tag(pd.concat([base, adj]), "post-hoc: occupation recorded (state-module subsample)")
            st, pooled = state_heterogeneity(d, c, "stunted", sp)
            res["state_estimates_5"] = tag(st, "post-hoc: State-specific within-PSU estimates")
            res["state_pooled_5"] = tag(pooled, "post-hoc: random-effects pooling of State estimates")
    for k, v in res.items():
        v.to_csv(os.path.join(OUT, f"review_{k}.csv"), index=False)
        print(f"\n== {k}")
        cols = [x for x in ("analysis", "outcome", "exposure", "pp", "lo", "hi", "p", "or", "delta_for_zero",
                            "beta_star_delta1_pp", "pooled_re_pp", "i2_pct", "states", "states_same_sign", "n") if x in v.columns]
        print(v[cols].round(3).to_string(index=False)[:3500])


if __name__ == "__main__":
    main()
