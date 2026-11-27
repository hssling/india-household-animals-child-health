"""Run the pre-specified analysis (docs/ANALYSIS_PLAN.md v1.0) and write aggregate results to outputs/.

Usage:  python -I src/02_analysis.py --round 5      (primary)
        python -I src/02_analysis.py --round 4      (replication)

Each section is labelled with the plan section it implements. Coding clarifications that the plan did not
spell out are listed in docs/DEVIATIONS.md.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from fe_lpm import fit_fe_lpm  # noqa: E402
from survey_glm import fit_survey_logistic  # noqa: E402

OUT = os.path.join(REPO, "outputs")
PRIMARY_OUTCOMES = ["diarrhoea", "antibiotic", "stunted", "anaemic"]
SECONDARY_OUTCOMES = ["wasted", "bloody"]
SPECIES = ["cattle", "goat_sheep", "poultry", "pig"]


# ---------------------------------------------------------------- plan section 6: covariates
def covariates(d: pd.DataFrame) -> pd.DataFrame:
    c = pd.DataFrame(index=d.index)
    c["age_cat"] = pd.cut(d["age_m"], [-1, 5, 11, 23, 35, 47, 59], labels=["0-5", "6-11", "12-23", "24-35", "36-47", "48-59"])
    c["female"] = (d["b4"] == 2).astype(float)
    c["bord"] = pd.cut(d["bord"], [0, 1, 2, 3, 99], labels=["1", "2", "3", "4+"])
    c["mage"] = pd.cut(d["v012"], [0, 19, 24, 29, 34, 99], labels=["<20", "20-24", "25-29", "30-34", "35+"])
    c["educ"] = d["v106"].map({0: "none", 1: "primary", 2: "secondary", 3: "higher"})
    c["religion"] = d["v130"].map({1: "hindu", 2: "muslim", 3: "christian"}).fillna("other").where(d["v130"].notna())
    c["caste"] = d["s116"].map({1: "sc", 2: "st", 3: "obc", 4: "other"})          # 8 = don't know -> missing
    c["u5"] = pd.cut(d["v137"], [-1, 1, 2, 99], labels=["1", "2", "3+"])
    w = d["v113"]
    c["water_improved"] = np.where(w.isin([32, 42, 43, 96]), 0.0, 1.0)
    c.loc[w.isna() | (w == 97), "water_improved"] = np.nan
    t = d["v116"]
    c["toilet"] = np.select([t.isin([11, 12, 13, 15, 21, 22, 41]), t == 31], ["improved", "none"], "unimproved")
    c.loc[t.isna() | (t == 97), "toilet"] = np.nan
    f = d["v161"]
    c["fuel"] = np.select([f.isin([1, 2, 3, 4, 95]), f == 11], ["clean", "dung"], "solid")
    c.loc[f.isna() | (f == 97), "fuel"] = np.nan
    s = d["v465"]
    c["stool"] = np.select([s.isin([1, 2, 5]), s.isin([3, 4, 9, 96])], ["safe", "unsafe"], "not_asked")
    c["land"] = (d["hv244"] == 1).astype(float).where(d["hv244"].isin([0, 1]))
    c["month"] = d["v006"].astype("Int64").astype(str)
    return c


def asset_index(d: pd.DataFrame) -> pd.Series:
    """Plan section 6: first principal component of goods and housing materials, household level."""
    items = [x for x in d.columns if x.startswith("asset_")]
    hh = d.groupby("hh_id")[items].first()
    X = pd.DataFrame(index=hh.index)
    for c in items:
        v = hh[c]
        if c in ("asset_hv213", "asset_hv214", "asset_hv215"):
            X[c] = ((v >= 31) & (v <= 37)).astype(float).where(v.notna() & (v != 97))
        else:
            X[c] = (v == 1).astype(float).where(v.isin([0, 1]))
    X = X.dropna()
    Z = (X - X.mean()) / X.std(ddof=0)
    _, _, vt = np.linalg.svd(Z.to_numpy(), full_matrices=False)
    pc1 = Z.to_numpy() @ vt[0]
    if np.corrcoef(pc1, X.sum(axis=1))[0, 1] < 0:          # orient so that more assets = higher score
        pc1 = -pc1
    q = pd.qcut(pd.Series(pc1, index=X.index), 5, labels=["Q1", "Q2", "Q3", "Q4", "Q5"])
    return d["hh_id"].map(q)


def design(c: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    parts = []
    for col in cols:
        v = c[col]
        if v.dtype.kind in "fi" and set(pd.unique(v.dropna())) <= {0.0, 1.0}:
            parts.append(v.rename(col).astype(float))
        else:
            parts.append(pd.get_dummies(v.astype("string"), prefix=col, drop_first=True, dtype=float))
    return pd.concat(parts, axis=1)


ADJ = ["age_cat", "female", "bord", "mage", "educ", "religion", "caste", "u5", "water_improved", "toilet",
       "fuel", "stool", "land", "assetq", "month"]


def bh(p: pd.Series) -> pd.Series:
    p = p.astype(float)
    order = p.sort_values()
    m = len(order)
    q = (order * m / np.arange(1, m + 1))[::-1].cummin()[::-1].clip(upper=1)
    return q.reindex(p.index)


def fe_models(d, c, outcomes, exposures_sets, adj, label, cluster_col="psu", sample=None):
    rows = []
    for y in outcomes:
        for expo in exposures_sets:
            cols = expo + adj
            frame = pd.concat([d[[y, "psu", "district", "weight"]], d[expo], c[[a for a in adj]]], axis=1)
            if sample is not None:
                frame = frame[sample]
            frame = frame.dropna()
            X = design(frame.drop(columns=[y, "psu", "district", "weight"]).assign(**{e: frame[e] for e in expo}), cols)
            fit = fit_fe_lpm(frame[y].to_numpy(), X, groups=frame["psu"].to_numpy(),
                             clusters=frame[cluster_col].to_numpy(), weights=frame["weight"].to_numpy())
            for e in expo:
                if e in fit.params.index:
                    rows.append({"analysis": label, "outcome": y, "exposure": e, "model": "+".join(expo),
                                 "pp": 100 * fit.params[e], "lo": 100 * fit.lo[e], "hi": 100 * fit.hi[e],
                                 "p": fit.p[e], "n": fit.n, "clusters": fit.n_clusters,
                                 "prevalence": 100 * np.average(frame[y], weights=frame["weight"])})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--round", type=int, choices=[4, 5], default=5)
    a = ap.parse_args()
    d = pd.read_parquet(os.path.join(REPO, "data_local", f"nfhs{a.round}_children.parquet"))
    d = d[(d["age_m"] >= 0) & (d["age_m"] <= 59)].copy()
    d["district"] = d["district"].fillna(-1)
    c = covariates(d)
    c["assetq"] = asset_index(d)
    species = [s for s in SPECIES if d[s].notna().any()]
    tag = f"nfhs{a.round}"
    os.makedirs(OUT, exist_ok=True)

    # descriptive counts (plan section 10, Table 1 inputs)
    desc = []
    for e in ["any_livestock"] + species:
        for v in (0, 1):
            s = d[e] == v
            row = {"exposure": e, "level": v, "children": int(s.sum()),
                   "weighted_pct": 100 * d.loc[s, "weight"].sum() / d.loc[d[e].notna(), "weight"].sum()}
            for y in PRIMARY_OUTCOMES + SECONDARY_OUTCOMES:
                m = s & d[y].notna()
                row[y] = 100 * np.average(d.loc[m, y], weights=d.loc[m, "weight"]) if m.any() else np.nan
            desc.append(row)
    pd.DataFrame(desc).to_csv(os.path.join(OUT, f"{tag}_descriptive.csv"), index=False)

    # plan section 7: primary within-cluster models, 5 exposures x 4 outcomes
    prim = fe_models(d, c, PRIMARY_OUTCOMES, [["any_livestock"], species], ADJ, "primary")
    prim["q"] = bh(prim["p"])
    prim.to_csv(os.path.join(OUT, f"{tag}_primary.csv"), index=False)
    print(prim[["outcome", "exposure", "pp", "lo", "hi", "p", "q", "n"]].round(3).to_string(index=False))

    # secondary outcomes
    fe_models(d, c, SECONDARY_OUTCOMES, [["any_livestock"], species], ADJ, "secondary_outcomes") \
        .to_csv(os.path.join(OUT, f"{tag}_secondary_outcomes.csv"), index=False)

    # plan section 7: effect modification (improved toilet; age 6-23 months)
    em = []
    for mod_name, mod in (("improved_toilet", (c["toilet"] == "improved").astype(float).where(c["toilet"].notna())),
                          ("age_6_23", ((d["age_m"] >= 6) & (d["age_m"] <= 23)).astype(float))):
        for y in PRIMARY_OUTCOMES:
            for expo in (["any_livestock"], species):
                frame = pd.concat([d[[y, "psu", "weight"]], d[expo], c[ADJ]], axis=1).assign(_mod=mod).dropna()
                X = design(frame.drop(columns=[y, "psu", "weight"]), expo + ADJ + ["_mod"])
                for e in expo:
                    X[f"{e}_x_mod"] = frame[e] * frame["_mod"]
                fit = fit_fe_lpm(frame[y].to_numpy(), X, frame["psu"].to_numpy(), frame["psu"].to_numpy(),
                                 frame["weight"].to_numpy())
                for e in expo:
                    k = f"{e}_x_mod"
                    if k in fit.params.index:
                        em.append({"modifier": mod_name, "outcome": y, "exposure": e, "interaction_pp": 100 * fit.params[k],
                                   "lo": 100 * fit.lo[k], "hi": 100 * fit.hi[k], "p": fit.p[k], "n": fit.n})
    pd.DataFrame(em).to_csv(os.path.join(OUT, f"{tag}_effect_modification.csv"), index=False)

    if a.round == 5:
        # plan section 7: secondary logistic model with state fixed effects (odds ratios)
        lg = []
        for y in PRIMARY_OUTCOMES:
            for expo in (["any_livestock"], species):
                frame = pd.concat([d[[y, "psu", "weight", "state"]], d[expo], c[ADJ]], axis=1).dropna()
                frame["state"] = frame["state"].astype(int).astype(str)
                terms = " + ".join(expo + [f"C({x})" if frame[x].dtype.kind not in "fi" else x for x in ADJ] + ["C(state)"])
                for x in ADJ:
                    if frame[x].dtype.kind not in "fi":
                        frame[x] = frame[x].astype(str)
                fit = fit_survey_logistic(f"{y} ~ {terms}", frame, frame["weight"].to_numpy(), frame["psu"].to_numpy())
                for e in expo:
                    lg.append({"outcome": y, "exposure": e, "or": np.exp(fit.params[e]), "lo": np.exp(fit.conf_low[e]),
                               "hi": np.exp(fit.conf_high[e]), "p": fit.p_values[e], "n": fit.n_observations,
                               "prevalence": 100 * np.average(frame[y], weights=frame["weight"])})
        lg = pd.DataFrame(lg)

        # plan section 8.3: E-values from the logistic estimates (RR ~ sqrt(OR) when the outcome is common)
        def evalue(rr):
            rr = rr if rr >= 1 else 1 / rr
            return rr + np.sqrt(rr * (rr - 1))
        rr = np.where(lg["prevalence"] > 15, np.sqrt(lg["or"]), lg["or"])
        rr_lo = np.where(lg["prevalence"] > 15, np.sqrt(lg["lo"]), lg["lo"])
        rr_hi = np.where(lg["prevalence"] > 15, np.sqrt(lg["hi"]), lg["hi"])
        lg["evalue_point"] = [evalue(r) for r in rr]
        lg["evalue_ci"] = [1.0 if (lo <= 1 <= hi) else evalue(lo if r > 1 else hi) for r, lo, hi in zip(rr, rr_lo, rr_hi)]
        lg.to_csv(os.path.join(OUT, f"{tag}_logistic_state_fe.csv"), index=False)

        # plan section 7: animal-source food, 6-23 months, youngest child living with the mother
        y623 = (d["age_m"] >= 6) & (d["age_m"] <= 23) & (d["midx"] == 1) & (d["b9"] == 0)
        flesh = d[["diet_v414a", "diet_v414g", "diet_v414m", "diet_v414n", "diet_v414t"]].eq(1).any(axis=1).astype(float)
        dairy = d[["diet_v411", "diet_v414p", "diet_v414v"]].eq(1).any(axis=1).astype(float)
        dvalid = d[["diet_v414g", "diet_v411"]].isin([0, 1]).all(axis=1)
        c2 = c.assign(flesh_eggs=flesh.where(dvalid), dairy=dairy.where(dvalid))
        med = []
        for y in ("stunted", "anaemic"):
            base = fe_models(d, c2, [y], [species], ADJ, "asf_base", sample=y623 & dvalid)
            plus = fe_models(d, c2, [y], [species], ADJ + ["flesh_eggs", "dairy"], "asf_adjusted", sample=y623 & dvalid)
            med += [base, plus]
        pd.concat(med).to_csv(os.path.join(OUT, f"{tag}_animal_source_food.csv"), index=False)

        # plan section 8: sensitivity analyses
        sens = []
        adj_off = [x for x in ADJ if x != "assetq"]
        c_off = c.assign(official_wealth=d["v190"].astype("Int64").astype(str))
        sens.append(fe_models(d, c_off, PRIMARY_OUTCOMES, [["any_livestock"], species], adj_off + ["official_wealth"], "s1_official_wealth_index"))
        # s2 negative control: land without animals vs neither, among households without animals
        no_anim = d["any_livestock"] == 0
        d2 = d.assign(land_no_animals=c["land"])
        sens.append(fe_models(d2, c, PRIMARY_OUTCOMES, [["land_no_animals"]], [x for x in ADJ if x != "land"], "s2_negative_control_land", sample=no_anim))
        sens.append(fe_models(d, c, PRIMARY_OUTCOMES, [["any_livestock"], species], ADJ, "s4_rural_only", sample=d["rural"] == 1))
        sens.append(fe_models(d, c, PRIMARY_OUTCOMES, [["any_livestock"], species], ADJ, "s5_cluster_district", cluster_col="district"))
        pd.concat(sens).to_csv(os.path.join(OUT, f"{tag}_sensitivity.csv"), index=False)
    print(f"\nwrote outputs/{tag}_*.csv")


if __name__ == "__main__":
    main()
