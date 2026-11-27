"""Descriptive inputs for the paper: sample flow, Table 1 and the State-level data for Figure 2.

Usage:  python -I src/03_descriptives.py
Writes aggregate tables only (outputs/nfhs5_flow.csv, nfhs5_table1.csv, nfhs5_state_poultry_antibiotic.csv,
nfhs4_flow.csv).
"""
import importlib.util
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(REPO, "outputs")

spec = importlib.util.spec_from_file_location("analysis", os.path.join(HERE, "02_analysis.py"))
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)


def wpct(x, w):
    m = x.notna()
    return 100 * np.average(x[m], weights=w[m]) if m.any() else np.nan


def flow(round_):
    d = pd.read_parquet(os.path.join(REPO, "data_local", f"nfhs{round_}_children.parquet"))
    d = d[(d["age_m"] >= 0) & (d["age_m"] <= 59)].copy()
    c = A.covariates(d)
    c["assetq"] = A.asset_index(d)
    expo = ["any_livestock", "cattle", "goat_sheep", "poultry"] + (["pig"] if d["pig"].notna().any() else [])
    rows = [("Living children aged 0-59 months", len(d))]
    for y in A.PRIMARY_OUTCOMES:
        has_y = d[y].notna()
        full = has_y & d[expo].notna().all(axis=1) & c[A.ADJ].notna().all(axis=1)
        rows.append((f"{y}: outcome recorded", int(has_y.sum())))
        rows.append((f"{y}: analysed (complete exposures and covariates)", int(full.sum())))
    pd.DataFrame(rows, columns=["step", "children"]).to_csv(os.path.join(OUT, f"nfhs{round_}_flow.csv"), index=False)
    return d, c


def main():
    d, c = flow(5)
    flow(4)
    w = d["weight"]

    # Table 1: weighted percentages by any livestock (NFHS-5)
    rows = []
    def add(label, series_by_level):
        r = {"characteristic": label}
        for lvl, name in ((0, "no_animals"), (1, "animals")):
            r[name] = series_by_level(d["any_livestock"] == lvl)
        r["all"] = series_by_level(d["any_livestock"].notna())
        rows.append(r)
    n = {lvl: int((d["any_livestock"] == lvl).sum()) for lvl in (0, 1)}
    rows.append({"characteristic": "Children (unweighted n)", "no_animals": n[0], "animals": n[1],
                 "all": int(d["any_livestock"].notna().sum())})
    for col, lab in (("cattle", "Owns cows, bulls or buffaloes"), ("goat_sheep", "Owns goats or sheep"),
                     ("poultry", "Owns chickens or ducks"), ("pig", "Owns pigs")):
        add(lab, lambda s, col=col: wpct(d.loc[s, col], w[s]))
    add("Rural residence", lambda s: wpct(d.loc[s, "rural"].astype(float), w[s]))
    for lab, cond in (("Mother: no schooling", c["educ"] == "none"), ("Mother: secondary or higher", c["educ"].isin(["secondary", "higher"])),
                      ("Scheduled caste or tribe", c["caste"].isin(["sc", "st"])),
                      ("Asset index, lowest quintile", c["assetq"] == "Q1"), ("Asset index, highest quintile", c["assetq"] == "Q5"),
                      ("Improved drinking water", c["water_improved"] == 1), ("Improved toilet facility", c["toilet"] == "improved"),
                      ("No toilet facility", c["toilet"] == "none"), ("Cooks with animal dung", c["fuel"] == "dung"),
                      ("Owns agricultural land", c["land"] == 1)):
        valid = cond.notna() & (c.loc[:, "educ"].notna() if lab.startswith("Mother") else True)
        x = cond.astype(float).where(valid)
        add(lab, lambda s, x=x: wpct(x[s], w[s]))
    for y, lab in (("diarrhoea", "Diarrhoea in past 2 weeks"), ("antibiotic", "Antibiotic in past 2 weeks"),
                   ("stunted", "Stunted"), ("anaemic", "Anaemic (6-59 months)")):
        add(lab, lambda s, y=y: wpct(d.loc[s, y], w[s]))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "nfhs5_table1.csv"), index=False)

    # Socioeconomic profile of households by species kept (children's weights), to describe who keeps which
    # animals; used to interpret the species-specific growth associations.
    prof = []
    kinds = {"Cattle or buffalo, no goats or sheep": (d["cattle"] == 1) & (d["goat_sheep"] == 0),
             "Goats or sheep, no cattle or buffalo": (d["goat_sheep"] == 1) & (d["cattle"] == 0),
             "Both": (d["cattle"] == 1) & (d["goat_sheep"] == 1),
             "Neither (may keep poultry or pigs)": (d["cattle"] == 0) & (d["goat_sheep"] == 0)}
    for lab, m in kinds.items():
        x = d[m]; cc = c[m]
        prof.append({"households_keeping": lab, "children": int(m.sum()),
                     "asset_lowest_quintile": wpct((cc["assetq"] == "Q1").astype(float).where(cc["assetq"].notna()), x["weight"]),
                     "asset_highest_quintile": wpct((cc["assetq"] == "Q5").astype(float).where(cc["assetq"].notna()), x["weight"]),
                     "owns_agricultural_land": wpct(cc["land"], x["weight"]),
                     "scheduled_caste_or_tribe": wpct(cc["caste"].isin(["sc", "st"]).astype(float).where(cc["caste"].notna()), x["weight"]),
                     "mother_no_schooling": wpct((cc["educ"] == "none").astype(float).where(cc["educ"].notna()), x["weight"]),
                     "stunted": wpct(x["stunted"], x["weight"])})
    pd.DataFrame(prof).to_csv(os.path.join(OUT, "nfhs5_species_profile.csv"), index=False)

    # Figure 2 data: State-level weighted poultry ownership and child antibiotic exposure
    g = d.dropna(subset=["poultry", "antibiotic"])
    st = g.groupby("state").apply(lambda x: pd.Series({
        "children": len(x), "poultry_pct": wpct(x["poultry"], x["weight"]),
        "antibiotic_pct": wpct(x["antibiotic"], x["weight"]), "goat_sheep_pct": wpct(x["goat_sheep"], x["weight"]),
        "cattle_pct": wpct(x["cattle"], x["weight"]), "stunted_pct": wpct(x["stunted"], x["weight"])}), include_groups=False)
    st.reset_index().to_csv(os.path.join(OUT, "nfhs5_state_poultry_antibiotic.csv"), index=False)
    print("wrote flow, Table 1 and State-level files")


if __name__ == "__main__":
    main()
