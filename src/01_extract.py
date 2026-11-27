"""Build the analysis extract from the NFHS children's (KR) and household (HR) recodes.

Reads only the columns the analysis plan needs, links children to households, and writes one
de-identified parquet file per survey round to data_local/ (never committed). Before writing, it checks
raw counts against the survey's own frequency (.FRQ) files and stops on any mismatch, so a wrong file,
a wrong link or a wrong recode cannot pass silently.

Usage:
  python -I src/01_extract.py --round 5 --zipdir "D:/Anemia research/NFHS 5 microdata" --workdir C:/nfhs_work
  python -I src/01_extract.py --round 4 --zipdir "D:/Anemia research/NFHS 4 microdata" --workdir C:/nfhs_work

The .DTA files are unpacked into --workdir (use a drive with space; the NFHS-5 household file is 5.2 GB)
and deleted afterwards unless --keep is given.
"""
import argparse
import os
import re
import sys
import zipfile

import pandas as pd
import pyreadstat

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(REPO, "data_local")

ROUND = {
    5: {"kr": ("IAKR7EDT.zip", "IAKR7EFL.DTA", "IAKR7EFL.FRQ"),
        "hr": ("IAHR7EDT.zip", "IAHR7EFL.DTA", "IAHR7EFL.FRQ"),
        "species": {"cattle": "sh67a", "goat_sheep": ["sh67d"], "poultry": "sh67f", "pig": "sh67e"}},
    4: {"kr": ("IAKR74DT.zip", "IAKR74FL.DTA", None),
        "hr": ("IAHR74DT.zip", "IAHR74FL.DTA", None),
        "species": {"cattle": "sh52a", "goat_sheep": ["sh52d", "sh52e"], "poultry": "sh52f", "pig": None}},
}

KR_COLS = ["v001", "v002", "v003", "bidx", "midx", "b5", "b4", "b9", "b19", "b3", "v008", "hw1", "bord",
           "v005", "v021", "v022", "v023", "v024", "v025", "sdist", "sdistri", "v006", "v007", "v012", "v106", "v130",
           "s116", "v137", "v113", "v116", "v161", "v465", "v190",
           "h11", "h11b", "h15", "h15b", "ml13n", "h22", "h31", "hw70", "hw72", "hw56",
           "v438", "v445", "v704", "v705", "v717", "v151",
           "v411", "v414a", "v414g", "v414m", "v414n", "v414t", "v414p", "v414v"]
# Household goods and housing materials for the rebuilt asset index (no animals, land, water, toilet or fuel).
ASSETS = ["hv206", "hv207", "hv208", "hv209", "hv210", "hv211", "hv212", "hv221", "hv243a", "hv243b",
          "hv243e", "hv213", "hv214", "hv215"]
HR_BASE = ["hv001", "hv002", "hv246", "hv244", "hv245", "hv025"]


def frq_count(zippath, frq, item, code):
    """Unweighted frequency of `code` for `item` in the survey's own .FRQ file."""
    t = zipfile.ZipFile(zippath).read(frq).decode("latin-1")
    m = re.search(r"Item %s:.*?(?=\n\s*Item |\Z)" % re.escape(item.upper()), t, re.S)
    if not m:
        raise SystemExit(f"{item} not in {frq}")
    for line in m.group(0).splitlines():
        mm = re.match(r"\s*(-?\d+)\s+.*?\s(\d+)\s+\d+\s", line)
        if mm and int(mm.group(1)) == code:
            return int(mm.group(2))
    raise SystemExit(f"code {code} of {item} not in {frq}")


def unpack(zipdir, zname, member, workdir):
    os.makedirs(workdir, exist_ok=True)
    dest = os.path.join(workdir, member)
    if not os.path.exists(dest):
        with zipfile.ZipFile(os.path.join(zipdir, zname)) as z, open(dest, "wb") as fh:
            with z.open(member) as src:
                while True:
                    buf = src.read(1 << 24)
                    if not buf:
                        break
                    fh.write(buf)
    return dest


def read(path, cols):
    _, meta = pyreadstat.read_dta(path, metadataonly=True)
    have = {c.lower(): c for c in meta.column_names}
    use = [have[c] for c in cols if c in have]
    missing = [c for c in cols if c not in have]
    df, _ = pyreadstat.read_dta(path, usecols=use, apply_value_formats=False)
    df.columns = [c.lower() for c in df.columns]
    return df, missing


def yes(s, code=1):
    """1 if equal to `code`, 0 if a valid other answer, missing for don't-know (8) or missing."""
    out = (s == code).astype("float")
    out[s.isna() | (s == 8)] = float("nan")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--round", type=int, choices=[4, 5], required=True)
    ap.add_argument("--zipdir", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    R = ROUND[a.round]
    species_cols = [c for v in R["species"].values() if v for c in ([v] if isinstance(v, str) else v)]

    kr_path = unpack(a.zipdir, R["kr"][0], R["kr"][1], a.workdir)
    hr_path = unpack(a.zipdir, R["hr"][0], R["hr"][1], a.workdir)
    kr, kr_missing = read(kr_path, KR_COLS)
    hr, hr_missing = read(hr_path, HR_BASE + species_cols + ASSETS)
    print(f"NFHS-{a.round}: KR {kr.shape}, HR {hr.shape}")
    if kr_missing or hr_missing:
        print("  columns absent (handled downstream):", kr_missing, hr_missing)

    # --- raw counts against the survey's frequency file (NFHS-5; NFHS-4 checked against NFHS-4 FRQ below)
    frq_zip_kr = os.path.join(a.zipdir, R["kr"][0])
    frq_name_kr = R["kr"][2] or next(n for n in zipfile.ZipFile(frq_zip_kr).namelist() if n.upper().endswith(".FRQ"))
    frq_zip_hr = os.path.join(a.zipdir, R["hr"][0])
    frq_name_hr = R["hr"][2] or next(n for n in zipfile.ZipFile(frq_zip_hr).namelist() if n.upper().endswith(".FRQ"))
    checks = [(kr, frq_zip_kr, frq_name_kr, "h11", 2), (kr, frq_zip_kr, frq_name_kr, "h15", 1),
              (kr, frq_zip_kr, frq_name_kr, "h11b", 1), (kr, frq_zip_kr, frq_name_kr, "ml13n", 1),
              (hr, frq_zip_hr, frq_name_hr, "hv246", 1), (hr, frq_zip_hr, frq_name_hr, R["species"]["poultry"], 1),
              (hr, frq_zip_hr, frq_name_hr, R["species"]["cattle"], 1)]
    bad = 0
    for df, zp, fn, item, code in checks:
        got, want = int((df[item] == code).sum()), frq_count(zp, fn, item, code)
        status = "ok" if got == want else "MISMATCH"
        bad += got != want
        print(f"  {item:6s} = {code}: extract {got:>8,}  frequency file {want:>8,}  {status}")
    if bad:
        raise SystemExit("raw counts disagree with the frequency file; stopping")

    # --- link children to households (1:1 on cluster and household number)
    hr = hr.rename(columns={"hv001": "v001", "hv002": "v002"})
    if hr.duplicated(["v001", "v002"]).any():
        raise SystemExit("household key not unique in HR")
    df = kr.merge(hr, on=["v001", "v002"], how="left", validate="m:1", indicator=True)
    unlinked = int((df["_merge"] != "both").sum())
    print(f"  children linked to a household: {len(df) - unlinked:,} of {len(df):,}")
    if unlinked:
        raise SystemExit(f"{unlinked} children did not link to a household")
    df = df.drop(columns="_merge")

    # --- derived variables, exactly as in docs/ANALYSIS_PLAN.md
    d = pd.DataFrame({"round": a.round, "psu": df["v021"], "strata": df["v022"], "state": df["v024"],
                      "district": df["sdist"] if "sdist" in df else df.get("sdistri"), "weight": df["v005"] / 1e6,
                      "rural": (df["v025"] == 2).astype(int)})
    d["alive"] = (df["b5"] == 1).astype(int)
    # Age in months: b19 where the recode has it (NFHS-5); otherwise interview date minus birth date in
    # century-month codes (NFHS-4), which is defined for every child, measured or not.
    age = df["b19"] if "b19" in df else (df["v008"] - df["b3"])
    d["age_m"] = age
    d["child_id"] = df["v001"].astype(str) + "-" + df["v002"].astype(str) + "-" + df["v003"].astype(str) + "-" + df["bidx"].astype(str)
    d["hh_id"] = df["v001"].astype(str) + "-" + df["v002"].astype(str)
    # exposures
    d["any_livestock"] = yes(df["hv246"])
    sp = R["species"]
    d["cattle"] = yes(df[sp["cattle"]])
    gs = [yes(df[c]) for c in sp["goat_sheep"]]
    d["goat_sheep"] = pd.concat(gs, axis=1).max(axis=1, skipna=True) if len(gs) > 1 else gs[0]
    d["poultry"] = yes(df[sp["poultry"]])
    d["pig"] = yes(df[sp["pig"]]) if sp["pig"] else float("nan")
    # outcomes
    d["diarrhoea"] = (df["h11"] == 2).astype(float)
    d.loc[df["h11"].isna() | (df["h11"] == 8), "diarrhoea"] = float("nan")
    abx_d = ((df["h15"] == 1) | (df["h15b"] == 1))
    abx_f = (df["ml13n"] == 1)
    d["antibiotic"] = (abx_d | abx_f).astype(float)
    d.loc[d["diarrhoea"].isna(), "antibiotic"] = float("nan")
    d["bloody"] = yes(df["h11b"])
    d["abx_diarrhoea"] = abx_d.astype(float).where(d["diarrhoea"] == 1)
    # Plan section 5 secondary outcomes: antibiotic among children with non-bloody diarrhoea, and among
    # children with fever or cough.
    d["abx_nonbloody_diarrhoea"] = abx_d.astype(float).where((d["diarrhoea"] == 1) & (df["h11b"] == 0))
    fc = (df["h22"] == 1) | (df["h31"].isin([1, 2]))
    d["fever_cough"] = fc.astype(float).where(df["h22"].isin([0, 1]) | df["h31"].isin([0, 1, 2]))
    d["abx_fever_cough"] = abx_f.astype(float).where(fc)
    # Additional covariates (post-hoc review analyses; see docs/DEVIATIONS.md)
    d["mheight_cm"] = (df["v438"] / 10).where((df["v438"] >= 1000) & (df["v438"] <= 2200))
    d["mbmi"] = (df["v445"] / 100).where((df["v445"] >= 1200) & (df["v445"] <= 6000))
    for c in ("v704", "v705", "v717", "v151", "hv245"):
        if c in df:
            d[c] = df[c]
    hw70 = df["hw70"].where(df["hw70"] < 9996)
    d["haz"] = hw70 / 100
    d["stunted"] = (hw70 < -200).astype(float).where(hw70.notna())
    hw72 = df["hw72"].where(df["hw72"] < 9996)
    d["wasted"] = (hw72 < -200).astype(float).where(hw72.notna())
    hb = df["hw56"].where((df["hw56"] >= 10) & (df["hw56"] < 250))
    d["hb"] = hb / 10
    d["anaemic"] = (hb < 110).astype(float).where(hb.notna() & (age >= 6))
    # covariates (codes are those of the NFHS recode; recoding categories happens in the analysis script)
    for c in ["b4", "bord", "v012", "v106", "v130", "s116", "v137", "v113", "v116", "v161", "v465", "v190",
              "v006", "v007", "hv244", "midx", "b9"]:
        if c in df:
            d[c] = df[c]
    for c in ASSETS:
        if c in df:
            d["asset_" + c] = df[c]
    for c in ["v411", "v414a", "v414g", "v414m", "v414n", "v414t", "v414p", "v414v"]:
        if c in df:
            d["diet_" + c] = df[c]

    d = d[d["alive"] == 1].drop(columns="alive")
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"nfhs{a.round}_children.parquet")
    d.to_parquet(path, index=False)
    print(f"  wrote {os.path.relpath(path, REPO)}: {len(d):,} living children, {d.shape[1]} columns")
    if not a.keep:
        for p in (kr_path, hr_path):
            os.remove(p)
        print("  removed unpacked .DTA files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
