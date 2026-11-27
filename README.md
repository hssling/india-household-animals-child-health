# Household animals, child diarrhoea, antibiotic use and growth in India

Code and derived outputs for a study by Siddalingaiah H S.

**Status:** v1.0.0, code and aggregate results accompanying the manuscript submitted to the Online Journal of Health and Allied Sciences. The analysis plan (`docs/ANALYSIS_PLAN.md`) was frozen on 27 November 2026 before any outcome was analysed and released as v0.1.0; changes and additions are logged in `docs/DEVIATIONS.md`.

**Archived:** analysis plan v0.1.0, doi:[10.5281/zenodo.24596396](https://doi.org/10.5281/zenodo.24596396) (concept DOI, always the latest version: [10.5281/zenodo.24596395](https://doi.org/10.5281/zenodo.24596395)).

## Data

NFHS-4 and NFHS-5 unit-level files from the DHS Program (registration required). Unit-level survey files are not redistributed; obtain them from the original provider.
Only code and aggregate outputs are kept here.

## Layout

| Folder | Content |
|---|---|
| `docs/` | Analysis plan and data notes |
| `src/` | Analysis code |
| `outputs/` | Aggregate results tables (no unit-level data) |
| `tests/` | Checks of the data build and of the results |

## How to reproduce

Requires Python 3 with pandas, numpy, scipy, statsmodels, pyreadstat and pyarrow, and the NFHS-4 and NFHS-5 children's (KR) and household (HR) recode files from the DHS Program.

| Step | Script | Output |
|---|---|---|
| 1 | `src/01_extract.py --round 5` and `--round 4` | `data_local/` extracts (not published); stops if counts differ from the survey frequency files |
| 2 | `src/02_analysis.py --round 5` and `--round 4` | Pre-specified primary, secondary, effect-modification and sensitivity results |
| 3 | `src/03_descriptives.py` | Sample flow, Table 1, species profile, State-level data for Figure 2 |
| 4 | `src/04_review_analyses.py` | Omitted pre-specified analyses and the labelled post-hoc robustness analyses |
| - | `tests/test_fe_lpm.py` | Verifies the fixed-effects estimator against statsmodels |

Result files in `outputs/` are prefixed `nfhs5_`, `nfhs4_` (plan analyses) and `review_` (analyses added after the primary results; see `docs/DEVIATIONS.md`).

## Licence

Code: MIT (`LICENSE`). Derived aggregate outputs: CC BY 4.0 (`LICENSE-DATA`).

## Citation

See `CITATION.cff`; cite the concept DOI 10.5281/zenodo.24596395.
