# Deviations from and clarifications of the analysis plan (v1.0)

Each entry gives the date, what changed, and why. Clarifications fill in coding detail the plan did not spell
out; they were decided before the related model was run.

| Date | Type | Item | Detail and reason |
|---|---|---|---|
| 2026-11-27 | Clarification | Household items coded 97 ("not a de jure resident") | Water, toilet and cooking fuel are not recorded for visiting mothers (11,071 NFHS-5 children); set to missing and excluded in complete-case analysis. |
| 2026-11-27 | Clarification | Stool disposal (V465) | Asked only about the mother's youngest child under 2. Other children form a "not asked" category rather than being excluded; "don't know" joins that category. |
| 2026-11-27 | Clarification | Category codes | Improved water: all codes except unprotected well or spring, surface water and other. Improved toilet: flush to sewer, septic tank or pit; flush, location unknown; ventilated improved pit; pit with slab; composting. "No facility" kept separate. Clean fuel: electricity, LPG, natural gas, biogas, and no food cooked. Finished floor, walls and roof: codes 31–37. |
| 2026-11-27 | Clarification | NFHS-4 child age | NFHS-4 has no `b19`; age in months is interview date minus birth date (V008 − B3), defined for every child. |
| 2026-11-27 | Clarification | Confidence intervals | Fixed-effects intervals use the t distribution with (clusters − 1) degrees of freedom. |
| 2026-11-27 | Omission corrected | Secondary outcomes (plan section 5) | Antibiotic among children with non-bloody diarrhoea, and among children with fever or cough, were pre-specified but not run in the first analysis; added at internal review (src/04_review_analyses.py). |
| 2026-11-27 | Omission corrected | Sensitivity analysis 6 (plan section 8) | Conditional logit with cluster fixed effects was pre-specified but not run in the first analysis; added at internal review. Unweighted, as conditional logit does not take survey weights. |
| 2026-11-27 | Post-hoc addition | Extended adjustment | Added after the primary results were seen, to address residual confounding: maternal height and body mass index, agricultural land area and sex of the household head. |
| 2026-11-27 | Post-hoc addition | Continuous outcomes | Height-for-age z-score and haemoglobin as continuous outcomes, so the stunting and anaemia results do not depend on cut-offs. |
| 2026-11-27 | Post-hoc addition | Occupation | Husband's and mother's occupation are recorded only in the state-module subsample (about 15% of children); the stunting model is compared with and without occupation in that subsample. |
| 2026-11-27 | Post-hoc addition | Coefficient stability | Oster (2019) bounds for the stunting associations within primary sampling units (short model: species and fixed effects; long model: primary model; R-max = 1.3 x long-model within R-squared). |
| 2026-11-27 | Post-hoc addition | State heterogeneity | State-specific within-PSU estimates of the stunting associations (States with at least 50 primary sampling units), pooled by random-effects meta-analysis with I-squared and prediction intervals. |
| 2026-11-27 | Authorship | Paper and repository | The paper is single-authored (Siddalingaiah H S). The frozen analysis plan v1.0 names two authors and is left unchanged as the dated record; repository and Zenodo metadata now list the single author. |
