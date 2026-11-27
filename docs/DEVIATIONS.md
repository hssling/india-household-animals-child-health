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
