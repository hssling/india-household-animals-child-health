# Analysis plan: household animals and children's diarrhoea, antibiotic exposure and growth in India

Version 1.0, frozen 27 November 2026, before any association between an exposure and an outcome was
estimated. Authors: Siddalingaiah H S, Preethi Ashok Masali (Shridevi Institute of Medical Sciences and
Research Hospital, Tumkur, India).

**What had been seen at freezing.** Only the published unweighted marginal frequencies in the survey's
frequency files: for example, in NFHS-5, 15,334 children with diarrhoea in the past two weeks, 2,346 given
antibiotic pills or syrups for it, and 51% of households owning livestock. No cross-tabulation of an
exposure with an outcome had been made.

Any departure from this plan will be reported in the paper as a deviation, with its reason.

## 1. Background and aim

Half of Indian households keep animals, and national schemes promote backyard poultry and dairy for income
and nutrition. Animals may support child growth through animal-source food and income, and harm it through
enteric infection and the antibiotic courses that follow. A scoping review of 911 Indian studies found only
4 that examined nutrition, zoonotic infection and antimicrobial resistance together.

Aim: among children under five in India, compare households within the same village or urban block, and
estimate the association of household animal ownership with diarrhoea, antibiotic exposure, stunting and
anaemia.

## 2. Data

- Primary: National Family Health Survey 2019–21 (NFHS-5), children's recode (KR) linked to the household
  recode (HR) by cluster and household number.
- Replication: NFHS-4 (2015–16), same recodes and the same models.
- Unit-level files are obtained from the DHS Program under its terms of use and are not redistributed.

## 3. Sample

Living children (B5 = 1) aged 0–59 months, born to interviewed women.
- Diarrhoea and antibiotic outcomes: all living children 0–59 months.
- Stunting: children with a valid height-for-age z-score (HW70 below 9996; flagged values excluded).
- Anaemia: children 6–59 months with a valid haemoglobin value (HW56).
- "Don't know" answers are set to missing. Analyses are complete-case, and the number excluded is reported
  for each outcome.

## 4. Exposures (household)

| Exposure | NFHS-5 | NFHS-4 |
|---|---|---|
| Any livestock, herds or farm animals | HV246 | HV246 |
| Cows, bulls or buffaloes | SH67A | SH52A |
| Goats or sheep | SH67D | SH52D or SH52E |
| Chickens or ducks | SH67F | SH52F |
| Pigs | SH67E | not available |

"Any livestock" is modelled on its own. The four species are entered together, so each is adjusted for the
others.

## 5. Outcomes (child)

Primary:
1. Diarrhoea in the past two weeks (H11 = 2).
2. Antibiotic exposure in the past two weeks: antibiotic pills or syrup (H15) or injection (H15B) for
   diarrhoea, or an antibiotic for fever or cough (ML13N). The denominator is all living children 0–59
   months; children without these illnesses count as unexposed. This is a population exposure, not a
   treatment rate among sick children.
3. Stunting: HW70 below −200 (z-score below −2).
4. Anaemia: HW56 below 110 (haemoglobin below 11.0 g/dl, altitude-adjusted).

Secondary: wasting (HW72 below −200); bloody diarrhoea (H11B = 1); among children with diarrhoea,
antibiotic for non-bloody diarrhoea; among children with fever or cough, an antibiotic.

## 6. Covariates

- Child: age in months (categories 0–5, 6–11, 12–23, 24–35, 36–47, 48–59), sex, birth order (1, 2, 3, 4+).
- Mother: education (none, primary, secondary, higher), age group.
- Household:
  - religion; caste group (scheduled caste, scheduled tribe, other backward class, other)
  - number of children under five
  - drinking water source (improved or not)
  - toilet facility (improved, unimproved, none)
  - cooking fuel (clean, solid non-dung, dung)
  - disposal of the youngest child's stools (safe or not)
  - agricultural land (HV244)
- Asset index (primary): the first principal component of household consumer goods and housing materials,
  **excluding animals, land, water, sanitation and cooking fuel**. Animals and land are left out because they
  are the exposure or close to it; water, sanitation and fuel are left out because they are adjusted for
  separately. Entered as quintiles.
- Season: month of interview.
- Rural or urban residence is absorbed by the cluster fixed effects.

## 7. Statistical analysis

**Primary model.**
- Linear probability model with cluster (primary sampling unit) fixed effects, so each comparison is between
  households in the same village or urban block.
- Survey weights (V005/1,000,000), and standard errors robust to clustering at the primary sampling unit.
- Results are percentage-point differences with 95% confidence intervals.

**Secondary model.** Survey-weighted logistic regression with state fixed effects and the same covariates,
with design-based standard errors; reported as odds ratios.

**Multiplicity.** Five exposures by four primary outcomes give 20 primary tests. Benjamini–Hochberg q
values are reported beside each P value. Findings are described as supported when q is below 0.05.

**Effect modification**, each tested by an interaction term in the primary model:
- improved toilet facility (yes or no);
- child age 6–23 months against 24–59 months.

**Secondary analysis: animal-source food (6–23 months, youngest child).** Add the child's consumption in
the past 24 hours to the growth models and report the change in the animal coefficient. Groups:
- flesh foods and eggs: V414A, V414G, V414M, V414N, V414T;
- dairy: V411, V414P, V414V.

Species-matched expectations: poultry operating through eggs and flesh foods, cattle through dairy. This is
descriptive and makes no causal mediation claim.

## 8. Sensitivity analyses (all pre-specified)

1. The survey's own wealth index (V190) in place of the rebuilt asset index.
2. Negative-control exposure: households with agricultural land and no animals, compared with households
   with neither. These share farming and income but not the animal pathway.
3. E-values for each primary association.
4. Rural households only.
5. Standard errors clustered at the district instead of the primary sampling unit.
6. Logistic regression with cluster fixed effects (conditional logit) for the binary outcomes.
7. Replication in NFHS-4 with identical definitions; pigs omitted.

## 9. Interpretation

The design is cross-sectional. Within-village comparison removes community-level confounding, but not
differences between households such as farming practice or food security. Results are reported as
associations. Antibiotic exposure is maternal report over two weeks and does not identify the drug.

## 10. Outputs

- Tables:
  - characteristics by animal ownership;
  - within-village associations for the 20 primary tests;
  - effect modification;
  - sensitivity analyses and the NFHS-4 replication.
- Figures:
  - forest plot of the primary associations;
  - State-level poultry ownership against child antibiotic exposure.

Reported following STROBE for cross-sectional studies.
