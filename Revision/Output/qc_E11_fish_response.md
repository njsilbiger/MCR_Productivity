# QC summary: Task 7.2 fish response models (E11)

Generated: 2026-10-02 15:28:15.542604

Estimand E11 (Benthos -> Herbivores / Corallivores, total effect), adjustment set {Time, Benthos_lag, Herb_lag} per Table S-DAG. Fit at transect level across all 6 backreef sites (LTER_1-6).

## A. Contemporaneous Benthos (ilr coordinates, site-year)

- Zero cells in the Coral/Algae/CCA/Other count matrix (120 site-years x 4 parts): 19. Replaced via `zCompositions::cmultRepl` before closing to proportions and taking ilr (`compositions::ilr()`, default sequential binary partition), matching 03b_dsep_test.R's Benthos-node construction.
- Same-YEAR (not lagged) composition, because E11 is the Benthos -> Herb/Corall total effect in the same year, not a lagged recovery effect.

## B. Herb_lag (adjustment variable)

- Site-mean Herbivore biomass (g/m^2, across the 4 fish transects), lagged +1 year, z-standardised (mean 16.05, sd 6.97). Identical construction to 04_benthic_gompertz.R.

## C. Model data

- `fish_resp_data`: 456 transect-year rows (6 sites, 19 years: 2007-2025).
- Rows dropped for missing `herb_lag_z` (first observed year per site, 2006) or unmatched benthic survey: 24.
- Outcome zero counts: herbivore_biomass = 0 of 456; corallivore_biomass = 17 of 456.

## D. Model results

- `herbivore_biomass`: `Gamma(link = "log")` (no zeros at transect level).
- `corallivore_biomass`: `hurdle_gamma(link = "log")` (17 zero-biomass transect-years, 3.7%).

### E11: Benthos -> Herbivore biomass

| Coefficient | Estimate | 95% CI |
|---|---|---|
| Intercept | 3.094 | [2.828, 3.354] |
| ilr1 | -0.087 | [-0.217, 0.033] |
| ilr2 | 0.112 | [0.045, 0.178] |
| ilr3 | 0.045 | [-0.082, 0.162] |
| herb_lag_z | 0.198 | [0.108, 0.285] |
| Year_c | 0.054 | [0.035, 0.073] |

### E11: Benthos -> Corallivore biomass

| Coefficient | Estimate | 95% CI |
|---|---|---|
| Intercept | 0.438 | [0.005, 0.868] |
| ilr1 | -0.211 | [-0.418, -0.014] |
| ilr2 | 0.043 | [-0.066, 0.146] |
| ilr3 | 0.170 | [-0.023, 0.365] |
| herb_lag_z | -0.047 | [-0.169, 0.072] |
| Year_c | 0.008 | [-0.019, 0.035] |

## E. Residual ACF per site

- Standardised response residuals, averaged to site-year then lag-1 Ljung-Box per site per model (12 site x model series tested; BH-corrected across all series).
- Significant at raw p<0.05: 1 of 12. Significant at BH-adjusted p<0.05: 0 of 12.

## F. Back-transformed marginal effects on Coral/Algae/CCA share

- ilr1-3 are coordinates of an abstract basis, not named-part effects. Re-expressed as the average marginal effect (AME) of a +1 percentage-point increase in one part's share (Coral, Algae or CCA), with the other three parts rescaled proportionally to close the simplex, on predicted herbivore/corallivore biomass -- same g-computation approach as the E1-E3 estimands in 04_benthic_gompertz.R.
- Caveat: observed mean share is Coral 20.4%, Algae 55.2%, CCA 2.2% (max 17.4%), so a +1-percentage-point step is a much larger RELATIVE perturbation for CCA than for Coral or Algae.

| Outcome | Part (+1 pp) | AME (mean) | 95% CI | P(AME > 0) | % change vs. mean prediction |
|---|---|---|---|---|---|
| Herbivore biomass (g/m2) | Coral | 0.0062 | [-0.2069, 0.2113] | 0.520 | 0.04% |
| Corallivore biomass (g/m2) | Coral | 0.0147 | [-0.0024, 0.0316] | 0.958 | 1.16% |
| Herbivore biomass (g/m2) | Algae | -0.0903 | [-0.1928, 0.0048] | 0.030 | -0.55% |
| Corallivore biomass (g/m2) | Algae | -0.0119 | [-0.0244, -0.0003] | 0.023 | -0.94% |
| Herbivore biomass (g/m2) | CCA | 1.9246 | [0.8234, 3.0244] | 1.000 | 11.69% |
| Corallivore biomass (g/m2) | CCA | -0.0238 | [-0.1278, 0.0803] | 0.327 | -1.87% |

