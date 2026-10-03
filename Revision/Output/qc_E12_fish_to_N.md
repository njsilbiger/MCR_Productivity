# QC summary: Task 7.3 fish -> N model (E12), LTER_1 Turbinaria tissue N

Generated: 2026-10-02 15:28:21.215161

Estimand E12 (Herb+Corall+OtherFish -> N, total effect). Adjustment set {Benthos, Time}, taken from the "Benthos -> N" DAG variant in `table_s_dag.csv` (its E12 row has a single minimal adjustment set, vs. 4 alternative sets under the "Rd -> N" variant, which require Rd/DayTemp/Flow/Season -- not available at this annual, LTER_1-only resolution). `ilr1-3` = contemporaneous Benthos; `Year_c` = Time.

## G1. Model data

- `e12_data`: 18 annual rows, LTER_1 only (2007-2024). `N_percent_mean` is the mean of `n_samples` Turbinaria tissue-N samples per year (range 5-10 samples/year).
- `fish_Nexcretion_z`: LTER_1's annual stoichiometric fish N-excretion estimate (Task 7.1's placeholder O:N = 20 conversion), z-standardised.

## G2. Collinearity among exposure and adjustment covariates

Pairwise Pearson correlations among `fish_Nexcretion_z`, `ilr1-3` (Benthos) and `Year_c` (Time), at n = 18 annual points:

| | fish_Nexcretion_z | ilr1 | ilr2 | ilr3 | Year_c |
|---|---|---|---|---|---|
| fish_Nexcretion_z | 1.00 | 0.80 | -0.11 | 0.53 | 0.68 |
| ilr1 | 0.80 | 1.00 | -0.27 | 0.71 | 0.91 |
| ilr2 | -0.11 | -0.27 | 1.00 | -0.35 | -0.55 |
| ilr3 | 0.53 | 0.71 | -0.35 | 1.00 | 0.63 |
| Year_c | 0.68 | 0.91 | -0.55 | 0.63 | 1.00 |

**`fish_Nexcretion_z` correlates r = 0.80 with `ilr1` and r = 0.68 with `Year_c`, and `ilr1` correlates r = 0.91 with `Year_c`** -- at n = 18, the exposure and two of its adjustment covariates are nearly collinear. The model below is fit as specified, but the individual coefficients (especially `fish_Nexcretion_z` vs. `ilr1`/`Year_c`) should be read as weakly identified from each other, not as cleanly separated effects.

## G3. Model results

- Primary: `log(N_percent_mean) ~ fish_Nexcretion_z + ilr1 + ilr2 + ilr3 + Year_c`, Gaussian likelihood on the log scale (= lognormal on `N_percent_mean`), n = 18.
- Robustness check: the same linear predictor with `Gamma(link = "log")` directly on `N_percent_mean` (not logged).

| Coefficient | Lognormal estimate | 95% CI | Gamma estimate | 95% CI |
|---|---|---|---|---|
| Intercept | -0.522 | [-1.174, 0.125] | -0.529 | [-1.177, 0.123] |
| fish_Nexcretion_z | 0.047 | [-0.048, 0.140] | 0.045 | [-0.047, 0.137] |
| ilr1 | 0.054 | [-0.334, 0.453] | 0.062 | [-0.333, 0.451] |
| ilr2 | -0.029 | [-0.108, 0.049] | -0.030 | [-0.108, 0.050] |
| ilr3 | -0.023 | [-0.155, 0.103] | -0.017 | [-0.149, 0.115] |
| Year_c | -0.048 | [-0.086, -0.010] | -0.049 | [-0.087, -0.010] |

- The two likelihoods agree closely on every coefficient -- the result below is not an artefact of the lognormal/Gamma choice.
- Posterior coefficient correlations (lognormal fit): b_ilr1 vs. b_Year_c = -0.86, b_fish_Nexcretion_z vs. b_ilr1 = -0.41, b_fish_Nexcretion_z vs. b_Year_c = 0.05 -- confirms the raw-data collinearity (G2) propagates into the posterior, widening and correlating the coefficient estimates rather than being resolved by the weakly informative priors.

**Result: `fish_Nexcretion_z` has a small, highly uncertain positive coefficient (lognormal: 0.047, 95% CI [-0.048, 0.140]) that does not exclude zero. `Year_c` has a negative coefficient whose 95% CI excludes zero (-0.048, 95% CI [-0.086, -0.010]), i.e. tissue N declines over the study period after adjusting for contemporaneous Benthos and fish N excretion -- but because `Year_c` and the Benthos `ilr1` coordinate are themselves correlated r = 0.91 at n = 18, this should be read as a covariate-of-time effect more than a cleanly isolated secular effect. No `ilr` coefficient's CI excludes zero.**

**Honesty note (plan Task 7.3): this is an n = 18 annual-point, single-site estimate with 5 regression coefficients and substantial collinearity among them (G2) -- it should be reported as suggestive at best, not as a precise or clearly causally separated estimate of the fish -> N effect.**

