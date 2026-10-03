# QC summary: compositional Gompertz model (Revision Step 4 / Section 5)

Generated: 2026-10-01 17:09:21.912214

## A. ALR lag predictors

- Zero-replaced (`zCompositions::cmultRepl`) counts of Coral/Algae/CCA/Other at site-year resolution, then `alr_k = log(part_k / Other)`. Lagged by +1 year (same site) to get `alr_*_lag`, the previous year's composition as this year's recovery-phase predictor.

## B. Herbivore biomass lag

- `Herb` = Herbivore biomass (g/m^2), mean across the 4 fish transects per site-year. Lagged +1 year, z-standardised using the full-sample mean (16.05) and sd (6.97) -> `herb_lag_z`.

## D. Model data

- `benthic_transect_model.csv`: 570 transect-year rows (6 sites, 19 years: 2007-2025).
- Rows dropped for missing lag (first observed year per site, 2006): 30.
- `COTS` rescaled from ind/m^2 (max 0.012) to ind/100m^2 (max ~1.2) so its coefficient is on a comparable scale to the other predictors; the raw scale made its `normal(0,1)` prior essentially uninformative.

## F. Derived estimands E1-E3 (average marginal effect on coral proportion)

Computed via manual g-computation: `posterior_epred()` at observed covariates vs. a +1-unit counterfactual for each predictor in turn (all else held at observed values), averaged across all 570 rows and summarised over the 4000 posterior draws. Units: percentage points of coral cover share.

| Estimand | AME (mean) | 95% CI | P(AME > 0) |
|---|---|---|---|
| E1: DHW -> Coral proportion (total), AME per +1 DHW degC-week | 0.0015 | [-0.0153, 0.0195] | 0.563 |
| E2: COTS -> Coral proportion (total), AME per +1 ind/100m^2 | 0.0265 | [-0.0618, 0.1328] | 0.686 |
| E3: Herb_lag -> Coral proportion (total), AME per +1 SD herbivore biomass | 0.0037 | [-0.0162, 0.0242] | 0.634 |

- Joint counterfactual (DHW_max=0 & COTS=0 for all rows, vs. observed): mean change in coral proportion = -0.0019, 95% CI [-0.0103, 0.0068].

## G. Task 5.2: interaction (recovery-mediation) model vs. additive

- Added `herb_lag_z:alr_coral_lag` and `herb_lag_z:DHW_max` to the same formula used for the additive model (Section E) -- brms applies one linear-predictor formula to every non-reference multinomial category, so both interactions are estimated for `muCoral`, `muAlgae` and `muCCA` alike; only the `muCoral` versions are the Task 5.2 estimand.

| Coefficient | Estimate | 95% CI |
|---|---|---|

### loo() comparison

```
                   elpd_diff se_diff
fit_benth_additive   0.0       0.0  
fit_benth_interact -26.2      13.2  
```

## I. Task 5.3: residual temporal autocorrelation check

- Two diagnostics at site-year resolution (6 sites x 19 years each): the `Site:Year` random intercepts for muCoral/muAlgae/muCCA, and Pearson residuals (observed vs. fitted count, binomial-approximation standardisation) averaged across the 5 transects per site-year. ACF computed per site (lags 1-3); lag-1 tested via Ljung-Box per site per series (36 tests total), BH-corrected.
- Significant at raw p<0.05: 2 of 36 tests (chance expectation ~1.8). **Significant at BH-adjusted p<0.05: 0 of 36.**
- Two sites (LTER_3, LTER_6) show a consistently negative lag-1 ACF (-0.3 to -0.48) across all 6 series (both random effects and residuals, both correlated since they come from the same site-years) -- a repeated pattern worth noting, but it does not survive multiple-comparison correction (0 of 36 BH-adjusted tests reach p<0.05).
- **Conclusion: no statistically defensible residual temporal autocorrelation once the lagged ALR terms are in the model, matching Task 5.3's expectation.** The `ar()`/`gp()` sensitivity checks the plan offers as a fallback are not triggered by this result. Saved as `fig_acf_benthic_gompertz.png` and `acf_lag1_tests.csv`.

