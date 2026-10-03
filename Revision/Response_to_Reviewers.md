# Response to Reviewers

*Draft. This file is assembled incrementally as each section of the Revision plan (`Reviewer_Response_Plan.md`) is completed. Sections marked **[PENDING]** have not been written yet and should not be read as representing a finished response — they are placeholders tracking what Section 12 of the plan still owes.*

For each comment we (i) concede the point where it is valid, (ii) describe what we changed, (iii) point to where it now lives in the manuscript/SI, and (iv) report the key new result.

---

## Comment 1 — The DAG/SEM conflation (causal interpretation of a jointly-estimated SEM)

**[PENDING — Section 4/9 of the Revision plan; not yet drafted.]**

---

## Comment 2 — The DAG is incomplete: fish are only sinks, and time/feedback is collapsed

We agree that the original model treated fish purely as a response variable and omitted their role as drivers of ecosystem respiration and nutrient supply, and that the coral–algae–herbivore relationship was represented as a single static arrow rather than a feedback unrolled in time. The revised, time-indexed DAG (new Fig. 3; `Reviewer_Response_Plan.md` Section 4) now includes `Benthos → Herb`, `Benthos → Corall`, `Herb/Corall/OtherFish → Rd`, and `Herb/Corall/OtherFish → N`, together with lagged `Benthos_lag` and `Herb_lag` nodes that unroll the feedback. Two of the twelve resulting causal estimands concern fish specifically, and are reported here.

### E11: does benthic composition drive herbivore and corallivore biomass?

*Estimand:* total effect of `Benthos → Herb` / `Benthos → Corall`. *Adjustment set* (backdoor criterion, Table S-DAG): `{Time, Benthos_lag, Herb_lag}`.

We modelled herbivore and corallivore biomass at the fish-transect level (4 transects × 6 backreef sites × 19 years, 2007–2025; n = 456 transect-years), regressing each on the contemporaneous benthic composition — expressed as three ilr (isometric log-ratio) coordinates of the Coral/Algae/CCA/Other composition, which respects the closure constraint that a proportional composition is not a set of independent covariates (Comment 3) — the site's herbivore biomass the previous year (`herb_lag_z`, standing in for `Herb_lag` and, via the already-fitted Benthos–Herb feedback, blocking the backdoor path through `Benthos_lag`), and centred year (`Year_c`, instrumenting `Time`). Herbivore biomass had no zero transect-years and was modelled with a Gamma likelihood (log link); corallivore biomass had a modest fraction of zeros (3.5%) and was modelled with a hurdle-Gamma likelihood. Both models included `(1 | Site) + (1 | Site:Year)` random intercepts and converged cleanly (R-hat = 1.00 for every parameter; one divergent transition out of 4,000 draws in each fit, below any threshold of concern).

| Outcome | ilr1 | ilr2 | ilr3 | Herb_lag (`herb_lag_z`) | Time (`Year_c`) |
|---|---|---|---|---|---|
| Herbivore biomass (log scale) | −0.087 [−0.217, 0.033] | **0.112 [0.045, 0.178]** | 0.045 [−0.082, 0.162] | **0.198 [0.108, 0.285]** | **0.054 [0.035, 0.073]** |
| Corallivore biomass (log scale) | **−0.211 [−0.418, −0.014]** | 0.043 [−0.066, 0.146] | 0.170 [−0.023, 0.365] | −0.047 [−0.169, 0.072] | 0.008 [−0.019, 0.035] |

(Posterior mean [95% credible interval]; bold = interval excludes zero.)

Because the `ilr` coefficients are coordinates of an abstract sequential-binary-partition basis rather than effects of any named benthic part, we back-transformed them into average marginal effects (AMEs) on herbivore/corallivore biomass of a +1-percentage-point increase in one part's share (Coral, Algae, or CCA), with the other three parts rescaled proportionally to close the simplex — the same g-computation approach used for the benthic Gompertz model's E1–E3 estimands.

| Benthic part (+1 pp share) | AME on herbivore biomass (g/m²) | AME on corallivore biomass (g/m²) |
|---|---|---|
| Coral | +0.006 [−0.21, 0.21] | **+0.015 [0.00, 0.03] (+1.2%)** |
| Algae | −0.090 [−0.19, 0.005] | **−0.012 [−0.02, 0.00] (−0.9%)** |
| CCA | **+1.92 [0.82, 3.02] (+11.7%)** | −0.024 [−0.13, 0.08] |

Only two of the six part×outcome AMEs have a 95% interval that excludes zero: herbivore biomass rises with CCA share, and corallivore biomass rises with coral share — both ecologically plausible (CCA as grazing substrate; coral as corallivore food) — alongside the autoregressive effect of last year's own herbivore biomass and a residual secular increase in herbivore biomass over the study period. We note explicitly that a +1-percentage-point step is a small relative perturbation for Coral (mean share 20%) and Algae (mean share 55%) but close to a 50% relative increase for CCA (mean share only 2.2%, maximum 17%), so the CCA estimand should be read as a larger relative manipulation than the other two, not treated as an equivalent "unit" effect. A residual-autocorrelation check (lag-1 ACF per site per model, Ljung-Box, Benjamini-Hochberg corrected across all 12 series) found no site/model combination significant after correction.

### E12: do fish contribute to ecosystem nutrient supply?

*Estimand:* total effect of fish (herbivores + corallivores + other fish, combined) on macroalgal tissue nitrogen. *Adjustment set* (Table S-DAG, "Benthos → N" DAG variant): `{Benthos, Time}`.

Tissue %N in *Turbinaria ornata* is measured at LTER_1 only, giving an annual panel of n = 18 years (2007–2024) — far too few points to support the richer adjustment sets available under the DAG variant that keeps `Rd → N` (those require day/season-level Rd, DayTemp, Flow, which do not exist at this annual, single-site resolution). We therefore used the simpler DAG variant, for which `{Benthos, Time}` is the unique minimal adjustment set, and fit `log(N_percent) ~ fish_Nexcretion_z + ilr1 + ilr2 + ilr3 + Year_c` (lognormal), with `fish_Nexcretion_z` a standardised, stoichiometrically derived estimate of community fish nitrogen excretion (from the same bioenergetics pipeline used for the fish-respiration estimand below). A Gamma(log link) fit on the untransformed outcome gave essentially identical coefficients, so the result is not an artefact of the likelihood choice.

**We flag, rather than paper over, a collinearity problem specific to this estimand.** At n = 18, the fish nitrogen-excretion exposure correlates r = 0.80 with the `ilr1` (coral-dominance) coordinate of Benthos and r = 0.68 with centred year, and `ilr1` itself correlates r = 0.91 with year — three variables that are supposed to be entered as separate terms are, at this sample size, close to interchangeable proxies for the same secular trajectory. The fitted model reflects this: `fish_Nexcretion_z` has a small coefficient whose 95% interval includes zero (0.047 [−0.048, 0.140]), `Year_c` has a coefficient whose interval excludes zero (−0.048 [−0.086, −0.010], i.e. tissue N declining over the study period after nominal adjustment), and none of the three `ilr` coefficients is distinguishable from zero. Given the r = 0.91 correlation between `Year_c` and `ilr1`, we read the `Year_c` coefficient as "covaries with time" rather than as a secular effect cleanly separated from Benthos — and we are not asserting a fish → N effect, detectable or otherwise, from this fit. **We report E12 as a small-sample, collinearity-limited estimate (n = 18, single site, 5 coefficients), not as a precise or causally isolated effect, and will present it in the SI with this caveat attached rather than as a headline result in the main text.**

### Fish respiration: a direct, independent check on "fish contribute to ecosystem respiration"

**[PENDING summary text — the underlying bioenergetics pipeline (Task 7.1, `05_fish.R`) and the year-by-year comparison against measured Rd at LTER_1 are built and produce `fish_vs_measured_Rd.csv`/`fig_fish_pct_of_Rd.png`; the headline percentage-of-Rd number and its interpretation still need to be pulled into this letter.]**

### Decline vs. recovery, separated in time

**[PENDING — Section 5's compositional Gompertz model, E1–E3, is fitted (`04_benthic_gompertz.R`) and will be summarised here alongside E11/E12 once this section is finalised.]**

---

## Comment 3 — Benthic composition should use a compositional likelihood

**[PENDING summary text — the multinomial-on-point-counts model (`04_benthic_gompertz.R`) and the ilr-coordinate treatment of composition as a covariate (used above for E11 and in `06_metabolism.R`) are both built; this section needs a consolidated write-up.]**

---

## Comment 4 — Twenty annual points, with autocorrelation ignored

**[PENDING — partially addressed throughout (transect/quadrat/hourly disaggregation; residual ACF checks reported per model, e.g. E11 above); needs a consolidated statement with the sample-size transparency table (plan Task 6.4).]**

---

## Comment 5 — Residuals used as predictors (`yearresid`)

We agree that regressing calendar year on temperature and using the residual as a predictor (`yearresid <- resid(lm(Year ~ Max_temp, ...))`) is the practice Freckleton (2002) warns against, and we have removed it. None of the Revision pipeline's models use `yearresid`; wherever the causal DAG identifies calendar year as a confounder proxy for unmeasured secular drivers, centred year (`Year_c`) is entered directly alongside the exposure of interest — in the benthic compositional Gompertz model (E1–E3), the hourly photosynthesis–irradiance model, and the E11/E12 fish models described above.

Because `Year_c` and the heat-stress covariate (`DHW_max`) both change over the study period, we checked whether they are too collinear for the heat-stress effect to be separable (Task 8.2), rather than assume it. At the raw-data level (126 site-years, all 6 sites), Year and `DHW_max` correlate r = 0.374 (95% CI 0.21–0.52) — a real but moderate association, reflecting that detectable heat stress is concentrated in specific recent event years (2019 onward) rather than trending smoothly across 2005–2025. In the one fitted model carrying both covariates simultaneously (the additive benthic Gompertz model), the *posterior* coefficient correlations between `DHW_max` and `Year_c` are modest and in fact of the opposite sign (−0.24 to −0.25) from the raw correlation, for all three benthic categories (Coral, Algae, CCA) — well short of the near-unidentifiability we found for the E12 fish→N model's exposure and adjustment covariates (Comment 2 above, r up to 0.91 at n = 18). We conclude Year and DHW_max are not so collinear that the heat-stress effect is unidentifiable in this model, but we state plainly that it is identified mainly from the event years rather than from smooth variation, which is an honest limitation of an 18–20-year series with a handful of thermal-stress years, not a precise estimate.

---

## Projections

**[PENDING — Section 11 of the plan; not started.]**
