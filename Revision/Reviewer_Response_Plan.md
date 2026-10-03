# Revision plan: response to the major statistical review

This plan sets out how to rebuild the causal and statistical analysis in `MCR_Productivity_Analysis.qmd` to address all five major comments. It is written so another model or an analyst can carry out each step. Every task lists its inputs, the exact modelling choices, its outputs and an acceptance check.

---

## 0. Strategy and framing

### 0.1 What the reviewer is asking for

| # | Comment | Root problem in current code | Core fix |
|---|---------|------------------------------|----------|
| 1 | DAG ≠ SEM; parameters estimated jointly are not causal | `brms_sem_full` fits all 8 equations at once, and `fig-sem-dag` presents those coefficients as causal effects | Use the Structural Causal Model (SCM) workflow: one estimand gives one adjustment set (backdoor criterion), which gives one model |
| 2 | DAG incomplete: fish→respiration is missing, and time/feedback is ignored | Fish only appear as sinks, and the coral↔algae feedback is collapsed into a single static arrow | Build a new time-indexed DAG that includes fish→Rd, fish→N, disturbances and lagged states; test it against the data |
| 3 | Benthos is compositional | log-z-scored coral and algae each get their own Gaussian likelihood | Use a compositional likelihood (multinomial on point counts, or Dirichlet/logistic-normal), with compositional covariates expressed as ilr/alr coordinates |
| 4 | 20 annual points; autocorrelation ignored | Everything goes through `Year_Averages`, and AR(1) is applied to coral only | Hierarchical models on the raw replication (quadrat/transect/site × year; day/hour for metabolism) plus explicit population dynamics (a Gompertz state model) |
| 5 | Residuals used as predictors | `yearresid <- resid(lm(Year ~ Max_temp))` | Drop it. Put the time-trend proxy into the model directly as a covariate when the DAG requires it (Freckleton 2002) |

### 0.2 Literature to anchor the revision

These are methodological precedents that fit the reviewer's comments. Check every citation against the actual paper before using it.

- **SCM / DAG-based estimation (Comment 1):**
  - Arif & MacNeil (2023) *Ecol. Monogr.* 93:e1554, "Applying the structural causal model framework for observational causal inference in ecology". This is the main template: one DAG, then the backdoor criterion, then a separate model for each causal query.
  - Arif & MacNeil (2022) *Ecol. Lett.* 25:1741–1745, "Predictive models aren't for causal inference". Use it to explain why the reduced SEM, which dropped paths whose CI crossed zero, is not valid causal practice.
  - Arif, Graham, Wilson & MacNeil (2022) *Ecosphere* 13:e3956, "Causal drivers of climate-mediated coral reef regime shifts". This is the closest coral-reef example of a DAG followed by per-effect models. Copy its layout: a DAG figure, a table of adjustment sets, then one model per effect.
  - Arif & MacNeil (2022) *Ecosphere* 13:e4009, "Utilizing causal diagrams across quasi-experimental approaches".
  - Byrnes & Dee (2025) *Ecol. Lett.*, "Causal inference with observational data and unobserved confounding variables". Covers fixed-effects and panel designs for unobserved confounders, and becomes relevant once the data are expanded to 6 sites.
  - Textor et al. (2016) *Int. J. Epidemiol.* (the `dagitty` package); Ankan, Wortel & Textor (2021) *Curr. Protoc.* (testing DAGs against data).
  - Shipley (2000, 2009). The d-separation test of DAG-implied conditional independences.
- **Population dynamics / disturbance–recovery (Comments 2 and 4):**
  - MacNeil et al. (2019) *Nat. Ecol. Evol.* 3:620–627, "Water quality mediates resilience on the Great Barrier Reef". The code is at https://github.com/mamacneil/GBR_Gompertz; `GBR_coral_model.ipynb` holds the model and `GBR_future_simulations.ipynb` the projections. It uses a hierarchical Bayesian Gompertz model of coral cover, in which disturbances cause losses and a covariate (water quality) modifies the intrinsic recovery rate and disturbance susceptibility. Projections are made by forward simulation of disturbance regimes. **Read the notebook and mirror its parameterisation.** The local analogue here is that herbivory and nutrients mediate coral recovery at Moorea.
  - Dennis et al. (2006) *Ecol. Monogr.* 76:323–341 (Gompertz state-space model with process and observation error).
  - Ives et al. (2003) *Ecol. Monogr.* 73:301–330 (MAR(1) multivariate Gompertz for interacting species).
- **Compositional data (Comment 3):**
  - Aitchison (1982) *JRSS-B* 44:139–177; Aitchison (1986) *The Statistical Analysis of Compositional Data*.
  - Egozcue et al. (2003) *Math. Geol.* (ilr transform).
  - van den Boogaart & Tolosana-Delgado (2013) *Analyzing Compositional Data with R* (`compositions`).
  - Pawlowsky-Glahn, Egozcue & Tolosana-Delgado (2015) *Modeling and Analysis of Compositional Data*.
  - Hron, Filzmoser & Thompson (2012) *J. Appl. Stat.* (compositional covariates in regression, pivot coordinates).
  - Palarea-Albaladejo & Martín-Fernández (2015) (`zCompositions`, zero replacement).
  - Douma & Weedon (2019) *Methods Ecol. Evol.* 10:1412–1430 (Dirichlet regression in ecology).
- **Residuals (Comment 5):** Freckleton (2002) *J. Anim. Ecol.* 71:542–545.
- **Fish metabolism/excretion (Comment 2):**
  - Barneche et al. (2014) *Ecol. Lett.* (scaling individual metabolism to reef-fish communities).
  - Schiettekatte et al. (2020) *Funct. Ecol.* (`fishflux` R package: individual respiration, consumption and N excretion from size, species and temperature).
  - Allgeier et al. (2014) *Glob. Change Biol.* (fish-mediated nutrient supply).
- **Moorea disturbance history:**
  - Kayal et al. (2012) *PLoS ONE* (the 2006–2010 *Acanthaster* outbreak).
  - Adam et al. (2011) *PLoS ONE* (herbivore response to that perturbation).
  - Adjeroud et al. (2009) *Coral Reefs*.
  - Holbrook et al. (2018) *Sci. Rep.* (herbivory and resilience).
  - Literature on Cyclone Oli (Feb 2010) and the 2019 Moorea bleaching event.
- **Robustness:**
  - Cinelli & Hazlett (2020) *JRSS-B* (`sensemakr`, sensitivity to omitted confounders).
  - Kallioinen et al. (2024) *Stat. Comput.* (`priorsense`, prior sensitivity).
  - Donner et al. (2005) *Glob. Change Biol.* (degree heating months from monthly SST, needed for the CMIP6 projections).

### 0.3 Things to be honest about in the response

- Metabolism is measured at **one site** on about 150 complete diel cycles across 18 years. Season coverage is uneven: 2010, 2020 and 2025 are summer-only, and 2008, 2011, 2013 and 2021 are winter-only. Benthic-to-metabolism effects therefore rest on between-year variation at one site. Temperature-to-metabolism effects can also use **within-year, day-to-day** temperature variation, which gives more leverage.
- The backreef temperature file contains **one logger, at LTER_2 (2 m)**, which is used for LTER_1. Heat stress therefore varies only across years, not across sites. With n ≈ 20 years and a few event years (2019 especially), the heat effect on coral will have wide uncertainty. Say so.
- The future projections in the current paper (CMIP6 SSPs to 2100) extrapolate far beyond the data. The plan rebuilds them as **conditional scenario simulations** (MacNeil 2019 style) and recommends moving them to the SI or toning them down.

### 0.4 Response-letter etiquette

Reviews are anonymous. Following one researcher's methodological approach is fine, but do not name or guess the reviewer in the response letter or the cover letter. Address the comments on their merits and cite the methods literature above.

---

## 1. Project set-up

**Task 1.1. Folder and file layout.** Leave the original qmd untouched so the old analysis can be reproduced. Create:

```
Revision/
  Reviewer_Response_Plan.md        (this file)
  R/
    00_packages.R
    01_data_disaggregated.R        # builds quadrat/transect/day-level tables
    02_disturbance_covariates.R    # DHW, COTS, cyclone
    03_dag.R                       # DAG, adjustment sets, implied CI tests
    04_benthic_composition.R       # compositional dynamics (Gompertz/MAR)
    05_fish.R                      # fish models + fish respiration/excretion
    06_metabolism.R                # PI-curve model with covariates on Pmax/Rd
    07_estimands.R                 # one model per causal estimand
    08_simulation_check.R          # DAG simulation: old SEM vs new approach
    09_projections.R               # forward scenario simulation
    10_sensitivity.R
  Output/                          # all new figures/tables/model objects (.rds)
  Response_to_Reviewers.md
  Methods_revised.md
```

**Task 1.2. Packages** (`00_packages.R`):

```r
library(tidyverse); library(here); library(lubridate)
library(brms); library(cmdstanr)       # use backend = "cmdstanr" for speed
library(posterior); library(tidybayes); library(loo); library(priorsense)
library(dagitty); library(ggdag)
library(compositions); library(zCompositions)
library(fishflux)                      # remotes::install_github("nschiett/fishflux") if not on CRAN
library(rerddap)                       # NOAA Coral Reef Watch DHW
library(sensemakr)
options(mc.cores = 4, brms.backend = "cmdstanr")
```

Seeds: use `set.seed(4817)` and pass `seed = 4817` to every `brm()` call.

---

## 2. Data disaggregation (Comment 4)

Stop working from `Year_Averages` (n = 20). Build the following analysis tables from the raw files in `Data/raw_data/`. The checks during planning established:

- Benthic file: 6 backreef sites (LTER_1–6) × ~20 years × 5 transects × 10 quadrats. `Percent_Cover` values are integers, about 99% are multiples of 4, and quadrat totals are usually 100.
- Fish file: 6 backreef sites × 4 transects × 2 swaths × 20 years, with individual-level `Total_Length`, `Count`, `Biomass`, `Taxonomy`, `Fine_Trophic`.
- Metabolism (`QC_PP`): hourly data, about 150 complete diel days across 18 years, with `Season`, `PAR`, `Temperature_mean`, `Flow_mean` available per hour.

**Task 2.1. Benthic composition at quadrat level** → `benthic_quad`.
- Keep the existing functional-group mapping, but use **four parts**: `Coral`, `Fleshy Macroalgae and Turf` (call it `Algae`), `Crustose Corallines` (`CCA`) and `Other/Sand`. Put everything not in the first three into `Other`.
- **QC first:**
  1. Confirm the number of points per quadrat from the EDI metadata for the MCR benthic dataset. The multiples of 4 suggest 25 points per quadrat, but check this.
  2. Tabulate the number of quadrats per site-year. The current code hard-codes `/5000`, which assumes 50 quadrats, so confirm that holds every year.
  3. Flag quadrats whose total ≠ 100.
- If points per quadrat = N_pts, create integer counts `n_k = round(Percent_Cover_k / 100 * N_pts)` and `n_total = sum(n_k)`. Keep `Site, Year, Transect, Quadrat`.
- Also build `benthic_site` (site × year proportions) and `benthic_transect` (transect × site × year) for plotting and lag construction.

**Task 2.2. Fish at transect level** → `fish_transect`.
- Use all 6 backreef sites. Keep the same shark (> 8000 g) and negative-biomass filters, but report how many rows each filter removes.
- Groups: `Herbivore` (current definition), `Corallivore`, `Planktivore`, `Invertivore`, `Piscivore`, `Other`. Compute biomass in g m⁻² per transect, adjusting the denominator for swath width (1 m vs 5 m swaths are surveyed for different size classes). Check the MCR fish protocol and do not just divide by 1200.
- Keep the individual-level table (`fish_ind`) for Task 7.1.

**Task 2.3. Metabolism at day and hour level** → `pp_hour`, `pp_day`.
- Start from `All_PP_data` after the existing QC (unit conversion, removal of bad dates, 24-hour complete days).
- `pp_day` (one row per `DielDate`) contains: `Year`, `Season`, mean daily water temperature, mean flow, total daily PAR, night-time R (mean PP when PAR = 0), daytime NP, and the UP/DN deployment ID (`UPDN`).
- Join the **benthic composition of LTER_1 in the same year** and **fish biomass of LTER_1 in the same year**. Task 6.3 handles the uncertainty this join introduces.

**Acceptance check for Section 2:** a short QC table in `Revision/Output/qc_disaggregation.md` listing row counts, number of site-years and the number of days per year × season.

---

## 3. Disturbance and heat-stress covariates (Comments 2 and 4)

**Task 3.1. Heat stress.** Replace "annual max daily temperature" with a bleaching-relevant metric.
- Primary: NOAA Coral Reef Watch 5 km **Degree Heating Weeks** (daily, 1985–present) for the Moorea pixel(s), downloaded via `rerddap` (CoastWatch ERDDAP, dataset "NOAA_DHW"; verify the dataset ID). Compute annual max DHW for the 12 months before each benthic survey, since surveys happen around April–May. The window should end at the survey date, not the calendar year.
- Secondary (sensitivity): the in situ logger record (LTER_2 backreef, plus the forereef loggers in the same file). Compute in situ DHW using the CRW MMM climatology for the pixel.
- For projections (Section 11), compute **degree heating months** from the bias-corrected CMIP6 monthly `tos` that already exists in the qmd (Donner et al. 2005).

**Task 3.2. Crown-of-thorns starfish.** Search EDI for the MCR LTER corallivorous invertebrate / *Acanthaster* survey dataset. If it exists, compute COTS density per site-year. Otherwise use a literature-based outbreak indicator (2006–2010; Kayal et al. 2012), coded as site-year intensity where sources allow.

**Task 3.3. Cyclone.** Use an indicator for Cyclone Oli (Feb 2010), which affected the 2010 survey. Check whether it is backreef-relevant at each site, because backreef impacts were much smaller than forereef. Optionally add a wave-energy proxy.

**Task 3.4. Secular-trend proxy.** Keep calendar `Year` as an observed proxy for unmeasured secular drivers (acidification, fishing, observer changes). It enters models **directly as a covariate**, never as a residual (Section 8).

Output: `disturbance_site_year.csv` with `Site, Year, DHW_max, COTS, Cyclone`.

---

## 4. A revised, time-indexed DAG (Comments 1 and 2)

**Task 4.1. Write the DAG in `dagitty`** (`03_dag.R`). The benthos is a **single compositional node** (`Benthos`), so the closure constraint is not drawn as fake causal arrows such as "Coral → Algae". The coral↔algae feedback is unrolled in time, which removes the cycle. Fish now affect respiration and nutrients.

```r
dag_mcr <- dagitty('dag {
  Time      [pos="0,0"]
  DHW       [exposure, pos="1,0"]
  COTS      [pos="1,1"]
  Cyclone   [pos="1,2"]
  Benthos_lag [pos="0,3"]
  Herb_lag  [pos="0,4"]
  N_lag     [pos="0,5"]
  Benthos   [pos="2,2"]
  Herb      [pos="3,3"]
  Corall    [pos="3,4"]
  OtherFish [pos="3,5"]
  N         [pos="4,5"]
  Season    [pos="3,0"]
  DayTemp   [pos="4,0"]
  Flow      [pos="4,1"]
  Rd        [outcome, pos="5,2"]
  Pmax      [pos="5,3"]

  Time -> DHW
  Time -> Benthos
  Time -> Herb
  Time -> OtherFish
  Time -> N
  Benthos_lag -> COTS
  DHW -> Benthos
  COTS -> Benthos
  Cyclone -> Benthos
  Benthos_lag -> Benthos
  Herb_lag -> Benthos
  N_lag -> Benthos
  Benthos_lag -> Herb
  Herb_lag -> Herb
  Benthos -> Herb
  Benthos -> Corall
  Benthos -> OtherFish
  Herb -> N
  OtherFish -> N
  Corall -> N
  DHW -> DayTemp
  Season -> DayTemp
  Season -> Flow
  Benthos -> Rd
  Herb -> Rd
  Corall -> Rd
  OtherFish -> Rd
  DayTemp -> Rd
  Flow -> Rd
  Benthos -> Pmax
  DayTemp -> Pmax
  Flow -> Pmax
  N -> Pmax
  Rd -> N
}')
```

Notes for whoever implements this:

- **`Rd → Pmax` is removed.** Respiration does not cause photosynthetic capacity. They covary because they share causes (benthic biomass, temperature). In the new framework that covariance is modelled as a residual correlation between the two parameters inside the PI model (Section 6), not as a causal path. This change is also a direct example of the DAG-vs-SEM distinction the reviewer raises.
- **`Rd → N`** is kept only if there is a mechanism to state, such as remineralisation raising tissue N. Otherwise replace it with `Benthos → N`. Decide this with the PI, write the justification into the methods, and draw the DAG both ways in the SI.
- **`Time`** is a proxy node for secular confounders. Its effects on Benthos, fish and N are exactly what `yearresid` was trying to soak up.
- Add DAG edges the authors believe in. Each *missing* arrow is a testable claim (Task 4.3).

**Task 4.2. Adjustment sets for every estimand.** Run `adjustmentSets(dag_mcr, exposure = X, outcome = Y, effect = "total")` and also `effect = "direct"` for the estimands listed below. Save the results as Table S-DAG (columns: estimand, effect type, minimal adjustment set, model in which it is estimated).

| ID | Estimand | Effect type | Expected adjustment (verify with dagitty) |
|----|----------|-------------|-------------------------------------------|
| E1 | DHW → Benthos (coral share) | total | {Time, Benthos_lag} or similar |
| E2 | COTS → Benthos | total | {Benthos_lag, Time} |
| E3 | Herb_lag → Benthos (recovery mediation) | total | {Benthos_lag, Time, ...} |
| E4 | Benthos → Rd | direct (not via fish) | {fish groups, DayTemp, Flow, Time?} |
| E5 | Benthos → Rd | total | {Time, DHW, ...} |
| E6 | Fish (all groups) → Rd | direct | {Benthos, DayTemp, Flow} |
| E7 | DayTemp → Rd (physiological) | direct | {Benthos, fish, Season / Flow} |
| E8 | DHW → Rd | total | {Time} |
| E9 | Benthos → Pmax | direct | {DayTemp, Flow, N} |
| E10 | DayTemp → Pmax | direct | {Benthos, Flow, N} |
| E11 | Benthos → Corallivores / Herbivores | total | {Time, Benthos_lag, Herb_lag} |
| E12 | Fish → N | total | {Benthos, Time} |

Important: each estimand gets **its own model with only its adjustment set** (Arif & MacNeil 2023). Do not read E5 off the model built for E4.

**Task 4.3. Test the DAG against data.** Use `impliedConditionalIndependencies(dag_mcr)` and `localTests()` (for continuous data use the `"cis"` / linear test on site-year data; for the compositional node, use ilr coordinates). Also run a Shipley d-sep test (Fisher's C) as a summary. Report any violated independences, and either add the missing arrow (with a rationale) or explain why it is not added. Put the results in Table S-dsep.

**Task 4.4. DAG figure.** Replace Figure 3 with a **pure DAG** figure that has *no coefficients on the arrows* (via `ggdag`, or keep the existing hand-drawn ggplot style but remove the β labels). Effect estimates go in a separate forest plot of estimands (new Figure 4). Retitle: Figure 3 is "Hypothesised causal structure", not "SEM".

---

## 5. Benthic composition dynamics: compositional Gompertz (Comments 2, 3 and 4)

This is the centrepiece and mirrors MacNeil et al. (2019). Coral dynamics are modelled as a **population process**: disturbances (heat, COTS, cyclone) cause losses, and recovery depends on intrinsic growth modulated by herbivory and the space occupied by algae. That separates the **decline phase (disturbance → coral → algae)** from the **recovery phase (herbivory → algae → coral)**, as the reviewer asks.

**Task 5.1. Likelihood: compositional.** Option A is preferred.

*Option A: multinomial on quadrat point counts* (handles zeros and respects closure exactly).

```r
benthic_quad$Y <- with(benthic_quad, cbind(Other, Coral, Algae, CCA))  # Other = reference
f_benth <- bf(
  Y | trials(n_total) ~ 1 + alr_coral_lag + alr_algae_lag + alr_cca_lag +
    DHW_max + COTS + Cyclone + herb_lag_z + Year_c +
    (1 | Site) + (1 | Site:Year) + (1 | Site:Year:Transect)
)
fit_benth <- brm(f_benth, family = multinomial(), data = benthic_quad,
                 prior = c(prior(normal(0, 1), class = "b", dpar = "muCoral"),
                           prior(normal(0, 1), class = "b", dpar = "muAlgae"),
                           prior(normal(0, 1), class = "b", dpar = "muCCA"),
                           prior(exponential(2), class = "sd", dpar = "muCoral"),
                           prior(exponential(2), class = "sd", dpar = "muAlgae"),
                           prior(exponential(2), class = "sd", dpar = "muCCA")),
                 chains = 4, iter = 3000, warmup = 1500, seed = 4817)
```

- The multinomial-logit linear predictors are **additive log-ratios** `log(p_k / p_Other)`. A regression of alr(composition_t) on alr(composition_{t-1}) is a **multivariate discrete-time Gompertz / MAR(1)** model (Ives et al. 2003) on compositions. The diagonal lag coefficients measure density dependence and the off-diagonals measure competition for space. State this explicitly in the methods.
- `alr_*_lag` are site-level alr values from the previous year's survey, computed from `benthic_site` after multiplicative zero replacement (`zCompositions::cmultRepl`) if any site-year part is zero.
- `(1 | Site:Year)` absorbs overdispersion and site-year process noise (a logistic-normal multinomial). `(1 | Site:Year:Transect)` handles transect clustering.
- If the multinomial is too slow, fit at **transect level** (aggregate counts over the 10 quadrats), which keeps the same likelihood with fewer rows.

*Option B (if point counts can't be reconstructed): Dirichlet* on transect-level proportions with zero replacement, `family = dirichlet()`, same formula (the brms Dirichlet also uses an alr/logit parameterisation with the first category as reference). Alternatively, `logistic_normal()` in brms.

*Option C (sensitivity): ilr coordinates + multivariate normal* on site-year compositions. This is Aitchison's original approach, and the reviewer will accept it because the transform respects the simplex. It differs from the current model, which is log-z on each part separately.

**Task 5.2. Recovery mediation, following MacNeil (2019).** Fit a second version where herbivory **modifies the recovery rate** rather than entering additively:
- Add the interaction `herb_lag_z:alr_coral_lag` (herbivory changes density dependence / recovery speed). Also add `herb_lag_z:DHW_max` (herbivory changes susceptibility). This is the analogue of water quality modifying recovery and disturbance impacts in the GBR paper.
- Compare the additive and interaction models with `loo()` (PSIS-LOO, leave-future-out preferred; see Task 10.4). Report both.

**Task 5.3. Temporal autocorrelation check.** After fitting, extract site-year random effects and Pearson residuals, compute the ACF per site and plot them. With the lagged state in the model, residual autocorrelation should be negligible. If not, add `ar(time = Year, gr = Site, p = 1)` (only available for some families). Otherwise, add a site-specific random walk via `gp(Year, by = Site, k = 5)` as a sensitivity check.

**Task 5.4. Derived quantities** (from posterior draws, `posterior_epred`):
- Expected coral proportion trajectory per site with and without each disturbance (counterfactual: set COTS = 0, DHW = baseline).
- Recovery rate: years to return to the pre-disturbance coral share, as a function of herbivore biomass.
- **Estimands E1–E3**: the average marginal effect of DHW, COTS and herbivory on coral **proportion** (not log-odds), using `marginaleffects::avg_comparisons()` or manual g-computation over posterior draws.

Note: E1–E3 come from this model only if its covariates match the dagitty adjustment set. If they don't, fit a separate model per estimand with the correct set.

---

## 6. Ecosystem metabolism: put covariates inside the PI model (Comments 1, 2 and 4)

The current workflow fits Pmax and Rd per year as fixed effects, takes the point estimates (`fixef(...)$Estimate`), and then uses them as data in the SEM. That **discards the PI-curve uncertainty**, even though the code comment claims it propagates. It also **confounds year with season**, because some years are summer-only and others winter-only.

**Task 6.1. Hierarchical PI model with Pmax and Rd as regression functions.** Fit at the **hourly** level, with log-parameterisation for positivity:

```r
f_pi <- bf(
  PP ~ (exp(la) * exp(lP) * PAR) / (exp(la) * PAR + exp(lP)) - exp(lR),
  la ~ 1 + (1 | Year),
  lP ~ 1 + ilr1 + ilr2 + ilr3 + invkT_c + log_flow_c + N_z + Season +
       (1 | p | Year) + (1 | q | Year:DielDate),
  lR ~ 1 + ilr1 + ilr2 + ilr3 + fish_resp_z + invkT_c + log_flow_c + Season +
       (1 | p | Year) + (1 | q | Year:DielDate),
  nl = TRUE
) + student()
```

- `ilr1..3` are **pivot coordinates** of the LTER_1 composition for that year (Hron et al. 2012). Put coral first in the pivot so that `ilr1 = sqrt(3/4) * log(Coral / gmean(Algae, CCA, Other))` reads as "coral relative to everything else". Re-pivot with algae first to get the algae effect. This is the standard way to use a composition as a **predictor** (Comment 3).
- `invkT_c = 1/(k*T) − 1/(k*T_ref)` with T in Kelvin and k = 8.617e-5 eV/K. The coefficient on `lR` is −E (activation energy). This connects to the Metabolic Theory of Ecology and gives a physiologically interpretable temperature effect estimated from **day-to-day** temperature variation.
- `(1 | p | Year)` and `(1 | q | Year:DielDate)` share IDs `p` and `q`, so **year- and day-level deviations in Pmax and Rd are correlated**. That is where the Pmax–Rd covariance belongs, rather than in a causal `Rd → Pmax` path.
- `Season` removes the summer/winter imbalance.
- Priors: `normal(log(100), 0.5)` on the lP intercept, `normal(log(80), 0.5)` on the lR intercept (check these against current estimates), `normal(0, 0.5)` on slopes, `exponential(2)` on SDs, `normal(0.65, 0.2)` on the `invkT_c` coefficient for lR (MTE prior; run a sensitivity analysis with `normal(0, 1)`).
- Check residual autocorrelation within days (hourly ACF). If it is strong, sensitivity options are thinning to every 2nd hour, or fitting a daily-integrated model (daily R and daily GPP as responses, `pp_day`).

**Task 6.2. Per-estimand versions (E4, E6, E7, E9, E10).** From the dagitty sets, create separate fits of the above with **only** the required covariates in `lP`/`lR`. Example: E7 (DayTemp → Rd, direct) needs {Benthos, fish, Flow, Season} and must exclude N. E5 (Benthos → Rd, total) must **exclude fish** (fish are mediators) and must include whatever the DAG says confounds Benthos and Rd (likely `Time`, `DHW`).

**Task 6.3. Propagate benthic and fish uncertainty.** Composition and fish biomass for each year are estimates, not known values. Two options:
- (Preferred, simple) Draw M = 20–50 posterior draws of LTER_1 composition (from the Section 5 model, `posterior_epred` at Site = LTER_1) and of fish respiration (Task 7.1), build M datasets, and fit with `brm_multiple()`. Posterior draws are pooled automatically.
- (Alternative) `me(ilr1, se_ilr1)` measurement-error terms.

**Task 6.4. Sample-size transparency.** In the methods, state: number of days, hours, years, and days per year × season (the table from Task 2.3). Report the **effective number of years** informing benthos→metabolism effects (18), so the reviewer sees we are not overclaiming.

---

## 7. Fish: inside the DAG, contributing to respiration and N (Comment 2)

**Task 7.1. Fish community respiration (mechanistic, independent of the in situ flux measurement).** Use `fish_ind` (sizes per individual) with `fishflux` (Schiettekatte et al. 2020) or Barneche et al. (2014) allometric scaling, I = i0 · M^α · e^(−E/kT). Compute **fish respiration (mmol O₂ m⁻² h⁻¹)** and **N excretion** per transect × site × year, using the year's mean water temperature.
- Compare fish respiration with the in situ measured Rd at LTER_1, as the percentage of Rd attributable to fish per year. This is a direct, quantitative answer to "fish definitely contribute to … ecosystem respiration". Whatever the result (likely a few percent to tens of percent), report it. Note that because Rd is measured as an open-flow, upstream-downstream flux rather than an enclosed incubation, it already integrates fish respiration occurring within the flow path between sensors — the bioenergetics estimate is an independent cross-check on the magnitude, not a wholly separate measurement.
- Use `fish_resp_z` (standardised) as the fish covariate in the Rd model (Task 6.1), rather than raw biomass.

**Task 7.2. Fish response models (E11).** Model herbivore and corallivore biomass at transect level across 6 sites:

```r
bf(herb_biomass ~ ilr1 + ilr2 + ilr3 + herb_lag_z + Year_c + (1 | Site) + (1 | Site:Year))
# family = hurdle_gamma() if zeros exist, else Gamma(link = "log")
```

Adjustment sets come from dagitty. Check residual ACF per site.

**Task 7.3. Fish → N (E12).** Turbinaria tissue N is LTER_1 only, n ≈ 19 years. Fit `log(N_percent) ~ fish_Nexcretion_z + <adjustment set>` with `Gamma` or lognormal, and be explicit that this is a small-sample estimate.

---

## 8. Remove the residual covariate (Comment 5)

**Task 8.1.** Delete `yearresid` everywhere. Where the DAG says `Time` confounds an exposure–outcome pair, include `Year_c` (centred year) **directly** in that model alongside the exposure, as Freckleton (2002) prescribes. A sensitivity version uses a low-flexibility smooth `s(Year, k = 4)`. Do not use more than k = 4, because a flexible smooth absorbs the disturbance signal.

**Task 8.2. Collinearity diagnostic.** Report the correlation between `Year` and `DHW_max` and the posterior correlation of their coefficients (`pairs(fit, variable = c("b_DHW_max", "b_Year_c"))`). If they are strongly correlated, say so. The heat-stress effect is then identified mainly from event years, which is an honest limitation.

---

## 9. Simulation check: show the new approach recovers known effects (supports Comments 1 and 4)

Arif & MacNeil (2023) use simulations from a known DAG. Do the same so the reviewer can see the problem with the old method, and the n = 20 limitation, quantitatively.

**Task 9.1.** Simulate 500 datasets from `dag_mcr`:
- Linear-Gaussian structural equations, with effect sizes set near the new posterior medians.
- Same sample structure as the real data: 20 years × 6 sites for the benthos, and 1 site × 150 days for metabolism.
- Benthos generated on the alr scale and closed with softmax.

**Task 9.2.** On each dataset, fit:
1. The old structure (8-equation simultaneous model on 20 annual means with `yearresid`).
2. The new per-estimand models. For speed, use `lm`/`glm` analogues in the simulation, not brms.

Report bias, RMSE and 95% interval coverage for E1, E4, E7 and E8.

**Task 9.3.** Repeat with a collapsed n = 20 annual dataset versus the disaggregated data, to show the precision gain. Output Figure S-sim.

---

## 10. Robustness and diagnostics

**10.1. Convergence:** Rhat < 1.01, bulk/tail ESS > 400 per chain, zero divergences. Report these in a table.

**10.2. Posterior predictive checks:**
- Benthos: `pp_check(type = "bars_grouped")` per category.
- Metabolism: observed vs predicted diel curves for 6 representative days.
- Fish: `pp_check(type = "dens_overlay")`.

**10.3. Prior sensitivity:** `priorsense::powerscale_sensitivity()` for each estimand model.

**10.4. Out-of-sample checks for the dynamic benthic model:** leave-future-out CV (Bürkner, Gabry & Vehtari 2020, *J. Stat. Comput. Simul.*). Predict year t+1 from data up to t, starting from 2012. This is also the best validation of the projection model.

**10.5. Unmeasured confounding:** for E1, E4 and E8, refit frequentist analogues and run `sensemakr`. Report robustness values: how strong an omitted confounder would have to be to explain away each effect.

**10.6. Alternative DAGs:** refit E4/E5/E8 under (a) `Rd → N` replaced by `Benthos → N`, and (b) no `Time` node. Show the estimates in an SI forest plot.

**10.7. Drop the "reduced SEM".** Removing the paths whose CIs spanned zero, then projecting with the pruned model, is data-driven model selection. Delete that section. Paths are decided a priori by the DAG.

---

## 11. Projections, rebuilt as scenario simulations (addresses "undermine … future projections")

Follow the structure of `GBR_future_simulations.ipynb` in MacNeil et al. (2019).

**Task 11.1.** For each SSP (1-2.6, 2-4.5, 5-8.5) and each CMIP6 model (existing `ann_tos_bc`), compute annual max DHM (or a DHW proxy) at Moorea to 2100.

**Task 11.2.** Disturbance regimes:
- Draw COTS outbreaks and cyclones as stochastic events with empirical recurrence rates (Moorea literature: COTS outbreaks roughly every 15–25 years, and a cyclone annual probability from the regional record).
- Bleaching impacts come from the fitted DHW coefficient (E1 model).

**Task 11.3.** Forward-simulate the benthic composition with the fitted compositional Gompertz:
- Start from the 2025 composition.
- Use 1000 posterior draws × CMIP6 models × disturbance realisations.
- Use process noise from the `Site:Year` SD.

**Task 11.4.** Push each simulated composition and temperature through the metabolism model (E5 total-effect version for Rd, and the corresponding Pmax model) by **g-computation**: predict Rd and Pmax with `posterior_epred`, setting covariates to simulated values. This is the causally coherent way to chain separate estimand models, because no product-of-coefficients mediation is needed.

**Task 11.5.** Present results as **conditional scenarios with full uncertainty**. Truncate to 2050 in the main text and put 2100 in the SI. Add a statement on the range of DHW observed vs projected (extrapolation flag).

---

## 12. Response to reviewers: skeleton (`Response_to_Reviewers.md`)

For each comment: (i) thank and concede the point where it is valid; (ii) what we changed; (iii) where it is in the manuscript; (iv) the key new result.

- **R1, DAG vs SEM.** "We agree, and have reframed the analysis under the Structural Causal Model framework (Pearl 2009; Arif & MacNeil 2023). The DAG (new Fig. 3) now encodes assumptions only. For each of 12 causal estimands we derived the minimal adjustment set via the backdoor criterion (Table S-DAG) and fitted a separate model conditional on that set. Path coefficients from the simultaneous model are no longer interpreted causally. The `Rd → Pmax` path, which reflected shared biomass rather than causation, is replaced by a residual covariance between year- and day-level Pmax and Rd within the photosynthesis–irradiance model. We tested the DAG's implied conditional independences (Table S-dsep)."
- **R2, incomplete DAG and time.**
  - The DAG now includes fish → Rd and fish → N, disturbances (DHW, COTS, Cyclone Oli) and lagged states, which unrolls the coral–algae–herbivore feedback in time.
  - Fish respiration estimated from bioenergetics accounts for X% of Rd (new result).
  - Decline and recovery are separated in a compositional Gompertz model in which disturbances drive losses and herbivory modulates recovery, following the population-dynamic approach of MacNeil et al. (2019).
- **R3, compositional data.** Benthic composition is now modelled with a multinomial-logit (or Dirichlet) likelihood on point counts (Aitchison 1986). Composition enters metabolism models as ilr pivot coordinates (Egozcue et al. 2003; Hron et al. 2012).
- **R4, 20 data points and autocorrelation.** The models now use 6 sites × ~20 years × transects/quadrats for benthos and fish, and ~150 diel cycles (~3,600 hourly observations) for metabolism, with hierarchical random effects. Autocorrelation is handled structurally through lagged states (Gompertz/MAR(1)), and the residual ACFs are shown (Fig. S-ACF). The simulation study (Fig. S-sim) quantifies bias and coverage under the real sample structure.
- **R5, residuals.** `yearresid` has been removed. Year enters directly as a covariate where the DAG identifies it as a confounder proxy (Freckleton 2002).
- **Projections.** These are rebuilt as scenario simulations from the dynamic model, with g-computation to metabolism, truncated to 2050 in the main text, with explicit extrapolation caveats.

---

## 13. Execution order and checkpoints

| Step | Script | Depends on | Checkpoint / deliverable |
|------|--------|-----------|--------------------------|
| 1 | 00, 01 | — | `qc_disaggregation.md`; points-per-quadrat confirmed |
| 2 | 02 | 1 | `disturbance_site_year.csv`; DHW time series plot |
| 3 | 03 | — | DAG figure, Table S-DAG, Table S-dsep (stop and review with PI before modelling) |
| 4 | 04 | 1, 2, 3 | Benthic Gompertz fit, PPCs, ACF, E1–E3 |
| 5 | 05 | 1 | Fish respiration/excretion table; % of Rd from fish; E11, E12 |
| 6 | 06 | 4, 5 | PI hierarchical fit; E4–E10 per-estimand fits |
| 7 | 07 | 3–6 | Forest plot of all estimands (new Fig. 4) |
| 8 | 08 | 3 | Simulation figure S-sim |
| 9 | 10 | 4–7 | Robustness tables |
| 10 | 09 | 4, 6 | Projection figures (main: to 2050; SI: to 2100) |
| 11 | — | all | `Methods_revised.md`, `Response_to_Reviewers.md` |

**Hard stops where a human decision is needed** (do not let the executing model guess):
1. The `Rd → N` vs `Benthos → N` mechanism (Task 4.1).
2. Any DAG edge added after the d-separation tests fail (Task 4.3).
3. Whether to include forereef/fringing habitats to increase leverage on heat stress (more data, but a habitat-specific DAG is needed).
4. Whether to keep projections in the main text.

**Compute notes:** the hourly PI model with covariates and correlated random effects will be slow. Develop it on a 25% subsample of days first, then run the full fit. Use `cmdstanr`, 4 chains, `iter = 3000`, `warmup = 1500`, `adapt_delta = 0.95`, and save every fit with `file = here("Revision","Output","fits","<name>")` so it is not refitted.

---

## 14. Testing log

Working rule for this project: **all new work lives under `Revision/`.** Nothing in the repository root (`MCR_Productivity_Analysis.qmd`, `Data/`, `Output/`) is read-write — original scripts, raw data, and original output are read-only inputs. Every script, derived data file, model object, and figure produced while executing this plan goes under `Revision/R/`, `Revision/Data/derived/`, or `Revision/Output/`. This section is updated after each step is executed, recording what was run, what it found, and what still needs a decision.

### Step 1 — Data disaggregation (Section 2 of the plan)

**Status:** done. **Script:** `Revision/R/00_packages.R`, `Revision/R/01_data_disaggregated.R`. **Outputs:** `Revision/Data/derived/{benthic_quad,benthic_transect,benthic_site,fish_transect,fish_ind,pp_hour,pp_day}.csv`, `Revision/Output/qc_disaggregation.md`. Run via `source("Revision/R/01_data_disaggregated.R")` from the repository root; reads only from `Data/raw_data/`, writes only under `Revision/`.

**What it does:** rebuilds the benthic, fish, and metabolism data at their native replication (quadrat/transect/site/hour/day) instead of the 20-row `Year_Averages` table, per Task 2.1–2.3.

**Key findings, to carry into later steps:**

1. **Points per quadrat is not constant across years — CONFIRMED 2026-10-01.** 25 points/quadrat in every year except **2020, which uses 50 points/quadrat**. Verified by two independent methods that agree exactly: (i) the smallest non-zero `Percent_Cover` increment per year, and (ii) the per-quadrat GCD of non-zero cover values restricted to quadrats with ≥4 taxa (so the GCD is diagnostic rather than a degenerate artefact of sparse composition) — this stricter check still gives 25 points/quadrat in 96–100% of diagnostic quadrats every year, and 50 points/quadrat in 96.4% of the 196 diagnostic quadrats in 2020, ruling out coincidence or a data-entry artefact. This was not documented anywhere in the original analysis. Likely explanation: a COVID-19-era change in photoquadrat image-analysis protocol — 2020 is also the one year metabolism sampling was summer-only (item 7), consistent with disrupted fieldwork that year. `benthic_quad.csv` uses a year-specific lookup (`n_pts_inferred`) as the multinomial `trials()` denominator in Task 5.1. **Resolved for pipeline purposes;** the one remaining recommendation (not a blocker) is to cross-check the MCR LTER data-package version history/EDI changelog for a documented 2020 protocol note before the Section 5 model results are finalized for submission.
2. **Quadrat design is otherwise fully consistent:** every Site × Year has exactly 50 quadrats (5 transects × 10 quadrats), for all 6 backreef sites across all 20 years (2006–2025) — no missing site-years. This gives 120 site-years / 6,000 quadrat-years of benthic data versus the 20 LTER_1-only annual points used in `Year_Averages`.
3. **Quadrat cover totals mostly close to 100%,** with small exceptions: 21 of 6,000 quadrats (0.35%) have Coral+Algae+CCA+Other ≠ 100% ± 0.01, concentrated in 2020 (17 quadrats) and a few in 2010, 2019, 2021. These are kept as-is with their observed `n_total` (not forced to the nominal count) so Task 5.1's `trials(n_total)` reflects what was actually recorded.
4. **Fish filters:** of 29,089 raw backreef fish records (all 6 sites), 63 rows were removed as shark outliers (Biomass > 8000 g) and 62 as the negative-biomass missing-data code, retaining 28,964 (99.57%). Matches the scale of exclusions in the original single-site processing, now quantified across all sites.
5. **Fish survey area derived from first principles, not hard-coded.** `Location` strings and the `Swath` column confirm two fixed swath widths (1 m and 5 m) are recorded for every Site × Year × Transect combination (verified: all 480 combinations have exactly 2 distinct swath values, every year 2006–2025). Combined with the MCR LTER fish protocol (four 50 m transects per site-habitat), area per transect = 50 × (1+5) = 300 m², and 4 transects/site reproduces the original code's hard-coded 1200 m²/site denominator exactly. The derivation is now explicit and kept at transect resolution in `fish_transect.csv`, rather than a magic number.
6. **Fish trophic grouping — CONFIRMED 2026-10-01.** Beyond the original Herbivore/Corallivore split, `fish_ind.csv`/`fish_transect.csv` add Planktivore, Invertivore, Piscivore, Omnivore, and Other, needed for the fish→Rd and fish→N estimands (Tasks 6.2, 7.1–7.3). `Fine_Trophic == "Omnivore"` appears under three different `Coarse_Trophic` labels (Planktivore, Primary Consumer, Secondary Consumer); checked these by `Family` composition rather than assuming they're interchangeable. Planktivore-coded Omnivore rows are overwhelmingly Pomacentridae (damselfish) and are grouped with Planktivore. Primary/Secondary-Consumer-coded Omnivore rows are taxonomically mixed (Pomacentridae, Pomacanthidae, Tetraodontidae, Balistidae, Zanclidae, Ostraciidae) with no single dominant family, and together represent ~6% of total retained biomass — too large and biologically distinct a group to fold into the genuinely rare "Other" catch-all (Fish Scale Consumer, Sediment Sucker, 2 unidentified/no-fish rows, <0.5% of rows). These are now broken out as their own **"Omnivore"** trophic group (7 groups total, up from 6). **Resolved;** no longer treated as a hard-stop item, though the PI should sanity-check the family-level breakdown in the updated QC report (`Revision/Output/qc_disaggregation.md`, Section B2) before the Section 7 fish models are fitted.
7. **Metabolism coverage confirmed as planned:** 29 raw hourly files, 157 diel days seen, 122 retained after the existing complete-24-hour filter (77.7%), spanning 18 years (2008–2025, with gaps). The year × season day-count table in the QC report reproduces and extends the imbalance already flagged in Section 0.3 (e.g., 2010, 2020, 2025 summer-only; 2008, 2011, 2013, 2021 winter-only), confirming `Season` must be a covariate in the Task 6.1 PI-curve model rather than left implicit in annual means.
8. **Sample-size uplift achieved:** benthic 20 → 120 site-years (6,000 quadrat-years); fish 20 → 480 transect-years (3,360 transect-year × trophic-group rows across 7 groups); metabolism 20 annual Pmax/Rd points → 2,928 hourly rows / 122 day-level rows. Metabolism stays single-site (LTER_1 only — the in situ Lagrangian (upstream-downstream) sensor arrays used to measure PP were only ever deployed at that site), so between-site replication is not available for the respiration/Pmax models; this limitation is unchanged from Section 0.3 and should stay explicit in the methods and response letter.
9. **Correction (post-hoc):** earlier drafts of this plan and the Step 1 script referred to the PP measurement as a "flume". It is not an enclosed chamber/flume incubation — it is an **in situ Lagrangian (upstream-downstream) flux measurement**: paired UP/DN sensor arrays record oxygen, temperature and flow as water moves across the reef, and net community production/respiration is derived from the upstream-to-downstream change combined with flow velocity (`UP_Oxy`, `DN_Oxy`, `UP_Velocity_mps`, `DN_Velocity_mps`, the `UPDN` deployment ID). This has been corrected throughout the plan and in `01_data_disaggregated.R`. It also means the measurement already integrates whatever benthos and fish sit within the flow path between the two sensors — relevant context for Task 7.1's comparison of bioenergetics-derived fish respiration against measured Rd, and for interpreting `Flow_mean` and `UPDN` as covariates rather than incidental metadata in the Section 6 PI-curve model.

**Open items resolved (2026-10-01):** both items previously flagged for PI confirmation were investigated with additional data-internal checks (not external metadata access) and resolved well enough to proceed:
- *2020 benthic points-per-quadrat:* confirmed via two independent, mutually agreeing methods (item 1 above). Residual recommendation (non-blocking): verify against the MCR LTER data-package changelog if/when that documentation is available, for the methods section's citation trail.
- *Omnivore trophic-group assignment:* confirmed via `Family`-level taxonomic composition (item 6 above) — the three `Coarse_Trophic`-coded Omnivore subsets are not interchangeable, and a dedicated "Omnivore" group is now used instead of folding the larger two subsets into "Other". Residual recommendation (non-blocking): have the PI sanity-check the family breakdown in `qc_disaggregation.md` Section B2 against field notes/expert judgement.

Neither item blocks downstream work; both are now implemented in `01_data_disaggregated.R` and reflected in the regenerated `benthic_quad.csv` and `fish_transect.csv`/`fish_ind.csv`.

**Not yet done:** Sections 3 (disturbance covariates: DHW, COTS, cyclone), 4 (DAG and adjustment sets), and everything downstream. `benthic_site`/`benthic_transect`/`fish_transect`/`pp_day` are ready to be joined to the Section 3 disturbance table once it exists.

### Step 2 — Disturbance and heat-stress covariates (Section 3 of the plan)

**Status:** done. **Script:** `Revision/R/02_disturbance_covariates.R`. **Outputs:** `Revision/Data/derived/disturbance_site_year.csv`, `Revision/Output/qc_disturbance.md`, plus cached intermediates under `Revision/Data/derived/_cache/` (raw satellite/logger pulls, so re-running the script does not re-hit the network). Run via `source("Revision/R/02_disturbance_covariates.R")` from the repository root.

**What it does:** builds a Site x Year (6 sites x 21 years = 126 rows) disturbance table with satellite DHW (primary), in-situ logger DHW (secondary/sensitivity), COTS density, and a Cyclone Oli indicator, per Task 3.1–3.4.

**Process note (network reliability):** both the EDI/PASTA API and large ERDDAP `griddap()` pulls were unreliable in this environment — EDI's web portal returned a bot-verification page and the PASTA REST API returned an explicit `403 ... Public Access ... not authorized` for `knb-lter-mcr.1039` even via direct `httr::GET()` with a browser user agent (not just the fetch tool); separately, a large multi-year `griddap()` pull caused the R session to become unresponsive and ultimately crash with no work saved. Resolution: (1) the COTS file was downloaded manually by the PI and placed in `Data/raw_data/`; (2) `02_disturbance_covariates.R` was restructured so every slow/network step (the two ERDDAP pulls, the LTER_2 logger extraction/aggregation) writes its result to `Revision/Data/derived/_cache/` immediately on success and is skipped on re-run if the cache file exists, so a crash mid-script loses at most one in-flight step. This pattern should be reused for the remaining network-dependent steps (e.g., any `fishflux` lookups in Task 7.1).

**Key findings, to carry into later steps:**

1. **Benthic survey timing is austral-summer (mostly January), not April–May.** An earlier draft of the plan assumed April–May survey timing; the actual median survey date is in January for 20 of 21 years (2005 is the one May survey). This is used as the anchor for a trailing 364-day DHW window per year rather than the calendar year. 2023 has an unusually wide within-"Year" date spread (2022-02-27 to 2023-01-24) — flagged for the PI, not resolved further since it doesn't block the covariate calculation.
2. **Per-site DHW coordinates — resolved with the PI as a documented simplification.** The 6 backreef sites span all 3 shores of Moorea (not one 5 km satellite pixel), but only LTER_1's coordinate is documented in this repo. After EDI/PASTA access and web search both failed to turn up authoritative LTER_2–6 coordinates, the PI chose to use a single island-wide DHW series (from the MCR LTER network centroid, -17.4909/-149.826) applied to all 6 sites, rather than guessing per-site coordinates. **Consequence carried into Section 4/5:** `DHW_max`, `DHW_logger_max`, and `Cyclone` are identical across sites for a given year by construction — only `COTS_density_m2` varies by site in this table. Between-site contrasts attributed to heat stress or the cyclone in the Section 5 Gompertz model will be confounded with between-site differences in starting benthic state, not independent exposure. State this limitation explicitly in the methods.
3. **Satellite DHW (primary, NOAA Coral Reef Watch `CRW_DHW`, daily, 5 km) pulled 2004-01-01 to 2025-03-01, 7726 days, 5 missing.** Sanity check passed: the 2020 survey window (ending Jan 2020) shows `DHW_max = 3.36`, consistent with the documented 2019 Moorea bleaching event (plan Section 0.3) — the window logic is capturing real heat-stress history, not an artefact. 2025 also shows an elevated value (2.30).
4. **In-situ logger DHW (secondary/sensitivity) is an approximation, not NOAA's validated product — flagged as such.** The LTER_2 backreef (2 m) logger is the only backreef logger (consistent with Step 1). NOAA's official MMM baseline raster isn't served on the ERDDAP instance used here, so MMM was proxied as the maximum calendar-month-mean of satellite SST at the same pixel over 1985–2025 (28.88°C, March) and applied to the logger's daily means with the standard CRW HotSpot/12-week accumulation rule. The resulting series agrees in direction with the satellite series (both peak in 2020, at 4.08°C-weeks logger-based vs 3.36 satellite-based) but should be used only as a sensitivity check, not a primary estimate.
5. **COTS data required manual download — automated EDI/PASTA access was blocked.** See the process note above. Raw file: `MCR_LTER_COTS_abundance_2005-2025_20250310.csv`, 1512 rows (6 sites x 3 habitats x 21 years x 4 transects), 504 backreef rows kept, matching the full 6x21x4 design exactly (no missing site-year-transects). Density computed as count / (4 transects x 250 m² belt-transect area, same 5x50 m geometry as the fish survey).
6. **COTS density is very low throughout (max 0.012 ind/m² at any site-year); the backreef data show both documented outbreaks only weakly, but a forereef cross-check (pulled for comparison, not used as a covariate) confirms both are real.** Backreef totals (6 sites): low single digits 2005–2007, a modest rise to 6 (2008) and 12 (2009), back to 0 by 2014 and staying at 0–1 through 2022, then a **pulse of 16 individuals in 2024** (concentrated at LTER_2 and LTER_3, 8 and 6 individuals respectively — checked at the individual-transect level, not a data-entry artefact; the PI independently confirmed this is a real event). Pulling the forereef series for the same years shows both pulses far more clearly — 108 individuals in 2008 and 70 in 2009 (the documented Kayal et al. 2012 outbreak), and 17 (2023) / 56 (2024) for a **second, previously undiscussed 2023–2024 outbreak** — confirming (a) the muted backreef signal in 2008–2009 is expected, since Moorea COTS outbreaks are predominantly a forereef phenomenon (Task 3.3), and (b) the 2024 backreef pulse is a real, synchronised event and not an artefact. **Action:** write up 2023–2024 as a second documented Moorea COTS outbreak alongside 2006–2010 in the methods/disturbance-history text (Section 0.2); the forereef comparison table is in `qc_disturbance.md` Section D1.
7. **Cyclone Oli's indicator attributes impact to the 2011 survey, not 2010 — a change from the original analysis's implicit framing.** Closest approach to Moorea was taken as 2010-02-16. The 2010 benthic survey (survey_date 2010-01-01) predates the cyclone; the first survey taken after it is 2011-01-17 (335 days later), which the 365-day rule flags as `Cyclone = 1`. This shifts which annual benthic transition the Section 5 Gompertz model will attribute to the cyclone, relative to how the original qmd's static per-year framing likely treated it. The indicator also does not yet distinguish shore-facing exposure (backreef wave exposure during Oli was much lower than forereef, per Section 0.2/3.3) — a wave-energy proxy was not attempted.

**Open items for the PI (non-blocking, flagged for review before Section 4/5 modelling):**
- Confirm whether the 2011 (not 2010) cyclone attribution is intended, or whether a different closest-approach date / window rule should be used.
- **Resolved (2026-10-01):** the 2024 COTS pulse (item 6) is confirmed real by the PI and corroborated by the forereef cross-check (much larger there, 56 vs. 16 individuals) — write up 2023–2024 as a second Moorea COTS outbreak in the disturbance-history section.
- Decide whether the island-wide-DHW simplification (item 2) is acceptable for submission, or whether per-site coordinates for LTER 2–6 should be sourced (e.g., from a GIS layer the PI has access to) and the DHW pull redone per-site.
- The in-situ logger DHW's locally-derived MMM proxy (item 4) should ideally be replaced with NOAA's official MMM baseline if it can be obtained, before being presented as more than an internal sensitivity check.

### Step 3 — Time-indexed DAG and adjustment sets (Section 4 of the plan, Tasks 4.1/4.2/4.4)

**Status:** done for Tasks 4.1, 4.2 and 4.4. Task 4.3 (testing the DAG against data via `impliedConditionalIndependencies()`/`localTests()`/a Shipley d-sep test) deliberately deferred to a follow-up step — see "Not yet done" below. **Script:** `Revision/R/03_dag.R`. **Outputs:** `Revision/Output/table_s_dag.csv` / `.md` (Table S-DAG), `Revision/Output/dag_mcr_rd_to_n.png` / `.pdf`, `Revision/Output/dag_mcr_benthos_to_n.png` / `.pdf`. Run via `source("Revision/R/03_dag.R")` from the repository root.

**Hard-stop decision (Section 13, item 1) — resolved with the PI 2026-10-01:** whether `N` (nitrogen) is caused by `Rd` (respiration-linked remineralisation) or directly by `Benthos`. The PI chose **"include both as alternatives"**: the DAG is built in two variants (`dag_mcr_rd_to_n`, `dag_mcr_benthos_to_n`), both are plotted, and both flow through Table S-DAG's adjustment sets. The choice is to be revisited after Task 4.3's d-separation test is run (not yet done), rather than picked now on priors alone.

**What it does:**
1. Defines both DAG variants in `dagitty`, sharing every edge in the plan's draft (Task 4.1) except the one differing `N`-mechanism edge. Confirmed both parse as valid acyclic graphs.
2. Computes minimal sufficient adjustment sets for all 13 estimand rows (E1–E10, E12, plus E11 split into E11a/E11b since dagitty defines adjustment sets per exposure-outcome pair, not per the plan's combined "Corallivores/Herbivores" outcome label) x both DAG variants x both `"total"`/`"direct"` effect types as specified per estimand (Task 4.2), using `adjustmentSets(..., type = "minimal")`. E6 and E12 ("Fish → Rd" and "Fish → N") use dagitty's generalised backdoor criterion for a joint vector exposure (`c("Herb","Corall","OtherFish")`) rather than picking one fish group arbitrarily.
3. Plots both DAG variants as coefficient-free causal-structure figures (Task 4.4, replacing the old "SEM" Figure 3 with its β labels) — iterated once after the first draft (circular nodes) clipped several labels (`OtherFish`, `DayTemp`, `Benthos_lag`); the final version uses small points with `ggrepel`-style external labels, which rendered all 17 nodes legibly.

**Key findings, to carry into Section 5/7 modelling:**

1. **Adjustment sets largely match the plan's rough expectations, with the dagitty-verified sets sometimes smaller or including a DHW/Season/Flow substitution the plan didn't anticipate.** E.g. E1 (`DHW -> Benthos`, total) needs only `{Time}` — the plan's guess of `{Time, Benthos_lag}` turned out to include an unnecessary covariate, since `Benthos_lag` doesn't open a backdoor path from `DHW` (it isn't a common cause of `DHW` and `Benthos`). Full sets are in `table_s_dag.csv`/`.md`.
2. **Several estimands have more than one valid minimal adjustment set** (up to 8, for E7/E10 in some variants) because `Season`, `Flow`, and `DHW` sit on an interchangeable chain (`Season -> DayTemp`, `Season -> Flow`, `DHW -> DayTemp`), so blocking any one of several nodes on that chain is sufficient. This is a genuine feature of the DAG, not an error, but it means the Section 7 modelling step will need to **choose one adjustment set per estimand** from the listed options (e.g. prefer the most parsimonious, or the one matching variables already in the disaggregated tables) rather than there being a single unambiguous answer.
3. **The two DAG variants produce different adjustment-set counts for the Pmax estimands (E9/E10) specifically**, because the `Rd -> N` variant opens an additional path `Benthos -> Rd -> N -> Pmax` that the `Benthos -> N` variant does not have — the `Rd -> N` variant's alternate minimal sets therefore include `Rd` itself as a required covariate where the `Benthos -> N` variant doesn't need it. This is a concrete, checkable consequence of the hard-stop decision, and should sharpen the Task 4.3 d-sep test's ability to discriminate between the two variants once it's run.

**Not yet done (at the time Step 3 was written):** Task 4.3 (testing both DAG variants against data) was deliberately left out of this step's scope. Running it requires several things Steps 1–2 did not build: ilr coordinates for the compositional `Benthos` node (from `benthic_site.csv`), a nitrogen (`N`) covariate (built in Step 2b below), and `Rd`/`Pmax` proxies ahead of the real Section 6 PI-curve fit (e.g. `R_mean`/`GP_mean` from `pp_day.csv`, clearly labelled as proxies, not the final modelled quantities). Section 5 (compositional Gompertz) and everything downstream also remains to be done.

### Step 3b — Task 4.3: testing both DAG variants against data

**Status:** done. **Script:** `Revision/R/03b_dsep_test.R`. **Outputs:** `Revision/Output/table_s_dsep.csv` / `.md` (one row per `impliedConditionalIndependencies(..., type="missing.edge")` test, both DAG variants), `Revision/Output/qc_dsep.md` (Shipley/Fisher's C global test plus full interpretation). Also writes `Revision/Data/derived/_cache/dsep_panel_full.csv` (120 rows, 6 sites) and `dsep_panel_lter1.csv` (18 rows, LTER_1 only) — the joint analysis datasets built for this test.

**What it does:** builds ilr coordinates for the compositional `Benthos` node (zero-replaced via `zCompositions::cmultRepl` on the underlying point counts, per Task 5.1's approach), site-year fish trophic-group summaries (`Herb`/`Corall`/`OtherFish`), and Rd/Pmax/DayTemp/Flow/Season PROXIES from `pp_day.csv` (explicitly labelled as proxies, not the Section 6 model's eventual output); joins everything into a 6-site panel (for tests not involving a metabolism node) and an LTER_1-only panel (for tests that do, since metabolism is single-site); tests both DAG variants via `dagitty::ciTest(..., type="cis.pillai")` (a canonical-correlation test that handles the compositional node's full 3-dimensional ilr vector natively, not reduced to one scalar); and computes a Shipley (2000) Fisher's C global statistic from the basis-set implications.

**Headline result: 25 of 440 testable missing-edge independencies are violated (BH-adjusted p<0.05), with a clear, interpretable, and highly consistent pattern across both DAG variants** (details and exact statistics in `qc_dsep.md`):

1. **`Benthos_lag`, `Herb_lag`, and `N_lag` are each strongly correlated with `Time`** (r = 0.67/0.56/-0.57, all p < 1e-9, n = 107–114, full 6-site panel). This is a genuine gap, not a surprising result: these nodes are literally last year's `Benthos`/`Herb`/`N`, which the DAG already connects to `Time`, but the DAG never draws the analogous `Time -> *_lag` edges — a lagged quantity does not stop being time-varying just because it is one year behind. **Candidate fix: add `Time -> Benthos_lag`, `Time -> Herb_lag`, `Time -> N_lag`.**
2. The three lag nodes are also pairwise correlated with each other (r = 0.32–0.43, p < 0.001) — plausibly fully explained by (1) via their shared `Time` parent; worth re-testing after adding the `Time ->` edges before concluding anything further is missing.
3. **`COTS` is correlated with `N` and `N_lag`** (r = 0.30–00.47, p < 0.002, n = 107) despite having no path to `N` and no connection to `Time` in the current DAG. Step 2's finding that COTS outbreaks cluster in specific years (2008–2009, 2023–2024) rather than tracking site-level lagged coral cover smoothly is consistent with **`Time -> COTS` also being a missing edge** (a temporal/regional outbreak-event signal, not purely local density-dependence).
4. `COTS`–`Pmax` correlation at LTER_1 (r = -0.70 to -0.81) — flagged with lower confidence (small n = 17–18, `Pmax` is a raw proxy, single site).
5. **`Pmax`–`Rd` remain correlated given their shared measured causes** (r = 0.76–0.78) — this is NOT a new problem. It is exactly what Task 4.1 anticipated when it deliberately removed the `Rd -> Pmax` edge and assigned their shared variance to a residual correlation inside the Section 6 PI-curve model instead of a DAG arrow. This result confirms that modelling choice was necessary, rather than indicating the DAG is misspecified.

**Important methodological caveat, reported rather than hidden:** the Shipley/Fisher's C global test rejects both DAG variants overwhelmingly (C ≈ 58, df = 14, p < 0.0001, identical for both variants). This should NOT be read as strong evidence the causal structure is wrong — `impliedConditionalIndependencies(type="basis.set")` tests each node against the *joint block* of all its non-descendant/non-parent nodes at once, and several of these blocks have 7–10 variables tested at n = 17–18 (the LTER_1-only metabolism panel): a near-saturated regression regime that will reject almost any model, true or false. 10 of 17 basis-set tests were not even computable given the available n. The pairwise, low-conditioning-dimension Table S-dsep tests above are far more informative and are the basis for the recommendations here; the Fisher's C result is a genuine data-size limitation (Section 0.3), not a flaw in the DAG, and should be reported in the methods as inconclusive given sample size rather than as a failed global test.

**Hard stop (Section 13, item 2) — resolved with the PI 2026-10-01:** add all four candidate edges and re-test. `Time -> Benthos_lag`, `Time -> Herb_lag`, `Time -> N_lag`, and `Time -> COTS` were added to both DAG variants in `03_dag.R` (documented with an inline comment and date), and `03_dag.R` / `03b_dsep_test.R` were re-run.

**Result of the re-test — a clean before/after comparison:**

| | Before (original DAG) | After (4 edges added) |
|---|---|---|
| Missing-edge violations (BH-adjusted p<0.05) | 25 / 440 | **15 / 440** |
| Shipley/Fisher's C (both variants identical) | C = 58.03, df = 14, p < 0.0001 | **C = 18.05, df = 14, p = 0.2046** |

The drop in Fisher's C (58 → 18, rejection → non-rejection) is a meaningful, consistent signal that the four added edges captured real missing structure, not noise — while still treating the *passing* global test as corroborating rather than decisive, given the basis-set test's small-n limitations documented above. Adding `Time -> Benthos_lag/Herb_lag/N_lag` **fully resolved** the `Benthos_lag`–`N_lag` violation but only **partially** resolved `Benthos_lag`–`Herb_lag` (still violated after conditioning on `Time`, r = 0.37, both variants) — a genuine residual relationship between lagged coral cover and lagged herbivore biomass not explained by the shared secular trend. Adding `Time -> COTS` did **not** resolve the `COTS`–`N`/`N_lag` violations (still significant conditioning on `Time`) — pointing to a real, unexplained COTS–nitrogen relationship. One literature-supported candidate mechanism is the nutrient-enrichment hypothesis for COTS outbreaks (e.g. Fabricius et al. 2010), i.e. a possible `N -> COTS` edge — **not added here**, since it would be a new substantive causal claim requiring its own literature check and PI sign-off, flagged as an open item for the next DAG review rather than acted on from a single correlation. The `COTS`–`Pmax` (LTER_1, low confidence) and `Pmax`–`Rd` (expected, per Task 4.1's design) findings are unaffected by this edit, as anticipated.

Updated `Revision/Output/table_s_dag.csv`/`.md`, `dag_mcr_*.png`/`.pdf`, `table_s_dsep.csv`/`.md`, and `qc_dsep.md` all reflect the post-edit DAG. Note the adjustment sets for E2 (`COTS -> Benthos`) and E3 (`Herb_lag -> Benthos`) changed as a direct, expected consequence: both now require `Time` in their adjustment set, since `Time` is newly a confounder on those paths.

### Step 2b — Nitrogen (N) covariate (prep for Task 4.3; not in the plan's original file list)

**Status:** done. **Script:** `Revision/R/02b_nutrient_covariate.R` (numbered 02b, not in Section 1.1's original file list — added because the plan didn't anticipate needing an `N` covariate until Task 4.3/Section 7's E9/E10/E12 estimands came into view). **Outputs:** `Revision/Data/derived/N_sample.csv` (814 individual Turbinaria tissue samples, 6 sites, 2007–2025), `Revision/Data/derived/N_site_year.csv` (108 rows, 6 sites x 2007–2024, mean/sd/n per site-year), `Revision/Output/qc_nutrient.md`.

**Decision (made with the PI 2026-10-01, see AskUser record):** build both a sample-level table and a site-year summary (mirroring `benthic_quad.csv`/`benthic_site.csv`), per Comment 4's "disaggregate, don't pre-average" principle, rather than a site-year mean only.

**Source decision and reasoning:**
1. **Macroalgal tissue %N (Turbinaria ornata, Backreef), not `WaterColumnN.csv`.** The water-column dissolved N+N file is single-site (LTER_1) and stops in 2018 — it misses the 2019 bleaching event, the 2020 DHW peak, and the entire 2023–2024 COTS outbreak documented in Step 2. The original qmd already used tissue %N as the primary nutrient proxy for this reason, with water-column N+N only as a validation check.
2. **Reproducing that validation check surfaces a correction to the original qmd's narrative.** The original comment describing this relationship calls it a "tight linear relationship." Reproducing the regression (`Nitrite_and_Nitrate ~ N_percent`, LTER_1, n = 12 overlapping years) gives Pearson r = 0.63 (95% CI 0.09–0.88), p = 0.028, R² = 0.40 — a real but only moderate relationship with a wide, data-limited CI (n = 12), not a "tight" one. Worth softening this language in the revised methods text.
3. **Turbinaria only, not pooled with Sargassum.** Checked explicitly: mean tissue %N is 1.14% for Sargassum vs. 0.68% for Turbinaria (~2 SD apart), and Sargassum sampling stopped after 2014 (protocol shift to Turbinaria-only) — pooling would conflate a genus effect with a time/protocol-era effect.
4. **All 6 sites, 2007–2024 core window.** Coverage is essentially complete (5–10 samples per site-year) with exactly one gap (LTER_5, 2007) in that window. 2005–2006 have no Turbinaria samples at all (Sargassum-only or no CHN sampling that early) and 2025 is nearly empty (5 samples, 1 site — the same lab-processing-lag pattern already seen elsewhere in this project for the current year); both are excluded from `N_site_year.csv`'s core window but retained as recorded in `N_sample.csv`.

**Open item for later modelling:** the one remaining gap (LTER_5, 2007) in `N_site_year.csv` needs an explicit missing-data strategy (e.g. a `Site:Year` random effect absorbing it, or leaving it `NA` and letting `brms` handle it) when this table is used in Task 4.3 or the Section 7 N-related estimand models — not resolved here.

### Step 4 — Compositional Gompertz model of benthic dynamics (Section 5, Tasks 5.1 and 5.4; Tasks 5.2/5.3 deferred)

**Status:** Task 5.1 (additive multinomial model) and the E1-E3 half of Task 5.4 (average marginal effects / counterfactual) done. Tasks 5.2 (interaction/recovery-mediation model + `loo()` comparison) and 5.3 (residual autocorrelation check) deliberately deferred to a follow-up step. **Script:** `Revision/R/04_benthic_gompertz.R`. **Outputs:** `Revision/Data/derived/benthic_transect_model.csv`, `Revision/Output/fits/benthic_gompertz_additive.rds` (cached `brm()` fit), `Revision/Output/estimands_E1_E3_coral_proportion.csv`, `Revision/Output/counterfactual_no_disturbance.csv`, `Revision/Output/qc_benthic_gompertz.md`.

**Resolution used (documented, not a silent shortcut):** fit at **transect level** (benthic_transect.csv, 570 rows = 6 sites x 19 years x 5 transects, after dropping each site's first year for the lag), not quadrat level (6000 rows) — Task 5.1's own named fallback ("if the multinomial is too slow, fit at transect level"), chosen up front given this session's history of long-running compute. The `(1|Site:Year:Transect)` term is dropped accordingly (no quadrat replication left to absorb at this resolution); `(1|Site)` and `(1|Site:Year)` are retained.

**Compute lesson learned, fixed for future steps:** the first fit ran with brms' default rstan backend and **no `mc.cores` set**, so the 4 chains ran sequentially — chain 1 took ~9 minutes, chain 3 anomalously took ~48 minutes (likely system contention), for a total wall time of **76 minutes** for a 570-row model. Added `options(mc.cores = min(4, parallel::detectCores()))` to `00_packages.R` so all chains run in parallel; the refit (after also fixing the COTS scaling issue below) completed in **~12 minutes**. This should be the default going forward for every `brm()` call in Sections 6-7.

**Data issue found and fixed before trusting any COTS coefficient:** `COTS` (density, ind/m^2, max 0.012) is tiny relative to its `normal(0,1)` prior and the other predictors' scales, which left its posterior almost exactly equal to the prior (CI roughly -1.9 to 1.9, uninformative). Rescaled to ind/100m^2 (max ~1.2) before fitting; the rescaled coefficient's CI narrowed substantially (e.g. `muCoral_COTS`: -0.42 to 1.91 with a mean pinned near the prior before rescaling, -0.42 to 0.63 after -- removing most of the scale-driven posterior inflation, though still not a credible non-zero effect).

**Model converged cleanly:** all Rhat <= 1.01, Bulk/Tail ESS >= ~1100 for every parameter across 4000 post-warmup draws; 2 divergent transitions out of 4000 (0.05%, minor, noted not chased). Headline coefficients (full detail in `qc_benthic_gompertz.md`): strong near-unit density dependence in coral share (`muCoral_alr_coral_lag` = 1.01 [0.90, 1.12]), a credible negative coral-algae competition term (`muCoral_alr_algae_lag` = -0.35 [-0.59, -0.13]), a declining secular trend in CCA share (`muCCA_Year_c` = -0.16 [-0.23, -0.09]), and a borderline-positive DHW-CCA association (`muCCA_DHW_max` = 0.29 [0.00, 0.57]).

**E1-E3 average marginal effects on coral proportion (Task 5.4) — all three are statistically null, reported honestly rather than overclaimed:**

| Estimand | AME (mean) | 95% CI | P(AME > 0) |
|---|---|---|---|
| E1: DHW → coral proportion, per +1 DHW °C-week | 0.0015 | [-0.0153, 0.0195] | 0.563 |
| E2: COTS → coral proportion, per +1 ind/100m² | 0.0265 | [-0.0618, 0.1328] | 0.686 |
| E3: Herb_lag → coral proportion, per +1 SD herbivore biomass | 0.0037 | [-0.0162, 0.0242] | 0.634 |

A joint "no-disturbance" counterfactual (DHW_max = 0 & COTS = 0 everywhere, vs. observed) gives a mean change of -0.0019 [-0.0103, 0.0068] — also null. **Interpretation:** at this 19-year, 6-site, transect-level resolution, none of DHW, COTS, or lagged herbivory show a detectable average effect on coral share; all three disturbance/recovery covariates are heavily zero-inflated (DHW_max = 0 in most site-years; COTS = 0 in the large majority), limiting power. The best-supported signal in the model is the lagged-composition (Gompertz) structure itself (persistence + competition), not the disturbance terms — state this plainly in the methods rather than letting the null E1-E3 results read as "no effect exists."

**Adjustment-set check (Task 5.4's note):** confirmed against `table_s_dag.csv` that this single additive model's covariates satisfy all three estimands' adjustment sets simultaneously (E1/E3 need `{Time}`; E2 needs `{Benthos_lag, Time}`; the model includes Time, Benthos_lag (as the three `alr_*_lag` terms), and Herb_lag together, none of which are descendants of DHW/COTS/Herb_lag in the DAG) — so one model validly answers E1-E3 without violating the "one model per estimand" principle.

**Not yet done:** Task 5.2 (herbivory-modifies-recovery interaction model: `herb_lag_z:alr_coral_lag`, `herb_lag_z:DHW_max`, plus a `loo()` comparison against this additive model) and Task 5.3 (residual autocorrelation check). Flagged as the next step rather than folded in here, consistent with this project's incremental-scope pattern.

### Step 4b — Task 5.2: recovery-mediation interaction model vs. additive

**Status:** done. **Script:** `Revision/R/04_benthic_gompertz.R` (Section G, appended). **Outputs:** `Revision/Output/fits/benthic_gompertz_interaction.rds` (cached fit; `benthic_gompertz_additive.rds` was also re-saved automatically once `loo` criteria were attached to both), updated `Revision/Output/qc_benthic_gompertz.md` (Sections G-H).

**What it does:** adds `herb_lag_z:alr_coral_lag` (herbivory modifies coral's density dependence/recovery speed) and `herb_lag_z:DHW_max` (herbivory modifies coral's heat-stress susceptibility) to the additive model's formula, refits, and compares the two models via `loo_compare()`. Because brms applies one linear-predictor formula across all non-reference multinomial categories, both interactions were estimated for `muCoral`, `muAlgae`, and `muCCA` alike; only the `muCoral` versions are Task 5.2's estimand of interest.

**Result: the interaction model is preferred by neither the parameter estimates nor (nominally) the LOO comparison.** Both new coefficients are null for coral: `muCoral_alr_coral_lag:herb_lag_z` = 0.01 [-0.14, 0.15], `muCoral_DHW_max:herb_lag_z` = -0.02 [-0.19, 0.14] (both well-converged: Rhat = 1.00, ESS > 2500). `loo_compare()` favours the simpler additive model (elpd_diff = -26.2, se_diff = 13.2 for the interaction model) — roughly 2 SE worse, consistent with adding two uninformative parameters.

**Important reliability caveat, reported rather than glossed over:** `loo()`'s Pareto-k diagnostics are poor for BOTH models — only ~36% of the 570 observations have k ≤ 0.7; 27-28% are "bad" and 36-37% are "very bad" for both fits. This is a known problem for multinomial/binomial-style models with large trial counts (`n_total` ≈ 125-250 here), where single highly-informative observations break PSIS-LOO's importance-sampling approximation regardless of whether the model is correctly specified. **The `loo_compare()` numeric result should be read as suggestive, not decisive** — a trustworthy comparison would need `kfold()` cross-validation (estimated ≈2 hours per model given the ~12-minute parallel fit time established in Step 4, not run here) or `reloo()` (infeasible: 363 of 570 points exceed k = 0.7). This is flagged as a follow-up rather than silently accepted or discarded.

**Bottom line for the response letter:** combining the (unreliable but directionally consistent) LOO result with the individually-null, well-converged interaction coefficients, the honest conclusion is that **this dataset gives no evidence that herbivory modifies coral's recovery rate or heat-stress susceptibility** — not that it definitely doesn't (the same zero-inflation/power limitations noted for E1-E3 in Step 4 apply here too), but there is no signal pointing toward Task 5.2's hypothesized mechanism either. Task 5.3 (residual autocorrelation check) remains the one undone piece of Section 5's core model-fitting tasks.

### Step 4c — Task 5.3: residual temporal autocorrelation check

**Status:** done. **Script:** `Revision/R/04_benthic_gompertz.R` (Section H/I, appended). **Outputs:** `Revision/Output/acf_lag1_tests.csv` (36 Ljung-Box tests: 6 sites x 6 series), `Revision/Output/fig_acf_benthic_gompertz.png`, updated `qc_benthic_gompertz.md`.

**What it does:** extracts the additive model's `Site:Year` random intercepts (muCoral/muAlgae/muCCA) and computes Pearson residuals (observed vs. fitted count, binomial-approximation standardisation, since `brms::residuals()` doesn't support "pearson" for multinomial/categorical families) averaged across the 5 transects per site-year; computes the ACF (lags 1-3) per site for both diagnostics; and runs a lag-1 Ljung-Box test per site per series (36 tests total) with a Benjamini-Hochberg correction across all 36.

**Result: no statistically defensible residual autocorrelation, matching Task 5.3's expectation that the lagged ALR terms would absorb it.** Only 2 of 36 raw tests reach p<0.05 (chance expectation under the null is ~1.8) and **0 of 36 survive BH correction**. Two sites (LTER_3, LTER_6) show a visually consistent negative lag-1 ACF (-0.3 to -0.48) across all 6 series — worth a passing mention since the pattern repeats across otherwise-distinct series at the same two sites — but this does not rise to statistical significance once corrected for the 36 comparisons performed, and the plot confirms every bar sits within the (uncorrected) white-noise band. The `ar()`/`gp()` sensitivity checks the plan offers as a fallback for this scenario are therefore not triggered.

**Section 5's core model-fitting tasks (5.1-5.4) are now complete.** Remaining for Section 5 only if revisited later: the fuller counterfactual coral-proportion trajectories and "years to recovery" curve from Task 5.4's first two bullets (the E1-E3 average marginal effects and joint counterfactual were already done in Step 4). Next up per the plan's execution order (Section 13): Section 6 (hierarchical PI-curve metabolism model).

### Step 6 — Hierarchical PI-curve metabolism model (Section 6, Task 6.1; Tasks 6.2-6.4 deferred)

**Status:** Task 6.1 done, WITHOUT the `fish_resp_z` covariate (Task 7.1, fish bioenergetics, not yet built — Section 13's execution order lists Step 6 as depending on Step 5/fish; this is a known, flagged gap, not an oversight). Tasks 6.2 (per-estimand refits), 6.3 (uncertainty propagation), 6.4 (formal sample-size writeup) deferred. **Script:** `Revision/R/06_metabolism.R`. **Outputs:** `Revision/Data/derived/pi_model_hourly.csv`, `Revision/Output/fits/pi_model_dev.rds` (25% of days), `Revision/Output/fits/pi_model_full.rds` (all 117 days), `Revision/Output/qc_metabolism.md`.

**What it does:** builds LTER_1's coral-first pivot (ilr) coordinates of benthic composition per year (Hron et al. 2012), joins day-to-day temperature (converted to an MTE-style `invkT_c` term), flow, nitrogen, and season onto the hourly flux data, and fits the plan's nonlinear PI-curve model (`PP ~ (exp(la)*exp(lP)*PAR)/(exp(la)*PAR+exp(lP)) - exp(lR)`, Student-t family, correlated year- and day-level random effects shared between the `lP` and `lR` nonlinear parameters). Per the plan's own compute guidance, developed first on a 25%-of-days, Year-stratified subsample (624 rows) before running the full fit (2808 rows, 117 days, 17 years).

**A likely sign error in the plan's own prior was found and corrected before the full fit, with PI sign-off.** The plan states the coefficient on `invkT_c` for `lR` is `-E` (activation energy), but then specifies a prior of `normal(+0.65, 0.2)` for that coefficient — inconsistent with its own formula, since activation energy is conventionally positive (~0.65 eV) for a process (respiration) that speeds up with warming, which requires the *coefficient itself* to be negative. Confirmed both algebraically (`ln B(T) = ln B(T_ref) - E · invkT_c`) and empirically: the dev fit under the plan's literal `+0.65` prior produced a posterior of `lR_invkT_c` = 0.38 [0.03, 0.74] — positive, implying respiration *decreases* with warming, the biologically backwards direction. Corrected to `normal(-0.65, 0.2)` (chosen over keeping the literal prior or running all three variants); the dev fit under the corrected prior gave `lR_invkT_c` = -0.72 [-1.08, -0.37], now biologically sensible.

**Both fits converged cleanly** (dev: all R-hat = 1.00, ~7 minutes; full: all R-hat ≤ 1.01, ~52 minutes — compile time plus 4.5× more rows than the dev fit).

**Headline results from the full fit (2808 hourly rows, 117 days, 17 years), full detail in `qc_metabolism.md`:**

1. **Respiration and photosynthetic capacity respond to day-to-day temperature in opposite directions.** `lR_invkT_c` = -1.17 [-1.50, -0.83] (respiration increases with temperature, implied activation energy ≈1.17 eV, same direction as but larger than the canonical ~0.65 eV MTE value). `lP_invkT_c` = 1.93 [1.33, 2.56] (photosynthetic capacity *decreases* with temperature — consistent with thermal impairment of the photosynthetic apparatus). Framed together this is a "rising respiratory cost, falling photosynthetic capacity" pattern under warming, worth presenting as a headline physiological finding rather than a background covariate.
2. **The Pmax-Rd covariance now lives exactly where Task 4.1 said it should.** `cor(lP_Intercept, lR_Intercept)` at the Year:DielDate level = 0.48 [0.31, 0.63] — clearly credible and positive, confirming the decision to model this relationship as a residual correlation inside the PI model rather than a causal `Rd → Pmax` DAG edge.
3. **The original year/season confound was real, not hypothetical.** `lP_SeasonWinter` = 0.45 [0.27, 0.64], `lR_SeasonWinter` = 0.27 [0.16, 0.38], both credible.
4. **Flow is strongly associated with both Pmax and Rd** (`log_flow_c`: 0.89 [0.76, 1.02] for lP, 0.52 [0.44, 0.60] for lR), consistent with known boundary-layer effects on in-situ flux measurements.
5. **Coral cover (ilr1) is associated with higher respiration** (`lR_ilr1` = 0.32 [0.13, 0.50]); no credible compositional effect on Pmax was detected, and nitrogen's effect on Pmax (`lP_N_z` = -0.20 [-0.42, 0.09]) is not credible.

**Not yet done:** Task 6.2 (per-estimand covariate-restricted refits for E4/E6/E7/E9/E10), Task 6.3 (propagating benthic-composition and fish-biomass posterior uncertainty via `brm_multiple()`), Task 6.4 (formal sample-size-transparency writeup — raw numbers already recorded above), and refitting once Section 7's `fish_resp_z` exists. Section 7 (fish bioenergetics) is the natural next step, both because it unblocks this model's remaining gap and because it's independently required for estimands E11/E12.

### Step 6b — Task 6.2: per-estimand screening refits (E4, E6, E7, E9/E10)

**Status:** done, as SCREENING fits (25%-of-days dev subsample, not the full 2808-row dataset — see compute note below). **Script:** `Revision/R/06_metabolism.R` (Section E, appended). **Outputs:** `Revision/Output/fits/pi_E4_dev.rds`, `pi_E6_dev.rds`, `pi_E7_dev.rds`, `pi_E9_E10_dev.rds`, updated `Revision/Output/qc_metabolism.md` (Sections G-H).

**Adjustment sets used (from `table_s_dag.csv`, post Step 3b's edge additions):** for each of these five estimands, the smallest minimal adjustment set happens to be identical across both DAG variants (`Rd -> N` vs `Benthos -> N`), so one model per estimand is reported rather than two. E9 and E10 resolve to the exact same covariate list ({Benthos, DayTemp, Flow, N}) once their respective exposures are folded in, so they share a single fitted model (the same legitimate coincidence already seen for E1-E3 in Section 5). New covariates needed beyond what Task 6.1 already had: `Herb_z`/`Corall_z`/`OtherFish_z` (LTER_1 fish biomass by trophic group, year-level, z-standardised — these are the DAG's `Herb`/`Corall`/`OtherFish` nodes, **not** Task 7.1's bioenergetics-derived `fish_resp_z`, which remains unbuilt) and `DHW_z` (z-standardised `DHW_max` at LTER_1).

**Compute decision:** all four fits ran on the same 624-row, Year-stratified 25%-of-days subsample already validated for Task 6.1's dev fit — the full Task 6.1 fit alone took ~52 minutes, and four more full-data fits would cost several hours of sequential compute. All four converged cleanly (R-hat ≈ 1.00, zero divergent transitions) in a combined ~28 minutes. **These are reported as screening results; scaling any of them to the full dataset is a flagged follow-up, not done automatically.**

**Results, with the DAG-relevant (exposure) coefficient for each estimand:**

| Estimand | Exposure coefficient(s) | Estimate [95% CrI] | Credible? |
|---|---|---|---|
| E4: Benthos → Rd (direct, not via fish) | `lR_ilr1/ilr2/ilr3` | 0.26 [-0.15, 0.66]; -0.11 [-0.44, 0.23]; -0.04 [-0.31, 0.22] | No — all include 0 |
| E6: Fish → Rd (direct) | `lR_Herb_z/Corall_z/OtherFish_z` | -0.12 [-0.36, 0.12]; 0.07 [-0.13, 0.27]; 0.06 [-0.13, 0.26] | No — all include 0 |
| E7: DayTemp → Rd (direct) | `lR_invkT_c` | -0.68 [-1.05, -0.31] | **Yes** |
| E9: Benthos → Pmax (direct) | `lP_ilr1/ilr2/ilr3` | 0.04 [-0.36, 0.42]; 0.08 [-0.27, 0.47]; 0.08 [-0.23, 0.42] | No — all include 0 |
| E10: DayTemp → Pmax (direct) | `lP_invkT_c` | 0.62 [-0.19, 1.47] | No — includes 0 |

**Two findings worth flagging before any of these are treated as final:**

1. **E7 robustly replicates Task 6.1's headline temperature-respiration finding** using its own, much smaller DAG-derived adjustment set ({DHW, Flow} only, no Benthos/fish/N/Season needed) — `lR_invkT_c` = -0.68 [-1.05, -0.31] here vs. -1.17 [-1.50, -0.83] in the full Task 6.1 model. Same sign, overlapping interpretation, smaller magnitude (plausibly just the smaller dev-subsample's wider uncertainty pulling the point estimate toward the prior, since the corrected prior for this coefficient is centred at -0.65).
2. **E10's Pmax-temperature result does NOT replicate Task 6.1's finding** — `lP_invkT_c` = 0.62 [-0.19, 1.47] here (not credible) vs. 1.93 [1.33, 2.56] in the full model (strongly credible). The DAG's minimal adjustment set for E9/E10 does not require conditioning on `Season`, unlike Task 6.1's full model, which does — since day-to-day temperature correlates with season, dropping `Season` plausibly inflates this coefficient's uncertainty even though the backdoor criterion says it isn't required for identification. **This discrepancy should be resolved by refitting E9/E10 on the full dataset before reporting either number as the estimand's answer** — it is flagged as the single highest-priority item among the "not yet done" list below, since it directly affects whether Task 6.1's headline "Pmax declines with warming" finding holds up once properly estimand-matched.
3. **E4 and E6 (fish-related) are both null in this screening fit**, consistent with (not contradicting) Task 6.1's inconclusive `N_z`/fish findings, but the small dev subsample limits how much weight to put on this either way — also flagged for a full-data refit before concluding fish's direct contribution to Rd is genuinely negligible.

**Not yet done:** full-dataset versions of all four per-estimand models (E9/E10 highest priority, per finding 2 above); Task 6.3 (propagating benthic-composition and fish-biomass posterior uncertainty); Task 6.4 (formal sample-size-transparency writeup); and `fish_resp_z`/Section 7 integration.

### Step 5 — Task 7.1: fish community respiration (mechanistic), N excretion, and fish_resp_z

**Status:** Task 7.1 done (respiration, Rd comparison, `fish_resp_z` built and tested). N excretion built only as a placeholder. Tasks 7.2 (E11, fish response models) and 7.3 (E12, fish → N) deferred. **Script:** `Revision/R/05_fish.R`. **Outputs:** `Revision/Data/derived/fish_respiration_transect.csv`, `fish_respiration_site_year.csv`, `fish_vs_measured_Rd.csv`, `Revision/Output/fig_fish_pct_of_Rd.png`, `Revision/Output/qc_fish.md`, plus two new cached sensitivity fits (`pi_E6_fishresp_dev.rds`, `pi_E4_fishresp_dev.rds`).

**Method decision:** `fishflux` (Schiettekatte et al. 2020, GitHub-only, needs a full species-trait reference database) was not installed — a GitHub package pull with a large reference dataset was judged a real risk given this project's repeated network-reliability issues (EDI/PASTA access failures in Steps 2/2b), and was not attempted blind. Used instead, as the plan's own stated alternative: the general teleost mass-scaling equation of **Clarke & Johnston (1999, *J. Anim. Ecol.* 68:893-905)**, `ln(Rb) = 0.80·ln(M) − 5.43` (Rb in mmol O₂/h, M in g wet mass), extended with a standard Arrhenius/MTE temperature correction (activation energy 0.6 eV, the midpoint of the 0.4-0.8 eV range reported for fish respiration across Barneche et al. 2014, Barneche & Allen 2018, and Brown et al. 2004). **This is explicitly an approximation using general-teleost literature parameters, not a reproduction of Barneche et al. (2014)'s species/trophic-group-specific fitted coefficients** — flagged for replacement with `fishflux`'s precise values before final submission if warranted. The 20°C reference temperature assumed for Clarke & Johnston's intercept is likewise an assumption, not independently confirmed against their paper.

**Headline result (the direct, quantitative answer Task 7.1 asks for) — smaller than the plan anticipated, reported as found:** mechanistic fish respiration accounts for a mean of **0.28% (range 0.07%–0.56%)** of measured ecosystem respiration at LTER_1 across the 18 years both series are available (2008–2025) — not the "few percent to tens of percent" the plan's draft guessed. Even generously correcting for field vs. standard metabolic rate (Barneche's documented activity-scope ratios of 1.2–3.2×), the estimate would still only reach ~0.2–1.7%. This is reported honestly rather than adjusted toward the a priori expectation; it is consistent with (not contradicted by) Rd being an open-flow flux that already integrates whatever fish respiration occurs within the sensor flow path, per Task 7.1's own framing — fish are evidently a small fraction of total community respiration at this backreef site, mechanistically estimated.

**`fish_resp_z` was built and tested directly against the "give E4/E6 more power" premise — honest finding: it does not clearly help, at screening scale.** Refitting E6 (Fish → Rd, direct) with `fish_resp_z` substituted for the biomass trio gives `lR_fish_resp_z` = -0.061 [-0.274, 0.150], no more credible or precise than the original biomass-based null result. Refitting E4 (Benthos → Rd, direct) with `fish_resp_z` added alongside the existing biomass adjustment leaves the Benthos coefficients essentially unchanged and `fish_resp_z` itself comes back with a wide, uninformative CI (-0.118 [-0.730, 0.485]). This makes sense rather than being a coding error: `fish_resp_z` is computed *from* the same biomass/count data as the trophic-group biomass covariates (via the allometric transform) — it is a re-weighting of the same information, not an independent measurement, so there is no a priori reason it should sharpen inference in a regression that already has biomass in it. **`fish_resp_z`'s real contribution is as a physically interpretable quantity directly comparable to measured Rd** (the headline 0.28%-of-Rd result above), not as a regression covariate with more statistical power than biomass — the premise behind this step's request doesn't hold up once tested, and is reported as such rather than reframed as a success.

**N excretion** was computed only as an illustrative placeholder (O:N atomic ratio = 20, not trophic-group-specific) — order-of-magnitude only, pending a proper excretion model in Task 7.3.

**Not yet done:** Task 7.2 (herbivore/corallivore biomass response models, E11, across all 6 sites with a residual-ACF check) and Task 7.3 (fish → N, E12, using the N excretion estimate above once it's properly trophic-group-specific). `fishflux` installation remains a recommended follow-up if more precise, species-specific respiration estimates are needed for final submission.

### Step 5b — Task 7.1 revised: species-specific fish respiration via fishflux/Barneche & Allen (2018)

**Status:** done, superseding Step 5's generic Clarke & Johnston (1999) approximation at the PI's request after `fishflux` and `rfishbase` were installed. **Script:** `Revision/R/05_fish.R` (fully rewritten). **Outputs:** same files as Step 5, now species-specific, plus `Revision/Data/derived/_cache/fishbase_{ecology,morphometrics,popgrowth,species}.csv` and `fishflux_family_metabolism.csv`.

**Initial blocker, then resolved:** `fishflux`'s own per-species wrapper functions (`trophic_level()`, `aspect_ratio()`, `growth_params()`) depend on `rfishbase`, which first hit a sandbox-level block on its remote data path (`huggingface.co/datasets/cboettig/fishbase`) — a permission boundary, not a flaky network issue, and not something attempted to route around. Reinstalling `rfishbase` alone did not fix it; the actual fix was upgrading `duckdbfs` from 0.1.0 (missing an exported `duckdb_config` function that the installed `rfishbase` 5.0.3 required) to 0.1.2. Once that dependency conflict was resolved, bulk (not per-species) calls to the underlying `rfishbase::ecology()`/`morphometrics()`/`popgrowth()` functions worked and were fast (each a single call for all 238 species, not a per-species loop — `fishflux`'s own wrappers loop one species at a time via these same functions, which would have been impractically slow).

**Method:** family- and trophic-level-specific resting metabolic rate following **Barneche & Allen (2018, *Ecology Letters*, doi:10.1111/ele.12947)** "Model 2", via `fishflux::metabolism()` + `metabolic_rate()`. Trophic level, caudal-fin aspect ratio, and von Bertalanffy growth parameters (K, Loo) come from bulk FishBase queries with a species → genus → family → global-mean fallback hierarchy (extending `fishflux`'s own species → genus-only fallback, since every individual record needs a value). Max body size and daily growth rate combine FishBase's K/Loo with this project's own survey-fitted length-weight relationship (FishBase's own `Weight` field is too sparse for reef fish to rely on). Activity scope (f=2) and respiratory quotient (RQ=0.8, for converting metabolic carbon loss to O₂ consumption) remain literature-typical constants, not species-specific — flagged as the main remaining approximations.

**Two real bugs caught and fixed during development** (documented in the script header so they aren't silently reintroduced): (1) `Total_Length` in `fish_ind.csv` is in millimetres while FishBase's `Loo` is in centimetres — combining them directly silently corrupted `m_max`/`growth_g_day` (wrong magnitude, no error), caught by sanity-checking a known species' size (a reef shark). (2) The global length-weight fallback coefficients were extracted via `lw_fit_global["a"]`/`["b"]`, but R's `coef()` preserves full term names, so this silently returned `NA` for 26 of 28,964 individual records — caught by checking for unexpected `NA` counts rather than assuming a clean run meant correctness.

**Result: the species-specific estimate is smaller than the generic approximation (roughly half) but tells the same story, which is itself a useful robustness check:**

| | Generic (Step 5) | Species-specific (Step 5b) |
|---|---|---|
| Mean % of Rd from fish, LTER_1 2008-2025 | 0.28% | **0.13%** |
| Range | 0.07%-0.56% | 0.03%-0.23% |
| `fish_resp_z` credible in E6 (Fish → Rd, direct)? | No | No |
| `fish_resp_z` credible in E4 (Benthos → Rd, direct) when added? | No | No |

The more rigorous, genuinely taxon-resolved model **confirms, rather than overturns, Step 5's two headline findings**: fish respiration is a small fraction of measured ecosystem respiration at this backreef site (now estimated at ~0.1%, even smaller than the already-modest generic estimate), and the mechanistic respiration estimate does not sharpen inference over raw trophic-group biomass in either E4 or E6. Cross-method agreement between a crude cross-species approximation and a genuinely species-specific model, arriving at the same qualitative conclusion, is stronger evidence than either method alone.

**Not yet done:** Task 7.2 (E11, fish response models) and Task 7.3 (E12, fish → N) remain deferred, as in Step 5. N excretion remains a non-species-specific stoichiometric placeholder (O:N = 20).
