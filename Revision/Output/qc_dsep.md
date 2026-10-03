# QC summary: Task 4.3 d-separation test (Revision Step 3b)

Generated: 2026-10-01 14:57:35.930536

## A. Benthos ilr coordinates

- Zero cells in the Coral/Algae/CCA/Other count matrix (120 site-years x 4 parts): 19 of 480. Replaced via `zCompositions::cmultRepl` (multiplicative simple replacement for compositional count data), matching the plan's Task 5.1 approach, before closing to proportions and taking ilr.
- ilr basis: `compositions::ilr()` default sequential binary partition on (Coral, Algae, CCA, Other). Three coordinates (`ilr1`,`ilr2`,`ilr3`) jointly represent the Benthos node -- tested as a 3-dimensional vector via `cis.pillai` (canonical correlation), not reduced to one scalar.

## B. Fish trophic groups (site-year)

- `Herb` = Herbivore biomass (g/m^2, mean across the 4 transects per site-year). `Corall` = Corallivore biomass. `OtherFish` = Planktivore + Invertivore + Piscivore + Omnivore + Other biomass, summed per transect then averaged across transects -- the DAG's catch-all fish node.
- Rows: 120 site-years.

## C. Nitrogen node

- `N` = Turbinaria tissue %N, site-year mean (from Step 2b). 108 rows.

## D. Metabolism PROXIES (LTER_1 only -- NOT the Section 6 PI-curve model)

- `Rd_proxy` = -mean(daily night-time R_mean); `Pmax_proxy` = mean(daily GP_mean); `DayTemp_proxy`/`Flow_proxy` = annual means of `pp_day.csv`'s daily means; `Season_frac_summer` = fraction of that year's complete diel days classified Summer. 18 annual rows, LTER_1 only, 2008-2025 (with gaps -- see Step 1's QC for the season-coverage imbalance by year).

## F. Analysis panels

- `panel_full`: 120 site-year rows (6 sites), used for implied independencies not involving a metabolism node.
- `panel_lter1`: 18 annual rows, LTER_1 only, used only for independencies involving Rd_proxy/Pmax_proxy/DayTemp_proxy/Flow_proxy/Season_frac_summer. This is a small sample -- treat any test run on this panel as low-powered, per Section 0.3 of the plan.

## G. D-separation test results (Task 4.3)

- Table S-dsep (missing-edge tests): 442 rows total, 440 testable, 15 violated at BH-adjusted p<0.05, 2 not testable (insufficient df).

- Shipley/Fisher's C [Rd -> N]: C = 18.05, df = 14 (k = 7 basis-set tests used), global p = 0.2046
- Shipley/Fisher's C [Benthos -> N]: C = 18.05, df = 14 (k = 7 basis-set tests used), global p = 0.2046

**Fisher's C note:** the Shipley/Fisher's C global test was run TWICE: once against the original DAG (before this step's edits) and once after adding the four candidate edges below. Before: C = 58.03, df = 14, p < 0.0001 (both variants) -- an overwhelming rejection. After: **C = 18.05, df = 14, p = 0.2046 (both variants)** -- no longer significant. This drop (58 -> 18) is a meaningful, consistent signal that the four added edges captured real missing structure, not noise. That said, `impliedConditionalIndependencies(type="basis.set")` tests each node against the JOINT block of ALL its non-descendant/non-parent nodes at once -- several of these blocks have 7-10 variables tested at n=17-18 (the LTER_1-only metabolism panel), a near-saturated regime, and 10 of 17 basis-set tests were not even computable. Treat the *passing* p=0.21 result as corroborating, not decisive, evidence; the pairwise Table S-dsep tests below remain the primary basis for the findings.

### What was found, decided, and the net effect of the DAG update

Four edges were added to BOTH DAG variants in `03_dag.R` after review with the PI (hard stop, Section 13 item 2): `Time -> Benthos_lag`, `Time -> Herb_lag`, `Time -> N_lag`, `Time -> COTS`. Rationale and effect:

1. **`Benthos_lag`, `Herb_lag` and `N_lag` were each marginally correlated with `Time`** (r = 0.67, 0.56, -0.57; all p < 1e-9, n = 107-114, `panel_full`) in the original DAG -- expected, since these lag nodes are literally last year's `Benthos`/`Herb`/`N`, which the DAG already connects to `Time`, but the lagged versions had no such edge. **Adding `Time -> Benthos_lag/Herb_lag/N_lag` fully resolved the `Benthos_lag`-`N_lag` violation** (p_adj_BH 0.003 -> 0.96) but **only partially resolved `Benthos_lag`-`Herb_lag`**, which remains violated after conditioning on `Time` (r = 0.37, p_adj_BH = 0.03-0.04, both variants) -- a genuine residual relationship between lagged coral cover and lagged herbivore biomass beyond the shared secular trend, not resolved by this edit. Candidate explanations: unmeasured site-level habitat quality, or this pair may be better treated like Rd/Pmax (finding 4 below) -- two correlated initial conditions handled via a residual correlation in whatever model uses them, rather than a new causal arrow. Left open rather than forcing a specific direction.
2. **`COTS` was correlated with `N`/`N_lag`** (no path existed in the original DAG, and `COTS` had no connection to `Time` at all), consistent with Step 2's finding that COTS outbreaks cluster in specific years (2008-2009, 2023-2024) rather than purely tracking lagged local coral cover. **Adding `Time -> COTS` did NOT resolve this** -- `COTS`-`N`/`N_lag` remain violated even conditioning on `Time` (r = 0.30-0.50, p_adj_BH < 0.05, both variants). This points to a real, unexplained relationship between COTS density and nitrogen, beyond a shared time trend. One literature-supported candidate mechanism is the nutrient-enrichment hypothesis for COTS outbreaks (elevated nutrients favouring larval survival; see e.g. Fabricius et al. 2010) -- i.e. a possible `N -> COTS` edge. **Not added here**: this would be a substantive new causal claim requiring its own literature check and PI sign-off, not something to infer from one significant correlation. Flagged as an open item for the next DAG review, not a hard requirement before Section 5/7 modelling proceeds.
3. **`COTS` is correlated with `Pmax` at LTER_1** (r = -0.70 to -0.81, p_adj_BH < 0.05, n = 17-18) -- unaffected by this step's edits (COTS has no path to Pmax in either DAG variant, with or without `Time -> COTS`). Lower confidence than (1)-(2): single site, small n, and `Pmax` here is a raw daily-mean PROXY, not the real Section 6 model output. Worth re-checking once real Pmax estimates exist rather than treated as a finding now.
4. **`Pmax` and `Rd` remain correlated given their shared measured causes** (r = 0.76-0.78, p_adj_BH < 0.05, n = 17-18) -- unaffected by, and unrelated to, this step's edits. This is NOT a new problem: it is exactly what Task 4.1 anticipated when it deliberately removed the `Rd -> Pmax` edge and assigned their shared variance to a residual correlation inside the Section 6 PI-curve model instead of a DAG arrow. This result confirms that choice was necessary.

**Net effect of the four added edges:** violated missing-edge tests dropped from 25/440 to 15/440 (BH-adjusted p<0.05), and the Shipley Fisher's C global test went from strongly rejecting both DAG variants (p<0.0001) to not rejecting them (p=0.20). Two genuine open items remain (`Benthos_lag`-`Herb_lag`, `COTS`-`N`) that were not resolved by this edit and were deliberately left as flagged findings rather than prompting further ad hoc edges without a stated mechanism.

