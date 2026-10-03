# QC summary: data disaggregation (Revision Step 1)

Generated: 2026-10-01 12:30:30.210834

This report documents the checks run while building disaggregated benthic, fish, and metabolism tables from the raw MCR LTER files, and the assumptions those tables rely on. See Section 2 of `Revision/Reviewer_Response_Plan.md` for the rationale.

## A. Benthic cover

- Raw backreef rows (taxon x quadrat x year): 63713

### A1. Points per quadrat (two independent checks)

Method 1: smallest non-zero Percent_Cover increment per year (coarse, whole-year signal). Method 2: per-quadrat GCD of non-zero cover values, restricted to quadrats with >= 4 taxa so the GCD is diagnostic rather than a degenerate artefact of sparse composition, then taking the modal value across quadrats each year.

| Year | Method 1: N pts/quadrat | Method 2: modal N pts/quadrat (reliable quadrats, % agreeing) |
|---|---|---|
| 2006 | 25 | 25 (n=59, 100.0%) |
| 2007 | 25 | 25 (n=54, 100.0%) |
| 2008 | 25 | 25 (n=119, 100.0%) |
| 2009 | 25 | 25 (n=122, 100.0%) |
| 2010 | 25 | 25 (n=136, 100.0%) |
| 2011 | 25 | 25 (n=151, 100.0%) |
| 2012 | 25 | 25 (n=183, 100.0%) |
| 2013 | 25 | 25 (n=163, 100.0%) |
| 2014 | 25 | 25 (n=192, 99.5%) |
| 2015 | 25 | 25 (n=182, 100.0%) |
| 2016 | 25 | 25 (n=212, 100.0%) |
| 2017 | 25 | 25 (n=132, 99.2%) |
| 2018 | 25 | 25 (n=99, 100.0%) |
| 2019 | 25 | 25 (n=101, 100.0%) |
| 2020 | 50 | 50 (n=196, 96.4%) |
| 2021 | 25 | 25 (n=117, 100.0%) |
| 2022 | 25 | 25 (n=109, 100.0%) |
| 2023 | 25 | 25 (n=99, 100.0%) |
| 2024 | 25 | 25 (n=103, 100.0%) |
| 2025 | 25 | 25 (n=106, 100.0%) |

**CONFIRMED (2026-10-01):** both independent methods agree exactly: every year uses 25 points/quadrat (4% increments) *except* **2020, which uses 50 points/quadrat** (2% increments). Method 2, restricted to the 2635 quadrats with enough taxa to be diagnostic, agrees with the modal value at 96-100% of quadrats in every single year including 2020 (96.4%), ruling out a data-entry coincidence. This is a previously undocumented change in point-count density for 2020 only; the most likely explanation is a COVID-19-era change in photoquadrat image-analysis protocol (2020 metabolism sampling was also summer-only that year, consistent with disrupted fieldwork — see Section C). `benthic_quad.csv` uses this year-specific lookup (`n_pts_inferred`) as the `trials()` denominator for the Section 5 multinomial model. Recommended before submission: cross-check against the MCR LTER data-package version history/EDI metadata changelog for a documented 2020 protocol note, but the within-data evidence here is sufficient to proceed.

### A2. Quadrat design consistency

- Every Site x Year combination has exactly 50 quadrats (5 transects x 10 quadrats): **TRUE**
- All 6 backreef sites (LTER_1-6) are present every year, giving ~120 site-years of benthic data versus the 20 annual points used in `Year_Averages`.

### A3. Functional-group recoding and closure check

- Parts: Coral, Algae (fleshy macroalgae + turf, same taxon list as the original analysis), CCA (Crustose Corallines, broken out separately), Other (everything else).
- Quadrats whose Coral+Algae+CCA+Other does not sum to 100% (+/-0.01): 21 of 6000 (0.35%)
- By year:

| Year | Quadrats with total != 100% |
|---|---|
| 2010 | 1 |
| 2019 | 1 |
| 2020 | 17 |
| 2021 | 2 |

These are retained in `benthic_quad.csv` with their observed `n_total` (not forced to the nominal point count), so the multinomial model in Section 5 of the plan can use `trials(n_total)` rather than assuming a fixed denominator.

### A4. Resulting tables

- `benthic_quad.csv`: 6000 rows (quadrat x year), 6 sites x 20 years.
- `benthic_transect.csv`: 600 rows (transect x year).
- `benthic_site.csv`: 120 rows (site x year) — compare to the 20 rows in the original `Year_Averages` (LTER_1 only).

## B. Fish biomass

- Raw backreef fish rows (all 6 sites): 29089
- Rows removed as large-shark outliers (Biomass > 8000 g): 63
- Rows removed as missing-value code (Biomass < 0): 62
- Rows retained: 28964 (99.57% of raw)

### B1. Swath / area derivation

- Distinct swath widths recorded (m): 1, 5
- Every Site x Year x Transect combination records both swath widths: **TRUE** (design is constant across the time series)
- Area surveyed per transect = 50 m x (1 + 5) m = 300 m^2.
- With 4 transects/site, this reproduces the original code's hard-coded denominator of 1200 m^2/site (1200). It is now computed from the data rather than hard-coded, and is kept at transect resolution (300 m^2 per transect) so biomass density can be estimated per transect rather than only pooled to the site level.

### B2. Trophic grouping (confirmed 2026-10-01 — see plan Section 14)

Fine_Trophic == "Omnivore" is taxonomically distinct across its three Coarse_Trophic labels (checked via Family composition): the Planktivore-coded rows are overwhelmingly Pomacentridae and are grouped with Planktivore; the Primary/Secondary-Consumer-coded rows are a taxonomically mixed, ~6%-of-biomass group broken out as its own "Omnivore" category rather than folded into the much smaller, rarer "Other" catch-all (Fish Scale Consumer, Sediment Sucker, 2 unidentified/no-fish rows).

| trophic_group | Fine_Trophic values included |
|---|---|
| Corallivore | Corallivore |
| Herbivore | Browser, Brusher, Concealed Cropper, Cropper, Excavator, Herbivore/Detritivore, Scraper |
| Invertivore | Benthic Invertebrate Consumer |
| Omnivore | Omnivore |
| Other | Fish Scale Consumer, Sediment Sucker |
| Piscivore | Piscivore |
| Planktivore | Omnivore, Planktivore_exclusively |

### B3. Resulting tables

- `fish_transect.csv`: 480 transect-years x 7 trophic groups.
- `fish_ind.csv`: 28964 individual-count records retained (for Task 7.1 bioenergetics calculations).

## C. Ecosystem metabolism (in situ Lagrangian flux data)

- Raw hourly PP files read: 29
- Raw hourly rows: 3375
- Diel days seen in raw data: 157
- Diel days with all 24 hourly measurements (kept): 122 (77.7%)
- Resulting `pp_hour.csv`: 2928 hourly rows across 18 years.

### C1. Days per year x season

Confirms the season imbalance flagged in Section 0.3 of the plan (some years summer-only, some winter-only) — this is why `Season` must be a covariate in the PI-curve model (Section 6).

| Year | Winter | Summer |
|---|---|---|
| 2008 | 4 | 0 |
| 2009 | 2 | 3 |
| 2010 | 0 | 4 |
| 2011 | 1 | 0 |
| 2012 | 7 | 6 |
| 2013 | 6 | 0 |
| 2014 | 3 | 4 |
| 2015 | 5 | 2 |
| 2016 | 4 | 6 |
| 2017 | 4 | 5 |
| 2018 | 6 | 6 |
| 2019 | 6 | 6 |
| 2020 | 0 | 5 |
| 2021 | 3 | 0 |
| 2022 | 2 | 3 |
| 2023 | 3 | 3 |
| 2024 | 4 | 4 |
| 2025 | 0 | 5 |

- `pp_day.csv`: 122 day-level rows (one per complete diel cycle x deployment).

## D. Summary: disaggregated vs. original sample sizes

| Component | Original (`Year_Averages`) | Disaggregated |
|---|---|---|
| Benthic cover | 20 annual points (LTER_1 only) | 6000 quadrat-years (6 sites) / 120 site-years |
| Fish biomass | 20 annual points (LTER_1 only) | 3360 transect-year x trophic-group rows (6 sites) |
| Metabolism | 20 annual Pmax/Rd estimates (1 site) | 2928 hourly rows / 122 day-level rows (1 site, 122 complete diel cycles) |

Metabolism remains single-site, consistent with Section 0.3 of the plan: benthos-to-metabolism effects still rest on between-year variation at one site, but temperature-to-metabolism effects can now use within-year, day-to-day temperature variation once the hourly/daily tables are used instead of annual means.
