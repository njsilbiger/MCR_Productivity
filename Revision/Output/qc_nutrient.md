# QC summary: nutrient (N) covariate (Revision Step 2b)

Generated: 2026-10-01 14:44:06.023719

## A. Two candidate sources

| | WaterColumnN.csv (dissolved N+N) | Macroalgal CHN file (tissue %N) |
|---|---|---|
| Sites | LTER 1 Backreef Water Column | LTER_1, LTER_2, LTER_3, LTER_4, LTER_5, LTER_6 |
| Years | 2005-2018 | 2005-2025 |
| Rows | 166 | 3104 |

`WaterColumnN.csv` is single-site (LTER_1) and stops in 2018 -- it misses the 2019 bleaching event, the 2020 DHW peak, and the entire 2023-2024 COTS outbreak documented in Step 2's QC. **Decision: use macroalgal tissue %N as the covariate; keep water-column N+N only as a documentation/validation cross-check (Section B), reproducing the original qmd's approach.**

## B. Validation: does tissue %N track dissolved N+N at LTER_1?

Reproduces the original qmd's check (`modN <- lm(Nitrite_and_Nitrate ~ N_percent)`), restricted to the years both series overlap (2007-2018, n = 12 years).

- Pearson r = 0.630 (95% CI 0.087 to 0.884), p = 0.0282
- lm(Nitrite_and_Nitrate ~ N_percent): slope = 0.788, R^2 = 0.397

## C. Genus check: Sargassum vs Turbinaria tissue %N

| Genus | mean %N | sd %N | n samples | years |
|---|---|---|---|---|
| Sargassum | 1.144 | 0.244 | 165 | 2005-2014 |
| Turbinaria | 0.675 | 0.225 | 814 | 2007-2025 |

Sargassum tissue %N is ~70% higher than Turbinaria's on average (a ~2 SD difference), and Sargassum sampling stops after 2014 while Turbinaria runs 2007-2025 -- a genus effect confounded with a protocol-era effect. **Decision: Turbinaria only**, not pooled, matching the original qmd.

## D. Resulting tables

- `N_sample.csv`: 814 individual Turbinaria tissue samples (Backreef, all 6 sites, 2007-2025, as recorded).
- `N_site_year.csv`: 108 rows (6 sites x 18 years, 2007-2024). Missing site-years (no Turbinaria samples that year): 1.

| Site | Year | n samples |
|---|---|---|
| LTER_5 | 2007 | 0 |

Per-sample-count distribution (non-missing site-years): min = 5, median = 9, max = 10.

2005-2006 (no Turbinaria samples at all that early -- Sargassum-only or no CHN sampling) and 2025 (5 samples, 1 site -- processing lag, same pattern as other 2025-incomplete data in this project) are excluded from `N_site_year.csv`'s core 2007-2024 window but retained as recorded in `N_sample.csv`. Downstream models using `N_site_year.csv` will need an explicit missing-data strategy for the one remaining gap (LTER_5, 2007) within the core window -- left unresolved here, not imputed.

