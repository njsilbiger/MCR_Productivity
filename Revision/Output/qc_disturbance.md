# QC summary: disturbance covariates (Revision Step 2)

Generated: 2026-10-01 14:30:42.102456

## A. Survey-date window per year

Benthic survey dates are mostly early-to-mid January (austral summer, post-peak-heat-stress season), not April-May as an earlier draft of the plan assumed -- the only exception is 2005 (surveyed in late May). The `survey_date` (median survey date per year) anchors a trailing 364-day DHW window for that year. Note 2023 has a wide date range (2022-02-27 to 2023-01-24) -- flagged for the PI; not resolved further here since it does not block the DHW calculation (median still falls in Jan 2023).

| Year | survey_date (median) | min_date | max_date |
|---|---|---|---|
| 2005 | 2005-05-24 | 2005-05-20 | 2005-05-30 |
| 2006 | 2006-01-14 | 2006-01-10 | 2006-01-17 |
| 2007 | 2007-01-15 | 2007-01-08 | 2007-01-24 |
| 2008 | 2008-01-08 | 2008-01-08 | 2008-01-08 |
| 2009 | 2009-01-11 | 2009-01-11 | 2009-01-11 |
| 2010 | 2010-01-01 | 2010-01-01 | 2010-01-01 |
| 2011 | 2011-01-17 | 2011-01-10 | 2011-01-23 |
| 2012 | 2012-01-16 | 2012-01-10 | 2012-01-24 |
| 2013 | 2013-01-15 | 2013-01-08 | 2013-01-21 |
| 2014 | 2014-01-15 | 2014-01-07 | 2014-01-24 |
| 2015 | 2015-01-12 | 2015-01-04 | 2015-01-27 |
| 2016 | 2016-01-09 | 2016-01-03 | 2016-01-17 |
| 2017 | 2017-01-10 | 2017-01-03 | 2017-01-17 |
| 2018 | 2018-01-09 | 2018-01-03 | 2018-01-15 |
| 2019 | 2019-01-08 | 2019-01-04 | 2019-01-14 |
| 2020 | 2020-01-13 | 2020-01-04 | 2020-01-19 |
| 2021 | 2021-02-04 | 2021-01-31 | 2021-02-08 |
| 2022 | 2022-02-21 | 2022-02-13 | 2022-02-28 |
| 2023 | 2023-01-18 | 2022-02-27 | 2023-01-24 |
| 2024 | 2024-01-16 | 2024-01-07 | 2024-01-24 |
| 2025 | 2025-01-17 | 2025-01-06 | 2025-01-23 |

## B. Primary heat stress: satellite DHW (NOAA Coral Reef Watch)

- Source: NOAA Coral Reef Watch daily 5 km `CRW_DHW` (degree heating weeks), via ERDDAP (`coastwatch.pfeg.noaa.gov/erddap`, dataset `NOAA_DHW`), pulled for 2004-01-01 to 2025-03-01.
- Single island-wide pixel used for all 6 sites (see code comment for rationale) -- nearest grid node to the MCR LTER network centroid (-17.4909, -149.826) is -17.475, -149.825.
- Missing days in the pulled series: 5 of 7731.
- `DHW_max` per year is the maximum `CRW_DHW` in the 364 days up to and including that year's median benthic-survey date (Task 3.1's 12-month, survey-date-anchored window).

| Year | survey_date | DHW_max (deg C-weeks) | days in window |
|---|---|---|---|
| 2005 | 2005-05-24 | 0.00 | 365 |
| 2006 | 2006-01-14 | 0.00 | 365 |
| 2007 | 2007-01-15 | 0.00 | 365 |
| 2008 | 2008-01-08 | 0.32 | 365 |
| 2009 | 2009-01-11 | 0.00 | 365 |
| 2010 | 2010-01-01 | 0.00 | 365 |
| 2011 | 2011-01-17 | 0.00 | 365 |
| 2012 | 2012-01-16 | 0.00 | 365 |
| 2013 | 2013-01-15 | 0.46 | 365 |
| 2014 | 2014-01-15 | 0.46 | 365 |
| 2015 | 2015-01-12 | 0.00 | 365 |
| 2016 | 2016-01-09 | 0.32 | 365 |
| 2017 | 2017-01-10 | 0.66 | 365 |
| 2018 | 2018-01-09 | 0.00 | 365 |
| 2019 | 2019-01-08 | 0.00 | 365 |
| 2020 | 2020-01-13 | 3.36 | 365 |
| 2021 | 2021-02-04 | 0.29 | 365 |
| 2022 | 2022-02-21 | 0.00 | 365 |
| 2023 | 2023-01-18 | 0.00 | 365 |
| 2024 | 2024-01-16 | 0.00 | 365 |
| 2025 | 2025-01-17 | 2.30 | 360 |

The 2020 survey (window ending Jan 2020) shows DHW_max = 3.36, consistent with the known 2019 Moorea bleaching event (Section 0.3 of the plan) -- a sanity check that the window logic is capturing real heat-stress history rather than an artefact.

## C. Secondary / sensitivity heat stress: in-situ logger DHW

- Logger: LTER_2 backreef, 2 m depth (the only backreef logger in the temperature file; applied here to all 6 sites, same limitation as the satellite series). Daily coverage: 7066 days, 2005-05-31 to 2025-07-21 (291 missing days).
- **This is an approximation, not NOAA's official in-situ DHW product.** MMM (maximum monthly mean) is proxied as the maximum calendar-month mean of satellite `sea_surface_temperature` at the same Moorea pixel across the full 1985-2025 record (28.88 deg C, month 3), because NOAA's official MMM baseline raster is not served on this ERDDAP instance. HotSpot = max(0, logger_temp - MMM); DHW = trailing 12-week (84-day) sum of HotSpot/7 for days with HotSpot >= 1 deg C, following the standard CRW accumulation rule applied to the logger series instead of satellite SST.
- Annual `DHW_logger_max` uses the same trailing-364-day-to-survey-date window as the satellite metric, for comparability.

| Year | DHW_logger_max (deg C-weeks) | logger days in window |
|---|---|---|
| 2005 | NA | 0 |
| 2006 | 0.00 | 229 |
| 2007 | 0.00 | 252 |
| 2008 | 0.00 | 204 |
| 2009 | 0.00 | 365 |
| 2010 | 0.00 | 365 |
| 2011 | 0.00 | 365 |
| 2012 | 0.00 | 365 |
| 2013 | 0.00 | 365 |
| 2014 | 0.00 | 364 |
| 2015 | 0.00 | 365 |
| 2016 | 0.00 | 365 |
| 2017 | 1.07 | 365 |
| 2018 | 0.00 | 342 |
| 2019 | 0.00 | 365 |
| 2020 | 4.08 | 365 |
| 2021 | 0.92 | 365 |
| 2022 | 0.00 | 365 |
| 2023 | 0.00 | 365 |
| 2024 | 0.00 | 365 |
| 2025 | 1.71 | 365 |

Treat this strictly as a sensitivity check (Task 3.1): it uses a locally-derived MMM proxy, not NOAA's validated baseline, and -- like the satellite series -- is a single island-wide/single-logger value applied to all 6 sites, so it adds local temperature variability but not real between-site heat-stress contrast.

## D. Crown-of-thorns starfish (COTS)

- Source: `MCR_LTER_COTS_abundance_2005-2025_20250310.csv` (manually downloaded by the PI from EDI; automated access failed -- see testing log). Citation: Moorea Coral Reef LTER and A. Brooks. 2026. MCR LTER: Coral Reef: Long-term Population Dynamics of *Acanthaster planci*, ongoing since 2005 ver 13. Environmental Data Initiative. https://doi.org/10.6073/pasta/0601a5aa24c8f35fda90a99b4f1a50bd.
- Raw rows (all habitats, all sites): 1512. Backreef rows kept: 504.
- Design confirmed: 6 sites x 21 years (2005-2025) x 4 backreef transects = 504 expected rows; observed 504.
- Density = total COTS counted / (n transects x 250 m^2), where 250 m^2 is the 5 x 50 m belt-transect area used for the MCR COTS survey (same geometry as the annual fish survey's wide swath).
- Site-year COTS density range: 0.0000 to 0.0120 ind/m^2.

| Year | total backreef COTS (6 sites) | max site density (ind/m^2) |
|---|---|---|
| 2005 | 5 | 0.0030 |
| 2006 | 3 | 0.0010 |
| 2007 | 2 | 0.0010 |
| 2008 | 6 | 0.0060 |
| 2009 | 12 | 0.0120 |
| 2010 | 5 | 0.0020 |
| 2011 | 3 | 0.0020 |
| 2012 | 2 | 0.0010 |
| 2013 | 7 | 0.0050 |
| 2014 | 0 | 0.0000 |
| 2015 | 0 | 0.0000 |
| 2016 | 0 | 0.0000 |
| 2017 | 0 | 0.0000 |
| 2018 | 1 | 0.0010 |
| 2019 | 0 | 0.0000 |
| 2020 | 0 | 0.0000 |
| 2021 | 0 | 0.0000 |
| 2022 | 0 | 0.0000 |
| 2023 | 1 | 0.0010 |
| 2024 | 16 | 0.0080 |
| 2025 | 0 | 0.0000 |

### D1. Cross-habitat check (forereef, not used as a covariate)

Pulled for comparison only, to check the backreef pattern against the published (largely forereef) outbreak history -- not joined into `disturbance_site_year.csv`.

| Year | Total COTS: backreef (6 sites) | Total COTS: forereef (6 sites) |
|---|---|---|
| 2005 | 5 | 5 |
| 2006 | 3 | 0 |
| 2007 | 2 | 13 |
| 2008 | 6 | 108 |
| 2009 | 12 | 70 |
| 2010 | 5 | 5 |
| 2011 | 3 | 0 |
| 2012 | 2 | 0 |
| 2013 | 7 | 0 |
| 2014 | 0 | 0 |
| 2015 | 0 | 0 |
| 2016 | 0 | 0 |
| 2017 | 0 | 0 |
| 2018 | 1 | 0 |
| 2019 | 0 | 0 |
| 2020 | 0 | 0 |
| 2021 | 0 | 0 |
| 2022 | 0 | 2 |
| 2023 | 1 | 17 |
| 2024 | 16 | 56 |
| 2025 | 0 | 2 |

**Confirms two real outbreak pulses, both much larger on the forereef:** (1) the documented 2006-2010 outbreak (Kayal et al. 2012) is clearly visible on the forereef (108 individuals in 2008, 70 in 2009) and only weakly visible on the backreef (6 in 2008, 12 in 2009) -- consistent with COTS outbreaks being a largely forereef phenomenon at Moorea, so the muted backreef signal is expected, not a data problem. (2) A second, previously undiscussed outbreak pulse appears in 2023-2024: forereef totals rise to 17 (2023) and 56 (2024), with the backreef pulse (1 in 2023, 16 in 2024; concentrated at LTER_2 and LTER_3) tracking the same timing at smaller magnitude. The 2024 backreef pulse (plan testing log, Step 2 item 6) is therefore corroborated by the forereef data as a real, synchronised outbreak event -- not an artefact -- and should be written up as a second documented Moorea COTS outbreak (2023-2024) alongside 2006-2010 in the methods/disturbance-history text.

## E. Cyclone Oli (Feb 2010)

- Closest approach to Moorea taken as 2010-02-16. `Cyclone = 1` for the first survey taken within 365 days after that date, else 0.
- Because the 2010 survey predates Oli (survey_date = 2010-01-01) and the 2011 survey (2011-01-17) falls 335 days after it, the indicator flags **2011**, not 2010, as the cyclone-affected survey year. Flag this for the PI: it changes which annual benthic transition the Section 5 Gompertz model should attribute to the cyclone, relative to the original (uncyclone-indicator) analysis.
- Backreef wave exposure during Oli was much lower than forereef exposure (Section 0.2/3.3 of the plan); this binary indicator does not yet distinguish shore-facing exposure across the 6 sites.

| Year | Cyclone |
|---|---|
| 2005 | 0 |
| 2006 | 0 |
| 2007 | 0 |
| 2008 | 0 |
| 2009 | 0 |
| 2010 | 0 |
| 2011 | 1 |
| 2012 | 0 |
| 2013 | 0 |
| 2014 | 0 |
| 2015 | 0 |
| 2016 | 0 |
| 2017 | 0 |
| 2018 | 0 |
| 2019 | 0 |
| 2020 | 0 |
| 2021 | 0 |
| 2022 | 0 |
| 2023 | 0 |
| 2024 | 0 |
| 2025 | 0 |

## F. Assembled table

- `disturbance_site_year.csv`: 126 rows (6 sites x 21 years).
- Columns: Site, Year, DHW_max (satellite, island-wide), DHW_logger_max (in-situ sensitivity, island-wide), COTS_count_total, COTS_density_m2 (site-year), Cyclone (0/1, island-wide).
- DHW_max and Cyclone are identical across all 6 sites for a given year by construction (Sections B, C, E) -- only COTS varies by site. This is an explicit limitation to carry into Section 4's DAG/adjustment-set work and the Section 5 methods text: between-site contrasts in the Gompertz model for DHW/Cyclone effects come entirely from between-site variation in starting benthic composition and COTS, not from independent exposure.

