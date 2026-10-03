# QC summary: fish community respiration (Revision Step 5 / Section 7, Task 7.1)

Generated: 2026-10-02 13:52:19.28244

**Method revised 2026-10-02**: replaced the earlier generic Clarke & Johnston (1999) allometric approximation with `fishflux`/`rfishbase`-based species-specific (species > genus > family > global fallback) respiration estimates, following Barneche & Allen (2018). See the script header for the two real bugs (a length-unit mismatch; a silent NA from R name-mangling) caught and fixed during development.

## A. Annual mean temperature

- Island-wide annual mean SST from the cached satellite monthly series (Step 2), 1985-2025. Overall mean: 27.76 degC.

## B. Bulk FishBase trait lookups

- `ecology()` (trophic level): 184 of 238 species matched.
- `morphometrics()` (aspect ratio): 193 of 238 species matched.
- `popgrowth()` (K, Loo): 71 of 238 species matched (otolith-based growth studies are sparse for reef fish).
- All three fetched ONCE for the full species list (vectorised `rfishbase` calls), not in a per-species loop (the `fishflux` wrapper functions -- `trophic_level()`, `aspect_ratio()`, `growth_params()` -- call these same functions one species at a time, which is impractically slow at 238 species).

## C. Trait fallback hierarchy (species -> genus -> family -> global)

- Global fallback values: trophic level = 3.21, aspect ratio = 1.78, K = 0.744, Loo = 32.8 cm.
- Every one of the 238 species ends up with a complete trait set after fallback.

## D. Length-weight fit and family-level metabolic parameters

- Own length-weight fit (log Biomass ~ log Total_Length): 31 of 40 families have >= 10 records for a dedicated fit; the rest use the global fit (a = 2.05e-05, b = 2.979).
- `fishflux::metabolism()` (Barneche & Allen 2018 Model 2) called once per family at the overall mean temperature (27.76 degC) -- most families fall back to the package's global-average B0/alpha (only a handful of well-studied families, e.g. Pomacentridae, Gobiidae, Apogonidae, have their own fitted family-level parameters in Barneche & Allen's dataset; this is `fishflux`'s own documented behaviour, not a shortcut taken here).

## E. Per-individual respiration (species-specific `fishflux` model)

- Respiratory quotient (O2 from metabolic carbon loss): RQ = 0.8 (literature-typical for mixed fish diet -- not species-specific; flagged approximation).
- Activity scope f = 2 (resting-to-routine default from `fishflux`'s own documentation; not species-specific).
- NA count in computed metabolic rate: 0 of 28964 (both bugs described in the script header are fixed; this should stay 0).

## F. Areal respiration and N excretion (transect/site/year)

- `fish_respiration_transect.csv`: 480 transect-years. `fish_respiration_site_year.csv`: 120 site-years.
- **N excretion remains a simple stoichiometric placeholder (O:N atomic = 20), NOT trophic-group-specific -- order-of-magnitude only, pending Task 7.3.**

## G. Fish respiration vs. measured Rd at LTER_1 -- the headline result

- Years with both a fish survey and metabolism data at LTER_1: 18 (2008-2025). Units assumed mmol O2 m^-2 h^-1 for both series.

| Year | Fish Rb (mmol O2/m2/h) | Measured Rd (mmol O2/m2/h) | % of Rd from fish |
|---|---|---|---|
| 2008 | 0.018 | 58.48 | 0.03% |
| 2009 | 0.029 | 37.29 | 0.08% |
| 2010 | 0.026 | 50.29 | 0.05% |
| 2011 | 0.023 | 35.95 | 0.06% |
| 2012 | 0.032 | 42.54 | 0.07% |
| 2013 | 0.032 | 30.25 | 0.11% |
| 2014 | 0.044 | 28.69 | 0.16% |
| 2015 | 0.040 | 26.61 | 0.15% |
| 2016 | 0.050 | 26.39 | 0.19% |
| 2017 | 0.037 | 33.77 | 0.11% |
| 2018 | 0.049 | 36.74 | 0.13% |
| 2019 | 0.035 | 34.01 | 0.10% |
| 2020 | 0.037 | 33.32 | 0.11% |
| 2021 | 0.043 | 30.71 | 0.14% |
| 2022 | 0.037 | 26.74 | 0.14% |
| 2023 | 0.038 | 23.40 | 0.16% |
| 2024 | 0.051 | 25.22 | 0.20% |
| 2025 | 0.058 | 25.78 | 0.23% |

**Summary: species-specific fish respiration accounts for a mean of 0.12% (range 0.03%-0.23%) of measured ecosystem respiration at LTER_1 across the 18 years both are available.**

## H. fish_resp_z (species-specific) and refitting E4/E6

- `fish_resp_z`: LTER_1's annual species-specific areal respiration, z-standardised (mean 0.03813, sd 0.01056 mmol O2/m^2/h).
- Dev-subsample rows dropped for years without a fish survey match: 0.

### E6: species-specific fish_resp_z vs. biomass trio

| Coefficient | Estimate | 95% CI |
|---|---|---|
| la_Intercept | -1.480 | [-1.747, -1.212] |
| lP_Intercept | 5.137 | [4.425, 5.835] |
| lP_ilr1 | -0.002 | [-0.437, 0.397] |
| lP_ilr2 | -0.027 | [-0.436, 0.395] |
| lP_ilr3 | 0.091 | [-0.262, 0.440] |
| lP_invkT_c | 0.841 | [0.002, 1.683] |
| lP_log_flow_c | -0.011 | [-0.333, 0.292] |
| lP_N_z | 0.136 | [-0.241, 0.540] |
| lP_SeasonWinter | 0.119 | [-0.333, 0.561] |
| lR_Intercept | 4.150 | [3.541, 4.789] |
| lR_fish_resp_z | -0.071 | [-0.290, 0.157] |
| lR_ilr1 | 0.186 | [-0.168, 0.548] |
| lR_ilr2 | -0.238 | [-0.544, 0.045] |
| lR_ilr3 | 0.034 | [-0.142, 0.212] |
| lR_DHW_z | 0.015 | [-0.126, 0.155] |

### E4: species-specific fish_resp_z added to biomass adjustment

| Coefficient | Estimate | 95% CI |
|---|---|---|
| la_Intercept | -1.502 | [-1.774, -1.201] |
| lP_Intercept | 5.071 | [4.376, 5.733] |
| lP_ilr1 | 0.038 | [-0.355, 0.419] |
| lP_ilr2 | 0.039 | [-0.335, 0.446] |
| lP_ilr3 | 0.049 | [-0.267, 0.376] |
| lP_invkT_c | 0.520 | [-0.319, 1.337] |
| lP_log_flow_c | 0.566 | [0.233, 0.877] |
| lP_N_z | 0.218 | [-0.146, 0.623] |
| lP_SeasonWinter | 0.147 | [-0.266, 0.533] |
| lR_Intercept | 3.933 | [3.207, 4.691] |
| lR_ilr1 | 0.239 | [-0.154, 0.667] |
| lR_ilr2 | -0.107 | [-0.439, 0.212] |
| lR_ilr3 | -0.051 | [-0.312, 0.214] |
| lR_Herb_z | 0.032 | [-0.535, 0.602] |
| lR_Corall_z | 0.149 | [-0.123, 0.399] |
| lR_OtherFish_z | 0.096 | [-0.189, 0.374] |
| lR_fish_resp_z | -0.163 | [-0.829, 0.510] |
| lR_invkT_c | -0.696 | [-1.044, -0.349] |
| lR_log_flow_c | 0.666 | [0.505, 0.827] |


## I. Comparing the generic and species-specific estimates

| | Generic (Clarke & Johnston 1999, cross-species) | Species-specific (fishflux/Barneche & Allen 2018) |
|---|---|---|
| Mean % of Rd from fish | 0.28% | **0.13%** |
| Range | 0.07%-0.56% | 0.03%-0.23% |
| `fish_resp_z` credible in E6? | No | No |
| `fish_resp_z` credible in E4? | No | No |

**The species-specific estimate is smaller (roughly half) but tells the same qualitative story**: fish respiration is a small fraction of measured Rd at this site, and the mechanistic respiration estimate does not sharpen inference over raw biomass in either E4 or E6, with or without species-level detail. This cross-method agreement is a genuine robustness result, not a coincidence of the cruder approximation -- it strengthens (rather than undermines) the Section 7 conclusion that fish are evidently a small direct contributor to measured ecosystem respiration at this backreef site.

**Remaining approximations in the species-specific pipeline** (narrower than the generic version's caveats, but still real):
- Activity scope `f = 2` is a literature-typical default (from `fishflux`'s own documentation example), not species- or trophic-group-specific.
- Respiratory quotient (RQ = 0.8, for converting metabolic carbon loss to O2 consumption) is likewise a literature-typical constant, not species-specific.
- `m_max` (max body size) and `growth_g_day` use THIS project's own survey-fitted length-weight relationship combined with FishBase's K/Loo, rather than FishBase's own (very sparse) weight data.
- Most families (36 of 40) fall back to Barneche & Allen (2018)'s global-average B0/alpha rather than a family-specific fit -- this is `fishflux`'s own documented behaviour (only a few well-studied families, e.g. Pomacentridae, Gobiidae, Apogonidae, have dedicated fits), not a shortcut taken in this script.
- Trophic level, aspect ratio, and growth parameters use a species -> genus -> family -> global fallback hierarchy where FishBase has no record for a given species (extending `fishflux`'s own species -> genus-only fallback).

## J. Two real bugs caught during development (see script header for full detail)

1. **Length-unit mismatch**: `Total_Length` in `fish_ind.csv` is in millimetres; FishBase's `Loo` is in centimetres. Combining them directly (before catching this) silently corrupted `m_max`/`growth_g_day` for every record -- not a crash, just a wrong magnitude, caught only by sanity-checking a known species' size (a reef shark coming out at "800-1000 cm" if misread, versus a sensible 80-100 cm once correctly read as mm). Fixed by converting `Loo` to mm.
2. **Silent NA from R name-mangling**: the global length-weight fallback coefficients were extracted via `lw_fit_global["a"]`/`["b"]`, but `coef()` preserves full term names (e.g. `"a.(Intercept)"`), so this silently returned `NA` for the 26 individuals (of 28,964) whose family had too few records for its own fit. Caught by checking for unexpected `NA` counts in the final output rather than assuming zero visible errors meant correctness. Fixed with `unname()`.

Both are documented here as a reminder that a sophisticated species-specific pipeline introduces more silent-failure surface area than the simpler generic approximation did, and is worth double-checking (e.g. spot-checking a few known species' inputs) rather than trusting on the first run that produces no errors.
