# QC summary: PI-curve metabolism model (Revision Step 6 / Section 6)

Generated: 2026-10-02 12:10:20.669946

## A. LTER_1 benthic composition: pivot coordinates

- Zero cells in LTER_1's Coral/Algae/CCA/Other count matrix (20 years x 4 parts): 0 (zero-replacement is effectively a no-op here, but run for consistency with Sections 4-5).
- Pivot coordinates (Hron et al. 2012), coral first: `ilr1` = coral vs. geometric mean of (algae, CCA, other); `ilr2` = algae vs. geometric mean of (CCA, other); `ilr3` = CCA vs. other. Computed from LTER_1's site-level composition for each of the 20 years 2006-2025.

## B. Hourly model data

- `pi_model_hourly.csv`: 2808 hourly rows (of 2928 before dropping rows with missing covariates) -- 120 rows dropped, all from 2025 (`N_site_year.csv` does not yet cover 2025 -- the same lab-processing lag already documented for this project's current year elsewhere).
- `invkT_c` centred on the sample mean temperature (28.62 degC), not an external MTE reference temperature -- the lR coefficient is interpretable as the activation-energy-scaled deviation from the average-condition respiration rate in this dataset, not an absolute physiological reference.
- `N_z` uses `N_site_year.csv`'s LTER_1 series (mean 0.685, sd 0.147).

## C. Development fit: 25% of days (stratified by Year)

- Dev subsample: 624 hourly rows, 26 of 117 days (every year represented).

## D. Full fit (117 days, 2808 hourly rows)

- Same specification as the dev fit, confirmed on the full data. `invkT_c` coefficient prior for lR corrected to `normal(-0.65, 0.2)` (see code comment) after the dev fit on the plan's literal `normal(0.65, 0.2)` produced a biologically backwards posterior (respiration decreasing with warming); decided with the PI 2026-10-01.

## G. Task 6.2 covariates: fish biomass by trophic group, DHW (LTER_1, year-level)

- `Herb_z`/`Corall_z`/`OtherFish_z`: mean biomass (g/m^2) across the 4 fish transects at LTER_1 per year, z-standardised (needed for E4/E6, the DAG's `Herb`/`Corall`/`OtherFish` nodes -- NOT the bioenergetics `fish_resp_z` of Task 7.1, which remains unbuilt).
- `DHW_z`: z-standardised `DHW_max` at LTER_1 (needed for E6/E7).

## H. Task 6.2 per-estimand screening fits (25%-of-days dev subsample)

All four fits below use `pi_dev_data` (624 hourly rows, the same stratified 25%-of-days subsample validated in Section C1), NOT the full 2808-row dataset -- the full Task 6.1 fit alone took ~52 minutes, and four more full fits would cost several hours of sequential compute. **These are screening/development results; scaling any of them to the full dataset is a flagged follow-up.**

### E4: Benthos -> Rd (direct); adjustment {Corall, DayTemp, Flow, Herb, OtherFish}

| Coefficient | Estimate | 95% CI |
|---|---|---|
| la_Intercept | -1.502 | [-1.768, -1.198] |
| lP_Intercept | 5.065 | [4.407, 5.717] |
| lP_ilr1 | 0.039 | [-0.341, 0.409] |
| lP_ilr2 | 0.039 | [-0.345, 0.449] |
| lP_ilr3 | 0.044 | [-0.277, 0.383] |
| lP_invkT_c | 0.532 | [-0.311, 1.394] |
| lP_log_flow_c | 0.565 | [0.234, 0.882] |
| lP_N_z | 0.222 | [-0.162, 0.611] |
| lP_SeasonWinter | 0.148 | [-0.259, 0.537] |
| lR_Intercept | 3.939 | [3.201, 4.667] |
| lR_ilr1 | 0.259 | [-0.146, 0.660] |
| lR_ilr2 | -0.106 | [-0.443, 0.225] |
| lR_ilr3 | -0.043 | [-0.306, 0.220] |
| lR_Herb_z | -0.098 | [-0.348, 0.155] |
| lR_Corall_z | 0.118 | [-0.113, 0.344] |
| lR_OtherFish_z | 0.057 | [-0.166, 0.277] |
| lR_invkT_c | -0.702 | [-1.063, -0.353] |
| lR_log_flow_c | 0.664 | [0.509, 0.819] |

### E6: Fish (Herb+Corall+OtherFish) -> Rd (direct); adjustment {Benthos, DHW}

| Coefficient | Estimate | 95% CI |
|---|---|---|
| la_Intercept | -1.477 | [-1.765, -1.188] |
| lP_Intercept | 5.122 | [4.379, 5.843] |
| lP_ilr1 | -0.013 | [-0.436, 0.395] |
| lP_ilr2 | -0.029 | [-0.431, 0.388] |
| lP_ilr3 | 0.076 | [-0.274, 0.423] |
| lP_invkT_c | 0.847 | [0.044, 1.692] |
| lP_log_flow_c | -0.014 | [-0.322, 0.294] |
| lP_N_z | 0.167 | [-0.232, 0.583] |
| lP_SeasonWinter | 0.133 | [-0.317, 0.584] |
| lR_Intercept | 4.109 | [3.411, 4.805] |
| lR_Herb_z | -0.119 | [-0.356, 0.116] |
| lR_Corall_z | 0.073 | [-0.129, 0.266] |
| lR_OtherFish_z | 0.061 | [-0.133, 0.260] |
| lR_ilr1 | 0.171 | [-0.214, 0.557] |
| lR_ilr2 | -0.235 | [-0.551, 0.073] |
| lR_ilr3 | -0.025 | [-0.253, 0.205] |
| lR_DHW_z | 0.049 | [-0.097, 0.196] |

### E7: DayTemp -> Rd (direct); adjustment {DHW, Flow}

| Coefficient | Estimate | 95% CI |
|---|---|---|
| la_Intercept | -1.506 | [-1.782, -1.206] |
| lP_Intercept | 4.977 | [4.270, 5.646] |
| lP_ilr1 | -0.059 | [-0.444, 0.322] |
| lP_ilr2 | 0.052 | [-0.326, 0.439] |
| lP_ilr3 | 0.027 | [-0.286, 0.352] |
| lP_invkT_c | 0.524 | [-0.284, 1.341] |
| lP_log_flow_c | 0.580 | [0.250, 0.889] |
| lP_N_z | 0.188 | [-0.176, 0.572] |
| lP_SeasonWinter | 0.138 | [-0.267, 0.540] |
| lR_Intercept | 3.619 | [3.470, 3.773] |
| lR_invkT_c | -0.682 | [-1.047, -0.311] |
| lR_DHW_z | -0.014 | [-0.193, 0.165] |
| lR_log_flow_c | 0.682 | [0.525, 0.833] |

### E9/E10: Benthos -> Pmax / DayTemp -> Pmax (direct); adjustment {Benthos, DayTemp, Flow, N}

| Coefficient | Estimate | 95% CI |
|---|---|---|
| la_Intercept | -1.509 | [-1.785, -1.195] |
| lP_Intercept | 5.094 | [4.377, 5.768] |
| lP_ilr1 | 0.045 | [-0.360, 0.424] |
| lP_ilr2 | 0.085 | [-0.273, 0.469] |
| lP_ilr3 | 0.083 | [-0.229, 0.417] |
| lP_invkT_c | 0.621 | [-0.186, 1.468] |
| lP_log_flow_c | 0.517 | [0.216, 0.811] |
| lP_N_z | 0.177 | [-0.187, 0.564] |
| lR_Intercept | 4.098 | [3.579, 4.584] |
| lR_ilr1 | 0.351 | [0.099, 0.607] |
| lR_ilr2 | -0.196 | [-0.486, 0.120] |
| lR_ilr3 | 0.021 | [-0.169, 0.212] |
| lR_invkT_c | -0.721 | [-1.086, -0.359] |
| lR_log_flow_c | 0.674 | [0.507, 0.830] |
| lR_SeasonWinter | 0.161 | [-0.103, 0.434] |

