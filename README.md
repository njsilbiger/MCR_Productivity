# Moorea Coral Reef Productivity Analysis — Data README

**Project:** Long-term Ecosystem Dynamics at MCR LTER Site 1, Moorea, French Polynesia  
**PI:** Nyssa Silbiger  
**Study site:** MCR LTER Site 1 backreef (17.46°S, 149.78°W), Moorea, French Polynesia  
**Study period:** 2005–2025  
**Analysis script:** `MCR_Productivity_Analysis.qmd`

---

## Directory Structure

```
Data/
├── raw_data/               # Unmodified source data
│   ├── QC_PP/              # Quality-controlled primary production deployments (29 files)
│   ├── MCR_LTER_Annual_Fish_Survey_20260304.csv
│   ├── MCR_LTER_Annual_Survey_Benthic_Cover_20251009.csv
│   ├── MCR_LTER_Macroalgal_CHN_2005_to_2024_20250616.csv
│   ├── MCR_LTER02_BTM_Backreef_Forereef_20251114.csv
│   └── WaterColumnN.csv
├── cmip6_sst/              # CMIP6 NetCDF files organised by experiment
│   ├── historical/
│   ├── ssp1_2_6/
│   ├── ssp2_4_5/
│   └── ssp5_8_5/
├── Benthic_summary_algae.csv
├── fish_clean.csv
├── fish_summary.csv
├── InSituTemp.csv
├── MooreaForereef.csv
├── OISST_data.csv
├── SST_year.csv
├── Turb_mean.csv
├── Water_N_mean.csv
└── Year_Averages_PP.csv
```

> **Note on raw data filenames:** The eight-digit suffix on MCR LTER files (e.g., `20260304`) is the download date from the EDI Data Portal (https://edirepository.org), not the date of data collection.

---

## Raw Data Files

### `raw_data/MCR_LTER_Annual_Fish_Survey_20260304.csv`
Annual belt-transect fish surveys at MCR LTER sites.  
**Source:** MCR LTER, EDI Data Portal  
**Rows:** 115,737 | **Years:** 2006–2025

| Column | Description |
|--------|-------------|
| Year | Survey year |
| Date | Survey date |
| Start / End | Transect start and end times |
| Location | General location name |
| Site | LTER site code (e.g., LTER_1) |
| Habitat | Reef habitat type (Backreef, Forereef, Fringing) |
| Depth | Survey depth (m) |
| Transect | Transect number |
| Swath | Transect width (m) |
| Taxonomy | Species binomial |
| Family | Fish family |
| Count | Number of individuals observed |
| Total_Length | Fork length of observed fish (cm) |
| Length_Anomaly | Flag for unusual length measurements |
| Biomass | Estimated wet biomass (g), converted from length–weight relationships; negative values indicate missing data |
| Comment | Field notes |
| Coarse_Trophic | Broad trophic category |
| Fine_Trophic | Detailed trophic category (e.g., Herbivore/Detritivore, Corallivore, Excavator) |
| Cloud_Cover | Cloud cover at time of survey |
| Wind_Velocity | Wind speed at time of survey |
| Sea_State | Sea state code |
| Swell | Swell height/direction |
| Visibility | Underwater visibility (m) |
| Surge | Surge intensity |
| Diver | Observer ID |

**Analysis use:** Filtered to Site LTER_1, Backreef habitat. Three large sharks (Biomass > 8,000 g) removed as outliers. Fish grouped into Herbivores, Corallivores, and Other trophic categories. Total biomass standardised to g m⁻² using a survey area of 1,200 m² (4 transects × 300 m²).

---

### `raw_data/MCR_LTER_Annual_Survey_Benthic_Cover_20251009.csv`
Annual point-intercept photoquadrat benthic surveys at MCR LTER sites.  
**Source:** MCR LTER, EDI Data Portal  
**Rows:** 234,522 | **Years:** 2005–2025

| Column | Description |
|--------|-------------|
| Year | Survey year |
| Date | Survey date |
| Location | General location name |
| Site | LTER site code |
| Habitat | Reef habitat type |
| Depth | Survey depth (m) |
| Transect | Transect number |
| Quadrat | Quadrat number within transect |
| Taxonomy_Substrate_Functional_Group | Benthic organism or substrate category at each point |
| Percent_Cover | Percent cover within quadrat |

**Analysis use:** Filtered to Backreef habitat. Taxa collapsed into functional groups: Coral, Fleshy Macroalgae and Turf, Crustose Corallines, Sand, and Other. Mean percent cover per group calculated from 5,000 points per site per year (5 transects × 1,000 points).

---

### `raw_data/MCR_LTER_Macroalgal_CHN_2005_to_2024_20250616.csv`
Carbon, hydrogen, and nitrogen (CHN) elemental analysis of macroalgal tissue samples collected at MCR LTER sites.  
**Source:** MCR LTER, EDI Data Portal  
**Rows:** 3,099 | **Years:** 2005–2024

| Column | Description |
|--------|-------------|
| Year | Collection year |
| Site | LTER site code |
| Habitat | Reef habitat type |
| Genus | Macroalgal genus |
| Dry_Weight | Sample dry weight (g) |
| C | Carbon content (% dry weight) |
| H | Hydrogen content (% dry weight) |
| N | Nitrogen content (% dry weight) |
| CN_ratio | Carbon-to-nitrogen ratio (molar) |
| Comment | Sample notes |

**Analysis use:** Filtered to *Turbinaria ornata* from LTER 1 Backreef. Annual mean %N, %C, and C:N ratio calculated. Tissue %N used as a proxy for water column dissolved inorganic nitrogen 

---

### `raw_data/MCR_LTER02_BTM_Backreef_Forereef_20251114.csv`
Continuous bottom-mounted temperature logger data from MCR LTER backreef and forereef moorings.  
**Source:** MCR LTER, EDI Data Portal  
**Rows:** 7,423,591 | **Date range:** 2005-05-23 to 2025-07-21 | **Resolution:** 15 minutes

| Column | Description |
|--------|-------------|
| site | LTER site code |
| time_local | Timestamp in French Polynesia local time (UTC−10) |
| time_utc | Timestamp in UTC |
| reef_type_code | Habitat (Backreef or Forereef) |
| sensor_type | Logger model/type |
| sensor_depth_m | Deployment depth (m) |
| temperature_c | Water temperature (°C) |

**Analysis use:** Filtered to Backreef records. Aggregated to daily means, then to annual means and maxima for the SEM temperature predictor.

---

### `raw_data/WaterColumnN.csv`
Bimonthly water column nutrient samples collected during MCR LTER research cruises.  
**Source:** MCR LTER  
**Rows:** 166 | **Date range:** 2005-08-10 to 2018-08-14

| Column | Description |
|--------|-------------|
| Cruise | Cruise identifier |
| Location | Sampling location name |
| Type | Sample type |
| Date | Collection date (M/D/YYYY) |
| Time | Collection time |
| Latitude | Decimal degrees N |
| Longitude | Decimal degrees E |
| Bottom_Depth | Water column depth at station (m) |
| Sample_Depth | Sample collection depth (m) |
| Phosphate | Phosphate concentration (µmol L⁻¹) |
| Silicate | Silicate concentration (µmol L⁻¹) |
| Nitrite | Nitrite concentration (µmol L⁻¹) |
| Nitrite_and_Nitrate | Combined nitrite + nitrate concentration (µmol L⁻¹) |
| Comments | Field notes |

**Analysis use:** Aggregated to annual means by site. Direct measurements available 2005–2018 only. 

---

### `raw_data/QC_PP/` — Quality-Controlled Primary Production Deployments
Individual CSV files from paired upstream–downstream oxygen sensor deployments used to measure open-water ecosystem metabolism at LTER Site 1 backreef. Each file corresponds to one seasonal deployment (typically January and June).  
**Files:** 29 deployments | **Temporal span:** 2008 June – 2025 January  
**Filename format:** `YYYY_Mon_LTER1_UPDN##_PPHour_UPvsDN.csv`  
**Resolution:** Hourly

Typical columns in each file:

| Column | Description |
|--------|-------------|
| DateTime | Date and time of measurement |
| Season | Deployment season (Jan or Jun) |
| UPDN | Sensor pair identifier |
| PAR | Photosynthetically active radiation (µmol photons m⁻² s⁻¹) |
| PP | Net community production (mmol O₂ m⁻² hr⁻¹; g O₂ m⁻² hr⁻¹ for deployments before April 2014) |
| UP_Oxy / DN_Oxy | Upstream and downstream dissolved oxygen (same units as PP) |
| UP_Temp / DN_Temp | Upstream and downstream water temperature (°C) |
| UP_Velocity_mps / DN_Velocity_mps | Current velocity at upstream and downstream sensors (m s⁻¹) |

**Unit note:** Deployments from 2007–2014 recorded oxygen flux in g O₂ m⁻² hr⁻¹. The analysis script converts these to mmol O₂ m⁻² hr⁻¹ by dividing by 32 and multiplying by 1,000.

---

### `cmip6_sst/` — CMIP6 Sea Surface Temperature (NetCDF)
Monthly sea surface temperature (*tos*) downloaded from the Copernicus Climate Data Store (DOI: 10.24381/cds.c866074c) for five CMIP6 ScenarioMIP models: ACCESS-CM2, IPSL-CM6A-LR, MIROC6, MPI-ESM1-2-LR, and UKESM1-0-LL. Spatially subsetted to a 3° bounding box centred on Moorea (~16°S to ~19°S, ~148°E to ~151°E).

| Subfolder | Experiment | Years |
|-----------|------------|-------|
| `historical/` | Historical simulation | 1950–2014 |
| `ssp1_2_6/` | SSP1-2.6 (strong mitigation, ~1.5–2°C warming by 2100) | 2015–2100 |
| `ssp2_4_5/` | SSP2-4.5 (intermediate stabilisation) | 2015–2100 |
| `ssp5_8_5/` | SSP5-8.5 (high emissions, ~4–5°C warming by 2100) | 2015–2100 |

Files are in NetCDF format (`.nc`). Each model has one file per experiment. All models were bias-corrected to the observed MCR 2005–2025 annual maximum temperature mean before use in projections.

---

## Processed Data Files

These files are generated by `MCR_Productivity_Analysis.qmd` and are stored here to avoid re-running slow model fits. **Do not edit manually** — re-run the relevant code chunk to regenerate.

---

### `Benthic_summary_algae.csv`
Annual mean percent cover of benthic functional groups at all MCR LTER sites.  
**Rows:** 469 | **Years:** 2006–2025 | **Generated by:** chunk `process-benthic`

| Column | Type | Description |
|--------|------|-------------|
| Year | integer | Survey year |
| Site | character | LTER site name (e.g., "LTER 1") |
| name | character | Functional group (Coral, Fleshy Macroalgae and Turf, Crustose Corallines, Sand, Other) |
| total_cover | numeric | Sum of raw Percent_Cover values for this group, site, and year |
| mean_cover | numeric | Mean percent cover (%), calculated as total_cover / 5,000 × 100 |

---

### `fish_summary.csv`
Annual fish biomass by trophic group at LTER Site 1 backreef.  
**Rows:** 60 | **Years:** 2006–2025 | **Generated by:** chunk `process-fish`

| Column | Type | Description |
|--------|------|-------------|
| Year | integer | Survey year |
| trophic_new | character | Trophic group: Herbivores, Corallivore, or Other |
| fish_g_m2 | numeric | Total fish biomass (g m⁻²), standardised to 1,200 m² survey area |

**Trophic group definitions:**
- **Herbivores:** Fine_Trophic values Brusher, Browser, Excavator, Concealed Cropper, Cropper, Scraper, Herbivore/Detritivore
- **Corallivore:** Fine_Trophic value Corallivore
- **Other:** All remaining trophic categories

---

### `fish_clean.csv`
Annual total fish biomass across all trophic groups at LTER Site 1 backreef.  
**Rows:** 20 | **Years:** 2006–2025 | **Generated by:** `Scripts/Fishdata.R`

| Column | Type | Description |
|--------|------|-------------|
| Year | integer | Survey year |
| total_fish_g_m2 | numeric | Total fish biomass across all species (g m⁻²), standardised to 1,200 m² |


---

### `InSituTemp.csv`
Annual temperature statistics derived from the continuous backreef temperature loggers.  
**Rows:** 21 | **Years:** 2005–2025 | **Generated by:** chunk `process-temperature`

| Column | Type | Description |
|--------|------|-------------|
| Year | integer | Calendar year |
| Mean_temp | numeric | Annual mean of daily mean temperatures (°C) |
| Max_temp | numeric | Annual maximum of daily mean temperatures (°C) |

---

### `OISST_data.csv`
Daily satellite sea surface temperature for the Moorea grid cell downloaded from NOAA Optimum Interpolation SST v2.1 (OISST) via ERDDAP.  
**Rows:** 7,299 | **Date range:** 2006-01-01 to 2025-12-31 | **Generated by:** `Scripts/SST_calc.R`

| Column | Type | Description |
|--------|------|-------------|
| lon | numeric | Longitude of grid cell (decimal degrees; negative = West) |
| lat | numeric | Latitude of grid cell (decimal degrees; negative = South) |
| t | date | Date (YYYY-MM-DD) |
| temp | numeric | Daily sea surface temperature (°C) |

**Grid cell:** 149.78°W, 17.46°S — the nearest OISST 0.25° grid node to LTER Site 1 backreef.

---

### `SST_year.csv`
Annual SST statistics derived from the OISST daily data.  
**Rows:** 20 | **Years:** 2006–2025 | **Generated by:** `Scripts/SST_calc.R`

| Column | Type | Description |
|--------|------|-------------|
| Year | integer | Calendar year |
| mean_SST | numeric | Annual mean SST from OISST (°C) |
| max_SST | numeric | Annual maximum SST from OISST (°C) |
| Year_Benthic | integer | Year + 1 (aligns SST to following year's benthic response) |

---

### `Turb_mean.csv`
Annual mean tissue elemental composition of *Turbinaria ornata* at LTER Site 1 backreef.  
**Rows:** 18 | **Years:** 2007–2024 | **Generated by:** chunk `process-nutrients`

| Column | Type | Description |
|--------|------|-------------|
| Year | integer | Collection year |
| N_percent | numeric | Mean tissue nitrogen content (% dry weight) |
| C_percent | numeric | Mean tissue carbon content (% dry weight) |
| CN | numeric | Mean molar carbon-to-nitrogen ratio |

---

### `Water_N_mean.csv`
Annual mean water column nutrient concentrations at LTER Site 1.  
**Rows:** 14 | **Years:** 2005–2018 | **Generated by:** chunk `process-nutrients`

| Column | Type | Description |
|--------|------|-------------|
| Year | integer | Calendar year |
| Phosphate | numeric | Annual mean phosphate concentration (µmol L⁻¹) |
| Silicate | numeric | Annual mean silicate concentration (µmol L⁻¹) |
| Nitrite | numeric | Annual mean nitrite concentration (µmol L⁻¹) |
| Nitrite_and_Nitrate | numeric | Annual mean nitrite + nitrate concentration (µmol L⁻¹) |
| Site | character | Site identifier ("LTER 1") |

---

### `Year_Averages_PP.csv`
Annual ecosystem metabolism parameters derived from the PP deployments and the Bayesian photosynthesis–irradiance (PI) curve model.  
**Rows:** 18 | **Years:** 2008–2025 | **Generated by:** chunk `bayesian-pi-curve`

| Column | Type | Description |
|--------|------|-------------|
| Year | integer | Calendar year |
| NP_mean | numeric | Annual mean net community production (mmol O₂ m⁻² hr⁻¹) |
| NP_SE | numeric | Standard error of NP_mean |
| NP_max | numeric | Annual maximum net community production (mmol O₂ m⁻² hr⁻¹) |
| GP_mean | numeric | Annual mean gross production (mmol O₂ m⁻² hr⁻¹) |
| GP_SE | numeric | Standard error of GP_mean |
| R_mean | numeric | Annual mean ecosystem respiration from raw dark measurements (mmol O₂ m⁻² hr⁻¹; negative values) |
| R_SE | numeric | Standard error of R_mean |
| Temperature_mean | numeric | Annual mean water temperature averaged across upstream and downstream sensors (°C) |
| Flow_mean | numeric | Annual mean current velocity averaged across sensors (m s⁻¹) |
| PAR_mean | numeric | Annual mean PAR during daylight hours (µmol photons m⁻² s⁻¹) |
| Pmax | numeric | Maximum photosynthetic capacity from Bayesian PI curve — posterior median (mmol O₂ m⁻² hr⁻¹) |
| Rd | numeric | Ecosystem respiration from Bayesian PI curve — posterior median, expressed as a positive rate (mmol O₂ m⁻² hr⁻¹) |

**Note on Pmax and Rd:** These are model-derived estimates, not direct measurements. See the Methods section of `MCR_Productivity_Analysis.qmd` for the PI curve model specification.

---

### `MooreaForereef.csv`
Annual point-intercept photoquadrat benthic surveys at LTER Site 1 **forereef** — a companion dataset to the backreef benthic data in the raw_data folder.  
**Rows:** 44,007 | **Date range:** 2005-04 to 2025-04

| Column | Type | Description |
|--------|------|-------------|
| Date | character | Survey month in YYYY-MM format |
| Location | character | Location name |
| Site | character | LTER site code (LTER_1) |
| Habitat | character | Reef habitat (Forereef) |
| Depth | numeric | Survey depth (m) |
| Section_of_Transect | character | Transect section identifier |
| Quadrat_within_section | integer | Quadrat number within section |
| Quad40 | character | Quadrat 40 identifier |
| Taxonomy_Substrate_or_Functional_Group | character | Benthic organism or substrate at each point |
| Percent_Cover | numeric | Percent cover within quadrat |

---

## Data Flow Summary

```
raw_data/MCR_LTER02_BTM_Backreef_Forereef_20251114.csv ──► InSituTemp.csv
raw_data/MCR_LTER_Annual_Survey_Benthic_Cover_20251009.csv ──► Benthic_summary_algae.csv
raw_data/MCR_LTER_Annual_Fish_Survey_20260304.csv ──► fish_summary.csv
                                                  ──► fish_clean.csv
raw_data/MCR_LTER_Macroalgal_CHN_2005_to_2024_20250616.csv ──► Turb_mean.csv
raw_data/WaterColumnN.csv ──► Water_N_mean.csv
raw_data/QC_PP/*.csv ──► Year_Averages_PP.csv  (via Bayesian PI curve model)
ERDDAP (NOAA OISST) ──► OISST_data.csv ──► SST_year.csv
cmip6_sst/*.nc ──► (used directly in projection code; no intermediate file written)

All processed files ──► Year_Averages (in-memory) ──► Bayesian SEM ──► Projections
```

---

## General Notes

- All backreef analyses are restricted to **LTER Site 1** unless otherwise noted.
- Fish biomass units throughout are **grams wet weight per square metre (g m⁻²)**.
- Oxygen flux units throughout are **mmol O₂ m⁻² hr⁻¹** after unit conversion (see QC_PP note).
- Raw data files should not be edited. To update the analysis with a new EDI download, replace the raw file and update the filename reference in the relevant processing chunk of `MCR_Productivity_Analysis.qmd`.
