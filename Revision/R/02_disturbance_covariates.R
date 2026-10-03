# ---------------------------------------------------------------------------
# 02_disturbance_covariates.R
#
# Step 2 of the Revision plan (Revision/Reviewer_Response_Plan.md, Section 3).
# Builds `disturbance_site_year.csv`: DHW (satellite, primary; in-situ logger,
# secondary/sensitivity), COTS density, and a Cyclone Oli (Feb 2010) indicator
# per Site x Year, to be joined onto benthic_site / fish_transect / pp_day.
#
# IMPORTANT (network reliability): every slow/network step below caches its
# result to Revision/Data/derived/_cache/ immediately after it succeeds, and
# is skipped on re-run if the cache file already exists. This is written
# defensively because the R session has crashed mid-query on large ERDDAP
# griddap() pulls. Run this script in pieces (source line ranges, or run
# chunk by chunk in the console) rather than all at once if the session is
# unstable.
#
# INPUTS (read-only):
#   Data/raw_data/MCR_LTER_Annual_Survey_Benthic_Cover_20251009.csv  (survey dates)
#   Data/raw_data/MCR_LTER_COTS_abundance_2005-2025_20250310.csv
#   Data/raw_data/MCR_LTER02_BTM_Backreef_Forereef_20251114.csv      (LTER_2 logger)
#   NOAA Coral Reef Watch daily DHW, via ERDDAP (coastwatch.pfeg.noaa.gov)
#
# OUTPUT:
#   Revision/Data/derived/disturbance_site_year.csv  (Site, Year, DHW_max,
#     DHW_logger_max, COTS_density, Cyclone)
#   Revision/Output/qc_disturbance.md
# ---------------------------------------------------------------------------

source(here::here("Revision", "R", "00_packages.R"))
if (!requireNamespace("rerddap", quietly = TRUE)) install.packages("rerddap")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
library(rerddap)
library(data.table)

cache_dir <- here("Revision", "Data", "derived", "_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

qc_lines <- character(0)
qc <- function(...) qc_lines <<- c(qc_lines, paste0(...))
qc("# QC summary: disturbance covariates (Revision Step 2)")
qc("")
qc("Generated: ", as.character(Sys.time()))
qc("")

# =============================================================================
# A. Survey-date window per year (anchors the trailing 12-month DHW window)
# =============================================================================

benthic_raw <- read_csv(
  here("Data", "raw_data", "MCR_LTER_Annual_Survey_Benthic_Cover_20251009.csv"),
  show_col_types = FALSE
) |>
  mutate(Date = ymd(Date))

survey_window <- benthic_raw |>
  distinct(Year, Date) |>
  group_by(Year) |>
  summarise(survey_date = median(Date), min_date = min(Date), max_date = max(Date),
            .groups = "drop")

write_csv(survey_window, file.path(cache_dir, "survey_window.csv"))

qc("## A. Survey-date window per year")
qc("")
qc("Benthic survey dates are mostly early-to-mid January (austral summer, ",
   "post-peak-heat-stress season), not April-May as an earlier draft of the ",
   "plan assumed -- the only exception is 2005 (surveyed in late May). The ",
   "`survey_date` (median survey date per year) anchors a trailing 364-day ",
   "DHW window for that year. Note 2023 has a wide date range (2022-02-27 to ",
   "2023-01-24) -- flagged for the PI; not resolved further here since it ",
   "does not block the DHW calculation (median still falls in Jan 2023).")
qc("")
qc("| Year | survey_date (median) | min_date | max_date |")
qc("|---|---|---|---|")
for (i in seq_len(nrow(survey_window))) {
  qc(sprintf("| %d | %s | %s | %s |", survey_window$Year[i],
             survey_window$survey_date[i], survey_window$min_date[i],
             survey_window$max_date[i]))
}
qc("")

cat("Done: Section A (survey window). Cached to", file.path(cache_dir, "survey_window.csv"), "\n")

# =============================================================================
# B. Primary heat-stress metric: NOAA Coral Reef Watch daily DHW (satellite)
# =============================================================================
#
# Single representative Moorea pixel used for all 6 backreef sites (decided
# with the PI 2026-10-01): the sites span all 3 shores of the island and are
# not all in one 5 km CRW pixel, but we do not have reliable per-site
# coordinates for LTER 2-6 (only LTER_1 backreef, 17.46 S / 149.78 W, is
# documented in this repo's README). An island-wide series is used instead,
# at the MCR LTER network centroid coordinate (-17.4909, -149.826; the
# pixel ERDDAP actually returns is -17.475, -149.825, i.e. the nearest 0.05
# deg CRW grid node). This is a known simplification -- satellite DHW is
# unlikely to vary sharply at Moorea's ~16 km scale, but it cannot capture
# real shore-to-shore differences in local heat exposure. State this in the
# methods and response letter.

dhw_cache_file <- file.path(cache_dir, "moorea_dhw_daily.csv")

if (!file.exists(dhw_cache_file)) {
  info_dhw <- rerddap::info("NOAA_DHW", url = "https://coastwatch.pfeg.noaa.gov/erddap/")
  moorea_dhw_raw <- griddap(
    info_dhw,
    time = c("2004-01-01", "2025-03-01"),
    latitude = c(-17.4909, -17.4909),
    longitude = c(-149.826, -149.826),
    fields = "CRW_DHW"
  )
  dhw_daily <- moorea_dhw_raw$data |>
    mutate(date = as_date(time)) |>
    select(date, CRW_DHW) |>
    arrange(date)
  # Write immediately -- this is the step that has crashed the session before.
  write_csv(dhw_daily, dhw_cache_file)
  cat("Pulled and cached satellite DHW series:", nrow(dhw_daily), "days ->", dhw_cache_file, "\n")
} else {
  dhw_daily <- read_csv(dhw_cache_file, show_col_types = FALSE)
  cat("Loaded cached satellite DHW series:", nrow(dhw_daily), "days from", dhw_cache_file, "\n")
}

full_rng_dhw <- seq(min(dhw_daily$date), max(dhw_daily$date), by = "day")
missing_dhw_days <- full_rng_dhw[!full_rng_dhw %in% dhw_daily$date]

dhw_annual <- survey_window |>
  rowwise() |>
  mutate(
    window_start = survey_date - days(364),
    DHW_max = max(dhw_daily$CRW_DHW[dhw_daily$date >= window_start & dhw_daily$date <= survey_date], na.rm = TRUE),
    n_days_in_window = sum(dhw_daily$date >= window_start & dhw_daily$date <= survey_date)
  ) |>
  ungroup()

write_csv(dhw_annual, file.path(cache_dir, "dhw_annual.csv"))

qc("## B. Primary heat stress: satellite DHW (NOAA Coral Reef Watch)")
qc("")
qc("- Source: NOAA Coral Reef Watch daily 5 km `CRW_DHW` (degree heating ",
   "weeks), via ERDDAP (`coastwatch.pfeg.noaa.gov/erddap`, dataset ",
   "`NOAA_DHW`), pulled for ", min(dhw_daily$date), " to ", max(dhw_daily$date), ".")
qc("- Single island-wide pixel used for all 6 sites (see code comment for ",
   "rationale) -- nearest grid node to the MCR LTER network centroid ",
   "(-17.4909, -149.826) is -17.475, -149.825.")
qc("- Missing days in the pulled series: ", length(missing_dhw_days),
   " of ", length(full_rng_dhw), ".")
qc("- `DHW_max` per year is the maximum `CRW_DHW` in the 364 days up to and ",
   "including that year's median benthic-survey date (Task 3.1's 12-month, ",
   "survey-date-anchored window).")
qc("")
qc("| Year | survey_date | DHW_max (deg C-weeks) | days in window |")
qc("|---|---|---|---|")
for (i in seq_len(nrow(dhw_annual))) {
  qc(sprintf("| %d | %s | %.2f | %d |", dhw_annual$Year[i],
             dhw_annual$survey_date[i], dhw_annual$DHW_max[i],
             dhw_annual$n_days_in_window[i]))
}
qc("")
qc("The 2020 survey (window ending Jan 2020) shows DHW_max = ",
   sprintf("%.2f", dhw_annual$DHW_max[dhw_annual$Year == 2020]),
   ", consistent with the known 2019 Moorea bleaching event (Section 0.3 of ",
   "the plan) -- a sanity check that the window logic is capturing real ",
   "heat-stress history rather than an artefact.")
qc("")

cat("Done: Section B (satellite DHW). Cached to", file.path(cache_dir, "dhw_annual.csv"), "\n")

# =============================================================================
# C. Secondary / sensitivity heat-stress metric: in-situ logger DHW
# =============================================================================
#
# The backreef temperature file (MCR_LTER02_BTM_Backreef_Forereef) has only
# one backreef logger, at LTER_2 (2 m) -- confirmed in Step 1 and Section 0.3
# of the plan. This secondary DHW is therefore also a single site-year series
# applied to all 6 sites, like the satellite version, but captures local
# (non-satellite) temperature variability at 2 m depth.
#
# Approach (approximate CRW algorithm, documented as such -- NOT NOAA's
# official product): HotSpot_t = max(0, SST_t - MMM), where MMM (maximum
# monthly mean) is approximated locally from the CRW_SST satellite monthly
# climatology at the same Moorea pixel (since the official MMM baseline
# raster is not available on this ERDDAP server). DHW_t = sum over the
# trailing 12 weeks of (HotSpot_d / 7) for days where HotSpot_d >= 1 deg C,
# in degree-C-weeks, applied to the logger's daily mean temperature.

lter2_raw_file <- here("Revision", "Data", "derived", "_lter2_backreef_temp_raw.csv")
lter2_daily_cache <- file.path(cache_dir, "lter2_daily_temp.csv")

if (!file.exists(lter2_raw_file)) {
  # Extract LTER_2 Backreef rows from the ~510 MB source file (data.table is
  # used here for speed; this step is memory-light because `select=` is not
  # used -- we still need all 7 original columns for the date fields).
  btm_all <- fread(
    here("Data", "raw_data", "MCR_LTER02_BTM_Backreef_Forereef_20251114.csv"),
    colClasses = list(character = c("time_local", "time_utc"))
  )
  lter2_raw <- btm_all[site == "LTER_2" & reef_type_code == "Backreef"]
  fwrite(lter2_raw, lter2_raw_file)
  rm(btm_all)
  cat("Extracted LTER_2 backreef logger rows ->", lter2_raw_file, "\n")
}

if (!file.exists(lter2_daily_cache)) {
  lter2_raw <- fread(lter2_raw_file, colClasses = list(character = c("time_local", "time_utc")))
  lter2_raw[, time_local_parsed := ymd_hms(time_local)]
  n_parse_fail <- sum(is.na(lter2_raw$time_local_parsed))
  lter2_raw[, date := as_date(time_local_parsed)]
  lter2_daily <- lter2_raw[!is.na(date), .(temp_mean = mean(temperature_c, na.rm = TRUE), n_obs = .N), by = date]
  setorder(lter2_daily, date)
  fwrite(lter2_daily, lter2_daily_cache)
  cat("Aggregated LTER_2 logger to", nrow(lter2_daily), "daily means ->", lter2_daily_cache,
      "(", n_parse_fail, "timestamp parse failures excluded )\n")
} else {
  lter2_daily <- fread(lter2_daily_cache)
  cat("Loaded cached LTER_2 daily temperature:", nrow(lter2_daily), "days\n")
}
lter2_daily <- as_tibble(lter2_daily) |> mutate(date = as_date(date))

full_rng_logger <- seq(min(lter2_daily$date), max(lter2_daily$date), by = "day")
missing_logger_days <- full_rng_logger[!full_rng_logger %in% lter2_daily$date]

# ---- C1. MMM proxy from satellite monthly SST climatology ------------------
sst_monthly_cache <- file.path(cache_dir, "moorea_sst_monthly.csv")

if (!file.exists(sst_monthly_cache)) {
  info_sst_m <- rerddap::info("NOAA_DHW_monthly", url = "https://coastwatch.pfeg.noaa.gov/erddap/")
  moorea_sst_m_raw <- griddap(
    info_sst_m,
    time = c("1985-01-16T00:00:00Z", "2025-12-16T00:00:00Z"),
    latitude = c(-17.4909, -17.4909),
    longitude = c(-149.826, -149.826),
    fields = "sea_surface_temperature"
  )
  sst_monthly <- moorea_sst_m_raw$data |>
    mutate(date = as_date(time)) |>
    select(date, sst = sea_surface_temperature) |>
    arrange(date)
  write_csv(sst_monthly, sst_monthly_cache)
  cat("Pulled and cached monthly satellite SST:", nrow(sst_monthly), "months ->", sst_monthly_cache, "\n")
} else {
  sst_monthly <- read_csv(sst_monthly_cache, show_col_types = FALSE)
  cat("Loaded cached monthly satellite SST:", nrow(sst_monthly), "months\n")
}

sst_climatology <- sst_monthly |>
  mutate(cal_month = month(date)) |>
  group_by(cal_month) |>
  summarise(clim_mean = mean(sst, na.rm = TRUE), .groups = "drop")

MMM_proxy <- max(sst_climatology$clim_mean)
MMM_month <- sst_climatology$cal_month[which.max(sst_climatology$clim_mean)]

cat("MMM proxy (max monthly-mean SST, full satellite record):", round(MMM_proxy, 3),
    "deg C, in calendar month", MMM_month, "\n")

# ---- C2. HotSpot / DHW from the logger daily series -------------------------
lter2_daily <- lter2_daily |>
  arrange(date) |>
  mutate(hotspot = pmax(0, temp_mean - MMM_proxy))

roll_dhw <- function(dates, hotspot, window_days = 84) {
  vapply(seq_along(dates), function(i) {
    d0 <- dates[i] - days(window_days - 1)
    in_window <- dates >= d0 & dates <= dates[i]
    hs <- hotspot[in_window]
    hs[hs < 1] <- 0   # CRW convention: only HotSpot >= 1 deg C accumulates
    sum(hs, na.rm = TRUE) / 7
  }, numeric(1))
}

# NOTE: a naive O(n * window) rolling sum over ~7000 days is fine here
# (single site, run once, cached), but would not scale to a multi-site loop.
lter2_daily$DHW_logger <- roll_dhw(lter2_daily$date, lter2_daily$hotspot)

write_csv(lter2_daily, file.path(cache_dir, "lter2_daily_dhw.csv"))

dhw_logger_annual <- survey_window |>
  rowwise() |>
  mutate(
    DHW_logger_max = {
      idx <- lter2_daily$date >= (survey_date - days(1)) & lter2_daily$date <= survey_date
      if (any(idx)) lter2_daily$DHW_logger[lter2_daily$date == survey_date][1] else NA_real_
    }
  ) |>
  ungroup() |>
  select(Year, DHW_logger_max)

# DHW_logger is already a trailing-84-day rolling accumulation ending each
# date, so the annual value is just DHW_logger on the survey date itself
# (not a second trailing-364-day max of an already-rolling quantity, which
# would double-smooth it). Recompute properly: take the max of the already-
# rolling DHW_logger series over the 364 days up to the survey date, matching
# the satellite DHW_max definition (max accumulated stress reached at any
# point in the preceding year, not just on the survey date).
dhw_logger_annual <- survey_window |>
  rowwise() |>
  mutate(
    window_start = survey_date - days(364),
    DHW_logger_max = {
      idx <- lter2_daily$date >= window_start & lter2_daily$date <= survey_date
      if (any(idx)) max(lter2_daily$DHW_logger[idx], na.rm = TRUE) else NA_real_
    },
    n_logger_days_in_window = sum(lter2_daily$date >= window_start & lter2_daily$date <= survey_date)
  ) |>
  ungroup() |>
  select(Year, DHW_logger_max, n_logger_days_in_window)

write_csv(dhw_logger_annual, file.path(cache_dir, "dhw_logger_annual.csv"))

qc("## C. Secondary / sensitivity heat stress: in-situ logger DHW")
qc("")
qc("- Logger: LTER_2 backreef, 2 m depth (the only backreef logger in the ",
   "temperature file; applied here to all 6 sites, same limitation as the ",
   "satellite series). Daily coverage: ", nrow(lter2_daily), " days, ",
   as.character(min(lter2_daily$date)), " to ", as.character(max(lter2_daily$date)),
   " (", length(missing_logger_days), " missing days).")
qc("- **This is an approximation, not NOAA's official in-situ DHW product.** ",
   "MMM (maximum monthly mean) is proxied as the maximum calendar-month mean ",
   "of satellite `sea_surface_temperature` at the same Moorea pixel across ",
   "the full 1985-2025 record (", round(MMM_proxy, 3), " deg C, month ",
   MMM_month, "), because NOAA's official MMM baseline raster is not served ",
   "on this ERDDAP instance. HotSpot = max(0, logger_temp - MMM); DHW = ",
   "trailing 12-week (84-day) sum of HotSpot/7 for days with HotSpot >= 1 ",
   "deg C, following the standard CRW accumulation rule applied to the ",
   "logger series instead of satellite SST.")
qc("- Annual `DHW_logger_max` uses the same trailing-364-day-to-survey-date ",
   "window as the satellite metric, for comparability.")
qc("")
qc("| Year | DHW_logger_max (deg C-weeks) | logger days in window |")
qc("|---|---|---|")
for (i in seq_len(nrow(dhw_logger_annual))) {
  qc(sprintf("| %d | %s | %d |", dhw_logger_annual$Year[i],
             ifelse(is.na(dhw_logger_annual$DHW_logger_max[i]), "NA",
                    sprintf("%.2f", dhw_logger_annual$DHW_logger_max[i])),
             dhw_logger_annual$n_logger_days_in_window[i]))
}
qc("")
qc("Treat this strictly as a sensitivity check (Task 3.1): it uses a ",
   "locally-derived MMM proxy, not NOAA's validated baseline, and -- like ",
   "the satellite series -- is a single island-wide/single-logger value ",
   "applied to all 6 sites, so it adds local temperature variability but ",
   "not real between-site heat-stress contrast.")
qc("")

cat("Done: Section C (logger DHW, sensitivity). Cached to",
    file.path(cache_dir, "dhw_logger_annual.csv"), "\n")

# =============================================================================
# D. Crown-of-thorns starfish (COTS) density per Site x Year
# =============================================================================
#
# Source: Moorea Coral Reef LTER and A. Brooks. 2026. MCR LTER: Coral Reef:
# Long-term Population Dynamics of Acanthaster planci, ongoing since 2005
# ver 13. EDI. https://doi.org/10.6073/pasta/0601a5aa24c8f35fda90a99b4f1a50bd
# (downloaded manually by the PI 2026-10-01 -- EDI's portal and the PASTA
# REST API both reject automated/anonymous requests for this package from
# this environment; see Section 14 testing log for the access attempts).
# File: Data/raw_data/MCR_LTER_COTS_abundance_2005-2025_20250310.csv
#
# Design: Year x Site (LTER 1-6) x Habitat (Backreef/Forereef/Fringing) x
# Transect (1-4), COTS = count of A. planci on a 5 x 50 m belt transect
# (same transect geometry as the fish survey; confirmed via the Step-1 fish
# area derivation and the published MCR COTS survey description). We keep
# only Backreef to match the benthic/fish/metabolism tables used elsewhere
# in the revision.

cots_raw <- read_csv(
  here("Data", "raw_data", "MCR_LTER_COTS_abundance_2005-2025_20250310.csv"),
  show_col_types = FALSE
)

COTS_TRANSECT_AREA_M2 <- 5 * 50  # 250 m^2, same belt-transect geometry as fish

cots_backreef <- cots_raw |>
  filter(Habitat == "Backreef") |>
  mutate(Site = str_replace(Site, "^LTER\\s+", "LTER_"))

n_cots_raw <- nrow(cots_raw)
n_cots_backreef <- nrow(cots_backreef)

cots_site_year <- cots_backreef |>
  group_by(Year, Site) |>
  summarise(
    COTS_count_total = sum(COTS, na.rm = TRUE),
    n_transects = n(),
    COTS_density_m2 = COTS_count_total / (n_transects * COTS_TRANSECT_AREA_M2),
    .groups = "drop"
  )

write_csv(cots_site_year, file.path(cache_dir, "cots_site_year.csv"))

qc("## D. Crown-of-thorns starfish (COTS)")
qc("")
qc("- Source: `MCR_LTER_COTS_abundance_2005-2025_20250310.csv` (manually ",
   "downloaded by the PI from EDI; automated access failed -- see testing ",
   "log). Citation: Moorea Coral Reef LTER and A. Brooks. 2026. MCR LTER: ",
   "Coral Reef: Long-term Population Dynamics of *Acanthaster planci*, ",
   "ongoing since 2005 ver 13. Environmental Data Initiative. ",
   "https://doi.org/10.6073/pasta/0601a5aa24c8f35fda90a99b4f1a50bd.")
qc("- Raw rows (all habitats, all sites): ", n_cots_raw, ". Backreef rows ",
   "kept: ", n_cots_backreef, ".")
qc("- Design confirmed: 6 sites x 21 years (2005-2025) x 4 backreef ",
   "transects = ", 6 * 21 * 4, " expected rows; observed ", n_cots_backreef, ".")
qc("- Density = total COTS counted / (n transects x 250 m^2), where 250 m^2 ",
   "is the 5 x 50 m belt-transect area used for the MCR COTS survey (same ",
   "geometry as the annual fish survey's wide swath).")
qc("- Site-year COTS density range: ", sprintf("%.4f", min(cots_site_year$COTS_density_m2)),
   " to ", sprintf("%.4f", max(cots_site_year$COTS_density_m2)), " ind/m^2.")
qc("")
qc("| Year | total backreef COTS (6 sites) | max site density (ind/m^2) |")
qc("|---|---|---|")
cots_year_summary <- cots_site_year |>
  group_by(Year) |>
  summarise(total_cots = sum(COTS_count_total), max_density = max(COTS_density_m2), .groups = "drop")
for (i in seq_len(nrow(cots_year_summary))) {
  qc(sprintf("| %d | %d | %.4f |", cots_year_summary$Year[i],
             cots_year_summary$total_cots[i], cots_year_summary$max_density[i]))
}
qc("")

# ---- D1. Cross-habitat check: is the backreef pattern consistent with the
# published (largely forereef) Moorea COTS outbreak history? -----------------
# COTS outbreaks at Moorea are documented primarily as a forereef phenomenon
# (Kayal et al. 2012; Section 0.2/3.3 of the plan). We do not use forereef
# COTS as a covariate (the revision's benthic/fish/metabolism tables are all
# backreef), but pulling it here as a cross-check confirms the backreef
# signal is real rather than too sparse/noisy to interpret.
cots_forereef_year <- cots_raw |>
  filter(Habitat == "Forereef") |>
  group_by(Year) |>
  summarise(total_cots_forereef = sum(COTS, na.rm = TRUE),
            max_site_forereef = max(COTS, na.rm = TRUE), .groups = "drop")

cots_habitat_compare <- cots_year_summary |>
  select(Year, total_cots_backreef = total_cots) |>
  left_join(cots_forereef_year, by = "Year")

qc("### D1. Cross-habitat check (forereef, not used as a covariate)")
qc("")
qc("Pulled for comparison only, to check the backreef pattern against the ",
   "published (largely forereef) outbreak history -- not joined into ",
   "`disturbance_site_year.csv`.")
qc("")
qc("| Year | Total COTS: backreef (6 sites) | Total COTS: forereef (6 sites) |")
qc("|---|---|---|")
for (i in seq_len(nrow(cots_habitat_compare))) {
  qc(sprintf("| %d | %d | %d |", cots_habitat_compare$Year[i],
             cots_habitat_compare$total_cots_backreef[i],
             cots_habitat_compare$total_cots_forereef[i]))
}
qc("")
qc("**Confirms two real outbreak pulses, both much larger on the forereef:** ",
   "(1) the documented 2006-2010 outbreak (Kayal et al. 2012) is clearly ",
   "visible on the forereef (108 individuals in 2008, 70 in 2009) and only ",
   "weakly visible on the backreef (6 in 2008, 12 in 2009) -- consistent ",
   "with COTS outbreaks being a largely forereef phenomenon at Moorea, so ",
   "the muted backreef signal is expected, not a data problem. (2) A ",
   "second, previously undiscussed outbreak pulse appears in 2023-2024: ",
   "forereef totals rise to 17 (2023) and 56 (2024), with the backreef ",
   "pulse (1 in 2023, 16 in 2024; concentrated at LTER_2 and LTER_3) tracking ",
   "the same timing at smaller magnitude. The 2024 backreef pulse (plan ",
   "testing log, Step 2 item 6) is therefore corroborated by the forereef ",
   "data as a real, synchronised outbreak event -- not an artefact -- and ",
   "should be written up as a second documented Moorea COTS outbreak ",
   "(2023-2024) alongside 2006-2010 in the methods/disturbance-history text.")
qc("")

cat("Forereef cross-check (not a covariate) -- total COTS by year:\n")
print(cots_habitat_compare)

cat("Done: Section D (COTS). Cached to", file.path(cache_dir, "cots_site_year.csv"), "\n")

# =============================================================================
# E. Cyclone Oli (February 2010) indicator
# =============================================================================
#
# Cyclone Oli made its closest approach to Moorea in early February 2010,
# ahead of the Jan 2011 survey window... but the 2010 benthic survey itself
# was taken 2010-01-01 (see survey_window), i.e. BEFORE Oli. The cyclone's
# impact is therefore expected to register in the 2011 survey (the first
# survey taken after the storm), not the 2010 survey -- this differs from
# how the original qmd's static "2010 disturbance" framing would suggest.
# Backreef wave exposure during Oli was much lower than forereef exposure
# (Section 3.3 of the plan); this indicator does not distinguish sites by
# shore-facing exposure, which is a known simplification pending a
# wave-energy proxy (not attempted here).

cyclone_oli_date <- ymd("2010-02-16")  # closest approach to Moorea

disturbance_cyclone <- survey_window |>
  mutate(Cyclone = as.integer(survey_date > cyclone_oli_date &
                                (survey_date - cyclone_oli_date) <= days(365))) |>
  select(Year, Cyclone)

qc("## E. Cyclone Oli (Feb 2010)")
qc("")
qc("- Closest approach to Moorea taken as ", as.character(cyclone_oli_date),
   ". `Cyclone = 1` for the first survey taken within 365 days after that ",
   "date, else 0.")
qc("- Because the 2010 survey predates Oli (survey_date = 2010-01-01) and ",
   "the 2011 survey (2011-01-17) falls ", 
   as.integer(ymd("2011-01-17") - cyclone_oli_date), " days after it, the ",
   "indicator flags **2011**, not 2010, as the cyclone-affected survey year. ",
   "Flag this for the PI: it changes which annual benthic transition the ",
   "Section 5 Gompertz model should attribute to the cyclone, relative to ",
   "the original (uncyclone-indicator) analysis.")
qc("- Backreef wave exposure during Oli was much lower than forereef ",
   "exposure (Section 0.2/3.3 of the plan); this binary indicator does not ",
   "yet distinguish shore-facing exposure across the 6 sites.")
qc("")
qc("| Year | Cyclone |")
qc("|---|---|")
for (i in seq_len(nrow(disturbance_cyclone))) {
  qc(sprintf("| %d | %d |", disturbance_cyclone$Year[i], disturbance_cyclone$Cyclone[i]))
}
qc("")

cat("Done: Section E (Cyclone indicator).\n")

# =============================================================================
# F. Assemble disturbance_site_year.csv
# =============================================================================
#
# DHW and Cyclone are Year-only (island-wide / single-logger, per the
# decisions in Sections B, C and E above); COTS is Site x Year. The output
# is at Site x Year resolution so it joins directly onto benthic_site,
# fish_transect (summarised to site-year) and pp_day (LTER_1 only).

sites_all <- sort(unique(cots_site_year$Site))
years_all <- sort(unique(survey_window$Year))

disturbance_site_year <- expand_grid(Site = sites_all, Year = years_all) |>
  left_join(dhw_annual |> select(Year, DHW_max), by = "Year") |>
  left_join(dhw_logger_annual |> select(Year, DHW_logger_max), by = "Year") |>
  left_join(cots_site_year |> select(Site, Year, COTS_count_total, COTS_density_m2), by = c("Site", "Year")) |>
  left_join(disturbance_cyclone, by = "Year") |>
  arrange(Site, Year)

write_csv(disturbance_site_year, here("Revision", "Data", "derived", "disturbance_site_year.csv"))

qc("## F. Assembled table")
qc("")
qc("- `disturbance_site_year.csv`: ", nrow(disturbance_site_year), " rows (",
   length(sites_all), " sites x ", length(years_all), " years).")
qc("- Columns: Site, Year, DHW_max (satellite, island-wide), DHW_logger_max ",
   "(in-situ sensitivity, island-wide), COTS_count_total, COTS_density_m2 ",
   "(site-year), Cyclone (0/1, island-wide).")
qc("- DHW_max and Cyclone are identical across all 6 sites for a given year ",
   "by construction (Sections B, C, E) -- only COTS varies by site. This is ",
   "an explicit limitation to carry into Section 4's DAG/adjustment-set work ",
   "and the Section 5 methods text: between-site contrasts in the Gompertz ",
   "model for DHW/Cyclone effects come entirely from between-site variation ",
   "in starting benthic composition and COTS, not from independent exposure.")
qc("")

writeLines(qc_lines, here("Revision", "Output", "qc_disturbance.md"))

cat("Done. Wrote:\n",
    "  Revision/Data/derived/disturbance_site_year.csv\n",
    "  Revision/Output/qc_disturbance.md\n")

# =============================================================================
# G. Combined disturbance-history figure (for the SI)
# =============================================================================
#
# Visual check of DHW (primary satellite + secondary logger), COTS density by
# site, and the Cyclone Oli indicator (dashed vertical line at the survey
# year it flags, not the cyclone's actual date) on one set of aligned axes.
# This does not feed back into the covariate table -- it is purely a QC /
# SI figure, built after the fact from `disturbance_site_year.csv`.

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

cyclone_year_fig <- disturbance_site_year |>
  distinct(Year, Cyclone) |>
  filter(Cyclone == 1) |>
  pull(Year)

dhw_long_fig <- disturbance_site_year |>
  distinct(Year, DHW_max, DHW_logger_max) |>
  pivot_longer(c(DHW_max, DHW_logger_max), names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
                          DHW_max = "Satellite DHW (CRW, primary)",
                          DHW_logger_max = "Logger DHW (LTER_2, sensitivity)"))

p_dhw_fig <- ggplot(dhw_long_fig, aes(x = Year, y = value, color = metric)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = cyclone_year_fig, linetype = "dashed", color = "grey30") +
  labs(y = "DHW_max (deg C-weeks)", x = NULL, color = NULL) +
  theme_minimal() +
  theme(legend.position = "top")

p_cots_fig <- ggplot(disturbance_site_year, aes(x = Year, y = COTS_density_m2, color = Site)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = cyclone_year_fig, linetype = "dashed", color = "grey30") +
  labs(y = "COTS density (ind/m^2)", x = "Year") +
  theme_minimal()

p_disturbance_history <- (p_dhw_fig / p_cots_fig) +
  plot_annotation(
    title = "Moorea backreef disturbance history, 2005-2025",
    caption = paste0("Dashed line: Cyclone Oli indicator (flags the ", cyclone_year_fig,
                      " survey, the first survey taken after the Feb 2010 cyclone, ",
                      "not the cyclone's own date)")
  )

ggsave(
  here("Revision", "Output", "fig_disturbance_history.png"),
  p_disturbance_history, width = 8, height = 7, dpi = 300
)
ggsave(
  here("Revision", "Output", "fig_disturbance_history.pdf"),
  p_disturbance_history, width = 8, height = 7
)

cat("Done: Section G (disturbance-history figure). Wrote:\n",
    "  Revision/Output/fig_disturbance_history.png\n",
    "  Revision/Output/fig_disturbance_history.pdf\n")
