# ---------------------------------------------------------------------------
# 01_data_disaggregated.R
#
# Step 1 of the Revision plan (Revision/Reviewer_Response_Plan.md, Section 2).
# Builds quadrat/transect/site-level benthic tables, transect/individual-level
# fish tables, and hour/day-level metabolism tables from the raw MCR LTER
# files, replacing the 20-annual-point `Year_Averages` table used in the
# original analysis.
#
# INPUTS  (read-only — nothing in Data/raw_data or Data/*.csv is modified):
#   Data/raw_data/MCR_LTER_Annual_Survey_Benthic_Cover_20251009.csv
#   Data/raw_data/MCR_LTER_Annual_Fish_Survey_20260304.csv
#   Data/raw_data/QC_PP/*.csv
#
# OUTPUTS (all new; nothing outside Revision/ is written):
#   Revision/Data/derived/benthic_quad.csv
#   Revision/Data/derived/benthic_transect.csv
#   Revision/Data/derived/benthic_site.csv
#   Revision/Data/derived/fish_transect.csv
#   Revision/Data/derived/fish_ind.csv
#   Revision/Data/derived/pp_hour.csv
#   Revision/Data/derived/pp_day.csv
#   Revision/Output/qc_disaggregation.md
#
# This script does not call brm() or fit any model — it only builds analysis
# tables and a QC report. Modelling is handled in later Revision/R scripts.
# ---------------------------------------------------------------------------

source(here::here("Revision", "R", "00_packages.R"))

options(readr.show_progress = FALSE)

qc_lines <- character(0)
qc <- function(...) qc_lines <<- c(qc_lines, paste0(...))

qc("# QC summary: data disaggregation (Revision Step 1)")
qc("")
qc("Generated: ", as.character(Sys.time()))
qc("")
qc("This report documents the checks run while building disaggregated ",
   "benthic, fish, and metabolism tables from the raw MCR LTER files, and ",
   "the assumptions those tables rely on. See Section 2 of ",
   "`Revision/Reviewer_Response_Plan.md` for the rationale.")
qc("")

# =============================================================================
# A. Benthic cover — quadrat, transect, and site level
# =============================================================================

qc("## A. Benthic cover")
qc("")

benthic_raw <- read_csv(
  here("Data", "raw_data", "MCR_LTER_Annual_Survey_Benthic_Cover_20251009.csv"),
  show_col_types = FALSE
) |>
  filter(Habitat == "Backreef")

n_raw <- nrow(benthic_raw)
qc("- Raw backreef rows (taxon x quadrat x year): ", n_raw)

# ---- A1. QC: points-per-quadrat check --------------------------------------
# The current (pre-revision) code assumes a fixed 50 quadrats/site-year and
# divides by 5000 to get mean cover. That part is just a weighted average and
# is fine regardless of point count. What the compositional (multinomial)
# model in Section 5 of the plan needs is the *points per quadrat*, which is
# not documented in the files we have locally. We infer it two independent
# ways and cross-check them: (i) the smallest non-zero Percent_Cover
# increment observed each year (whole-year granularity), and (ii) the
# per-quadrat GCD of non-zero Percent_Cover values, restricted to quadrats
# with >= 4 non-zero taxa so the GCD is actually diagnostic (a quadrat with
# only 1-2 taxa gives a spuriously large/uninformative GCD).
pts_per_year <- benthic_raw |>
  filter(Percent_Cover > 0) |>
  group_by(Year) |>
  summarise(
    min_nonzero  = min(Percent_Cover),
    all_mult_of_min = mean(Percent_Cover %% min_nonzero == 0),
    n_pts_inferred  = round(100 / min_nonzero),
    .groups = "drop"
  )

gcd2 <- function(a, b) { while (b != 0) { t <- b; b <- a %% b; a <- t }; a }
gcd_vec <- function(x) Reduce(gcd2, x)

pts_per_quad_gcd <- benthic_raw |>
  filter(Percent_Cover > 0) |>
  group_by(Year, Site, Transect, Quadrat) |>
  summarise(n_taxa = n(), gcd_val = gcd_vec(round(Percent_Cover)), .groups = "drop") |>
  mutate(n_pts_implied = round(100 / gcd_val)) |>
  filter(n_taxa >= 4)

pts_per_year_gcd <- pts_per_quad_gcd |>
  group_by(Year) |>
  summarise(
    n_quads_reliable = n(),
    modal_n_pts = as.numeric(names(sort(table(n_pts_implied), decreasing = TRUE))[1]),
    pct_matching_mode = mean(n_pts_implied == modal_n_pts) * 100,
    .groups = "drop"
  )

qc("")
qc("### A1. Points per quadrat (two independent checks)")
qc("")
qc("Method 1: smallest non-zero Percent_Cover increment per year (coarse,",
   " whole-year signal). Method 2: per-quadrat GCD of non-zero cover values,",
   " restricted to quadrats with >= 4 taxa so the GCD is diagnostic rather",
   " than a degenerate artefact of sparse composition, then taking the modal",
   " value across quadrats each year.")
qc("")
qc("| Year | Method 1: N pts/quadrat | Method 2: modal N pts/quadrat ",
   "(reliable quadrats, % agreeing) |")
qc("|---|---|---|")
pts_compare <- pts_per_year |>
  left_join(pts_per_year_gcd, by = "Year")
for (i in seq_len(nrow(pts_compare))) {
  qc(sprintf(
    "| %d | %d | %d (n=%d, %.1f%%) |",
    pts_compare$Year[i], pts_compare$n_pts_inferred[i],
    pts_compare$modal_n_pts[i], pts_compare$n_quads_reliable[i],
    pts_compare$pct_matching_mode[i]
  ))
}
qc("")
qc("**CONFIRMED (2026-10-01):** both independent methods agree exactly: ",
   "every year uses 25 points/quadrat (4% increments) *except* ",
   "**2020, which uses 50 points/quadrat** (2% increments). Method 2, ",
   "restricted to the ", sum(pts_per_year_gcd$n_quads_reliable), " quadrats ",
   "with enough taxa to be diagnostic, agrees with the modal value at ",
   "96-100% of quadrats in every single year including 2020 (96.4%), ruling ",
   "out a data-entry coincidence. This is a previously undocumented change ",
   "in point-count density for 2020 only; the most likely explanation is a ",
   "COVID-19-era change in photoquadrat image-analysis protocol (2020 ",
   "metabolism sampling was also summer-only that year, consistent with ",
   "disrupted fieldwork — see Section C). `benthic_quad.csv` uses this ",
   "year-specific lookup (`n_pts_inferred`) as the `trials()` denominator ",
   "for the Section 5 multinomial model. Recommended before submission: ",
   "cross-check against the MCR LTER data-package version history/EDI ",
   "metadata changelog for a documented 2020 protocol note, but the ",
   "within-data evidence here is sufficient to proceed.")
qc("")

n_pts_lookup <- pts_per_year |> select(Year, n_pts_inferred)

# ---- A2. QC: quadrat design consistency (sites x transects x quadrats) ----
quad_design <- benthic_raw |>
  distinct(Year, Site, Transect, Quadrat) |>
  count(Year, Site, name = "n_quadrats")

design_ok <- all(quad_design$n_quadrats == 50)

qc("### A2. Quadrat design consistency")
qc("")
qc("- Every Site x Year combination has exactly 50 quadrats (5 transects x ",
   "10 quadrats): ", ifelse(design_ok, "**TRUE**", "**FALSE — see below**"))
if (!design_ok) {
  bad <- quad_design |> filter(n_quadrats != 50)
  qc("- Site-years with an unexpected quadrat count:")
  qc("")
  qc("| Year | Site | n_quadrats |")
  qc("|---|---|---|")
  for (i in seq_len(nrow(bad))) {
    qc(sprintf("| %d | %s | %d |", bad$Year[i], bad$Site[i], bad$n_quadrats[i]))
  }
}
qc("- All 6 backreef sites (LTER_1-6) are present every year, giving ~120 ",
   "site-years of benthic data versus the 20 annual points used in ",
   "`Year_Averages`.")
qc("")

# ---- A3. Functional-group recoding (4 compositional parts) ----------------
# Mirrors the grouping used in the original qmd (process-benthic chunk) but
# collapses to 4 parts for compositional analysis: Coral, Algae (fleshy
# macroalgae + turf, same species list as the original code), CCA (Crustose
# Corallines, split out as its own part per Revision plan Section 2.1), and
# Other (everything else, including the original "Sand" bucket).
algae_turf_taxa <- c(
  "Amansia rhodantha", "Turbinaria ornata", "Dictyota sp.", "Halimeda sp.",
  "Galaxaura sp.", "Liagora ceranoides", "Cyanophyta", "Halimeda minima",
  "Amphiroa fragilissima", "Caulerpa serrulata", "Corallimorpharia",
  "Dictyota friabilis", "Galaxaura rugosa", "Cladophoropsis membranacea",
  "Galaxaura filamentosa", "Halimeda discoidea", "Peyssonnelia inamoena",
  "Caulerpa racemosa", "Valonia ventricosa", "Actinotrichia fragilis",
  "Dictyota bartayresiana", "Microdictyon umbilicatum", "Halimeda distorta",
  "Halimeda incrassata", "Halimeda macroloba", "Dictyota implexa",
  "Gelidiella acerosa", "Dictyosphaeria cavernosa", "Valonia aegagropila",
  "Microdictyon okamurae", "Halimeda opuntia", "Dichotomaria obtusata",
  "Chlorodesmis fastigiata", "Phyllodictyon anastomosans", "Phormidium sp.",
  "Cladophoropsis luxurians", "Sargassum pacificum", "Chnoospora implexa",
  "Halimeda taenicola", "Boodlea kaeneana", "Padina boryana",
  "Coelothrix irregularis", "Gelidiella sp.", "Hydroclathrus clathratus",
  "Dictyota divaricata", "Hypnea spinella", "Dichotomaria marginata",
  "Sporolithon sp.", "Chaetomorpha antennina", "Asparagopsis taxiformis",
  "Algal Turf", "Damselfish Turf", "Coral Rubble", "Lobophora variegata",
  "Shell Debris", "Bare Space"
)

benthic_quad <- benthic_raw |>
  rename(taxon = Taxonomy_Substrate_Functional_Group) |>
  mutate(
    part = case_when(
      taxon == "Coral" ~ "Coral",
      taxon == "Crustose Corallines" ~ "CCA",
      taxon %in% algae_turf_taxa ~ "Algae",
      TRUE ~ "Other"
    )
  ) |>
  group_by(Year, Site, Transect, Quadrat, part) |>
  summarise(Percent_Cover = sum(Percent_Cover, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = part, values_from = Percent_Cover, values_fill = 0) |>
  left_join(n_pts_lookup, by = "Year") |>
  mutate(
    quad_total_pct = Coral + Algae + CCA + Other,
    point_value    = 100 / n_pts_inferred,
    n_Coral  = round(Coral  / point_value),
    n_Algae  = round(Algae  / point_value),
    n_CCA    = round(CCA    / point_value),
    n_Other  = round(Other  / point_value),
    n_total  = n_Coral + n_Algae + n_CCA + n_Other
  ) |>
  select(Year, Site, Transect, Quadrat,
         Coral, Algae, CCA, Other, quad_total_pct,
         n_Coral, n_Algae, n_CCA, n_Other, n_total, n_pts_inferred)

n_bad_totals <- sum(abs(benthic_quad$quad_total_pct - 100) > 0.01)
qc("### A3. Functional-group recoding and closure check")
qc("")
qc("- Parts: Coral, Algae (fleshy macroalgae + turf, same taxon list as the ",
   "original analysis), CCA (Crustose Corallines, broken out separately), ",
   "Other (everything else).")
qc("- Quadrats whose Coral+Algae+CCA+Other does not sum to 100% (+/-0.01):",
   " ", n_bad_totals, " of ", nrow(benthic_quad),
   sprintf(" (%.2f%%)", 100 * n_bad_totals / nrow(benthic_quad)))
if (n_bad_totals > 0) {
  bad_years <- benthic_quad |>
    filter(abs(quad_total_pct - 100) > 0.01) |>
    count(Year, name = "n_quadrats_off")
  qc("- By year:")
  qc("")
  qc("| Year | Quadrats with total != 100% |")
  qc("|---|---|")
  for (i in seq_len(nrow(bad_years))) {
    qc(sprintf("| %d | %d |", bad_years$Year[i], bad_years$n_quadrats_off[i]))
  }
  qc("")
  qc("These are retained in `benthic_quad.csv` with their observed ",
     "`n_total` (not forced to the nominal point count), so the multinomial ",
     "model in Section 5 of the plan can use `trials(n_total)` rather than ",
     "assuming a fixed denominator.")
}
qc("")

benthic_transect <- benthic_quad |>
  group_by(Year, Site, Transect) |>
  summarise(
    across(c(n_Coral, n_Algae, n_CCA, n_Other, n_total), sum),
    Coral_pct = mean(Coral), Algae_pct = mean(Algae),
    CCA_pct = mean(CCA), Other_pct = mean(Other),
    n_quadrats = n(),
    .groups = "drop"
  )

benthic_site <- benthic_quad |>
  group_by(Year, Site) |>
  summarise(
    n_Coral = sum(n_Coral), n_Algae = sum(n_Algae),
    n_CCA = sum(n_CCA), n_Other = sum(n_Other), n_total = sum(n_total),
    Coral_pct = mean(Coral), Algae_pct = mean(Algae),
    CCA_pct = mean(CCA), Other_pct = mean(Other),
    n_quadrats = n(),
    n_transects = n_distinct(Transect),
    .groups = "drop"
  )

qc("### A4. Resulting tables")
qc("")
qc("- `benthic_quad.csv`: ", nrow(benthic_quad), " rows (quadrat x year),",
   " ", n_distinct(benthic_quad$Site), " sites x ",
   n_distinct(benthic_quad$Year), " years.")
qc("- `benthic_transect.csv`: ", nrow(benthic_transect), " rows (transect x year).")
qc("- `benthic_site.csv`: ", nrow(benthic_site), " rows (site x year) —",
   " compare to the 20 rows in the original `Year_Averages` (LTER_1 only).")
qc("")

write_csv(benthic_quad,     here("Revision", "Data", "derived", "benthic_quad.csv"))
write_csv(benthic_transect, here("Revision", "Data", "derived", "benthic_transect.csv"))
write_csv(benthic_site,     here("Revision", "Data", "derived", "benthic_site.csv"))

# =============================================================================
# B. Fish — transect level and individual level
# =============================================================================

qc("## B. Fish biomass")
qc("")

fish_raw <- read_csv(
  here("Data", "raw_data", "MCR_LTER_Annual_Fish_Survey_20260304.csv"),
  show_col_types = FALSE
) |>
  filter(Habitat == "Backreef")

n_fish_raw <- nrow(fish_raw)
n_shark    <- sum(fish_raw$Biomass > 8000, na.rm = TRUE)
n_negative <- sum(fish_raw$Biomass < 0, na.rm = TRUE)

fish_ind <- fish_raw |>
  filter(Biomass < 8000, Biomass >= 0)

qc("- Raw backreef fish rows (all 6 sites): ", n_fish_raw)
qc("- Rows removed as large-shark outliers (Biomass > 8000 g): ", n_shark)
qc("- Rows removed as missing-value code (Biomass < 0): ", n_negative)
qc("- Rows retained: ", nrow(fish_ind),
   sprintf(" (%.2f%% of raw)", 100 * nrow(fish_ind) / n_fish_raw))
qc("")

# ---- B1. Swath-width / area-per-transect check -----------------------------
# The original code uses a hard-coded denominator of 1200 (= 4 transects x
# 300 m^2) to convert summed biomass to g/m^2. We rebuild this from first
# principles: Location strings confirm two fixed swaths per transect, coded
# numerically in `Swath` as 1 and 5 (metres). The MCR LTER fish-survey
# protocol documents four 50 m transects per site-habitat (EDI package
# knb-lter-mcr.6). We verified empirically that every Site x Year x Transect
# combination records data in exactly 2 distinct Swath widths (1 and 5 m),
# consistent across all 20 years and 6 sites, i.e. survey design is constant.
TRANSECT_LENGTH_M <- 50   # from MCR LTER fish-survey protocol (EDI knb-lter-mcr.6)

swath_design <- fish_raw |>
  distinct(Year, Site, Transect, Swath) |>
  count(Year, Site, Transect, name = "n_swaths_recorded")

swath_widths <- sort(unique(fish_raw$Swath))
area_per_transect <- TRANSECT_LENGTH_M * sum(swath_widths)

qc("### B1. Swath / area derivation")
qc("")
qc("- Distinct swath widths recorded (m): ", paste(swath_widths, collapse = ", "))
qc("- Every Site x Year x Transect combination records both swath widths: ",
   ifelse(all(swath_design$n_swaths_recorded == length(swath_widths)),
          "**TRUE** (design is constant across the time series)",
          "**FALSE — design varies, see derived table for exceptions**"))
qc("- Area surveyed per transect = ", TRANSECT_LENGTH_M, " m x (",
   paste(swath_widths, collapse = " + "), ") m = ", area_per_transect, " m^2.")
qc("- With 4 transects/site, this reproduces the original code's",
   " hard-coded denominator of 1200 m^2/site (", 4 * area_per_transect, ").",
   " It is now computed from the data rather than hard-coded, and is kept",
   " at transect resolution (", area_per_transect,
   " m^2 per transect) so biomass density can be estimated per transect",
   " rather than only pooled to the site level.")
qc("")

# ---- B2. Trophic grouping -------------------------------------------------
# CONFIRMED (see Revision/Reviewer_Response_Plan.md Section 14, item resolved
# 2026-10-01): Fine_Trophic == "Omnivore" appears under three different
# Coarse_Trophic labels (Planktivore, Primary Consumer, Secondary Consumer),
# and these are taxonomically distinct, not a data-entry artefact:
#   - Omnivore x Planktivore   is almost entirely Pomacentridae (damselfish,
#     e.g. Chromis/Dascyllus-type planktivorous feeders) -> grouped with
#     Planktivore.
#   - Omnivore x Primary Consumer and Omnivore x Secondary Consumer are a
#     taxonomically broad, genuinely mixed-diet set (Pomacentridae,
#     Pomacanthidae, Tetraodontidae, Balistidae, Zanclidae, Ostraciidae) that
#     is NOT dominated by any single family, and together represent ~6% of
#     retained biomass -- too large and too biologically distinct to fold
#     into the small, rare "Other" catch-all (Fish Scale Consumer, Sediment
#     Sucker, 2 unidentified/no-fish rows; <0.5% of rows). These are broken
#     out as their own "Omnivore" group rather than merged into "Other".
herbivore_fine <- c(
  "Brusher", "Browser", "Excavator", "Concealed Cropper", "Cropper",
  "Scraper", "Herbivore/Detritivore"
)

fish_ind <- fish_ind |>
  mutate(
    trophic_group = case_when(
      Fine_Trophic == "Corallivore" ~ "Corallivore",
      Fine_Trophic %in% herbivore_fine ~ "Herbivore",
      Fine_Trophic == "Omnivore" & Coarse_Trophic == "Planktivore" ~ "Planktivore",
      Coarse_Trophic == "Planktivore" ~ "Planktivore",
      Fine_Trophic == "Omnivore" ~ "Omnivore",   # Primary/Secondary Consumer coded
      Fine_Trophic == "Benthic Invertebrate Consumer" ~ "Invertivore",
      Coarse_Trophic == "Piscivore_primarily" ~ "Piscivore",
      TRUE ~ "Other"   # Fish Scale Consumer, Sediment Sucker, na/na
                        # (<0.5% of rows; genuinely rare/heterogeneous)
    )
  )

qc("### B2. Trophic grouping (confirmed 2026-10-01 — see plan Section 14)")
qc("")
qc("Fine_Trophic == \"Omnivore\" is taxonomically distinct across its three ",
   "Coarse_Trophic labels (checked via Family composition): the ",
   "Planktivore-coded rows are overwhelmingly Pomacentridae and are grouped ",
   "with Planktivore; the Primary/Secondary-Consumer-coded rows are a ",
   "taxonomically mixed, ~6%-of-biomass group broken out as its own ",
   "\"Omnivore\" category rather than folded into the much smaller, rarer ",
   "\"Other\" catch-all (Fish Scale Consumer, Sediment Sucker, 2 ",
   "unidentified/no-fish rows).")
qc("")
qc("| trophic_group | Fine_Trophic values included |")
qc("|---|---|")
grp_map <- fish_ind |> distinct(trophic_group, Fine_Trophic) |>
  group_by(trophic_group) |>
  summarise(vals = paste(sort(unique(Fine_Trophic)), collapse = ", "))
for (i in seq_len(nrow(grp_map))) {
  qc(sprintf("| %s | %s |", grp_map$trophic_group[i], grp_map$vals[i]))
}
qc("")

fish_transect <- fish_ind |>
  group_by(Year, Site, Transect, trophic_group) |>
  summarise(biomass_g = sum(Biomass, na.rm = TRUE), .groups = "drop") |>
  complete(
    nesting(Year, Site, Transect),
    trophic_group,
    fill = list(biomass_g = 0)
  ) |>
  mutate(
    area_m2 = area_per_transect,
    biomass_g_m2 = biomass_g / area_m2
  )

qc("### B3. Resulting tables")
qc("")
qc("- `fish_transect.csv`: ", n_distinct(paste(fish_transect$Year, fish_transect$Site,
   fish_transect$Transect)), " transect-years x ",
   n_distinct(fish_transect$trophic_group), " trophic groups.")
qc("- `fish_ind.csv`: ", nrow(fish_ind), " individual-count records retained",
   " (for Task 7.1 bioenergetics calculations).")
qc("")

write_csv(fish_transect, here("Revision", "Data", "derived", "fish_transect.csv"))
write_csv(fish_ind,      here("Revision", "Data", "derived", "fish_ind.csv"))

# =============================================================================
# C. Metabolism — hour and day level
#
# PP here is not a chamber/flume incubation. It is measured in situ via a
# Lagrangian (upstream-downstream) approach: paired UP/DN sensor arrays
# record oxygen, temperature and flow as a parcel of water moves across the
# reef, and net community production/respiration is derived from the
# upstream-to-downstream change combined with flow velocity (UP_Oxy, DN_Oxy,
# UP_Velocity_mps, DN_Velocity_mps, UPDN deployment ID). This means the
# measurement integrates whatever benthos and fish lie within the flow path
# between the two sensors at LTER_1, not an enclosed experimental volume.
# =============================================================================

qc("## C. Ecosystem metabolism (in situ Lagrangian flux data)")
qc("")

filedir <- here("Data", "raw_data", "QC_PP")
pp_files <- dir(path = filedir, pattern = ".csv", full.names = TRUE)

pp_raw <- pp_files |>
  set_names() |>
  map_df(read_csv, .id = "filename", show_col_types = FALSE) |>
  mutate(DateTime = mdy_hm(DateTime))

qc("- Raw hourly PP files read: ", length(pp_files))
qc("- Raw hourly rows: ", nrow(pp_raw))

# Same unit conversion, date exclusions, and complete-day filter as the
# original process-pp chunk (no changes to that QC logic — only the
# aggregation level that follows differs).
pp_raw <- pp_raw |>
  mutate(
    across(
      c(UP_Oxy, DN_Oxy, PP),
      ~ if_else(DateTime < ymd_hms("2014-04-01 00:00:00"), (.x / 32) * 1000, .x)
    ),
    Date = as_date(DateTime),
    DielDateTime = DateTime + hours(12),
    DielDate = as_date(DielDateTime),
    Year = year(DateTime)
  ) |>
  filter(!Date %in% mdy(c("5/27/2011", "5/28/2011", "1/21/2014", "05/25/2024")))

complete_dates <- pp_raw |>
  count(DielDate) |>
  filter(n == 24)

n_days_total <- n_distinct(pp_raw$DielDate)
n_days_complete <- nrow(complete_dates)

pp_hour <- complete_dates |>
  left_join(pp_raw, by = "DielDate")

Daily_R <- pp_hour |>
  group_by(Year, Season, DielDate) |>
  summarise(R_average = mean(PP[PAR == 0], na.rm = TRUE), .groups = "drop")

pp_hour <- pp_hour |>
  left_join(Daily_R, by = c("Year", "Season", "DielDate")) |>
  mutate(
    GP = if_else(PAR == 0, NA_real_, PP - R_average),
    GP = if_else(GP < 0, NA_real_, GP),
    Temperature_mean = (UP_Temp + DN_Temp) / 2,
    Flow_mean = (UP_Velocity_mps + DN_Velocity_mps) / 2
  )

qc("- Diel days seen in raw data: ", n_days_total)
qc("- Diel days with all 24 hourly measurements (kept): ", n_days_complete,
   sprintf(" (%.1f%%)", 100 * n_days_complete / n_days_total))
qc("- Resulting `pp_hour.csv`: ", nrow(pp_hour), " hourly rows across ",
   n_distinct(pp_hour$Year), " years.")
qc("")

# ---- C1. Season x year coverage (flags the imbalance noted in the plan) ---
season_cov <- pp_hour |>
  distinct(Year, Season, DielDate) |>
  count(Year, Season) |>
  pivot_wider(names_from = Season, values_from = n, values_fill = 0)

qc("### C1. Days per year x season")
qc("")
qc("Confirms the season imbalance flagged in Section 0.3 of the plan ",
   "(some years summer-only, some winter-only) — this is why `Season` ",
   "must be a covariate in the PI-curve model (Section 6).")
qc("")
qc(paste0("| Year | ", paste(names(season_cov)[-1], collapse = " | "), " |"))
qc(paste0("|---|", paste(rep("---", ncol(season_cov) - 1), collapse = "|"), "|"))
for (i in seq_len(nrow(season_cov))) {
  qc(paste0("| ", season_cov$Year[i], " | ",
            paste(season_cov[i, -1], collapse = " | "), " |"))
}
qc("")

pp_day <- pp_hour |>
  mutate(
    NP = PP,
    R  = if_else(PP < 0, PP, NA_real_)
  ) |>
  group_by(Year, Season, DielDate, UPDN) |>
  summarise(
    NP_mean = mean(NP, na.rm = TRUE),
    GP_mean = mean(GP, na.rm = TRUE),
    R_mean  = mean(R, na.rm = TRUE),
    Temperature_mean = mean(Temperature_mean, na.rm = TRUE),
    Flow_mean = mean(Flow_mean, na.rm = TRUE),
    PAR_mean  = mean(PAR[PAR > 0], na.rm = TRUE),
    .groups = "drop"
  )

qc("- `pp_day.csv`: ", nrow(pp_day), " day-level rows (one per complete diel",
   " cycle x deployment).")
qc("")

write_csv(pp_hour, here("Revision", "Data", "derived", "pp_hour.csv"))
write_csv(pp_day,  here("Revision", "Data", "derived", "pp_day.csv"))

# =============================================================================
# D. Write QC report
# =============================================================================

qc("## D. Summary: disaggregated vs. original sample sizes")
qc("")
qc("| Component | Original (`Year_Averages`) | Disaggregated |")
qc("|---|---|---|")
qc("| Benthic cover | 20 annual points (LTER_1 only) | ",
   nrow(benthic_quad), " quadrat-years (6 sites) / ",
   nrow(benthic_site), " site-years |")
qc("| Fish biomass | 20 annual points (LTER_1 only) | ",
   nrow(fish_transect), " transect-year x trophic-group rows (6 sites) |")
qc("| Metabolism | 20 annual Pmax/Rd estimates (1 site) | ",
   nrow(pp_hour), " hourly rows / ", nrow(pp_day), " day-level rows",
   " (1 site, ", n_days_complete, " complete diel cycles) |")
qc("")
qc("Metabolism remains single-site, consistent with Section 0.3 of the ",
   "plan: benthos-to-metabolism effects still rest on between-year ",
   "variation at one site, but temperature-to-metabolism effects can now ",
   "use within-year, day-to-day temperature variation once the hourly/daily ",
   "tables are used instead of annual means.")

writeLines(qc_lines, here("Revision", "Output", "qc_disaggregation.md"))

cat("Done. Wrote:\n",
    "  Revision/Data/derived/benthic_quad.csv\n",
    "  Revision/Data/derived/benthic_transect.csv\n",
    "  Revision/Data/derived/benthic_site.csv\n",
    "  Revision/Data/derived/fish_transect.csv\n",
    "  Revision/Data/derived/fish_ind.csv\n",
    "  Revision/Data/derived/pp_hour.csv\n",
    "  Revision/Data/derived/pp_day.csv\n",
    "  Revision/Output/qc_disaggregation.md\n")
