# ---------------------------------------------------------------------------
# fish_respiration.R
#
# Standalone estimate of fish community respiration (mmol O2 m^-2 h^-1) for
# the MCR LTER backreef, for import into MCR_Productivity_Analysis.qmd.
#
# Method: individual-level bioenergetics following fishflux (Schiettekatte et
# al. 2020, Funct. Ecol.) with resting-metabolism parameters from Barneche &
# Allen (2018, Ecol. Lett., Model 2):
#
#   B0      = b0 * troph^0.81 * exp(Er * (1/kTc - 1/kT)) / (1 + Er/(Ei-Er) * exp(Ei * (1/kTopt - 1/kT)))
#   B_rest  = B0 * m_max^(alpha - 1) * m  +  Em * growth           [fishflux::metabolic_rate]
#   B_field = f * B_rest                                            [activity scope f]
#
# B0 is computed with a vectorised copy of fishflux's internal getB0 Stan
# function (verified identical to fishflux::metabolism()$f0_m), so that
# temperature varies by YEAR and trophic level varies by SPECIES. Calling
# fishflux::metabolism() once per family at a single temperature cannot do
# either.
#
# Corrections relative to Revision/R/05_fish.R:
#   1. B0: 05_fish.R passed metabolism()$b0_m (the raw normalisation constant,
#      before temperature and trophic-level adjustment) into metabolic_rate().
#      fishflux's own model_parameters() passes f0_m (the adjusted value).
#      Using b0_m underestimates resting metabolism ~4.7x at 28 degC.
#   2. Survey area: MCR fishes are counted on two swaths per 50 m transect --
#      5 m wide (250 m^2; mobile fishes) and 1 m wide (50 m^2; small,
#      site-attached fishes, mostly Pomacentridae/Cirrhitidae). Each record is
#      scaled by its own swath area rather than pooling all fish over 300 m^2.
#   3. Exclusions: unchanged from the original analysis -- survey rows with
#      Biomass >= 8000 g are dropped. This removes large transient sharks and
#      outlier schools (e.g. 750 Acanthurus triostegus, LTER_1 2008). Note the
#      threshold applies to ROW totals, so sharks recorded in rows < 8000 g
#      are retained.
#   4. Temperature: in-situ backreef logger (annual mean of daily means)
#      rather than island-wide satellite SST.
#   5. Taxonomy: MCR uses the revised Julidinae genera (e.g. Gomphosus for
#      Thalassoma, Hemitautoga for Halichoeres) that FishBase does not yet
#      recognise; these are crosswalked to FishBase names, and FishBase
#      synonyms are resolved to accepted names.
#   6. Loo: FishBase Loo is reported as SL, FL or TL; survey lengths are TL.
#      Loo is converted to TL (FishBase TLinfinity, else type-specific ratio).
#   7. Length-weight: species-level fits (which recover the MCR biomass
#      equations exactly), then family, then global.
#   8. Transects with zero retained fish are kept as zeros; records with the
#      missing-biomass code (-1) but a valid length are retained, with mass
#      from the length-weight fit.
#
# Remaining approximations (documented, not species-specific):
#   - Activity scope f = 2 (fishflux default). A resting-only (f = 1) column is
#     also written; nocturnal Rd may be closer to resting fish metabolism.
#   - RQ = 0.8 to convert respired carbon to O2 (see the RQ parameter below).
#   - Family-level alpha/b0 exist in fishflux only for a few families; others
#     use the across-family mean (fishflux's documented behaviour).
#
# Output: Data/fish_respiration.csv  (one row per Site x Year)
# ---------------------------------------------------------------------------

library(tidyverse)
library(here)
library(rfishbase)
library(fishflux)

# ---- Parameters -----------------------------------------------------------

HABITAT           <- "Backreef"
TRANSECT_LENGTH_M <- 50
ACTIVITY_SCOPE    <- 2
C_MOLAR_MASS      <- 12.011

# Respiratory quotient (RQ = mol CO2 produced / mol O2 consumed), used to
# convert fishflux's respired carbon (g C/day) to O2 consumption.
# Substrate-specific RQs are ~1.0 for carbohydrate, ~0.7 for lipid and ~0.8
# for protein. Fish fuel aerobic metabolism mainly with protein and lipid,
# so measured RQs in resting/routine fish are typically below 1, mostly
# ~0.8-0.9 (reviewed in Wang et al. 2021, Rev. Aquacult.,
# doi:10.1111/raq.12522; e.g. 0.82-0.88 in zebrafish, Ferreira et al. 2019,
# J. Exp. Biol.). RQ = 1 is also commonly assumed (e.g. after Kutty 1968,
# J. Fish. Res. Board Can. 25:1689-1728, doi:10.1139/f68-150), which
# would lower fish O2 estimates by 20%.
# The ideal value is whatever factor Barneche & Allen (2018) used to convert
# measured O2 consumption to carbon units for B0 (not verified here).
RQ <- 0.8

# Survey rows at or above this total biomass (g) are outliers (large sharks, outlier schools)
MAX_ROW_BIOMASS_G <- 8000

# MCR (revised Julidinae genera) -> FishBase accepted names
TAXON_CROSSWALK <- c(
  "Gomphosus amblycephalum"   = "Thalassoma amblycephalum",
  "Gomphosus hardwicke"       = "Thalassoma hardwicke",
  "Gomphosus lutescens"       = "Thalassoma lutescens",
  "Gomphosus purpureum"       = "Thalassoma purpureum",
  "Gomphosus quinquevittatum" = "Thalassoma quinquevittatum",
  "Gomphosus trilobatum"      = "Thalassoma trilobatum",
  "Allocoris gaimard"         = "Coris gaimard",
  "Hemicoris claudia"         = "Halichoeres claudia",
  "Hemicoris margaritaceus"   = "Halichoeres margaritaceus",
  "Hemitautoga hortulana"     = "Halichoeres hortulanus",
  "Hemitautoga trimaculata"   = "Halichoeres trimaculatus",
  "Hemiulis marginatus"       = "Halichoeres marginatus"
)

cache_dir <- here("Data", "fishbase_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Fish survey data ----------------------------------------------------

fish_raw <- read_csv(
  here("Data", "raw_data", "MCR_LTER_Annual_Fish_Survey_20260304.csv"),
  show_col_types = FALSE
) |>
  filter(Habitat == HABITAT)

stopifnot(all(fish_raw$Swath %in% c(1, 5)))

# Full survey design, taken before any filtering so empty transects stay as 0
survey_transects <- fish_raw |> distinct(Site, Year, Transect)

fish <- fish_raw |>
  filter(
    Biomass < MAX_ROW_BIOMASS_G,
    Count > 0,
    Total_Length > 0
  ) |>
  mutate(
    swath_area_m2 = TRANSECT_LENGTH_M * Swath,
    fb_name = coalesce(unname(TAXON_CROSSWALK[Taxonomy]), Taxonomy)
  )

message(sprintf(
  "Fish records: %d raw backreef, %d retained (%d rows >= %d g excluded, %d rows without a valid length/count).",
  nrow(fish_raw), nrow(fish),
  sum(fish_raw$Biomass >= MAX_ROW_BIOMASS_G), MAX_ROW_BIOMASS_G,
  sum(fish_raw$Biomass < MAX_ROW_BIOMASS_G) - nrow(fish)
))

# ---- 2. Length-weight relationships (species -> family -> global) -----------

fit_lw <- function(df) {
  if (n_distinct(df$Total_Length) < 2) return(tibble(lw_a = NA_real_, lw_b = NA_real_))
  cf <- unname(coef(lm(log(mass_g) ~ log(Total_Length), data = df)))
  tibble(lw_a = exp(cf[1]), lw_b = cf[2])
}

lw_data <- fish |> filter(Biomass > 0) |> mutate(mass_g = Biomass / Count)

lw_species <- lw_data |> group_by(Taxonomy) |> group_modify(~ fit_lw(.x)) |> ungroup() |>
  rename(lw_a_sp = lw_a, lw_b_sp = lw_b)
lw_family <- lw_data |> group_by(Family) |> group_modify(~ fit_lw(.x)) |> ungroup() |>
  rename(lw_a_fam = lw_a, lw_b_fam = lw_b)
lw_global <- fit_lw(lw_data)

fish <- fish |>
  left_join(lw_species, by = "Taxonomy") |>
  left_join(lw_family, by = "Family") |>
  mutate(
    # a and b must come from the same fit, so choose the level by availability of a
    lw_level = case_when(!is.na(lw_a_sp) ~ "species", !is.na(lw_a_fam) ~ "family", TRUE ~ "global"),
    lw_a = case_when(lw_level == "species" ~ lw_a_sp, lw_level == "family" ~ lw_a_fam, TRUE ~ lw_global$lw_a),
    lw_b = case_when(lw_level == "species" ~ lw_b_sp, lw_level == "family" ~ lw_b_fam, TRUE ~ lw_global$lw_b),
    # Biomass = -1 is the MCR missing-value code; recover mass from length
    mass_g = if_else(Biomass > 0, Biomass / Count, lw_a * Total_Length^lw_b)
  ) |>
  select(-ends_with("_sp"), -ends_with("_fam"))

# ---- 3. FishBase traits (cached, one vectorised call per table) -------------

# The queried species list is stored alongside each table so the cache is
# refreshed only when the species list changes
fetch_cached <- function(name, species, fetch_fun) {
  path  <- file.path(cache_dir, paste0(name, ".csv"))
  query <- file.path(cache_dir, paste0(name, "_query.txt"))
  if (file.exists(path) && file.exists(query) && all(species %in% readLines(query))) {
    return(read_csv(path, show_col_types = FALSE))
  }
  d <- fetch_fun(species)
  write_csv(d, path)
  writeLines(species, query)
  d
}

species_mcr <- sort(unique(fish$fb_name))

synonym_tbl <- fetch_cached("synonyms", species_mcr, function(sp) {
  rfishbase::synonyms(sp) |>
    select(query_species = synonym, Species, Status) |>
    mutate(across(everything(), as.character))
})

accepted_map <- synonym_tbl |>
  filter(!is.na(Species), Status %in% c("accepted name", "synonym")) |>
  mutate(rank = if_else(Status == "accepted name", 1L, 2L)) |>
  arrange(query_species, rank) |>
  distinct(query_species, .keep_all = TRUE) |>
  select(fb_name = query_species, fb_accepted = Species)

fish <- fish |>
  left_join(accepted_map, by = "fb_name") |>
  mutate(fb_accepted = coalesce(fb_accepted, fb_name), fb_genus = word(fb_accepted, 1))

species_fb <- sort(unique(fish$fb_accepted))

eco_tbl <- fetch_cached("ecology", species_fb, function(sp) {
  rfishbase::ecology(sp, fields = c("Species", "DietTroph", "FoodTroph")) |>
    mutate(query_species = Species)
})
morph_tbl <- fetch_cached("morphometrics", species_fb, function(sp) {
  rfishbase::morphometrics(sp) |>
    select(Species, AspectRatio) |>
    mutate(query_species = Species)
})
growth_tbl <- fetch_cached("popgrowth", species_fb, function(sp) {
  rfishbase::popgrowth(sp) |>
    select(Species, K, Loo, Type, TLinfinity) |>
    mutate(query_species = Species)
})

# Loo to total length: FishBase TLinfinity where available, otherwise a
# type-specific TL/Loo ratio estimated from records that have both
tl_ratio <- growth_tbl |>
  filter(!is.na(TLinfinity), !is.na(Loo), Loo > 0) |>
  group_by(Type) |>
  summarise(ratio = median(TLinfinity / Loo), .groups = "drop")

growth_tbl <- growth_tbl |>
  left_join(tl_ratio, by = "Type") |>
  mutate(
    Loo_TL_cm = case_when(
      !is.na(TLinfinity) ~ TLinfinity,
      Type == "TL" ~ Loo,
      !is.na(ratio) ~ Loo * ratio,
      TRUE ~ NA_real_
    )
  )

sp_genus_family <- fish |> distinct(Species = fb_accepted, Genus = fb_genus, Family)

traits_sp <- sp_genus_family |>
  left_join(
    eco_tbl |> group_by(Species) |>
      summarise(troph = mean(c(DietTroph, FoodTroph), na.rm = TRUE), .groups = "drop"),
    by = "Species"
  ) |>
  left_join(
    morph_tbl |> group_by(Species) |>
      summarise(asp = mean(AspectRatio, na.rm = TRUE), .groups = "drop"),
    by = "Species"
  ) |>
  left_join(
    growth_tbl |> filter(!is.na(K), !is.na(Loo_TL_cm)) |> group_by(Species) |>
      summarise(K = median(K), Loo_TL_cm = median(Loo_TL_cm), .groups = "drop"),
    by = "Species"
  ) |>
  mutate(across(c(troph, asp, K, Loo_TL_cm), ~ if_else(is.finite(.x), .x, NA_real_)))

# Fallback: species -> genus -> family -> global mean
fill_hierarchy <- function(d, var) {
  d |>
    group_by(Genus)  |> mutate(.g = mean(.data[[var]], na.rm = TRUE)) |>
    group_by(Family) |> mutate(.f = mean(.data[[var]], na.rm = TRUE)) |>
    ungroup() |>
    mutate(
      .level = case_when(
        !is.na(.data[[var]]) ~ "species", is.finite(.g) ~ "genus",
        is.finite(.f) ~ "family", TRUE ~ "global"
      ),
      "{var}_level" := .level,
      "{var}" := case_when(
        .level == "species" ~ .data[[var]], .level == "genus" ~ .g,
        .level == "family" ~ .f, TRUE ~ mean(.data[[var]], na.rm = TRUE)
      )
    ) |>
    select(-.g, -.f, -.level)
}

# K and Loo are filled together (same level) to keep growth pairs consistent
traits <- traits_sp |>
  fill_hierarchy("troph") |>
  fill_hierarchy("asp") |>
  mutate(KL = if_else(is.na(K) | is.na(Loo_TL_cm), NA_real_, 1)) |>
  fill_hierarchy("KL") |>
  group_by(Genus)  |> mutate(K_g = median(K, na.rm = TRUE), L_g = median(Loo_TL_cm, na.rm = TRUE)) |>
  group_by(Family) |> mutate(K_f = median(K, na.rm = TRUE), L_f = median(Loo_TL_cm, na.rm = TRUE)) |>
  ungroup() |>
  mutate(
    K_glob = median(K, na.rm = TRUE), L_glob = median(Loo_TL_cm, na.rm = TRUE),
    K         = case_when(KL_level == "species" ~ K, KL_level == "genus" ~ K_g,
                          KL_level == "family" ~ K_f, TRUE ~ K_glob),
    Loo_TL_cm = case_when(KL_level == "species" ~ Loo_TL_cm, KL_level == "genus" ~ L_g,
                          KL_level == "family" ~ L_f, TRUE ~ L_glob)
  ) |>
  select(Species, troph, troph_level, asp, asp_level, K, Loo_TL_cm, growth_level = KL_level) |>
  distinct(Species, .keep_all = TRUE)

stopifnot(!anyNA(traits$troph), !anyNA(traits$asp), !anyNA(traits$K), !anyNA(traits$Loo_TL_cm))

# ---- 4. In-situ backreef temperature (annual mean of daily means) -----------

temp_annual <- read_csv(
  here("Data", "raw_data", "MCR_LTER02_BTM_Backreef_Forereef_20251114.csv"),
  show_col_types = FALSE
) |>
  filter(reef_type_code == "Backreef") |>
  mutate(date = as.Date(time_local)) |>
  group_by(date) |>
  summarise(daily_temp = mean(temperature_c, na.rm = TRUE), .groups = "drop") |>
  mutate(Year = year(date)) |>
  group_by(Year) |>
  summarise(temp_c = mean(daily_temp, na.rm = TRUE), n_days = n(), .groups = "drop")

missing_temp_years <- setdiff(unique(fish$Year), temp_annual$Year)
if (length(missing_temp_years) > 0) {
  warning("No in-situ temperature for years ", paste(missing_temp_years, collapse = ", "),
          "; using the long-term mean.")
}

# ---- 5. Individual metabolic rates -----------------------------------------

# Vectorised copy of fishflux's getB0 (Barneche & Allen 2018, Model 2)
b0_adjust <- function(b0, troph, temp_c) {
  k <- 8.62e-5; Er <- 0.70; Ei <- 2.43; Tc <- 293.15; Topt <- 306.4
  Tk   <- temp_c + 273.15
  kt   <- 1 / (k * Tc) - 1 / (k * Tk)
  ktop <- 1 / (k * Topt) - 1 / (k * Tk)
  exp(log(b0) + 0.81 * log(troph) + Er * kt - log(1 + (Er / (Ei - Er)) * exp(Ei * ktop)))
}
stopifnot(abs(b0_adjust(0.00048083, 2, 27.5) - 0.001572901) < 1e-9)

met_pars <- fishflux::metabolic_parameters
met_all  <- met_pars |> filter(family == "all")
met_fam  <- fish |> distinct(Family) |>
  left_join(met_pars |> select(Family = family, alpha = a_m, b0 = b0_m), by = "Family") |>
  mutate(
    met_level = if_else(is.na(b0), "all-family mean", "family"),
    alpha = coalesce(alpha, met_all$a_m),
    b0    = coalesce(b0, met_all$b0_m)
  )

fish_calc <- fish |>
  left_join(traits, by = c("fb_accepted" = "Species")) |>
  left_join(met_fam, by = "Family") |>
  left_join(temp_annual |> select(Year, temp_c), by = "Year") |>
  mutate(
    temp_c = coalesce(temp_c, mean(temp_annual$temp_c)),
    Loo_mm = Loo_TL_cm * 10,                       # survey lengths are TL in mm
    m_max  = pmax(lw_a * Loo_mm^lw_b, mass_g),
    growth_g_day = pmax(0, lw_b * lw_a * Total_Length^(lw_b - 1) * K * (Loo_mm - Total_Length) / 365),
    B0 = b0_adjust(b0, troph, temp_c)
  )

mr_field <- fishflux::metabolic_rate(
  temp = fish_calc$temp_c, troph = fish_calc$troph, asp = fish_calc$asp,
  B0 = fish_calc$B0, m_max = fish_calc$m_max, m = fish_calc$mass_g, a = fish_calc$alpha,
  growth_g_day = fish_calc$growth_g_day, f = ACTIVITY_SCOPE
)

gC_to_mmolO2_h <- function(gC_day) gC_day / C_MOLAR_MASS / RQ * 1000 / 24

# Two respiration estimates are carried through to the output:
#   Rfish_resting -- resting metabolic rate (f = 1): maintenance metabolism of
#     an inactive, post-absorptive fish (Barneche & Allen 2018 B0 scaling),
#     plus the cost of somatic growth. A lower bound on fish O2 demand; the
#     closer analogue to nighttime respiration, when most diurnal reef fishes
#     shelter and are inactive.
#   Rfish_field -- field metabolic rate (f = ACTIVITY_SCOPE, default 2):
#     resting rate multiplied by the activity scope to include swimming,
#     foraging, digestion (specific dynamic action) and social activity. A
#     24-h average for free-living fish (Schiettekatte et al. 2020) and the
#     default estimate of total daily fish respiration.
# Because f is a single multiplier, Rfish_field = ACTIVITY_SCOPE * Rfish_resting
# for every record; the choice between them is about which fish activity
# state matches the respiration being compared (e.g. dark Rd vs daily R).

fish_calc <- fish_calc |>
  mutate(
    C_field_g_d   = mr_field$Total_metabolic_rate_C_g_d,
    C_rest_g_d    = C_field_g_d / ACTIVITY_SCOPE,
    growth_share  = mr_field$Cost_growth_j_g / mr_field$Resting_metabolic_rate_j_d,
    # Areal flux contributed by this record, scaled by its own swath area
    R_field_areal = gC_to_mmolO2_h(C_field_g_d) * Count / swath_area_m2,
    R_rest_areal  = gC_to_mmolO2_h(C_rest_g_d)  * Count / swath_area_m2,
    biomass_areal = mass_g * Count / swath_area_m2
  )

stopifnot(!anyNA(fish_calc$R_field_areal), all(fish_calc$R_field_areal >= 0))

# ---- 6. Aggregate: transect -> site-year ------------------------------------

fish_resp_transect <- fish_calc |>
  group_by(Site, Year, Transect) |>
  summarise(
    fish_biomass_g_m2         = sum(biomass_areal),
    Rfish_field_mmolO2_m2_h   = sum(R_field_areal),
    Rfish_resting_mmolO2_m2_h = sum(R_rest_areal),
    .groups = "drop"
  ) |>
  right_join(survey_transects, by = c("Site", "Year", "Transect")) |>
  mutate(across(c(fish_biomass_g_m2, starts_with("Rfish")), ~ replace_na(.x, 0)))

fish_respiration <- fish_resp_transect |>
  group_by(Site, Year) |>
  summarise(
    n_transects = n(),
    fish_biomass_g_m2            = mean(fish_biomass_g_m2),
    Rfish_field_se               = sd(Rfish_field_mmolO2_m2_h) / sqrt(n_transects),
    Rfish_field_mmolO2_m2_h      = mean(Rfish_field_mmolO2_m2_h),
    Rfish_resting_mmolO2_m2_h    = mean(Rfish_resting_mmolO2_m2_h),
    .groups = "drop"
  ) |>
  left_join(temp_annual |> select(Year, temp_c), by = "Year") |>
  relocate(Rfish_field_se, .after = Rfish_field_mmolO2_m2_h) |>
  arrange(Site, Year)

write_csv(fish_respiration, here("Data", "fish_respiration.csv"))

# ---- 7. Diagnostics ---------------------------------------------------------

trait_coverage <- fish_calc |>
  summarise(
    across(
      c(troph_level, asp_level, growth_level, lw_level, met_level),
      ~ 100 * sum(biomass_areal[.x %in% c("species", "family")]) / sum(biomass_areal)
    )
  )
message("Percent of biomass with species-level traits (family-level for LW / metabolic params):")
print(trait_coverage)

message(sprintf(
  "Growth cost as share of resting metabolism: median %.1f%%",
  100 * median(fish_calc$growth_share)
))

message("Fish respiration at LTER_1 (mmol O2 m-2 h-1):")
print(fish_respiration |> filter(Site == "LTER_1"), n = Inf)
