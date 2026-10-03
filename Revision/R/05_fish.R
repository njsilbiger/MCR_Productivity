# ---------------------------------------------------------------------------
# 05_fish.R
#
# Section 7 of the Revision plan (Revision/Reviewer_Response_Plan.md):
# fish community respiration (Task 7.1), mechanistic and independent of the
# in-situ flux measurement, to (a) answer "do fish contribute to ecosystem
# respiration, and how much" directly and quantitatively, and (b) build
# `fish_resp_z` to unblock the Rd model covariate Task 6.1 specified.
#
# METHOD (REVISED 2026-10-02): the PI installed `fishflux` (Schiettekatte et
# al. 2020) and `rfishbase`, replacing this script's earlier generic
# Clarke & Johnston (1999) approximation with genuine SPECIES- (falling
# back to genus- then family- then global-) SPECIFIC parameters:
#
#   - Resting metabolic rate follows Barneche & Allen (2018, Ecol. Lett.
#     doi:10.1111/ele.12947) "Model 2": family- and trophic-level-specific
#     B0/alpha via `fishflux::metabolism()`, combined with trophic level,
#     caudal-fin aspect ratio (activity proxy), max body size, and daily
#     growth rate via `fishflux::metabolic_rate()`.
#   - Trophic level (`rfishbase::ecology()`), aspect ratio
#     (`rfishbase::morphometrics()`), and von Bertalanffy growth parameters
#     K/Loo (`rfishbase::popgrowth()`) are pulled ONCE in bulk for all 238
#     species in `fish_ind.csv` (NOT per-species in a loop -- the
#     `fishflux` wrapper functions do that and are extremely slow; calling
#     the underlying vectorised `rfishbase` functions directly is not).
#   - Species -> genus -> family -> global-mean fallback is used wherever
#     FishBase has no record for a given species (fishflux's own wrappers
#     only fall back one level, species -> genus; extended here to family
#     and a global mean as a final backstop so every individual gets a
#     value).
#   - Max body size (`m_max`) and daily growth (`growth_g_day`) use a
#     length-weight relationship fit from THIS project's own survey data
#     (log Biomass ~ log Total_Length per family, falling back to a global
#     fit), combined with FishBase's Loo/K via the standard von Bertalanffy
#     growth-in-weight identity, rather than FishBase's own (very sparse,
#     Weight field is NA for most reef species) length-weight data.
#   - Activity scope `f = 2` (a literature-typical default, matching
#     `fishflux`'s own documentation example) -- NOT species-specific;
#     flagged as a remaining approximation.
#
# TWO REAL BUGS CAUGHT AND FIXED DURING DEVELOPMENT (documented so they are
# not silently reintroduced):
#   1. `Total_Length` in `fish_ind.csv` is in MILLIMETRES (confirmed via
#      Carcharhinus melanopterus, 800-1000 = 80-100cm, a sensible reef
#      shark size), but FishBase's `Loo` is in CENTIMETRES. Combining them
#      directly corrupted `m_max`/`growth_g_day` silently (wrong magnitude,
#      not always a visible error). Fixed by converting Loo to mm
#      (`Loo_mm <- Loo_final * 10`) before use.
#   2. The global length-weight fallback coefficients were extracted from
#      a named vector via `lw_fit_global["a"]`/`["b"]`, but R's `coef()`
#      preserves the original term names (e.g. "a.(Intercept)"), so the
#      lookup silently returned NA for every individual whose family fell
#      back to the global fit (26 of 28,964 records) -- caught by checking
#      for unexpected NAs in the output rather than assuming zero NAs meant
#      success. Fixed with `unname()`.
#
# OUTPUTS:
#   Revision/Data/derived/fish_respiration_transect.csv
#   Revision/Data/derived/fish_respiration_site_year.csv
#   Revision/Data/derived/fish_vs_measured_Rd.csv
#   Revision/Output/fig_fish_pct_of_Rd.png
#   Revision/Output/qc_fish.md
#   Revision/Data/derived/_cache/fishbase_{ecology,morphometrics,popgrowth,species}.csv
#   Revision/Data/derived/_cache/fishflux_family_metabolism.csv
# ---------------------------------------------------------------------------

source(here::here("Revision", "R", "00_packages.R"))
library(tidyverse)
if (!requireNamespace("fishflux", quietly = TRUE)) stop("fishflux not installed -- see plan Section 7 for install instructions")
if (!requireNamespace("rfishbase", quietly = TRUE)) stop("rfishbase not installed")
library(fishflux)
library(rfishbase)

qc_lines <- character(0)
qc <- function(...) qc_lines <<- c(qc_lines, paste0(...))
qc("# QC summary: fish community respiration (Revision Step 5 / Section 7, Task 7.1)")
qc("")
qc("Generated: ", as.character(Sys.time()))
qc("")
qc("**Method revised 2026-10-02**: replaced the earlier generic Clarke & ",
   "Johnston (1999) allometric approximation with `fishflux`/`rfishbase`-",
   "based species-specific (species > genus > family > global fallback) ",
   "respiration estimates, following Barneche & Allen (2018). See the ",
   "script header for the two real bugs (a length-unit mismatch; a silent ",
   "NA from R name-mangling) caught and fixed during development.")
qc("")

cache_dir <- here("Revision", "Data", "derived", "_cache")

# =============================================================================
# A. Annual mean temperature (island-wide satellite SST, all sites/years)
# =============================================================================

sst_monthly <- read_csv(file.path(cache_dir, "moorea_sst_monthly.csv"), show_col_types = FALSE)
sst_annual <- sst_monthly |> mutate(Year = year(date)) |> group_by(Year) |> summarise(T_mean_C = mean(sst, na.rm = TRUE), .groups = "drop")
T_OVERALL <- mean(sst_annual$T_mean_C)

qc("## A. Annual mean temperature")
qc("")
qc("- Island-wide annual mean SST from the cached satellite monthly series ",
   "(Step 2), ", min(sst_annual$Year), "-", max(sst_annual$Year),
   ". Overall mean: ", round(T_OVERALL, 2), " degC.")
qc("")

# =============================================================================
# B. Bulk FishBase trait lookups (ONE call per table, not per species)
# =============================================================================

fish_ind <- read_csv(here("Revision", "Data", "derived", "fish_ind.csv"), show_col_types = FALSE) |>
  mutate(Genus = word(Taxonomy, 1))
species_list <- sort(unique(fish_ind$Taxonomy))

fetch_or_load <- function(file, fetch_fun) {
  path <- file.path(cache_dir, file)
  if (file.exists(path)) return(read_csv(path, show_col_types = FALSE))
  d <- fetch_fun()
  write_csv(d, path)
  d
}

eco_bulk   <- fetch_or_load("fishbase_ecology.csv",      function() rfishbase::ecology(species_list))
morph_bulk <- fetch_or_load("fishbase_morphometrics.csv", function() rfishbase::morphometrics(species_list))
growth_bulk <- fetch_or_load("fishbase_popgrowth.csv",    function() rfishbase::popgrowth(species_list))

qc("## B. Bulk FishBase trait lookups")
qc("")
qc("- `ecology()` (trophic level): ", n_distinct(eco_bulk$Species[!is.na(eco_bulk$DietTroph) | !is.na(eco_bulk$FoodTroph)]),
   " of ", length(species_list), " species matched.")
qc("- `morphometrics()` (aspect ratio): ", n_distinct(morph_bulk$Species[!is.na(morph_bulk$AspectRatio)]),
   " of ", length(species_list), " species matched.")
qc("- `popgrowth()` (K, Loo): ", n_distinct(growth_bulk$Species[!is.na(growth_bulk$K) & !is.na(growth_bulk$Loo)]),
   " of ", length(species_list), " species matched (otolith-based growth studies are sparse for reef fish).")
qc("- All three fetched ONCE for the full species list (vectorised `rfishbase` calls), not in a per-species loop ",
   "(the `fishflux` wrapper functions -- `trophic_level()`, `aspect_ratio()`, `growth_params()` -- call these same ",
   "functions one species at a time, which is impractically slow at 238 species).")
qc("")

# =============================================================================
# C. Species -> genus -> family -> global trait fallback hierarchy
# =============================================================================

species_genus_map <- tibble(Species = species_list, Genus = word(species_list, 1))
sp_family_map <- fish_ind |> distinct(Taxonomy, Family) |> rename(Species = Taxonomy)

troph_sp <- eco_bulk |> mutate(troph = rowMeans(cbind(DietTroph, FoodTroph), na.rm = TRUE)) |>
  group_by(Species) |> summarise(troph = mean(troph, na.rm = TRUE), .groups = "drop") |> filter(is.finite(troph))
asp_sp <- morph_bulk |> group_by(Species) |> summarise(asp = mean(AspectRatio, na.rm = TRUE), .groups = "drop") |> filter(is.finite(asp))
growth_sp <- growth_bulk |> group_by(Species) |> summarise(K = mean(K, na.rm = TRUE), Loo = mean(Loo, na.rm = TRUE), .groups = "drop") |>
  filter(is.finite(K), is.finite(Loo))

troph_genus  <- troph_sp |> left_join(species_genus_map, by = "Species") |> group_by(Genus) |> summarise(troph_genus = mean(troph, na.rm = TRUE), .groups = "drop")
asp_genus    <- asp_sp   |> left_join(species_genus_map, by = "Species") |> group_by(Genus) |> summarise(asp_genus = mean(asp, na.rm = TRUE), .groups = "drop")
growth_genus <- growth_sp |> left_join(species_genus_map, by = "Species") |> group_by(Genus) |>
  summarise(K_genus = mean(K, na.rm = TRUE), Loo_genus = mean(Loo, na.rm = TRUE), .groups = "drop")

troph_family  <- troph_sp |> left_join(sp_family_map, by = "Species") |> group_by(Family) |> summarise(troph_family = mean(troph, na.rm = TRUE), .groups = "drop")
asp_family    <- asp_sp   |> left_join(sp_family_map, by = "Species") |> group_by(Family) |> summarise(asp_family = mean(asp, na.rm = TRUE), .groups = "drop")
growth_family <- growth_sp |> left_join(sp_family_map, by = "Species") |> group_by(Family) |>
  summarise(K_family = mean(K, na.rm = TRUE), Loo_family = mean(Loo, na.rm = TRUE), .groups = "drop")

troph_global <- mean(troph_sp$troph, na.rm = TRUE)
asp_global   <- mean(asp_sp$asp, na.rm = TRUE)
K_global     <- mean(growth_sp$K, na.rm = TRUE)
Loo_global   <- mean(growth_sp$Loo, na.rm = TRUE)

species_traits <- species_genus_map |>
  left_join(sp_family_map |> distinct(Species, Family), by = "Species") |>
  left_join(troph_sp, by = "Species") |> left_join(troph_genus, by = "Genus") |> left_join(troph_family, by = "Family") |>
  left_join(asp_sp, by = "Species") |> left_join(asp_genus, by = "Genus") |> left_join(asp_family, by = "Family") |>
  left_join(growth_sp, by = "Species") |> left_join(growth_genus, by = "Genus") |> left_join(growth_family, by = "Family") |>
  mutate(
    troph_level = coalesce(troph, troph_genus, troph_family, troph_global),
    asp_final   = coalesce(asp, asp_genus, asp_family, asp_global),
    K_final     = coalesce(K, K_genus, K_family, K_global),
    Loo_final   = coalesce(Loo, Loo_genus, Loo_family, Loo_global)
  ) |>
  dplyr::select(Species, Family, troph_level, asp_final, K_final, Loo_final)

qc("## C. Trait fallback hierarchy (species -> genus -> family -> global)")
qc("")
qc("- Global fallback values: trophic level = ", round(troph_global, 2),
   ", aspect ratio = ", round(asp_global, 2), ", K = ", round(K_global, 3),
   ", Loo = ", round(Loo_global, 1), " cm.")
qc("- Every one of the ", nrow(species_traits), " species ends up with a complete trait set after fallback.")
qc("")

# =============================================================================
# D. Own length-weight fit (per family, for m_max and growth_g_day) and
#    family-level B0/alpha via fishflux::metabolism()
# =============================================================================

fish_ind_lw <- fish_ind |> filter(Total_Length > 0, Biomass > 0, Count > 0) |> mutate(mass_g = Biomass / Count)

lw_fit <- fish_ind_lw |>
  group_by(Family) |>
  filter(n() >= 10) |>
  summarise(
    b_exp  = unname(coef(lm(log(mass_g) ~ log(Total_Length)))[2]),
    a_coef = unname(exp(coef(lm(log(mass_g) ~ log(Total_Length)))[1])),
    .groups = "drop"
  )
lw_global_lm <- lm(log(mass_g) ~ log(Total_Length), data = fish_ind_lw)
lw_fit_global <- c(b = unname(coef(lw_global_lm)[2]), a = unname(exp(coef(lw_global_lm)[1])))

fam_metab_cache <- file.path(cache_dir, "fishflux_family_metabolism.csv")
if (file.exists(fam_metab_cache)) {
  fam_metab <- read_csv(fam_metab_cache, show_col_types = FALSE)
} else {
  families_all <- unique(fish_ind$Family)
  fam_metab <- map_dfr(families_all, function(f) {
    troph_f <- mean(species_traits$troph_level[species_traits$Family == f], na.rm = TRUE)
    r <- fishflux::metabolism(family = f, temp = T_OVERALL, troph_m = troph_f)
    r$Family <- f
    r
  })
  write_csv(fam_metab, fam_metab_cache)
}

qc("## D. Length-weight fit and family-level metabolic parameters")
qc("")
qc("- Own length-weight fit (log Biomass ~ log Total_Length): ", nrow(lw_fit),
   " of ", n_distinct(fish_ind$Family), " families have >= 10 records for a dedicated fit; ",
   "the rest use the global fit (a = ", signif(lw_fit_global["a"], 3), ", b = ", round(lw_fit_global["b"], 3), ").")
qc("- `fishflux::metabolism()` (Barneche & Allen 2018 Model 2) called once per family at the overall mean ",
   "temperature (", round(T_OVERALL, 2), " degC) -- most families fall back to the package's global-average ",
   "B0/alpha (only a handful of well-studied families, e.g. Pomacentridae, Gobiidae, Apogonidae, have their own ",
   "fitted family-level parameters in Barneche & Allen's dataset; this is `fishflux`'s own documented behaviour, ",
   "not a shortcut taken here).")
qc("")

# =============================================================================
# E. Per-individual respiration (fishflux::metabolic_rate)
# =============================================================================

RQ <- 0.8  # respiratory quotient, literature-typical for a mixed fish diet -- flagged approximation
ACTIVITY_SCOPE <- 2  # f in metabolic_rate(); literature-typical default (fishflux's own example), not species-specific

fish_calc <- fish_ind |>
  filter(Count > 0, Biomass > 0, Total_Length > 0) |>
  left_join(species_traits |> dplyr::select(Species, troph_level, asp_final, K_final, Loo_final), by = c("Taxonomy" = "Species")) |>
  left_join(lw_fit, by = "Family") |>
  left_join(fam_metab |> dplyr::select(Family, B0 = b0_m, alpha = alpha_m), by = "Family") |>
  left_join(sst_annual, by = "Year") |>
  mutate(
    a_coef = coalesce(a_coef, lw_fit_global["a"]),
    b_exp  = coalesce(b_exp,  lw_fit_global["b"]),
    mass_g = Biomass / Count,
    Loo_mm = Loo_final * 10,  # Loo (FishBase) is in cm; Total_Length is in mm -- see header bug #1
    m_max  = pmax(a_coef * Loo_mm^b_exp, mass_g),
    growth_g_day = pmax(0, b_exp * a_coef * Total_Length^(b_exp - 1) * K_final * (Loo_mm - Total_Length) / 365),
    f_activity = ACTIVITY_SCOPE
  )

mr <- fishflux::metabolic_rate(
  temp = fish_calc$T_mean_C, troph = fish_calc$troph_level, asp = fish_calc$asp_final,
  B0 = fish_calc$B0, m_max = fish_calc$m_max, m = fish_calc$mass_g, a = fish_calc$alpha,
  growth_g_day = fish_calc$growth_g_day, f = fish_calc$f_activity
)
fish_calc$Cm_gC_day <- mr$Total_metabolic_rate_C_g_d

n_na_cm <- sum(is.na(fish_calc$Cm_gC_day))
stopifnot(n_na_cm == 0)  # both known bugs are fixed above; any new NA here needs investigation, not silent dropping

fish_resp_ind <- fish_calc |>
  mutate(
    mol_O2_day = (Cm_gC_day / 12.011) / RQ,     # gC/day -> mol CO2/day -> mol O2/day via RQ
    Rb_per_ind = mol_O2_day * 1000 / 24,         # -> mmol O2/h per individual
    Rb_total   = Rb_per_ind * Count
  )

qc("## E. Per-individual respiration (species-specific `fishflux` model)")
qc("")
qc("- Respiratory quotient (O2 from metabolic carbon loss): RQ = ", RQ,
   " (literature-typical for mixed fish diet -- not species-specific; flagged approximation).")
qc("- Activity scope f = ", ACTIVITY_SCOPE, " (resting-to-routine default from `fishflux`'s own documentation; not species-specific).")
qc("- NA count in computed metabolic rate: ", n_na_cm, " of ", nrow(fish_calc), " (both bugs described in the script header are fixed; this should stay 0).")
qc("")

cat("Rb_per_ind summary (mmol O2/h per individual), species-specific model:\n")
print(summary(fish_resp_ind$Rb_per_ind))

# =============================================================================
# F. Aggregate to transect/site/year; convert to areal flux
# =============================================================================

AREA_PER_TRANSECT_M2 <- 300  # from Step 1 (50 m x (1+5) m swaths)
ON_RATIO <- 20                # placeholder atomic O:N (see note below), not trophic-group-specific

fish_resp_transect <- fish_resp_ind |>
  group_by(Year, Site, Transect) |>
  summarise(Rb_total_mmol_h = sum(Rb_total, na.rm = TRUE), .groups = "drop") |>
  mutate(
    Rb_areal_mmolO2_m2_h = Rb_total_mmol_h / AREA_PER_TRANSECT_M2,
    N_excretion_mmolN_m2_h = Rb_areal_mmolO2_m2_h / ON_RATIO
  )
write_csv(fish_resp_transect, here("Revision", "Data", "derived", "fish_respiration_transect.csv"))

fish_resp_site_year <- fish_resp_transect |>
  group_by(Year, Site) |>
  summarise(
    Rb_areal_mmolO2_m2_h = mean(Rb_areal_mmolO2_m2_h),
    N_excretion_mmolN_m2_h = mean(N_excretion_mmolN_m2_h),
    n_transects = n(),
    .groups = "drop"
  )
write_csv(fish_resp_site_year, here("Revision", "Data", "derived", "fish_respiration_site_year.csv"))

qc("## F. Areal respiration and N excretion (transect/site/year)")
qc("")
qc("- `fish_respiration_transect.csv`: ", nrow(fish_resp_transect), " transect-years. ",
   "`fish_respiration_site_year.csv`: ", nrow(fish_resp_site_year), " site-years.")
qc("- **N excretion remains a simple stoichiometric placeholder (O:N atomic = ", ON_RATIO,
   "), NOT trophic-group-specific -- order-of-magnitude only, pending Task 7.3.**")
qc("")

cat("Areal respiration (mmol O2/m2/h) summary, all site-years, species-specific model:\n")
print(summary(fish_resp_site_year$Rb_areal_mmolO2_m2_h))

# =============================================================================
# G. LTER_1 series: compare to measured Rd (the headline Task 7.1 result)
# =============================================================================

fish_resp_lter1 <- fish_resp_site_year |> filter(Site == "LTER_1") |> arrange(Year)

pp_day <- read_csv(here("Revision", "Data", "derived", "pp_day.csv"), show_col_types = FALSE)
rd_annual_lter1 <- pp_day |>
  mutate(Rd_measured = -R_mean) |>
  group_by(Year) |>
  summarise(Rd_measured_mmolO2_m2_h = mean(Rd_measured, na.rm = TRUE), n_days = n(), .groups = "drop")

fish_vs_rd <- fish_resp_lter1 |>
  inner_join(rd_annual_lter1, by = "Year") |>
  mutate(pct_of_Rd = 100 * Rb_areal_mmolO2_m2_h / Rd_measured_mmolO2_m2_h)
write_csv(fish_vs_rd, here("Revision", "Data", "derived", "fish_vs_measured_Rd.csv"))

qc("## G. Fish respiration vs. measured Rd at LTER_1 -- the headline result")
qc("")
qc("- Years with both a fish survey and metabolism data at LTER_1: ", nrow(fish_vs_rd),
   " (", min(fish_vs_rd$Year), "-", max(fish_vs_rd$Year), "). Units assumed mmol O2 m^-2 h^-1 for both series.")
qc("")
qc("| Year | Fish Rb (mmol O2/m2/h) | Measured Rd (mmol O2/m2/h) | % of Rd from fish |")
qc("|---|---|---|---|")
for (i in seq_len(nrow(fish_vs_rd))) {
  r <- fish_vs_rd[i, ]
  qc(sprintf("| %d | %.3f | %.2f | %.2f%% |", r$Year, r$Rb_areal_mmolO2_m2_h, r$Rd_measured_mmolO2_m2_h, r$pct_of_Rd))
}
qc("")
qc(sprintf("**Summary: species-specific fish respiration accounts for a mean of %.2f%% (range %.2f%%-%.2f%%) of measured ecosystem respiration at LTER_1 across the %d years both are available.**",
           mean(fish_vs_rd$pct_of_Rd), min(fish_vs_rd$pct_of_Rd), max(fish_vs_rd$pct_of_Rd), nrow(fish_vs_rd)))
qc("")

cat("\nFish respiration as % of measured Rd at LTER_1 (species-specific model):\n")
print(fish_vs_rd |> dplyr::select(Year, Rb_areal_mmolO2_m2_h, Rd_measured_mmolO2_m2_h, pct_of_Rd))

p_fish_pct <- ggplot(fish_vs_rd, aes(x = Year, y = pct_of_Rd)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Year", y = "Fish respiration, % of measured Rd",
       title = "Species-specific fish respiration (fishflux) as a share of measured ecosystem respiration (LTER_1)") +
  theme_minimal()
ggsave(here("Revision", "Output", "fig_fish_pct_of_Rd.png"), p_fish_pct, width = 8, height = 5, dpi = 300)

writeLines(qc_lines, here("Revision", "Output", "qc_fish.md"))
cat("\nDone (through Section G). Wrote:\n",
    "  Revision/Data/derived/fish_respiration_transect.csv\n",
    "  Revision/Data/derived/fish_respiration_site_year.csv\n",
    "  Revision/Data/derived/fish_vs_measured_Rd.csv\n",
    "  Revision/Output/fig_fish_pct_of_Rd.png\n",
    "  Revision/Output/qc_fish.md\n")

# =============================================================================
# H. fish_resp_z (species-specific) and refitting E4/E6
# =============================================================================

if (!requireNamespace("brms", quietly = TRUE)) install.packages("brms")
library(brms)

fish_resp_mean <- mean(fish_resp_lter1$Rb_areal_mmolO2_m2_h)
fish_resp_sd   <- sd(fish_resp_lter1$Rb_areal_mmolO2_m2_h)
fish_resp_z_lter1 <- fish_resp_lter1 |>
  mutate(fish_resp_z = (Rb_areal_mmolO2_m2_h - fish_resp_mean) / fish_resp_sd) |>
  dplyr::select(Year, fish_resp_z)

pi_dev_data <- read_csv(here("Revision", "Data", "derived", "pi_model_hourly.csv"), show_col_types = FALSE)

fish_transect_05 <- read_csv(here("Revision", "Data", "derived", "fish_transect.csv"), show_col_types = FALSE)
fish_lter1_year_05 <- fish_transect_05 |>
  filter(Site == "LTER_1") |>
  mutate(dag_group = case_when(trophic_group == "Herbivore" ~ "Herb", trophic_group == "Corallivore" ~ "Corall", TRUE ~ "OtherFish")) |>
  group_by(Year, Transect, dag_group) |> summarise(biomass_g_m2 = sum(biomass_g_m2), .groups = "drop") |>
  group_by(Year, dag_group) |> summarise(biomass_g_m2 = mean(biomass_g_m2), .groups = "drop") |>
  pivot_wider(names_from = dag_group, values_from = biomass_g_m2, values_fill = 0) |>
  mutate(Herb_z = as.numeric(scale(Herb)), Corall_z = as.numeric(scale(Corall)), OtherFish_z = as.numeric(scale(OtherFish))) |>
  dplyr::select(Year, Herb_z, Corall_z, OtherFish_z)
dhw_lter1_year_05 <- read_csv(here("Revision", "Data", "derived", "disturbance_site_year.csv"), show_col_types = FALSE) |>
  filter(Site == "LTER_1") |> mutate(DHW_z = as.numeric(scale(DHW_max))) |> dplyr::select(Year, DHW_z)

pi_dev_data <- pi_dev_data |> left_join(fish_lter1_year_05, by = "Year") |> left_join(dhw_lter1_year_05, by = "Year")

set.seed(4817)
days_by_year <- pi_dev_data |> distinct(Year, DielDate)
dev_days <- days_by_year |> group_by(Year) |> slice_sample(prop = 0.25) |> ungroup()
missing_years <- setdiff(unique(days_by_year$Year), unique(dev_days$Year))
if (length(missing_years) > 0) {
  dev_days <- bind_rows(dev_days, days_by_year |> filter(Year %in% missing_years) |> group_by(Year) |> slice_sample(n = 1) |> ungroup())
}
pi_dev_data <- pi_dev_data |> semi_join(dev_days, by = c("Year", "DielDate")) |> left_join(fish_resp_z_lter1, by = "Year")
n_missing_fish_resp_z <- sum(is.na(pi_dev_data$fish_resp_z))
pi_dev_data <- pi_dev_data |> filter(!is.na(fish_resp_z))

qc("## H. fish_resp_z (species-specific) and refitting E4/E6")
qc("")
qc("- `fish_resp_z`: LTER_1's annual species-specific areal respiration, z-standardised (mean ",
   round(fish_resp_mean, 5), ", sd ", round(fish_resp_sd, 5), " mmol O2/m^2/h).")
qc("- Dev-subsample rows dropped for years without a fish survey match: ", n_missing_fish_resp_z, ".")
qc("")

priors_la2 <- c(prior(normal(0, 1), nlpar = "la", class = "b"), prior(exponential(2), nlpar = "la", class = "sd"))
lP_full <- "lP ~ 1 + ilr1 + ilr2 + ilr3 + invkT_c + log_flow_c + N_z + Season + (1 | p | Year) + (1 | q | Year:DielDate)"

fit_refit <- function(name, lP_rhs, lR_rhs, extra_priors) {
  f <- bf(as.formula("PP ~ (exp(la) * exp(lP) * PAR) / (exp(la) * PAR + exp(lP)) - exp(lR)"),
          as.formula("la ~ 1 + (1 | Year)"), as.formula(lP_rhs), as.formula(lR_rhs), nl = TRUE)
  p <- c(priors_la2, extra_priors)
  brm(f, data = pi_dev_data, family = student(), prior = p,
      chains = 4, iter = 2000, warmup = 1000, seed = 4817,
      control = list(adapt_delta = 0.95),
      file = here("Revision", "Output", "fits", paste0("pi_", name, "_dev")),
      file_refit = "on_change")
}

fit_E6_fishresp <- fit_refit(
  "E6_fishresp_v2",
  lP_full,
  "lR ~ 1 + fish_resp_z + ilr1 + ilr2 + ilr3 + DHW_z + (1 | p | Year) + (1 | q | Year:DielDate)",
  c(prior(normal(log(80), 0.5), nlpar = "lR", coef = "Intercept"), prior(normal(0, 0.5), nlpar = "lR", class = "b"),
    prior(exponential(2), nlpar = "lR", class = "sd"),
    prior(normal(log(100), 0.5), nlpar = "lP", coef = "Intercept"), prior(normal(0, 0.5), nlpar = "lP", class = "b"),
    prior(exponential(2), nlpar = "lP", class = "sd"))
)

fit_E4_fishresp <- fit_refit(
  "E4_fishresp_v2",
  lP_full,
  "lR ~ 1 + ilr1 + ilr2 + ilr3 + Herb_z + Corall_z + OtherFish_z + fish_resp_z + invkT_c + log_flow_c + (1 | p | Year) + (1 | q | Year:DielDate)",
  c(prior(normal(log(80), 0.5), nlpar = "lR", coef = "Intercept"), prior(normal(0, 0.5), nlpar = "lR", class = "b"),
    prior(normal(-0.65, 0.2), nlpar = "lR", coef = "invkT_c"), prior(exponential(2), nlpar = "lR", class = "sd"),
    prior(normal(log(100), 0.5), nlpar = "lP", coef = "Intercept"), prior(normal(0, 0.5), nlpar = "lP", class = "b"),
    prior(exponential(2), nlpar = "lP", class = "sd"))
)

cat("\n=== E6 with species-specific fish_resp_z ===\n"); print(fixef(fit_E6_fishresp))
cat("\n=== E4 with species-specific fish_resp_z added ===\n"); print(fixef(fit_E4_fishresp))

write_result_table <- function(fit, name) {
  fx <- fixef(fit)
  qc(paste0("### ", name)); qc("")
  qc("| Coefficient | Estimate | 95% CI |"); qc("|---|---|---|")
  for (i in seq_len(nrow(fx))) qc(sprintf("| %s | %.3f | [%.3f, %.3f] |", rownames(fx)[i], fx[i,"Estimate"], fx[i,"Q2.5"], fx[i,"Q97.5"]))
  qc("")
}
write_result_table(fit_E6_fishresp, "E6: species-specific fish_resp_z vs. biomass trio")
write_result_table(fit_E4_fishresp, "E4: species-specific fish_resp_z added to biomass adjustment")

writeLines(qc_lines, here("Revision", "Output", "qc_fish.md"))
cat("\nDone (Section H). qc_fish.md updated.\n")
