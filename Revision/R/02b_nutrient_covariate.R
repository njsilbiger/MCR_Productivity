# ---------------------------------------------------------------------------
# 02b_nutrient_covariate.R
#
# Supplementary to Step 2 of the Revision plan -- not in the plan's original
# file list (Section 1.1), because the plan did not anticipate needing a
# nitrogen (N) covariate until Task 4.3 (testing the DAG against data) and
# Section 7 (estimands E9/E10 "Benthos/DayTemp -> Pmax", E12 "Fish -> N")
# came into view. Added here, numbered 02b, because it is a covariate-
# building step like 02_disturbance_covariates.R, and both are needed before
# Task 4.3 can run.
#
# Builds the N covariate from macroalgal tissue nitrogen (Comment 4: keep
# raw replication, don't pre-average, where a later step might want it) --
# decided with the PI 2026-10-01 to build BOTH a sample-level table and a
# site-year summary, mirroring benthic_quad.csv/benthic_site.csv.
#
# DECISION (documented in the Revision plan testing log, Step "2b"):
#   - Source: macroalgal tissue %N (Turbinaria ornata, Backreef), NOT
#     Data/raw_data/WaterColumnN.csv. The water-column dissolved N+N file is
#     single-site (LTER_1 only) and stops in 2018 -- it misses the 2019
#     bleaching event, the 2020 DHW peak, and the entire 2023-2024 COTS
#     outbreak documented in Step 2. The original qmd already treated
#     Turbinaria %N as the primary nutrient proxy and used water-column N+N
#     only to validate it via a linear regression (reproduced below as a
#     QC/documentation step, not as the covariate itself).
#   - Genus: Turbinaria ONLY, not pooled with Sargassum. Mean tissue N% is
#     ~1.14% for Sargassum vs ~0.68% for Turbinaria (a ~2 SD genus effect,
#     checked explicitly below) and Sargassum sampling stopped after 2014
#     (protocol shift to Turbinaria-only) -- pooling genera would conflate a
#     species effect with a time/site effect.
#   - Years: restricted to 2007-2024 for the site-year summary. 2005-2006
#     have no Turbinaria samples at all (Sargassum-only or nothing that
#     early); 2025 is nearly empty (5 samples, 1 site) -- consistent with the
#     lab-processing lag already seen elsewhere in this project for the
#     current year. The sample-level table keeps 2005-2025 as recorded;
#     missingness is documented, not imputed.
#   - Resolution: all 6 backreef sites (upgrading the original LTER_1-only
#     analysis), matching Steps 1-2.
#
# OUTPUTS:
#   Revision/Data/derived/N_sample.csv     (one row per Turbinaria tissue sample)
#   Revision/Data/derived/N_site_year.csv  (Site x Year mean/sd/n, 2007-2024)
#   Revision/Output/qc_nutrient.md
# ---------------------------------------------------------------------------

source(here::here("Revision", "R", "00_packages.R"))

qc_lines <- character(0)
qc <- function(...) qc_lines <<- c(qc_lines, paste0(...))
qc("# QC summary: nutrient (N) covariate (Revision Step 2b)")
qc("")
qc("Generated: ", as.character(Sys.time()))
qc("")

# =============================================================================
# A. Load and characterise both candidate sources
# =============================================================================

chn_raw <- read_csv(
  here("Data", "raw_data", "MCR_LTER_Macroalgal_CHN_2005_to_2024_20250616.csv"),
  show_col_types = FALSE
)

water_raw <- read_csv(here("Data", "raw_data", "WaterColumnN.csv"), show_col_types = FALSE) |>
  mutate(Date = mdy(Date), Year = year(Date))

qc("## A. Two candidate sources")
qc("")
qc("| | WaterColumnN.csv (dissolved N+N) | Macroalgal CHN file (tissue %N) |")
qc("|---|---|---|")
qc(sprintf("| Sites | %s | %s |",
           paste(unique(water_raw$Location), collapse = "; "),
           paste(sort(unique(chn_raw$Site)), collapse = ", ")))
qc(sprintf("| Years | %d-%d | %d-%d |",
           min(water_raw$Year), max(water_raw$Year),
           min(chn_raw$Year), max(chn_raw$Year)))
qc(sprintf("| Rows | %d | %d |", nrow(water_raw), nrow(chn_raw)))
qc("")
qc("`WaterColumnN.csv` is single-site (LTER_1) and stops in 2018 -- it ",
   "misses the 2019 bleaching event, the 2020 DHW peak, and the entire ",
   "2023-2024 COTS outbreak documented in Step 2's QC. **Decision: use ",
   "macroalgal tissue %N as the covariate; keep water-column N+N only as a ",
   "documentation/validation cross-check (Section B), reproducing the ",
   "original qmd's approach.**")
qc("")

# =============================================================================
# B. Validation cross-check: tissue %N vs dissolved N+N (reproduces the
#    original qmd's justification for using tissue %N as the proxy)
# =============================================================================

turb_lter1 <- chn_raw |>
  filter(Habitat == "Backreef", Genus == "Turbinaria", Site == "LTER_1", !is.na(N)) |>
  group_by(Year) |>
  summarise(N_percent = mean(N, na.rm = TRUE), .groups = "drop")

water_lter1 <- water_raw |>
  group_by(Year) |>
  summarise(Nitrite_and_Nitrate = mean(Nitrite_and_Nitrate, na.rm = TRUE), .groups = "drop")

nutrient_validation <- turb_lter1 |> inner_join(water_lter1, by = "Year")
mod_n_validation <- lm(Nitrite_and_Nitrate ~ N_percent, data = nutrient_validation)
cor_validation <- cor.test(nutrient_validation$Nitrite_and_Nitrate, nutrient_validation$N_percent)

qc("## B. Validation: does tissue %N track dissolved N+N at LTER_1?")
qc("")
qc("Reproduces the original qmd's check (`modN <- lm(Nitrite_and_Nitrate ~ N_percent)`), ",
   "restricted to the years both series overlap (", min(nutrient_validation$Year),
   "-", max(nutrient_validation$Year), ", n = ", nrow(nutrient_validation), " years).")
qc("")
qc(sprintf("- Pearson r = %.3f (95%% CI %.3f to %.3f), p = %.4f",
           cor_validation$estimate, cor_validation$conf.int[1], cor_validation$conf.int[2],
           cor_validation$p.value))
qc(sprintf("- lm(Nitrite_and_Nitrate ~ N_percent): slope = %.3f, R^2 = %.3f",
           coef(mod_n_validation)[2], summary(mod_n_validation)$r.squared))
qc("")

cat("Validation (LTER_1, tissue %N vs dissolved N+N):\n")
print(cor_validation)

# =============================================================================
# C. Genus check: is it safe to pool Sargassum and Turbinaria?
# =============================================================================

genus_compare <- chn_raw |>
  filter(Habitat == "Backreef", !is.na(N)) |>
  group_by(Genus) |>
  summarise(mean_N = mean(N), sd_N = sd(N), n = n(),
            year_range = paste(range(Year), collapse = "-"), .groups = "drop")

genus_by_year <- chn_raw |>
  filter(Habitat == "Backreef", !is.na(N)) |>
  count(Year, Genus) |>
  pivot_wider(names_from = Genus, values_from = n, values_fill = 0)

qc("## C. Genus check: Sargassum vs Turbinaria tissue %N")
qc("")
qc("| Genus | mean %N | sd %N | n samples | years |")
qc("|---|---|---|---|---|")
for (i in seq_len(nrow(genus_compare))) {
  qc(sprintf("| %s | %.3f | %.3f | %d | %s |", genus_compare$Genus[i],
             genus_compare$mean_N[i], genus_compare$sd_N[i], genus_compare$n[i],
             genus_compare$year_range[i]))
}
qc("")
qc("Sargassum tissue %N is ~70% higher than Turbinaria's on average (a ~2 SD ",
   "difference), and Sargassum sampling stops after 2014 while Turbinaria ",
   "runs 2007-2025 -- a genus effect confounded with a protocol-era effect. ",
   "**Decision: Turbinaria only**, not pooled, matching the original qmd.")
qc("")

cat("Genus comparison:\n")
print(genus_compare)

# =============================================================================
# D. Build N_sample.csv (sample-level, all years as recorded) and
#    N_site_year.csv (Site x Year summary, 2007-2024)
# =============================================================================

N_sample <- chn_raw |>
  filter(Habitat == "Backreef", Genus == "Turbinaria", !is.na(N)) |>
  select(Year, Site, Habitat, Genus, Dry_Weight, C, H, N, CN_ratio)

write_csv(N_sample, here("Revision", "Data", "derived", "N_sample.csv"))

N_site_year_full <- N_sample |>
  group_by(Year, Site) |>
  summarise(
    N_percent_mean = mean(N, na.rm = TRUE),
    N_percent_sd   = sd(N, na.rm = TRUE),
    CN_ratio_mean  = mean(CN_ratio, na.rm = TRUE),
    n_samples      = n(),
    .groups = "drop"
  )

sites_all_n <- sort(unique(N_sample$Site))
years_core <- 2007:2024

N_site_year <- expand_grid(Site = sites_all_n, Year = years_core) |>
  left_join(N_site_year_full, by = c("Site", "Year")) |>
  arrange(Site, Year)

write_csv(N_site_year, here("Revision", "Data", "derived", "N_site_year.csv"))

n_missing_site_years <- sum(is.na(N_site_year$N_percent_mean))

qc("## D. Resulting tables")
qc("")
qc("- `N_sample.csv`: ", nrow(N_sample), " individual Turbinaria tissue samples ",
   "(Backreef, all 6 sites, ", min(N_sample$Year), "-", max(N_sample$Year), ", as recorded).")
qc("- `N_site_year.csv`: ", nrow(N_site_year), " rows (", length(sites_all_n),
   " sites x ", length(years_core), " years, ", min(years_core), "-", max(years_core),
   "). Missing site-years (no Turbinaria samples that year): ", n_missing_site_years, ".")
qc("")
missing_site_years <- N_site_year |> filter(is.na(N_percent_mean)) |> select(Site, Year)
if (nrow(missing_site_years) > 0) {
  qc("| Site | Year | n samples |")
  qc("|---|---|---|")
  for (i in seq_len(nrow(missing_site_years))) {
    qc(sprintf("| %s | %d | 0 |", missing_site_years$Site[i], missing_site_years$Year[i]))
  }
  qc("")
}
qc("Per-sample-count distribution (non-missing site-years): min = ",
   min(N_site_year_full$n_samples), ", median = ", median(N_site_year_full$n_samples),
   ", max = ", max(N_site_year_full$n_samples), ".")
qc("")
qc("2005-2006 (no Turbinaria samples at all that early -- Sargassum-only or ",
   "no CHN sampling) and 2025 (5 samples, 1 site -- processing lag, same ",
   "pattern as other 2025-incomplete data in this project) are excluded ",
   "from `N_site_year.csv`'s core 2007-2024 window but retained as recorded ",
   "in `N_sample.csv`. Downstream models using `N_site_year.csv` will need ",
   "an explicit missing-data strategy for the one remaining gap (LTER_5, ",
   "2007) within the core window -- left unresolved here, not imputed.")
qc("")

writeLines(qc_lines, here("Revision", "Output", "qc_nutrient.md"))

cat("Done. Wrote:\n",
    "  Revision/Data/derived/N_sample.csv\n",
    "  Revision/Data/derived/N_site_year.csv\n",
    "  Revision/Output/qc_nutrient.md\n")
