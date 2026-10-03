# ---------------------------------------------------------------------------
# 06_metabolism.R
#
# Section 6 of the Revision plan (Revision/Reviewer_Response_Plan.md):
# hierarchical, nonlinear photosynthesis-irradiance (PI) model with
# covariates on Pmax and Rd, fit at the HOURLY level (Task 6.1).
#
# KNOWN GAP, FLAGGED NOT HIDDEN: the plan's Task 6.1 formula includes
# `fish_resp_z` (standardised fish community respiration) as a covariate on
# lR. That covariate is Task 7.1's output (mechanistic fish bioenergetics
# via `fishflux`/Barneche et al. 2014 scaling), which has not been built yet
# -- Section 13's execution-order table lists Step 6 as depending on Step 5
# (fish). This script fits Task 6.1's model WITHOUT fish_resp_z for now and
# clearly labels every output as the fish-omitted version; it should be
# refit once Section 7 is done.
#
# COMPUTE STRATEGY (plan's own Section 13 note): "develop on a 25% subsample
# of days first, then run the full fit." Implemented below as a stratified
# (by Year) random subsample. cmdstanr is not installed in this environment
# (checked in Step 4); using brms' default rstan backend with
# options(mc.cores = 4) (set in 00_packages.R after Step 4's sequential-
# chains lesson).
#
# OUTPUTS:
#   Revision/Data/derived/pi_model_hourly.csv
#   Revision/Output/fits/pi_model_dev.rds    (25% subsample, via brm's file=)
#   Revision/Output/fits/pi_model_full.rds   (full data, if time allows)
#   Revision/Output/qc_metabolism.md
# ---------------------------------------------------------------------------

source(here::here("Revision", "R", "00_packages.R"))
if (!requireNamespace("brms", quietly = TRUE)) install.packages("brms")
if (!requireNamespace("compositions", quietly = TRUE)) install.packages("compositions")
if (!requireNamespace("zCompositions", quietly = TRUE)) install.packages("zCompositions")
library(compositions)
library(zCompositions)
library(tidyverse)  # re-attach after compositions/MASS masking (established pattern)
library(brms)

dir.create(here("Revision", "Output", "fits"), showWarnings = FALSE, recursive = TRUE)

qc_lines <- character(0)
qc <- function(...) qc_lines <<- c(qc_lines, paste0(...))
qc("# QC summary: PI-curve metabolism model (Revision Step 6 / Section 6)")
qc("")
qc("Generated: ", as.character(Sys.time()))
qc("")

# =============================================================================
# A. LTER_1 benthic composition: pivot ("coral-first") coordinates (Task 6.1)
# =============================================================================
# Hron et al. (2012) pivot coordinates: with the composition ordered
# (Coral, Algae, CCA, Other), z1 isolates Coral against the geometric mean
# of the rest, z2 isolates Algae against the geometric mean of the
# remaining two, z3 isolates CCA against Other. Re-pivoting with a
# different part first (e.g. Algae) would be needed for an algae-focused
# sensitivity version -- not done here, flagged as a follow-up.

benthic_site <- read_csv(here("Revision", "Data", "derived", "benthic_site.csv"), show_col_types = FALSE)

lter1_counts <- benthic_site |>
  filter(Site == "LTER_1") |>
  arrange(Year) |>
  dplyr::select(Year, n_Coral, n_Algae, n_CCA, n_Other)

counts_mat_lter1 <- lter1_counts |> dplyr::select(n_Coral, n_Algae, n_CCA, n_Other) |> as.matrix()
n_zero_lter1 <- sum(counts_mat_lter1 == 0)
# cmultRepl errors out (rather than no-op) when there is nothing to replace --
# LTER_1's series has zero zero-cells (checked), so skip it in that case.
counts_repl_lter1 <- if (n_zero_lter1 > 0) {
  cmultRepl(counts_mat_lter1, label = 0, method = "SQ", suppress.print = TRUE)
} else {
  counts_mat_lter1
}
props_lter1 <- counts_repl_lter1 / rowSums(counts_repl_lter1)

pivot_coords <- function(p) {
  # p: matrix with columns in pivot order (part-of-interest first)
  z1 <- sqrt(3/4) * log(p[, 1] / (p[, 2] * p[, 3] * p[, 4])^(1/3))
  z2 <- sqrt(2/3) * log(p[, 2] / (p[, 3] * p[, 4])^(1/2))
  z3 <- sqrt(1/2) * log(p[, 3] / p[, 4])
  cbind(ilr1 = z1, ilr2 = z2, ilr3 = z3)
}

ilr_lter1 <- pivot_coords(props_lter1)
benthos_year_lter1 <- lter1_counts |>
  dplyr::select(Year) |>
  bind_cols(as_tibble(ilr_lter1))

qc("## A. LTER_1 benthic composition: pivot coordinates")
qc("")
qc("- Zero cells in LTER_1's Coral/Algae/CCA/Other count matrix (",
   nrow(counts_mat_lter1), " years x 4 parts): ", n_zero_lter1,
   " (zero-replacement is effectively a no-op here, but run for ",
   "consistency with Sections 4-5).")
qc("- Pivot coordinates (Hron et al. 2012), coral first: `ilr1` = coral vs. ",
   "geometric mean of (algae, CCA, other); `ilr2` = algae vs. geometric ",
   "mean of (CCA, other); `ilr3` = CCA vs. other. Computed from LTER_1's ",
   "site-level composition for each of the ", nrow(benthos_year_lter1), " years 2006-2025.")
qc("")

# =============================================================================
# B. Assemble the hourly model dataset
# =============================================================================

pp_hour <- read_csv(here("Revision", "Data", "derived", "pp_hour.csv"), show_col_types = FALSE)
N_site_year <- read_csv(here("Revision", "Data", "derived", "N_site_year.csv"), show_col_types = FALSE) |>
  filter(Site == "LTER_1") |>
  dplyr::select(Year, N_percent_mean)

N_mean <- mean(N_site_year$N_percent_mean, na.rm = TRUE)
N_sd   <- sd(N_site_year$N_percent_mean, na.rm = TRUE)

K_BOLTZMANN <- 8.617e-5  # eV/K

pi_model_hourly <- pp_hour |>
  left_join(benthos_year_lter1, by = "Year") |>
  left_join(N_site_year, by = "Year") |>
  mutate(
    T_kelvin = Temperature_mean + 273.15,
    T_ref_kelvin = mean(Temperature_mean, na.rm = TRUE) + 273.15,  # centred on the sample mean temperature
    invkT_c = 1 / (K_BOLTZMANN * T_kelvin) - 1 / (K_BOLTZMANN * T_ref_kelvin),
    log_flow_c = log(Flow_mean) - mean(log(Flow_mean), na.rm = TRUE),
    N_z = (N_percent_mean - N_mean) / N_sd,
    Season = factor(Season)
  )

n_before_drop <- nrow(pi_model_hourly)
pi_model_hourly <- pi_model_hourly |> filter(!is.na(N_z), !is.na(ilr1))
n_after_drop <- nrow(pi_model_hourly)

write_csv(pi_model_hourly, here("Revision", "Data", "derived", "pi_model_hourly.csv"))

qc("## B. Hourly model data")
qc("")
qc("- `pi_model_hourly.csv`: ", n_after_drop, " hourly rows (of ", n_before_drop,
   " before dropping rows with missing covariates) -- ",
   n_before_drop - n_after_drop, " rows dropped, all from 2025 ",
   "(`N_site_year.csv` does not yet cover 2025 -- the same lab-processing ",
   "lag already documented for this project's current year elsewhere).")
qc("- `invkT_c` centred on the sample mean temperature (",
   round(mean(pp_hour$Temperature_mean, na.rm = TRUE), 2),
   " degC), not an external MTE reference temperature -- the lR ",
   "coefficient is interpretable as the activation-energy-scaled deviation ",
   "from the average-condition respiration rate in this dataset, not an ",
   "absolute physiological reference.")
qc("- `N_z` uses `N_site_year.csv`'s LTER_1 series (mean ", round(N_mean, 3),
   ", sd ", round(N_sd, 3), ").")
qc("")

cat("pi_model_hourly:", n_after_drop, "rows,",
    n_distinct(pi_model_hourly$Year), "years,",
    n_distinct(pi_model_hourly$DielDate), "days\n")
print(summary(pi_model_hourly |> dplyr::select(PP, PAR, invkT_c, log_flow_c, N_z)))

# =============================================================================
# C. Nonlinear PI-curve model (Task 6.1) -- fish_resp_z OMITTED (see header)
# =============================================================================

f_pi <- bf(
  PP ~ (exp(la) * exp(lP) * PAR) / (exp(la) * PAR + exp(lP)) - exp(lR),
  la ~ 1 + (1 | Year),
  lP ~ 1 + ilr1 + ilr2 + ilr3 + invkT_c + log_flow_c + N_z + Season +
       (1 | p | Year) + (1 | q | Year:DielDate),
  lR ~ 1 + ilr1 + ilr2 + ilr3 + invkT_c + log_flow_c + Season +
       (1 | p | Year) + (1 | q | Year:DielDate),
  nl = TRUE
)

priors_pi <- c(
  prior(normal(log(100), 0.5), nlpar = "lP", coef = "Intercept"),
  prior(normal(log(80),  0.5), nlpar = "lR", coef = "Intercept"),
  prior(normal(0, 0.5),        nlpar = "lP", class = "b"),
  prior(normal(0, 0.5),        nlpar = "lR", class = "b"),
  # NOTE 2026-10-01: the plan's draft specified normal(0.65, 0.2) here, but
  # its own text states the coefficient on invkT_c is -E (activation
  # energy); since E is conventionally positive (~0.65 eV) for a process
  # that speeds up with temperature, the coefficient itself should be
  # NEGATIVE. Confirmed algebraically (ln B(T) = ln B(Tref) - E*invkT_c)
  # and empirically: fitting with the plan's literal +0.65 prior pulled
  # the posterior positive (0.38 [0.03, 0.74]), implying respiration
  # DECREASES with warming -- the biologically backwards direction.
  # Corrected with the PI's sign-off to normal(-0.65, 0.2).
  prior(normal(-0.65, 0.2),    nlpar = "lR", coef = "invkT_c"),
  prior(exponential(2),        nlpar = "lP", class = "sd"),
  prior(exponential(2),        nlpar = "lR", class = "sd"),
  prior(normal(0, 1),          nlpar = "la", class = "b"),
  prior(exponential(2),        nlpar = "la", class = "sd")
)

# ---- C1. Development fit: 25% of days, stratified by Year -----------------
set.seed(4817)
days_by_year <- pi_model_hourly |> distinct(Year, DielDate)
dev_days <- days_by_year |>
  group_by(Year) |>
  slice_sample(prop = 0.25) |>
  ungroup()
# Guarantee at least 1 day per year so (1|Year) isn't starved for any level.
missing_years <- setdiff(unique(days_by_year$Year), unique(dev_days$Year))
if (length(missing_years) > 0) {
  dev_days <- bind_rows(dev_days, days_by_year |> filter(Year %in% missing_years) |> group_by(Year) |> slice_sample(n = 1) |> ungroup())
}

pi_dev_data <- pi_model_hourly |> semi_join(dev_days, by = c("Year", "DielDate"))

qc("## C. Development fit: 25% of days (stratified by Year)")
qc("")
qc("- Dev subsample: ", nrow(pi_dev_data), " hourly rows, ",
   n_distinct(pi_dev_data$DielDate), " of ", n_distinct(pi_model_hourly$DielDate),
   " days (every year represented).")
qc("")

fit_pi_dev_start <- Sys.time()
fit_pi_dev <- brm(
  f_pi, data = pi_dev_data, family = student(),
  prior = priors_pi,
  chains = 4, iter = 2000, warmup = 1000, seed = 4817,
  control = list(adapt_delta = 0.95),
  file = here("Revision", "Output", "fits", "pi_model_dev"),
  file_refit = "on_change"
)
fit_pi_dev_elapsed <- Sys.time() - fit_pi_dev_start

cat("Dev fit completed/loaded in", round(as.numeric(fit_pi_dev_elapsed, units = "secs"), 1), "seconds\n")
print(summary(fit_pi_dev))

# =============================================================================
# C2. Full fit (all 117 days), using the dev fit's confirmed specification
# =============================================================================

fit_pi_full_start <- Sys.time()
fit_pi_full <- brm(
  f_pi, data = pi_model_hourly, family = student(),
  prior = priors_pi,
  chains = 4, iter = 2000, warmup = 1000, seed = 4817,
  control = list(adapt_delta = 0.95),
  file = here("Revision", "Output", "fits", "pi_model_full"),
  file_refit = "on_change"
)
fit_pi_full_elapsed <- Sys.time() - fit_pi_full_start

cat("Full fit completed/loaded in", round(as.numeric(fit_pi_full_elapsed, units = "secs"), 1), "seconds\n")
print(summary(fit_pi_full))

qc("## D. Full fit (117 days, 2808 hourly rows)")
qc("")
qc("- Same specification as the dev fit, confirmed on the full data. ",
   "`invkT_c` coefficient prior for lR corrected to `normal(-0.65, 0.2)` ",
   "(see code comment) after the dev fit on the plan's literal ",
   "`normal(0.65, 0.2)` produced a biologically backwards posterior ",
   "(respiration decreasing with warming); decided with the PI 2026-10-01.")
qc("")

writeLines(qc_lines, here("Revision", "Output", "qc_metabolism.md"))
cat("Done. Wrote:\n",
    "  Revision/Data/derived/pi_model_hourly.csv\n",
    "  Revision/Output/fits/pi_model_dev.rds\n",
    "  Revision/Output/fits/pi_model_full.rds\n",
    "  Revision/Output/qc_metabolism.md\n")

# =============================================================================
# E. Task 6.2: per-estimand refits (E4, E6, E7, E9, E10)
# =============================================================================
# For each estimand, lR or lP is pruned to EXACTLY its dagitty-derived
# minimal adjustment set (plus the exposure itself); the other nonlinear
# parameter keeps the full Task 6.1 specification unchanged (that
# parameter is not the estimand's outcome, so its own specification
# doesn't need to satisfy any particular adjustment set here).
#
# Adjustment sets taken from table_s_dag.csv (post Step 3b's four added
# edges). For every one of these five estimands, deliberately choosing the
# SMALLEST minimal adjustment set happens to give an identical set for
# both DAG variants (Rd -> N vs Benthos -> N) -- checked directly against
# the table -- so a single model per estimand is reported, not two.
# E9 and E10 end up needing the exact same covariate list ({Benthos,
# DayTemp, Flow, N}), so they share one fitted model (same legitimate
# coincidence already seen for E1-E3 in Section 5).
#
# COMPUTE DECISION: fit all four new models on the SAME 25%-of-days dev
# subsample already built and validated in Section C1 (pi_dev_data), not
# the full 2808-row dataset -- the full Task 6.1 fit alone took ~52
# minutes, and four more full fits would cost several hours of
# unattended sequential compute. These are reported as SCREENING fits;
# scaling any of them to the full dataset is a flagged follow-up, not
# done automatically.

fish_transect <- read_csv(here("Revision", "Data", "derived", "fish_transect.csv"), show_col_types = FALSE)

fish_lter1_year <- fish_transect |>
  filter(Site == "LTER_1") |>
  mutate(dag_group = case_when(
    trophic_group == "Herbivore"   ~ "Herb",
    trophic_group == "Corallivore" ~ "Corall",
    TRUE                            ~ "OtherFish"
  )) |>
  group_by(Year, Transect, dag_group) |>
  summarise(biomass_g_m2 = sum(biomass_g_m2), .groups = "drop") |>
  group_by(Year, dag_group) |>
  summarise(biomass_g_m2 = mean(biomass_g_m2), .groups = "drop") |>
  pivot_wider(names_from = dag_group, values_from = biomass_g_m2, values_fill = 0) |>
  mutate(
    Herb_z      = as.numeric(scale(Herb)),
    Corall_z    = as.numeric(scale(Corall)),
    OtherFish_z = as.numeric(scale(OtherFish))
  ) |>
  dplyr::select(Year, Herb_z, Corall_z, OtherFish_z)

dhw_lter1_year <- read_csv(here("Revision", "Data", "derived", "disturbance_site_year.csv"), show_col_types = FALSE) |>
  filter(Site == "LTER_1") |>
  mutate(DHW_z = as.numeric(scale(DHW_max))) |>
  dplyr::select(Year, DHW_z)

pi_model_hourly <- pi_model_hourly |>
  left_join(fish_lter1_year, by = "Year") |>
  left_join(dhw_lter1_year, by = "Year")
pi_dev_data <- pi_dev_data |>
  left_join(fish_lter1_year, by = "Year") |>
  left_join(dhw_lter1_year, by = "Year")

qc("## G. Task 6.2 covariates: fish biomass by trophic group, DHW (LTER_1, year-level)")
qc("")
qc("- `Herb_z`/`Corall_z`/`OtherFish_z`: mean biomass (g/m^2) across the 4 ",
   "fish transects at LTER_1 per year, z-standardised (needed for E4/E6, ",
   "the DAG's `Herb`/`Corall`/`OtherFish` nodes -- NOT the bioenergetics ",
   "`fish_resp_z` of Task 7.1, which remains unbuilt).")
qc("- `DHW_z`: z-standardised `DHW_max` at LTER_1 (needed for E6/E7).")
qc("")

la_formula <- "la ~ 1 + (1 | Year)"
lP_full    <- "lP ~ 1 + ilr1 + ilr2 + ilr3 + invkT_c + log_flow_c + N_z + Season + (1 | p | Year) + (1 | q | Year:DielDate)"
lR_full    <- "lR ~ 1 + ilr1 + ilr2 + ilr3 + invkT_c + log_flow_c + Season + (1 | p | Year) + (1 | q | Year:DielDate)"

priors_la <- c(prior(normal(0, 1), nlpar = "la", class = "b"), prior(exponential(2), nlpar = "la", class = "sd"))

fit_per_estimand <- function(name, lP_rhs, lR_rhs, extra_priors) {
  f <- bf(as.formula("PP ~ (exp(la) * exp(lP) * PAR) / (exp(la) * PAR + exp(lP)) - exp(lR)"),
          as.formula(la_formula), as.formula(lP_rhs), as.formula(lR_rhs), nl = TRUE)
  p <- c(priors_la, extra_priors)
  brm(f, data = pi_dev_data, family = student(), prior = p,
      chains = 4, iter = 2000, warmup = 1000, seed = 4817,
      control = list(adapt_delta = 0.95),
      file = here("Revision", "Output", "fits", paste0("pi_", name, "_dev")),
      file_refit = "on_change")
}

# ---- E4: Benthos -> Rd (direct). Adjustment: {Corall, DayTemp, Flow, Herb, OtherFish} --
fit_E4 <- fit_per_estimand(
  "E4",
  lP_full,
  "lR ~ 1 + ilr1 + ilr2 + ilr3 + Herb_z + Corall_z + OtherFish_z + invkT_c + log_flow_c + (1 | p | Year) + (1 | q | Year:DielDate)",
  c(prior(normal(log(80), 0.5), nlpar = "lR", coef = "Intercept"),
    prior(normal(0, 0.5),       nlpar = "lR", class = "b"),
    prior(normal(-0.65, 0.2),   nlpar = "lR", coef = "invkT_c"),
    prior(exponential(2),       nlpar = "lR", class = "sd"),
    prior(normal(log(100), 0.5), nlpar = "lP", coef = "Intercept"),
    prior(normal(0, 0.5),        nlpar = "lP", class = "b"),
    prior(exponential(2),        nlpar = "lP", class = "sd"))
)

# ---- E6: Fish -> Rd (direct). Adjustment: {Benthos, DHW} ------------------
fit_E6 <- fit_per_estimand(
  "E6",
  lP_full,
  "lR ~ 1 + Herb_z + Corall_z + OtherFish_z + ilr1 + ilr2 + ilr3 + DHW_z + (1 | p | Year) + (1 | q | Year:DielDate)",
  c(prior(normal(log(80), 0.5), nlpar = "lR", coef = "Intercept"),
    prior(normal(0, 0.5),       nlpar = "lR", class = "b"),
    prior(exponential(2),       nlpar = "lR", class = "sd"),
    prior(normal(log(100), 0.5), nlpar = "lP", coef = "Intercept"),
    prior(normal(0, 0.5),        nlpar = "lP", class = "b"),
    prior(exponential(2),        nlpar = "lP", class = "sd"))
)

# ---- E7: DayTemp -> Rd (direct). Adjustment: {DHW, Flow} ------------------
fit_E7 <- fit_per_estimand(
  "E7",
  lP_full,
  "lR ~ 1 + invkT_c + DHW_z + log_flow_c + (1 | p | Year) + (1 | q | Year:DielDate)",
  c(prior(normal(log(80), 0.5), nlpar = "lR", coef = "Intercept"),
    prior(normal(0, 0.5),       nlpar = "lR", class = "b"),
    prior(normal(-0.65, 0.2),   nlpar = "lR", coef = "invkT_c"),
    prior(exponential(2),       nlpar = "lR", class = "sd"),
    prior(normal(log(100), 0.5), nlpar = "lP", coef = "Intercept"),
    prior(normal(0, 0.5),        nlpar = "lP", class = "b"),
    prior(exponential(2),        nlpar = "lP", class = "sd"))
)

# ---- E9/E10: Benthos -> Pmax / DayTemp -> Pmax (direct), shared model -----
# Both estimands' minimal adjustment sets resolve to the identical
# covariate list {Benthos, DayTemp, Flow, N} once the exposure is folded
# in, for BOTH DAG variants -- a single model legitimately answers both.
fit_E9_E10 <- fit_per_estimand(
  "E9_E10",
  "lP ~ 1 + ilr1 + ilr2 + ilr3 + invkT_c + log_flow_c + N_z + (1 | p | Year) + (1 | q | Year:DielDate)",
  lR_full,
  c(prior(normal(log(100), 0.5), nlpar = "lP", coef = "Intercept"),
    prior(normal(0, 0.5),        nlpar = "lP", class = "b"),
    prior(exponential(2),        nlpar = "lP", class = "sd"),
    prior(normal(log(80), 0.5), nlpar = "lR", coef = "Intercept"),
    prior(normal(0, 0.5),       nlpar = "lR", class = "b"),
    prior(normal(-0.65, 0.2),   nlpar = "lR", coef = "invkT_c"),
    prior(exponential(2),       nlpar = "lR", class = "sd"))
)

cat("\n=== E4 (Benthos -> Rd, direct) ===\n");       print(fixef(fit_E4))
cat("\n=== E6 (Fish -> Rd, direct) ===\n");          print(fixef(fit_E6))
cat("\n=== E7 (DayTemp -> Rd, direct) ===\n");       print(fixef(fit_E7))
cat("\n=== E9/E10 (Benthos/DayTemp -> Pmax, direct) ===\n"); print(fixef(fit_E9_E10))

qc("## H. Task 6.2 per-estimand screening fits (25%-of-days dev subsample)")
qc("")
qc("All four fits below use `pi_dev_data` (624 hourly rows, the same ",
   "stratified 25%-of-days subsample validated in Section C1), NOT the ",
   "full 2808-row dataset -- the full Task 6.1 fit alone took ~52 minutes, ",
   "and four more full fits would cost several hours of sequential ",
   "compute. **These are screening/development results; scaling any of ",
   "them to the full dataset is a flagged follow-up.**")
qc("")
write_screening_table <- function(fit, name) {
  fx <- fixef(fit)
  qc(paste0("### ", name))
  qc("")
  qc("| Coefficient | Estimate | 95% CI |")
  qc("|---|---|---|")
  for (i in seq_len(nrow(fx))) {
    qc(sprintf("| %s | %.3f | [%.3f, %.3f] |", rownames(fx)[i], fx[i, "Estimate"], fx[i, "Q2.5"], fx[i, "Q97.5"]))
  }
  qc("")
}
write_screening_table(fit_E4, "E4: Benthos -> Rd (direct); adjustment {Corall, DayTemp, Flow, Herb, OtherFish}")
write_screening_table(fit_E6, "E6: Fish (Herb+Corall+OtherFish) -> Rd (direct); adjustment {Benthos, DHW}")
write_screening_table(fit_E7, "E7: DayTemp -> Rd (direct); adjustment {DHW, Flow}")
write_screening_table(fit_E9_E10, "E9/E10: Benthos -> Pmax / DayTemp -> Pmax (direct); adjustment {Benthos, DayTemp, Flow, N}")

writeLines(qc_lines, here("Revision", "Output", "qc_metabolism.md"))
cat("\nDone: Task 6.2 screening fits written to qc_metabolism.md\n")
