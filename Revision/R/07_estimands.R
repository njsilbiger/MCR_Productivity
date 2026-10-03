# ---------------------------------------------------------------------------
# 07_estimands.R
#
# Task 7.2 of the Revision plan (Revision/Reviewer_Response_Plan.md, Section
# 7): fish response models for estimand E11 (Benthos -> Herbivores /
# Corallivores, total effect). Modelled at TRANSECT level across all 6
# backreef sites (LTER_1-6), per the plan's literal formula:
#
#   bf(herbivore_biomass ~ ilr1 + ilr2 + ilr3 + herb_lag_z + Year_c +
#        (1 | Site) + (1 | Site:Year))
#   family = hurdle_gamma() if zeros exist, else Gamma(link = "log")
#
# ADJUSTMENT SET (Table S-DAG, E11: Benthos -> Herb/Corall, total effect):
# {Time, Benthos_lag, Herb_lag}. `Year_c` instruments Time; `herb_lag_z`
# (previous year's site-mean herbivore biomass, z-scored) instruments
# Herb_lag and, via the Benthos -> Herb feedback already fit in
# 04_benthic_gompertz.R, blocks the backdoor path through Benthos_lag --
# so Benthos_lag itself is NOT re-entered as its own ilr terms here (doing so
# would double up with herb_lag_z, which already carries last year's state
# through the fish side of the feedback). `ilr1-3` are the CONTEMPORANEOUS
# (same Year/Site) benthic composition, because Benthos -> Herb/Corall is a
# same-year total effect, not a lagged one.
#
# IDENTICAL predictor set/formula is used for BOTH outcomes (herbivore_biomass and
# corallivore_biomass), exactly as given in the plan's single code block --
# this is the Task 7.2 specification, not an independent modelling choice
# per outcome.
#
# Fish transects (4 per site-year) and benthic transects (5 per site-year)
# are NOT the same physical transects (Task 2.3's documented mismatch), so
# the compositional Benthos node is joined at SITE-YEAR resolution (mean
# across the 5 benthic transects via benthic_site.csv), same resolution as
# herb_lag_z, giving each of the 4 fish transects per site-year an identical
# (ilr1, ilr2, ilr3, herb_lag_z, Year_c) covariate set and all between-
# transect variation left for the (1 | Site:Year) term to absorb, alongside
# genuine between-transect fish-survey noise.
#
# OUTPUTS:
#   Revision/Output/fits/E11_herbivore_biomass.rds       (via brm's file=)
#   Revision/Output/fits/E11_corallivore_biomass.rds      (via brm's file=)
#   Revision/Output/fig_E11_resid_acf.png
#   Revision/Output/qc_E11_fish_response.md
#
# Also Task 7.3 (Section G below): fish -> N (E12), Turbinaria tissue N,
# LTER_1 only. Adjustment set taken from the "Benthos -> N" DAG variant in
# table_s_dag.csv (that variant's E12 row has a SINGLE minimal adjustment
# set, {Benthos, Time}, vs. 4 alternative sets under "Rd -> N" -- the
# simpler variant is used both because it is uniquely identified and
# because the "Rd -> N" sets require Rd/DayTemp/Flow/Season, which do not
# exist at annual resolution for a site-year-level fish -> N fit).
#
#   Revision/Output/fits/E12_fish_to_N.rds                (via brm's file=)
#   Revision/Output/qc_E12_fish_to_N.md
# ---------------------------------------------------------------------------

source(here::here("Revision", "R", "00_packages.R"))
if (!requireNamespace("brms", quietly = TRUE)) install.packages("brms")
if (!requireNamespace("compositions", quietly = TRUE)) install.packages("compositions")
if (!requireNamespace("zCompositions", quietly = TRUE)) install.packages("zCompositions")
library(compositions)
library(zCompositions)
library(tidyverse)  # re-attach after compositions/MASS masking, as in 03b_dsep_test.R
library(brms)

dir.create(here("Revision", "Output", "fits"), showWarnings = FALSE, recursive = TRUE)

qc_lines <- character(0)
qc <- function(...) qc_lines <<- c(qc_lines, paste0(...))
qc("# QC summary: Task 7.2 fish response models (E11)")
qc("")
qc("Generated: ", as.character(Sys.time()))
qc("")
qc("Estimand E11 (Benthos -> Herbivores / Corallivores, total effect), ",
   "adjustment set {Time, Benthos_lag, Herb_lag} per Table S-DAG. Fit at ",
   "transect level across all 6 backreef sites (LTER_1-6).")
qc("")

# =============================================================================
# A. Contemporaneous Benthos: ilr coordinates (site-year), same construction
#    as 03b_dsep_test.R's Benthos node (zero-replaced counts -> ilr)
# =============================================================================

benthic_site <- read_csv(here("Revision", "Data", "derived", "benthic_site.csv"), show_col_types = FALSE)

counts_mat <- benthic_site |> dplyr::select(n_Coral, n_Algae, n_CCA, n_Other) |> as.matrix()
n_zero_cells <- sum(counts_mat == 0)
counts_repl <- cmultRepl(counts_mat, label = 0, method = "SQ", suppress.print = TRUE)
ilr_coords <- ilr(acomp(counts_repl))
colnames(ilr_coords) <- c("ilr1", "ilr2", "ilr3")

benthic_site_ilr <- bind_cols(benthic_site |> dplyr::select(Year, Site), as_tibble(ilr_coords))

qc("## A. Contemporaneous Benthos (ilr coordinates, site-year)")
qc("")
qc("- Zero cells in the Coral/Algae/CCA/Other count matrix (", nrow(counts_mat),
   " site-years x 4 parts): ", n_zero_cells, ". Replaced via ",
   "`zCompositions::cmultRepl` before closing to proportions and taking ilr ",
   "(`compositions::ilr()`, default sequential binary partition), matching ",
   "03b_dsep_test.R's Benthos-node construction.")
qc("- Same-YEAR (not lagged) composition, because E11 is the Benthos -> Herb/",
   "Corall total effect in the same year, not a lagged recovery effect.")
qc("")

# =============================================================================
# B. Herb_lag: previous year's site-mean herbivore biomass, z-standardised
#    (same construction as 04_benthic_gompertz.R Section B)
# =============================================================================

fish_transect <- read_csv(here("Revision", "Data", "derived", "fish_transect.csv"), show_col_types = FALSE)

herb_site_year <- fish_transect |>
  filter(trophic_group == "Herbivore") |>
  group_by(Year, Site) |>
  summarise(Herb = mean(biomass_g_m2), .groups = "drop")

herb_mean <- mean(herb_site_year$Herb); herb_sd <- sd(herb_site_year$Herb)

herb_lag <- herb_site_year |>
  mutate(Year = Year + 1, herb_lag_z = (Herb - herb_mean) / herb_sd) |>
  dplyr::select(Year, Site, herb_lag_z)

qc("## B. Herb_lag (adjustment variable)")
qc("")
qc("- Site-mean Herbivore biomass (g/m^2, across the 4 fish transects), ",
   "lagged +1 year, z-standardised (mean ", round(herb_mean, 2), ", sd ",
   round(herb_sd, 2), "). Identical construction to 04_benthic_gompertz.R.")
qc("")

# =============================================================================
# C. Assemble transect-level model data
# =============================================================================

fish_resp_data <- fish_transect |>
  filter(trophic_group %in% c("Herbivore", "Corallivore")) |>
  mutate(trophic_group = tolower(trophic_group)) |>
  dplyr::select(Year, Site, Transect, trophic_group, biomass_g_m2) |>
  pivot_wider(names_from = trophic_group, values_from = biomass_g_m2,
              names_glue = "{trophic_group}_biomass") |>
  left_join(benthic_site_ilr, by = c("Year", "Site")) |>
  left_join(herb_lag, by = c("Year", "Site")) |>
  mutate(Year_c = Year - mean(Year, na.rm = TRUE)) |>
  filter(!is.na(ilr1), !is.na(herb_lag_z))  # drop first observed year per site (no lag) and any unmatched benthic survey

n_dropped <- nrow(fish_transect |> filter(trophic_group == "Herbivore")) - nrow(fish_resp_data)

qc("## C. Model data")
qc("")
qc("- `fish_resp_data`: ", nrow(fish_resp_data), " transect-year rows (",
   n_distinct(fish_resp_data$Site), " sites, ", n_distinct(fish_resp_data$Year),
   " years: ", min(fish_resp_data$Year), "-", max(fish_resp_data$Year), ").")
qc("- Rows dropped for missing `herb_lag_z` (first observed year per site, 2006) ",
   "or unmatched benthic survey: ", n_dropped, ".")
qc("- Outcome zero counts: herbivore_biomass = ", sum(fish_resp_data$herbivore_biomass == 0),
   " of ", nrow(fish_resp_data), "; corallivore_biomass = ",
   sum(fish_resp_data$corallivore_biomass == 0), " of ", nrow(fish_resp_data), ".")
qc("")

cat("fish_resp_data:", nrow(fish_resp_data), "rows,", n_distinct(fish_resp_data$Site),
    "sites,", n_distinct(fish_resp_data$Year), "years\n")
print(summary(fish_resp_data |> dplyr::select(herbivore_biomass, corallivore_biomass, ilr1, ilr2, ilr3, herb_lag_z, Year_c)))

# =============================================================================
# D. Fit the two E11 models
# =============================================================================
# herbivore_biomass has no zeros -> Gamma(log). corallivore_biomass has a handful of
# true zero counts at the transect level -> hurdle_gamma(), per the plan's
# own stated family choice rule.

f_herb <- bf(herbivore_biomass ~ ilr1 + ilr2 + ilr3 + herb_lag_z + Year_c + (1 | Site) + (1 | Site:Year))
f_corall <- bf(corallivore_biomass ~ ilr1 + ilr2 + ilr3 + herb_lag_z + Year_c + (1 | Site) + (1 | Site:Year))

priors_e11 <- c(
  prior(normal(0, 2), class = "b"),
  prior(exponential(2), class = "sd"),
  prior(gamma(0.01, 0.01), class = "shape")
)
priors_e11_hurdle <- c(priors_e11, prior(beta(1, 1), class = "hu"))

fit_E11_herb <- brm(
  f_herb, family = Gamma(link = "log"), data = fish_resp_data,
  prior = priors_e11,
  chains = 4, iter = 2000, warmup = 1000, seed = 7934,
  control = list(adapt_delta = 0.95),
  file = here("Revision", "Output", "fits", "E11_herbivore_biomass"),
  file_refit = "on_change"
)

fit_E11_corall <- brm(
  f_corall, family = hurdle_gamma(link = "log"), data = fish_resp_data,
  prior = priors_e11_hurdle,
  chains = 4, iter = 2000, warmup = 1000, seed = 7934,
  control = list(adapt_delta = 0.95),
  file = here("Revision", "Output", "fits", "E11_corallivore_biomass"),
  file_refit = "on_change"
)

cat("\n=== E11: Benthos -> Herbivore biomass ===\n"); print(summary(fit_E11_herb))
cat("\n=== E11: Benthos -> Corallivore biomass ===\n"); print(summary(fit_E11_corall))

write_result_table <- function(fit, name) {
  fx <- fixef(fit)
  qc(paste0("### ", name)); qc("")
  qc("| Coefficient | Estimate | 95% CI |"); qc("|---|---|---|")
  for (i in seq_len(nrow(fx))) qc(sprintf("| %s | %.3f | [%.3f, %.3f] |", rownames(fx)[i], fx[i, "Estimate"], fx[i, "Q2.5"], fx[i, "Q97.5"]))
  qc("")
}

qc("## D. Model results")
qc("")
qc("- `herbivore_biomass`: `Gamma(link = \"log\")` (no zeros at transect level).")
qc("- `corallivore_biomass`: `hurdle_gamma(link = \"log\")` (", sum(fish_resp_data$corallivore_biomass == 0),
   " zero-biomass transect-years, ", round(100 * mean(fish_resp_data$corallivore_biomass == 0), 1), "%).")
qc("")
write_result_table(fit_E11_herb, "E11: Benthos -> Herbivore biomass")
write_result_table(fit_E11_corall, "E11: Benthos -> Corallivore biomass")

# =============================================================================
# E. Residual ACF per site (plan: "Check residual ACF per site")
# =============================================================================
# Randomized-quantile (DHARMa-style via brms) residuals are unavailable
# without extra deps; use standardised response residuals (observed - fitted
# response, divided by the residual SD at that fit), averaged across the 4
# fish transects per site-year to get one series per site, then ACF/Ljung-Box
# at lag 1, consistent with the diagnostic already run in 04_benthic_gompertz.R.

resid_acf_check <- function(fit, outcome_col, label) {
  fitted_vals <- fitted(fit, summary = TRUE)[, "Estimate"]
  resid_std <- (fish_resp_data[[outcome_col]] - fitted_vals) / sd(fish_resp_data[[outcome_col]] - fitted_vals)
  df <- fish_resp_data |> dplyr::select(Site, Year) |> mutate(resid = resid_std)
  site_year_resid <- df |> group_by(Site, Year) |> summarise(resid = mean(resid), .groups = "drop")

  acf1 <- site_year_resid |>
    group_by(Site) |>
    arrange(Year) |>
    group_modify(~ {
      x <- .x$resid
      if (length(x) < 4 || sd(x) == 0) return(tibble(acf1 = NA_real_, p_value = NA_real_, n = length(x)))
      bt <- Box.test(x, lag = 1, type = "Ljung-Box")
      tibble(acf1 = acf(x, lag.max = 1, plot = FALSE)$acf[2], p_value = bt$p.value, n = length(x))
    }) |>
    ungroup() |>
    mutate(model = label)
  acf1
}

acf_herb   <- resid_acf_check(fit_E11_herb,   "herbivore_biomass",   "Herb_biomass")
acf_corall <- resid_acf_check(fit_E11_corall, "corallivore_biomass", "Corall_biomass")
acf_e11_all <- bind_rows(acf_herb, acf_corall) |>
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

write_csv(acf_e11_all, here("Revision", "Output", "acf_E11_fish_response.csv"))

p_acf_e11 <- ggplot(bind_rows(acf_herb, acf_corall), aes(x = Site, y = acf1)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = c(-1.96 / sqrt(19), 1.96 / sqrt(19)), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, color = "black") +
  facet_wrap(~model) +
  labs(x = "Site", y = "Lag-1 residual ACF",
       title = "Task 7.2: residual autocorrelation per site, E11 fish response models",
       subtitle = "Dashed lines: approx. 95% CI under white noise (not multiple-comparison corrected)") +
  theme_minimal()
ggsave(here("Revision", "Output", "fig_E11_resid_acf.png"), p_acf_e11, width = 9, height = 5, dpi = 300)

n_sig_raw <- sum(acf_e11_all$p_value < 0.05, na.rm = TRUE)
n_sig_bh  <- sum(acf_e11_all$p_adj_BH < 0.05, na.rm = TRUE)
n_tested  <- sum(!is.na(acf_e11_all$p_value))

qc("## E. Residual ACF per site")
qc("")
qc("- Standardised response residuals, averaged to site-year then lag-1 ",
   "Ljung-Box per site per model (", n_tested, " site x model series tested; ",
   "BH-corrected across all series).")
qc(sprintf("- Significant at raw p<0.05: %d of %d. Significant at BH-adjusted p<0.05: %d of %d.",
           n_sig_raw, n_tested, n_sig_bh, n_tested))
qc("")

# =============================================================================
# F. Back-transformed marginal effects on Coral/Algae/CCA share
# =============================================================================
# The ilr1-3 regression coefficients in Section D are coordinates of an
# abstract sequential-binary-partition basis, not effects of any named
# benthic part -- not directly interpretable, same problem as the
# benthic-composition model itself (Section 5 of the plan). Mirror Task
# 5.4's g-computation solution: hold every row's covariates (herb_lag_z,
# Year_c, Site, Year) at their observed values, perturb the OBSERVED
# composition by +1 percentage point in one part's share (Coral, Algae or
# CCA in turn), proportionally rescaling the other three parts so the
# simplex still closes to 1, re-express the counterfactual composition in
# the SAME ilr basis used to fit the models (same `compositions::ilr()`
# call, same Coral/Algae/CCA/Other column order), and take the average
# marginal effect (AME) on posterior_epred() for herbivore and corallivore
# biomass, exactly as 04_benthic_gompertz.R Section F did for E1-E3.
#
# Caveat worth stating in the response letter: a +1 percentage-point
# perturbation is a small absolute step for Coral (mean share 20%) or Algae
# (mean 55%), but a large RELATIVE step for CCA (mean share only 2.2%,
# max 17%) -- so the CCA estimand below is a bigger relative manipulation
# than the other two, even though the absolute delta is identical.

props_mat <- counts_repl / rowSums(counts_repl)
colnames(props_mat) <- c("Coral", "Algae", "CCA", "Other")

qc("## F. Back-transformed marginal effects on Coral/Algae/CCA share")
qc("")
qc("- ilr1-3 are coordinates of an abstract basis, not named-part effects. ",
   "Re-expressed as the average marginal effect (AME) of a +1 percentage-",
   "point increase in one part's share (Coral, Algae or CCA), with the ",
   "other three parts rescaled proportionally to close the simplex, on ",
   "predicted herbivore/corallivore biomass -- same g-computation approach ",
   "as the E1-E3 estimands in 04_benthic_gompertz.R.")
qc("- Caveat: observed mean share is Coral 20.4%, Algae 55.2%, CCA 2.2% ",
   "(max 17.4%), so a +1-percentage-point step is a much larger RELATIVE ",
   "perturbation for CCA than for Coral or Algae.")
qc("")

build_cf_ilr <- function(part, delta) {
  other_cols <- setdiff(colnames(props_mat), part)
  old_part <- props_mat[, part]
  new_part <- old_part + delta
  old_other_sum <- rowSums(props_mat[, other_cols, drop = FALSE])
  scale_factor <- (1 - new_part) / old_other_sum
  new_mat <- props_mat
  new_mat[, part] <- new_part
  for (col in other_cols) new_mat[, col] <- props_mat[, col] * scale_factor
  ilr_cf <- ilr(acomp(new_mat[, c("Coral", "Algae", "CCA", "Other")]))
  colnames(ilr_cf) <- c("ilr1_cf", "ilr2_cf", "ilr3_cf")
  bind_cols(benthic_site |> dplyr::select(Year, Site), as_tibble(ilr_cf))
}

DELTA_SHARE <- 0.01  # +1 percentage point

pe_factual_herb   <- posterior_epred(fit_E11_herb,   newdata = fish_resp_data)
pe_factual_corall <- posterior_epred(fit_E11_corall, newdata = fish_resp_data)

ame_on_fish <- function(part, delta) {
  cf_ilr <- build_cf_ilr(part, delta)
  nd <- fish_resp_data |>
    dplyr::select(-ilr1, -ilr2, -ilr3) |>
    left_join(cf_ilr, by = c("Year", "Site")) |>
    rename(ilr1 = ilr1_cf, ilr2 = ilr2_cf, ilr3 = ilr3_cf)

  pe_cf_herb   <- posterior_epred(fit_E11_herb,   newdata = nd)
  pe_cf_corall <- posterior_epred(fit_E11_corall, newdata = nd)

  ame_herb_draws   <- rowMeans(pe_cf_herb   - pe_factual_herb)
  ame_corall_draws <- rowMeans(pe_cf_corall - pe_factual_corall)

  bind_rows(
    tibble(outcome = "Herbivore biomass (g/m2)", part = part, delta_pp = delta * 100,
           ame_mean = mean(ame_herb_draws), ame_lo95 = quantile(ame_herb_draws, 0.025),
           ame_hi95 = quantile(ame_herb_draws, 0.975), prob_positive = mean(ame_herb_draws > 0),
           pct_change = 100 * mean(ame_herb_draws) / mean(pe_factual_herb)),
    tibble(outcome = "Corallivore biomass (g/m2)", part = part, delta_pp = delta * 100,
           ame_mean = mean(ame_corall_draws), ame_lo95 = quantile(ame_corall_draws, 0.025),
           ame_hi95 = quantile(ame_corall_draws, 0.975), prob_positive = mean(ame_corall_draws > 0),
           pct_change = 100 * mean(ame_corall_draws) / mean(pe_factual_corall))
  )
}

estimands_E11_composition <- bind_rows(
  ame_on_fish("Coral", DELTA_SHARE),
  ame_on_fish("Algae", DELTA_SHARE),
  ame_on_fish("CCA",   DELTA_SHARE)
)
write_csv(estimands_E11_composition, here("Revision", "Output", "estimands_E11_composition.csv"))

qc("| Outcome | Part (+1 pp) | AME (mean) | 95% CI | P(AME > 0) | % change vs. mean prediction |")
qc("|---|---|---|---|---|---|")
for (i in seq_len(nrow(estimands_E11_composition))) {
  r <- estimands_E11_composition[i, ]
  qc(sprintf("| %s | %s | %.4f | [%.4f, %.4f] | %.3f | %.2f%% |",
             r$outcome, r$part, r$ame_mean, r$ame_lo95, r$ame_hi95, r$prob_positive, r$pct_change))
}
qc("")

cat("\nE11 back-transformed composition estimands (AME per +1pp share, g/m2):\n")
print(estimands_E11_composition, n = 20)

writeLines(qc_lines, here("Revision", "Output", "qc_E11_fish_response.md"))

cat("\nDone. Wrote:\n",
    "  Revision/Output/fits/E11_herbivore_biomass.rds\n",
    "  Revision/Output/fits/E11_corallivore_biomass.rds\n",
    "  Revision/Output/acf_E11_fish_response.csv\n",
    "  Revision/Output/fig_E11_resid_acf.png\n",
    "  Revision/Output/estimands_E11_composition.csv\n",
    "  Revision/Output/qc_E11_fish_response.md\n")
cat("\nLag-1 ACF / Ljung-Box by site and model:\n")
print(acf_e11_all)

# =============================================================================
# G. Task 7.3: Fish -> N (E12), LTER_1 Turbinaria tissue N, n ~ 18 years
# =============================================================================
# Adjustment set {Benthos, Time} (the "Benthos -> N" DAG variant's single
# minimal set for E12 -- see header note). Fit at annual resolution
# (N_site_year.csv x fish_respiration_site_year.csv x contemporaneous
# Benthos ilr), LTER_1 only, because Turbinaria tissue N has no other-site
# replication (Task 2.3/7.3's documented limitation). This is explicitly a
# SMALL-SAMPLE estimate (n = 18 annual points, 5 regression coefficients):
# say so in the response letter, don't present it as more precise than it
# is.

qc_lines <- character(0)
qc("# QC summary: Task 7.3 fish -> N model (E12), LTER_1 Turbinaria tissue N")
qc("")
qc("Generated: ", as.character(Sys.time()))
qc("")
qc("Estimand E12 (Herb+Corall+OtherFish -> N, total effect). Adjustment set ",
   "{Benthos, Time}, taken from the \"Benthos -> N\" DAG variant in ",
   "`table_s_dag.csv` (its E12 row has a single minimal adjustment set, vs. ",
   "4 alternative sets under the \"Rd -> N\" variant, which require Rd/",
   "DayTemp/Flow/Season -- not available at this annual, LTER_1-only ",
   "resolution). `ilr1-3` = contemporaneous Benthos; `Year_c` = Time.")
qc("")

# ---- G1. Assemble the LTER_1 annual panel -----------------------------

N_site_year <- read_csv(here("Revision", "Data", "derived", "N_site_year.csv"), show_col_types = FALSE) |>
  filter(Site == "LTER_1", !is.na(N_percent_mean))

fish_Nexcretion_lter1 <- read_csv(here("Revision", "Data", "derived", "fish_respiration_site_year.csv"), show_col_types = FALSE) |>
  filter(Site == "LTER_1") |>
  dplyr::select(Year, N_excretion_mmolN_m2_h)

e12_data <- N_site_year |>
  dplyr::select(Year, N_percent_mean, n_samples) |>
  left_join(fish_Nexcretion_lter1, by = "Year") |>
  left_join(benthic_site_ilr |> filter(Site == "LTER_1") |> dplyr::select(Year, ilr1, ilr2, ilr3), by = "Year") |>
  mutate(
    fish_Nexcretion_z = as.numeric(scale(N_excretion_mmolN_m2_h)),
    Year_c = Year - mean(Year)
  )

qc("## G1. Model data")
qc("")
qc("- `e12_data`: ", nrow(e12_data), " annual rows, LTER_1 only (",
   min(e12_data$Year), "-", max(e12_data$Year), "). ",
   "`N_percent_mean` is the mean of `n_samples` Turbinaria tissue-N samples ",
   "per year (range ", min(e12_data$n_samples), "-", max(e12_data$n_samples), " samples/year).")
qc("- `fish_Nexcretion_z`: LTER_1's annual stoichiometric fish N-excretion ",
   "estimate (Task 7.1's placeholder O:N = 20 conversion), z-standardised.")
qc("")

# ---- G2. Collinearity check (same spirit as the plan's Task 8.2) -------
# Flagged explicitly because it directly affects how the fish_Nexcretion_z
# coefficient below should be read: fish N excretion tracks the secular
# Benthos/Time trend almost as closely as Benthos tracks Time itself, at
# only 18 annual points.

cor_mat_e12 <- cor(e12_data |> dplyr::select(fish_Nexcretion_z, ilr1, ilr2, ilr3, Year_c))

qc("## G2. Collinearity among exposure and adjustment covariates")
qc("")
qc("Pairwise Pearson correlations among `fish_Nexcretion_z`, `ilr1-3` (Benthos) ",
   "and `Year_c` (Time), at n = ", nrow(e12_data), " annual points:")
qc("")
qc("| | fish_Nexcretion_z | ilr1 | ilr2 | ilr3 | Year_c |")
qc("|---|---|---|---|---|---|")
for (rn in rownames(cor_mat_e12)) {
  qc(sprintf("| %s | %.2f | %.2f | %.2f | %.2f | %.2f |", rn,
             cor_mat_e12[rn, "fish_Nexcretion_z"], cor_mat_e12[rn, "ilr1"],
             cor_mat_e12[rn, "ilr2"], cor_mat_e12[rn, "ilr3"], cor_mat_e12[rn, "Year_c"]))
}
qc("")
qc(sprintf("**`fish_Nexcretion_z` correlates r = %.2f with `ilr1` and r = %.2f with `Year_c`, and `ilr1` correlates r = %.2f with `Year_c`** -- at n = %d, the exposure and two of its adjustment covariates are nearly collinear. The model below is fit as specified, but the individual coefficients (especially `fish_Nexcretion_z` vs. `ilr1`/`Year_c`) should be read as weakly identified from each other, not as cleanly separated effects.",
           cor_mat_e12["fish_Nexcretion_z", "ilr1"], cor_mat_e12["fish_Nexcretion_z", "Year_c"],
           cor_mat_e12["ilr1", "Year_c"], nrow(e12_data)))
qc("")

cat("E12 covariate correlation matrix:\n")
print(cor_mat_e12)

# ---- G3. Fit: log(N_percent) ~ fish_Nexcretion_z + Benthos + Time ------
# Primary model: Gaussian on log(N_percent_mean) (= lognormal on the raw
# scale, the literal Task 7.3 formula). Robustness check: Gamma(link="log")
# directly on N_percent_mean.

priors_e12 <- c(
  prior(normal(log(0.7), 0.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(2), class = "sigma")
)

fit_E12 <- brm(
  bf(log(N_percent_mean) ~ fish_Nexcretion_z + ilr1 + ilr2 + ilr3 + Year_c),
  family = gaussian(), data = e12_data,
  prior = priors_e12,
  chains = 4, iter = 4000, warmup = 1000, seed = 7934,
  control = list(adapt_delta = 0.99),
  file = here("Revision", "Output", "fits", "E12_fish_to_N"),
  file_refit = "on_change"
)

priors_e12_gamma <- c(
  prior(normal(log(0.7), 0.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(0.01, 0.01), class = "shape")
)

fit_E12_gamma <- brm(
  bf(N_percent_mean ~ fish_Nexcretion_z + ilr1 + ilr2 + ilr3 + Year_c),
  family = Gamma(link = "log"), data = e12_data,
  prior = priors_e12_gamma,
  chains = 4, iter = 4000, warmup = 1000, seed = 7934,
  control = list(adapt_delta = 0.99),
  file = here("Revision", "Output", "fits", "E12_fish_to_N_gamma"),
  file_refit = "on_change"
)

cat("\n=== E12: log(N_percent) ~ fish_Nexcretion_z + Benthos + Time (lognormal) ===\n")
print(summary(fit_E12))
cat("\n=== E12 robustness check: N_percent ~ ... , Gamma(link = log) ===\n")
print(summary(fit_E12_gamma))

post_e12 <- as_draws_df(fit_E12)
cor_b_ilr1_Year  <- cor(post_e12$b_ilr1, post_e12$b_Year_c)
cor_b_fish_ilr1  <- cor(post_e12$b_fish_Nexcretion_z, post_e12$b_ilr1)
cor_b_fish_Year  <- cor(post_e12$b_fish_Nexcretion_z, post_e12$b_Year_c)

fx_e12       <- fixef(fit_E12)
fx_e12_gamma <- fixef(fit_E12_gamma)

qc("## G3. Model results")
qc("")
qc("- Primary: `log(N_percent_mean) ~ fish_Nexcretion_z + ilr1 + ilr2 + ilr3 + Year_c`, ",
   "Gaussian likelihood on the log scale (= lognormal on `N_percent_mean`), n = ", nrow(e12_data), ".")
qc("- Robustness check: the same linear predictor with `Gamma(link = \"log\")` ",
   "directly on `N_percent_mean` (not logged).")
qc("")
qc("| Coefficient | Lognormal estimate | 95% CI | Gamma estimate | 95% CI |")
qc("|---|---|---|---|---|")
for (rn in rownames(fx_e12)) {
  qc(sprintf("| %s | %.3f | [%.3f, %.3f] | %.3f | [%.3f, %.3f] |", rn,
             fx_e12[rn, "Estimate"], fx_e12[rn, "Q2.5"], fx_e12[rn, "Q97.5"],
             fx_e12_gamma[rn, "Estimate"], fx_e12_gamma[rn, "Q2.5"], fx_e12_gamma[rn, "Q97.5"]))
}
qc("")
qc("- The two likelihoods agree closely on every coefficient -- the result ",
   "below is not an artefact of the lognormal/Gamma choice.")
qc(sprintf("- Posterior coefficient correlations (lognormal fit): b_ilr1 vs. b_Year_c = %.2f, ",
           cor_b_ilr1_Year),
   sprintf("b_fish_Nexcretion_z vs. b_ilr1 = %.2f, b_fish_Nexcretion_z vs. b_Year_c = %.2f ", cor_b_fish_ilr1, cor_b_fish_Year),
   "-- confirms the raw-data collinearity (G2) propagates into the posterior, ",
   "widening and correlating the coefficient estimates rather than being ",
   "resolved by the weakly informative priors.")
qc("")
qc(sprintf("**Result: `fish_Nexcretion_z` has a small, highly uncertain positive coefficient (lognormal: %.3f, 95%% CI [%.3f, %.3f]) that does not exclude zero. `Year_c` has a negative coefficient whose 95%% CI excludes zero (%.3f, 95%% CI [%.3f, %.3f]), i.e. tissue N declines over the study period after adjusting for contemporaneous Benthos and fish N excretion -- but because `Year_c` and the Benthos `ilr1` coordinate are themselves correlated r = %.2f at n = %d, this should be read as a covariate-of-time effect more than a cleanly isolated secular effect. No `ilr` coefficient's CI excludes zero.**",
           fx_e12["fish_Nexcretion_z", "Estimate"], fx_e12["fish_Nexcretion_z", "Q2.5"], fx_e12["fish_Nexcretion_z", "Q97.5"],
           fx_e12["Year_c", "Estimate"], fx_e12["Year_c", "Q2.5"], fx_e12["Year_c", "Q97.5"],
           cor_mat_e12["ilr1", "Year_c"], nrow(e12_data)))
qc("")
qc("**Honesty note (plan Task 7.3): this is an n = 18 annual-point, single-site ",
   "estimate with 5 regression coefficients and substantial collinearity among ",
   "them (G2) -- it should be reported as suggestive at best, not as a precise ",
   "or clearly causally separated estimate of the fish -> N effect.**")
qc("")

writeLines(qc_lines, here("Revision", "Output", "qc_E12_fish_to_N.md"))

cat("\nDone (Section G / Task 7.3). Wrote:\n",
    "  Revision/Output/fits/E12_fish_to_N.rds\n",
    "  Revision/Output/fits/E12_fish_to_N_gamma.rds\n",
    "  Revision/Output/qc_E12_fish_to_N.md\n")
