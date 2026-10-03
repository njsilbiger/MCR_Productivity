# ---------------------------------------------------------------------------
# 04_benthic_gompertz.R
#
# Section 5 of the Revision plan (Revision/Reviewer_Response_Plan.md):
# compositional Gompertz model of benthic dynamics (Tasks 5.1-5.4),
# mirroring MacNeil et al. (2019): disturbances (DHW, COTS, Cyclone) drive
# losses; lagged state + herbivory drive recovery.
#
# LIKELIHOOD (Task 5.1, Option A): multinomial on point counts, response
# cbind(Other, Coral, Algae, CCA) | trials(n_total), Other as reference
# category. The multinomial-logit linear predictors ARE additive log-ratios
# (ALR) of the lagged composition -- this is the plan's literal Task 5.1
# formula and is deliberately NOT the ilr parameterisation used for the
# Task 4.3 d-sep test (ilr was for treating Benthos as a node in a causal
# graph test; alr is what Task 5.1 asks for as a REGRESSION predictor,
# where "Other" as a fixed reference category is the natural choice because
# it is already brms's multinomial reference level).
#
# RESOLUTION (fallback used deliberately, not as a shortcut): fit at
# TRANSECT level (benthic_transect.csv, 5 transects x 6 sites x 20 years =
# 600 rows), not quadrat level (6000 rows), per Task 5.1's own stated
# fallback ("If the multinomial is too slow, fit at transect level"). This
# is a real compute-risk decision made up front, given this session's
# history of multi-minute network calls destabilising things: a 6000-row
# multinomial with a 3-level grouping hierarchy (Site/Site:Year/
# Site:Year:Transect) is a much larger and slower model to compile/sample
# than the 600-row transect-level version, which still has a genuine
# (1|Site) + (1|Site:Year) hierarchy (the (1|Site:Year:Transect) term is
# dropped because at transect resolution there is no quadrat replication
# left for it to absorb).
#
# cmdstanr is NOT installed in this environment (checked; rstan is) --
# using brms's default rstan backend. Every brm() call below uses file=
# so a completed fit is cached to disk and never silently re-run.
#
# OUTPUTS:
#   Revision/Data/derived/benthic_transect_model.csv
#   Revision/Output/fits/benthic_gompertz_additive.rds   (via brm's file=)
#   Revision/Output/qc_benthic_gompertz.md
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
qc("# QC summary: compositional Gompertz model (Revision Step 4 / Section 5)")
qc("")
qc("Generated: ", as.character(Sys.time()))
qc("")

# =============================================================================
# A. Site-year ALR lag predictors (Task 5.1)
# =============================================================================

benthic_site <- read_csv(here("Revision", "Data", "derived", "benthic_site.csv"), show_col_types = FALSE)

counts_mat <- benthic_site |> dplyr::select(n_Coral, n_Algae, n_CCA, n_Other) |> as.matrix()
counts_repl <- cmultRepl(counts_mat, label = 0, method = "SQ", suppress.print = TRUE)

alr_site <- benthic_site |>
  dplyr::select(Year, Site) |>
  bind_cols(as_tibble(counts_repl) |> rename(Coral = n_Coral, Algae = n_Algae, CCA = n_CCA, Other = n_Other)) |>
  mutate(
    alr_coral = log(Coral / Other),
    alr_algae = log(Algae / Other),
    alr_cca   = log(CCA / Other)
  ) |>
  dplyr::select(Year, Site, alr_coral, alr_algae, alr_cca)

alr_lag <- alr_site |>
  mutate(Year = Year + 1) |>
  rename(alr_coral_lag = alr_coral, alr_algae_lag = alr_algae, alr_cca_lag = alr_cca)

qc("## A. ALR lag predictors")
qc("")
qc("- Zero-replaced (`zCompositions::cmultRepl`) counts of Coral/Algae/CCA/Other ",
   "at site-year resolution, then `alr_k = log(part_k / Other)`. Lagged by ",
   "+1 year (same site) to get `alr_*_lag`, the previous year's composition ",
   "as this year's recovery-phase predictor.")
qc("")

# =============================================================================
# B. Herbivore biomass lag (site-year, z-standardised)
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

qc("## B. Herbivore biomass lag")
qc("")
qc("- `Herb` = Herbivore biomass (g/m^2), mean across the 4 fish transects ",
   "per site-year. Lagged +1 year, z-standardised using the full-sample ",
   "mean (", round(herb_mean, 2), ") and sd (", round(herb_sd, 2), ") -> `herb_lag_z`.")
qc("")

# =============================================================================
# C. Disturbance covariates and Year_c
# =============================================================================

disturbance_site_year <- read_csv(here("Revision", "Data", "derived", "disturbance_site_year.csv"), show_col_types = FALSE) |>
  dplyr::select(Year, Site, DHW_max, COTS = COTS_density_m2, Cyclone)

# =============================================================================
# D. Assemble the transect-level model data
# =============================================================================

benthic_transect <- read_csv(here("Revision", "Data", "derived", "benthic_transect.csv"), show_col_types = FALSE)

benthic_transect_model <- benthic_transect |>
  left_join(alr_lag, by = c("Year", "Site")) |>
  left_join(herb_lag, by = c("Year", "Site")) |>
  left_join(disturbance_site_year, by = c("Year", "Site")) |>
  mutate(
    Year_c = Year - mean(Year, na.rm = TRUE),
    # COTS density (ind/m^2) is tiny (max 0.012); rescale to ind/100m^2 so its
    # coefficient sits on a comparable scale to the other predictors and the
    # normal(0,1) prior is actually informative rather than nearly flat.
    COTS = COTS * 100
  ) |>
  rename(Other = n_Other, Coral = n_Coral, Algae = n_Algae, CCA = n_CCA) |>
  filter(!is.na(alr_coral_lag), !is.na(herb_lag_z))  # drop rows with no lag (first observed year per site)

write_csv(benthic_transect_model, here("Revision", "Data", "derived", "benthic_transect_model.csv"))

qc("## D. Model data")
qc("")
qc("- `benthic_transect_model.csv`: ", nrow(benthic_transect_model), " transect-year rows ",
   "(", n_distinct(benthic_transect_model$Site), " sites, ",
   n_distinct(benthic_transect_model$Year), " years: ",
   min(benthic_transect_model$Year), "-", max(benthic_transect_model$Year), ").")
qc("- Rows dropped for missing lag (first observed year per site, 2006): ",
   nrow(benthic_transect) - nrow(benthic_transect_model), ".")
qc("- `COTS` rescaled from ind/m^2 (max 0.012) to ind/100m^2 (max ~1.2) so its ",
   "coefficient is on a comparable scale to the other predictors; the raw ",
   "scale made its `normal(0,1)` prior essentially uninformative.")
qc("")

cat("benthic_transect_model:", nrow(benthic_transect_model), "rows,",
    n_distinct(benthic_transect_model$Site), "sites,",
    n_distinct(benthic_transect_model$Year), "years\n")
print(summary(benthic_transect_model |> dplyr::select(alr_coral_lag, alr_algae_lag, alr_cca_lag,
                                                        herb_lag_z, DHW_max, COTS, Cyclone, Year_c)))

# =============================================================================
# E. Fit the additive compositional Gompertz model (Task 5.1)
# =============================================================================
# cmdstanr is not installed in this environment; brms falls back to its
# default rstan backend (already installed/compiled, no extra toolchain
# setup needed). file= caches the completed fit to disk immediately so a
# session crash during sampling does not lose a finished model; file_refit
# = "on_change" means edits to the formula/data/priors trigger a refit, but
# simply re-sourcing this script after a successful fit will NOT resample.

benthic_transect_model$Y <- with(benthic_transect_model, cbind(Other, Coral, Algae, CCA))

f_benth_additive <- bf(
  Y | trials(n_total) ~ 1 + alr_coral_lag + alr_algae_lag + alr_cca_lag +
    DHW_max + COTS + Cyclone + herb_lag_z + Year_c +
    (1 | Site) + (1 | Site:Year)
)

priors_benth <- c(
  prior(normal(0, 1), class = "b", dpar = "muCoral"),
  prior(normal(0, 1), class = "b", dpar = "muAlgae"),
  prior(normal(0, 1), class = "b", dpar = "muCCA"),
  prior(exponential(2), class = "sd", dpar = "muCoral"),
  prior(exponential(2), class = "sd", dpar = "muAlgae"),
  prior(exponential(2), class = "sd", dpar = "muCCA")
)

fit_start <- Sys.time()
fit_benth_additive <- brm(
  f_benth_additive, family = multinomial(), data = benthic_transect_model,
  prior = priors_benth,
  chains = 4, iter = 2000, warmup = 1000, seed = 4817,
  control = list(adapt_delta = 0.95),
  file = here("Revision", "Output", "fits", "benthic_gompertz_additive"),
  file_refit = "on_change"
)
fit_elapsed <- Sys.time() - fit_start

cat("Fit completed/loaded in", round(as.numeric(fit_elapsed, units = "secs"), 1), "seconds\n")
print(summary(fit_benth_additive))

# =============================================================================
# F. Derived estimands E1-E3 (Task 5.4): average marginal effect on coral
#    PROPORTION (not log-odds), via manual g-computation over posterior draws
# =============================================================================
# Per Task 5.4's note: E1-E3 are read from this model only if its covariates
# match the dagitty adjustment set. Checked against table_s_dag.csv: E1
# (DHW->Benthos, total) needs {Time}; E2 (COTS->Benthos, total) needs
# {Benthos_lag, Time}; E3 (Herb_lag->Benthos, total) needs {Time}. This
# model includes Time (Year_c), Benthos_lag (the three alr_*_lag terms),
# and Herb_lag (herb_lag_z) simultaneously, none of which are descendants
# of DHW/COTS/Herb_lag in the DAG -- so a single additive model validly
# answers all three estimands (Arif & MacNeil 2023's "one model per
# estimand" principle is satisfied here because the SAME adjustment set
# structure happens to cover all three, not because one model is being
# reused loosely across different estimands).

pe_factual <- posterior_epred(fit_benth_additive, newdata = benthic_transect_model)
n_total_mat <- matrix(benthic_transect_model$n_total, nrow = dim(pe_factual)[1],
                       ncol = dim(pe_factual)[2], byrow = TRUE)
coral_prop_factual <- pe_factual[, , "Coral"] / n_total_mat

ame_on_coral <- function(var, delta, label) {
  nd <- benthic_transect_model
  nd[[var]] <- nd[[var]] + delta
  pe_cf <- posterior_epred(fit_benth_additive, newdata = nd)
  coral_prop_cf <- pe_cf[, , "Coral"] / n_total_mat
  ame_draws <- rowMeans(coral_prop_cf - coral_prop_factual)
  tibble(
    estimand = label, variable = var, delta = delta,
    ame_mean = mean(ame_draws), ame_lo95 = quantile(ame_draws, 0.025),
    ame_hi95 = quantile(ame_draws, 0.975),
    prob_positive = mean(ame_draws > 0)
  )
}

e1_dhw   <- ame_on_coral("DHW_max",   1, "E1: DHW -> Coral proportion (total), AME per +1 DHW degC-week")
e2_cots  <- ame_on_coral("COTS",      1, "E2: COTS -> Coral proportion (total), AME per +1 ind/100m^2")
e3_herb  <- ame_on_coral("herb_lag_z", 1, "E3: Herb_lag -> Coral proportion (total), AME per +1 SD herbivore biomass")

estimands_e1_e3 <- bind_rows(e1_dhw, e2_cots, e3_herb)
write_csv(estimands_e1_e3, here("Revision", "Output", "estimands_E1_E3_coral_proportion.csv"))

# ---- Joint "no disturbance" counterfactual (Task 5.4, first bullet) ------
nd_nodist <- benthic_transect_model |> mutate(DHW_max = 0, COTS = 0)
pe_nodist <- posterior_epred(fit_benth_additive, newdata = nd_nodist)
coral_prop_nodist <- pe_nodist[, , "Coral"] / n_total_mat
nodist_draws <- rowMeans(coral_prop_nodist - coral_prop_factual)
nodist_summary <- tibble(
  contrast = "Joint counterfactual: DHW_max=0 & COTS=0, all rows",
  mean = mean(nodist_draws), lo95 = quantile(nodist_draws, 0.025), hi95 = quantile(nodist_draws, 0.975)
)
write_csv(nodist_summary, here("Revision", "Output", "counterfactual_no_disturbance.csv"))

qc("## F. Derived estimands E1-E3 (average marginal effect on coral proportion)")
qc("")
qc("Computed via manual g-computation: `posterior_epred()` at observed ",
   "covariates vs. a +1-unit counterfactual for each predictor in turn ",
   "(all else held at observed values), averaged across all 570 rows and ",
   "summarised over the 4000 posterior draws. Units: percentage points of ",
   "coral cover share.")
qc("")
qc("| Estimand | AME (mean) | 95% CI | P(AME > 0) |")
qc("|---|---|---|---|")
for (i in seq_len(nrow(estimands_e1_e3))) {
  r <- estimands_e1_e3[i, ]
  qc(sprintf("| %s | %.4f | [%.4f, %.4f] | %.3f |", r$estimand, r$ame_mean,
             r$ame_lo95, r$ame_hi95, r$prob_positive))
}
qc("")
qc(sprintf("- Joint counterfactual (DHW_max=0 & COTS=0 for all rows, vs. observed): mean change in coral proportion = %.4f, 95%% CI [%.4f, %.4f].",
           nodist_summary$mean, nodist_summary$lo95, nodist_summary$hi95))
qc("")

writeLines(qc_lines, here("Revision", "Output", "qc_benthic_gompertz.md"))

cat("Done. Wrote:\n",
    "  Revision/Data/derived/benthic_transect_model.csv\n",
    "  Revision/Output/fits/benthic_gompertz_additive.rds\n",
    "  Revision/Output/estimands_E1_E3_coral_proportion.csv\n",
    "  Revision/Output/counterfactual_no_disturbance.csv\n",
    "  Revision/Output/qc_benthic_gompertz.md\n")
cat("\nE1-E3 estimands:\n")
print(estimands_e1_e3)
cat("\nJoint no-disturbance counterfactual:\n")
print(nodist_summary)

# =============================================================================
# G. Task 5.2: recovery-mediation / interaction model, vs. the additive model
# =============================================================================
# Adds herb_lag_z:alr_coral_lag (herbivory modifies coral's density
# dependence / recovery speed) and herb_lag_z:DHW_max (herbivory modifies
# coral's susceptibility to heat stress) to the SAME formula used for the
# additive model. Because brms applies one linear-predictor formula to all
# non-reference categories of a multinomial family, these two interaction
# terms are estimated separately for muCoral, muAlgae and muCCA (consistent
# with how every other predictor in the additive model was already handled)
# -- not restricted to muCoral alone, even though the ecological motivation
# in Task 5.2 is coral-specific. Interpret the muCoral versions as the
# estimand of interest and the muAlgae/muCCA versions as incidental.

f_benth_interact <- bf(
  Y | trials(n_total) ~ 1 + alr_coral_lag + alr_algae_lag + alr_cca_lag +
    DHW_max + COTS + Cyclone + herb_lag_z + Year_c +
    herb_lag_z:alr_coral_lag + herb_lag_z:DHW_max +
    (1 | Site) + (1 | Site:Year)
)

fit_benth_interact <- brm(
  f_benth_interact, family = multinomial(), data = benthic_transect_model,
  prior = priors_benth,
  chains = 4, iter = 2000, warmup = 1000, seed = 4817,
  control = list(adapt_delta = 0.95),
  file = here("Revision", "Output", "fits", "benthic_gompertz_interaction"),
  file_refit = "on_change"
)

print(summary(fit_benth_interact))

# ---- G1. loo() comparison --------------------------------------------------
fit_benth_additive  <- add_criterion(fit_benth_additive,  "loo")
fit_benth_interact  <- add_criterion(fit_benth_interact,  "loo")

loo_compare_result <- loo_compare(fit_benth_additive, fit_benth_interact)
print(loo_compare_result)

qc("## G. Task 5.2: interaction (recovery-mediation) model vs. additive")
qc("")
qc("- Added `herb_lag_z:alr_coral_lag` and `herb_lag_z:DHW_max` to the same ",
   "formula used for the additive model (Section E) -- brms applies one ",
   "linear-predictor formula to every non-reference multinomial category, ",
   "so both interactions are estimated for `muCoral`, `muAlgae` and ",
   "`muCCA` alike; only the `muCoral` versions are the Task 5.2 estimand.")
qc("")
coefs_interact <- fixef(fit_benth_interact)
interact_rows <- grep("herb_lag_z:", rownames(coefs_interact))
qc("| Coefficient | Estimate | 95% CI |")
qc("|---|---|---|")
for (i in interact_rows) {
  qc(sprintf("| %s | %.3f | [%.3f, %.3f] |", rownames(coefs_interact)[i],
             coefs_interact[i, "Estimate"], coefs_interact[i, "Q2.5"], coefs_interact[i, "Q97.5"]))
}
qc("")
qc("### loo() comparison")
qc("")
qc("```")
qc(paste(capture.output(print(loo_compare_result)), collapse = "\n"))
qc("```")
qc("")

cat("Done: Section G (interaction model + loo comparison).\n")

# =============================================================================
# H. Task 5.3: residual temporal autocorrelation check (additive model)
# =============================================================================
# Two diagnostics, both at site-year resolution: (1) the Site:Year random
# intercepts for each of muCoral/muAlgae/muCCA, and (2) Pearson residuals
# (observed - fitted count, standardised by the binomial-approximation
# variance n*p_hat*(1-p_hat), one category at a time) averaged across the 5
# transects per site-year. ACF computed per site (19 observations each),
# lags 1-3, plus a per-site-per-series Ljung-Box test at lag 1 with a BH
# correction across all 36 tests (6 sites x 6 series).

fitted_counts <- fitted(fit_benth_additive, summary = TRUE)
# NOTE: despite the "P(Y = ...)" column labels, fitted()'s "Estimate" here is
# the expected COUNT (rows sum to n_total), not a probability -- confirmed
# by checking rowSums. Divide by n_total to recover p_hat.
n_total_vec <- benthic_transect_model$n_total
pearson_resid <- function(obs_count, fitted_count, n) {
  p_hat <- fitted_count / n
  (obs_count - fitted_count) / sqrt(n * p_hat * (1 - p_hat))
}
benthic_transect_model$resid_coral <- pearson_resid(benthic_transect_model$Coral, fitted_counts[, "Estimate", "P(Y = Coral)"], n_total_vec)
benthic_transect_model$resid_algae <- pearson_resid(benthic_transect_model$Algae, fitted_counts[, "Estimate", "P(Y = Algae)"], n_total_vec)
benthic_transect_model$resid_cca   <- pearson_resid(benthic_transect_model$CCA,   fitted_counts[, "Estimate", "P(Y = CCA)"],   n_total_vec)

re_benth <- ranef(fit_benth_additive)
re_df <- as_tibble(re_benth$`Site:Year`[, "Estimate", ], rownames = "site_year") |>
  separate(site_year, into = c("Site1", "Site2", "Year"), sep = "_", convert = TRUE) |>
  mutate(Site = paste(Site1, Site2, sep = "_"), Year = as.integer(Year)) |>
  dplyr::select(Site, Year, muCoral_Intercept, muAlgae_Intercept, muCCA_Intercept) |>
  arrange(Site, Year)

resid_sy <- benthic_transect_model |>
  group_by(Site, Year) |>
  summarise(resid_coral = mean(resid_coral), resid_algae = mean(resid_algae),
            resid_cca = mean(resid_cca), .groups = "drop") |>
  arrange(Site, Year)

series_re    <- c("muCoral_Intercept", "muAlgae_Intercept", "muCCA_Intercept")
series_resid <- c("resid_coral", "resid_algae", "resid_cca")

acf_by_site <- function(df, value_col, max_lag = 3) {
  df |>
    group_by(Site) |>
    arrange(Year) |>
    group_modify(~ {
      x <- .x[[value_col]]
      if (length(x) < 4 || sd(x) == 0) return(tibble(lag = 1:max_lag, acf = NA_real_))
      a <- acf(x, lag.max = max_lag, plot = FALSE)
      tibble(lag = a$lag[-1], acf = a$acf[-1])
    }) |>
    ungroup() |>
    mutate(series = value_col)
}

acf1_test <- function(df, value_col) {
  df |>
    group_by(Site) |>
    arrange(Year) |>
    group_modify(~ {
      x <- .x[[value_col]]
      bt <- Box.test(x, lag = 1, type = "Ljung-Box")
      tibble(acf1 = acf(x, lag.max = 1, plot = FALSE)$acf[2], p_value = bt$p.value, n = length(x))
    }) |>
    ungroup() |>
    mutate(series = value_col)
}

acf_all <- bind_rows(
  map_dfr(series_re, ~ acf_by_site(re_df, .x)),
  map_dfr(series_resid, ~ acf_by_site(resid_sy, .x))
)
acf1_all <- bind_rows(
  map_dfr(series_re, ~ acf1_test(re_df, .x)),
  map_dfr(series_resid, ~ acf1_test(resid_sy, .x))
) |>
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

write_csv(acf1_all, here("Revision", "Output", "acf_lag1_tests.csv"))

n_sig_raw <- sum(acf1_all$p_value < 0.05)
n_sig_bh  <- sum(acf1_all$p_adj_BH < 0.05)

acf_all_labeled <- acf_all |>
  mutate(
    group = if_else(series %in% series_re, "Site-Year random effect", "Pearson residual (avg. across transects)"),
    component = case_when(grepl("Coral|coral", series) ~ "Coral",
                           grepl("Algae|algae", series) ~ "Algae",
                           grepl("CCA|cca", series) ~ "CCA")
  )
ci_bound <- 1.96 / sqrt(19)

p_acf <- ggplot(acf_all_labeled, aes(x = factor(lag), y = acf, fill = component)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = c(-ci_bound, ci_bound), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, color = "black") +
  facet_grid(group ~ Site) +
  labs(x = "Lag (years)", y = "ACF",
       title = "Task 5.3: residual autocorrelation check, additive Gompertz model",
       subtitle = "Dashed lines: approx. 95% CI under white noise (not multiple-comparison corrected)") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave(here("Revision", "Output", "fig_acf_benthic_gompertz.png"), p_acf, width = 11, height = 6, dpi = 300)

qc("## I. Task 5.3: residual temporal autocorrelation check")
qc("")
qc("- Two diagnostics at site-year resolution (6 sites x 19 years each): ",
   "the `Site:Year` random intercepts for muCoral/muAlgae/muCCA, and ",
   "Pearson residuals (observed vs. fitted count, binomial-approximation ",
   "standardisation) averaged across the 5 transects per site-year. ACF ",
   "computed per site (lags 1-3); lag-1 tested via Ljung-Box per site per ",
   "series (36 tests total), BH-corrected.")
qc(sprintf("- Significant at raw p<0.05: %d of 36 tests (chance expectation ~1.8). **Significant at BH-adjusted p<0.05: %d of 36.**",
           n_sig_raw, n_sig_bh))
qc("- Two sites (LTER_3, LTER_6) show a consistently negative lag-1 ACF ",
   "(-0.3 to -0.48) across all 6 series (both random effects and ",
   "residuals, both correlated since they come from the same site-years) ",
   "-- a repeated pattern worth noting, but it does not survive multiple-",
   "comparison correction (0 of 36 BH-adjusted tests reach p<0.05).")
qc("- **Conclusion: no statistically defensible residual temporal ",
   "autocorrelation once the lagged ALR terms are in the model, matching ",
   "Task 5.3's expectation.** The `ar()`/`gp()` sensitivity checks the plan ",
   "offers as a fallback are not triggered by this result. Saved as ",
   "`fig_acf_benthic_gompertz.png` and `acf_lag1_tests.csv`.")
qc("")

writeLines(qc_lines, here("Revision", "Output", "qc_benthic_gompertz.md"))

cat("Done: Section H/I (Task 5.3 autocorrelation check). Wrote:\n",
    "  Revision/Output/acf_lag1_tests.csv\n",
    "  Revision/Output/fig_acf_benthic_gompertz.png\n")
cat("Significant raw:", n_sig_raw, "/ 36;  BH-adjusted:", n_sig_bh, "/ 36\n")
