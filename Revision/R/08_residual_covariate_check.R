# ---------------------------------------------------------------------------
# 08_residual_covariate_check.R
#
# Section 8 of the Revision plan (Revision/Reviewer_Response_Plan.md):
# remove the residual covariate (Comment 5).
#
# Task 8.1 ("Delete yearresid everywhere... include Year_c directly"): there
# is nothing to delete FROM the Revision/R/*.R scripts -- `yearresid`
# (`resid(lm(Year ~ Max_temp, ...))`) only ever existed in the original
# MCR_Productivity_Analysis.qmd. Every Revision script that needed a Time
# confounder (04_benthic_gompertz.R, 06_metabolism.R, 07_estimands.R's E11
# and E12 models) already enters calendar year directly as `Year_c`
# (centred year), per Freckleton (2002) -- confirmed below by grepping the
# Revision scripts for the string "yearresid" (zero matches). This script's
# job is therefore Task 8.2 only: the collinearity diagnostic between the
# Time proxy and the heat-stress covariate it was meant to guard against
# conflating.
#
# Task 8.2 ("Report the correlation between Year and DHW_max and the
# posterior correlation of their coefficients"): computed at two levels --
# (a) the raw Year-vs-DHW_max correlation in disturbance_site_year.csv
# (n = 126 site-years, all 6 sites), and (b) the POSTERIOR coefficient
# correlation between DHW_max and Year_c in the one already-fitted model
# that carries both as covariates: the additive compositional Gompertz
# model (04_benthic_gompertz.R, benthic_gompertz_additive.rds). No new
# model is fit here; this reuses that cached fit.
#
# OUTPUTS:
#   Revision/Output/collinearity_year_dhw.csv
#   Revision/Output/qc_residual_covariate.md
# ---------------------------------------------------------------------------

source(here::here("Revision", "R", "00_packages.R"))
if (!requireNamespace("brms", quietly = TRUE)) install.packages("brms")
library(brms)

qc_lines <- character(0)
qc <- function(...) qc_lines <<- c(qc_lines, paste0(...))
qc("# QC summary: Task 8.1-8.2, removing the residual covariate (Revision Section 8)")
qc("")
qc("Generated: ", as.character(Sys.time()))
qc("")

# =============================================================================
# A. Task 8.1: confirm yearresid is absent from every Revision script
# =============================================================================

# Exclude this script itself: its own comments necessarily discuss
# "yearresid" (the thing being audited for), which would otherwise register
# as a false-positive hit against itself.
r_scripts <- setdiff(
  list.files(here("Revision", "R"), pattern = "\\.R$", full.names = TRUE),
  here("Revision", "R", "08_residual_covariate_check.R")
)
yearresid_hits <- unlist(lapply(r_scripts, function(f) {
  lines <- readLines(f, warn = FALSE)
  hit_lines <- grep("yearresid", lines, ignore.case = TRUE)
  if (length(hit_lines) > 0) paste0(basename(f), ":", hit_lines) else character(0)
}))

qc("## A. Task 8.1: `yearresid` audit")
qc("")
qc("- Searched every script in `Revision/R/` (", length(r_scripts), " files) for ",
   "the string `yearresid` (case-insensitive): ", length(yearresid_hits), " matches.")
qc("- `yearresid` (`resid(lm(Year ~ Max_temp, data = Year_Averages))`, ",
   "`MCR_Productivity_Analysis.qmd` lines 991/1431) exists only in the ",
   "ORIGINAL, untouched qmd -- never ported into the Revision pipeline. ",
   "Every Revision model that needed a Time confounder already enters ",
   "`Year_c` (centred calendar year) directly: `04_benthic_gompertz.R` ",
   "(Task 5.1's additive/interaction Gompertz fits), `06_metabolism.R` ",
   "(the hourly PI model), and `07_estimands.R` (E11's fish-response ",
   "models and E12's fish-to-N model).")
qc("")
stopifnot(length(yearresid_hits) == 0)

cat("yearresid matches in Revision/R/*.R:", length(yearresid_hits), "(expect 0)\n")

# =============================================================================
# B. Task 8.2: raw Year vs. DHW_max collinearity (all 6 sites, 126 site-years)
# =============================================================================

disturbance_site_year <- read_csv(here("Revision", "Data", "derived", "disturbance_site_year.csv"), show_col_types = FALSE)
cor_year_dhw_raw <- cor(disturbance_site_year$Year, disturbance_site_year$DHW_max)
cor_test_raw <- cor.test(disturbance_site_year$Year, disturbance_site_year$DHW_max)

qc("## B. Task 8.2: raw Year vs. DHW_max correlation")
qc("")
qc(sprintf("- `disturbance_site_year.csv` (n = %d site-years, 6 sites): Pearson r = %.3f (95%% CI %.3f to %.3f), p = %.4f.",
           nrow(disturbance_site_year), cor_test_raw$estimate, cor_test_raw$conf.int[1], cor_test_raw$conf.int[2], cor_test_raw$p.value))
qc("- Moderate, not severe, at the raw-data level: DHW_max is concentrated ",
   "in a handful of recent heat-stress years (notably 2019 and after), not ",
   "spread smoothly across the full 2005-2025 span, so the correlation with ",
   "linear Year is well below the near-collinearity seen for the E12 fish-",
   "to-N model's exposure/adjustment set (Task 7.3, r up to 0.91).")
qc("")

cat(sprintf("Raw correlation, Year vs DHW_max (n = %d): %.3f\n", nrow(disturbance_site_year), cor_year_dhw_raw))

# =============================================================================
# C. Task 8.2: posterior coefficient correlation, DHW_max vs. Year_c
#    (benthic_gompertz_additive.rds -- the one fitted model with both)
# =============================================================================

fit_path <- here("Revision", "Output", "fits", "benthic_gompertz_additive.rds")
if (!file.exists(fit_path)) stop("benthic_gompertz_additive.rds not found -- run 04_benthic_gompertz.R first")
fit_benth_additive <- readRDS(fit_path)

post <- as_draws_df(fit_benth_additive)
dpars <- c("muCoral", "muAlgae", "muCCA")

coef_cor <- tibble(
  dpar = dpars,
  cor_DHW_Year = vapply(dpars, function(d) {
    cor(post[[paste0("b_", d, "_DHW_max")]], post[[paste0("b_", d, "_Year_c")]])
  }, numeric(1))
)

fx <- fixef(fit_benth_additive)
fx_rows <- fx[grepl("DHW_max|Year_c", rownames(fx)), , drop = FALSE]

write_csv(
  bind_rows(
    tibble(check = "raw_data", variable = "Year_vs_DHW_max", n = nrow(disturbance_site_year),
           cor = cor_year_dhw_raw, p_value = cor_test_raw$p.value),
    tibble(check = "posterior_coef", variable = paste0(coef_cor$dpar, "_DHW_max_vs_Year_c"),
           n = NA_integer_, cor = coef_cor$cor_DHW_Year, p_value = NA_real_)
  ),
  here("Revision", "Output", "collinearity_year_dhw.csv")
)

qc("## C. Task 8.2: posterior coefficient correlation (benthic Gompertz model)")
qc("")
qc("- Uses the already-fitted `benthic_gompertz_additive.rds` ",
   "(`04_benthic_gompertz.R`, Task 5.1), the one model in the Revision ",
   "pipeline that carries BOTH `DHW_max` and `Year_c` as covariates. No ",
   "new model fit for this diagnostic.")
qc("")
qc("| Multinomial category | Coefficient | Estimate | 95% CI |")
qc("|---|---|---|---|")
for (rn in rownames(fx_rows)) {
  qc(sprintf("| %s | %s | %.3f | [%.3f, %.3f] |",
             sub("_(DHW_max|Year_c)$", "", rn), sub("^.*_(DHW_max|Year_c)$", "\\1", rn),
             fx_rows[rn, "Estimate"], fx_rows[rn, "Q2.5"], fx_rows[rn, "Q97.5"]))
}
qc("")
qc("| Multinomial category | Posterior cor(b_DHW_max, b_Year_c) |")
qc("|---|---|")
for (i in seq_len(nrow(coef_cor))) {
  qc(sprintf("| %s | %.3f |", coef_cor$dpar[i], coef_cor$cor_DHW_Year[i]))
}
qc("")
qc(sprintf("- Posterior coefficient correlations are modest (%.2f to %.2f in magnitude) and NEGATIVE -- the opposite sign from the raw-data correlation (+%.3f). This is unsurprising given there are 8 covariates and two random-intercept levels in this model sharing the available signal; it does not indicate the severe, near-unidentifiable collinearity seen in the E12 fish-to-N model (Task 7.3), where the exposure and two adjustment covariates shared r > 0.8 at n = 18.",
           min(abs(coef_cor$cor_DHW_Year)), max(abs(coef_cor$cor_DHW_Year)), cor_year_dhw_raw))
qc("- `DHW_max`'s 95% CI excludes zero only for CCA (muCCA_DHW_max, positive); `Year_c`'s 95% CI excludes zero for muAlgae (positive) and muCCA (negative). Coral's response to both DHW_max and Year_c has a 95% CI that includes zero in this additive model.")
qc("")
qc("**Conclusion (Task 8.2): Year and DHW_max are not so collinear in this ",
   "dataset that the heat-stress effect is unidentifiable, but DHW_max IS ",
   "concentrated in specific event years (2019 especially) rather than ",
   "varying smoothly -- so, as the plan anticipates, any heat-stress effect ",
   "is identified mainly from those event years, which should be stated as ",
   "an honest limitation rather than implied precision.**")
qc("")

writeLines(qc_lines, here("Revision", "Output", "qc_residual_covariate.md"))

cat("\nPosterior coefficient correlations (b_DHW_max, b_Year_c):\n")
print(coef_cor)
cat("\nDone. Wrote:\n",
    "  Revision/Output/collinearity_year_dhw.csv\n",
    "  Revision/Output/qc_residual_covariate.md\n")
