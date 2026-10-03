# QC summary: Task 8.1-8.2, removing the residual covariate (Revision Section 8)

Generated: 2026-10-02 15:41:44.98259

## A. Task 8.1: `yearresid` audit

- Searched every script in `Revision/R/` (10 files) for the string `yearresid` (case-insensitive): 0 matches.
- `yearresid` (`resid(lm(Year ~ Max_temp, data = Year_Averages))`, `MCR_Productivity_Analysis.qmd` lines 991/1431) exists only in the ORIGINAL, untouched qmd -- never ported into the Revision pipeline. Every Revision model that needed a Time confounder already enters `Year_c` (centred calendar year) directly: `04_benthic_gompertz.R` (Task 5.1's additive/interaction Gompertz fits), `06_metabolism.R` (the hourly PI model), and `07_estimands.R` (E11's fish-response models and E12's fish-to-N model).

## B. Task 8.2: raw Year vs. DHW_max correlation

- `disturbance_site_year.csv` (n = 126 site-years, 6 sites): Pearson r = 0.374 (95% CI 0.213 to 0.515), p = 0.0000.
- Moderate, not severe, at the raw-data level: DHW_max is concentrated in a handful of recent heat-stress years (notably 2019 and after), not spread smoothly across the full 2005-2025 span, so the correlation with linear Year is well below the near-collinearity seen for the E12 fish-to-N model's exposure/adjustment set (Task 7.3, r up to 0.91).

## C. Task 8.2: posterior coefficient correlation (benthic Gompertz model)

- Uses the already-fitted `benthic_gompertz_additive.rds` (`04_benthic_gompertz.R`, Task 5.1), the one model in the Revision pipeline that carries BOTH `DHW_max` and `Year_c` as covariates. No new model fit for this diagnostic.

| Multinomial category | Coefficient | Estimate | 95% CI |
|---|---|---|---|
| muCoral | DHW_max | 0.025 | [-0.082, 0.130] |
| muCoral | Year_c | 0.017 | [-0.009, 0.043] |
| muAlgae | DHW_max | 0.011 | [-0.077, 0.101] |
| muAlgae | Year_c | 0.025 | [0.001, 0.048] |
| muCCA | DHW_max | 0.289 | [0.000, 0.571] |
| muCCA | Year_c | -0.158 | [-0.229, -0.089] |

| Multinomial category | Posterior cor(b_DHW_max, b_Year_c) |
|---|---|
| muCoral | -0.237 |
| muAlgae | -0.237 |
| muCCA | -0.245 |

- Posterior coefficient correlations are modest (0.24 to 0.25 in magnitude) and NEGATIVE -- the opposite sign from the raw-data correlation (+0.374). This is unsurprising given there are 8 covariates and two random-intercept levels in this model sharing the available signal; it does not indicate the severe, near-unidentifiable collinearity seen in the E12 fish-to-N model (Task 7.3), where the exposure and two adjustment covariates shared r > 0.8 at n = 18.
- `DHW_max`'s 95% CI excludes zero only for CCA (muCCA_DHW_max, positive); `Year_c`'s 95% CI excludes zero for muAlgae (positive) and muCCA (negative). Coral's response to both DHW_max and Year_c has a 95% CI that includes zero in this additive model.

**Conclusion (Task 8.2): Year and DHW_max are not so collinear in this dataset that the heat-stress effect is unidentifiable, but DHW_max IS concentrated in specific event years (2019 especially) rather than varying smoothly -- so, as the plan anticipates, any heat-stress effect is identified mainly from those event years, which should be stated as an honest limitation rather than implied precision.**

