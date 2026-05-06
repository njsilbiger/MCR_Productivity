# Reviewer Notes — SEM & CMIP6 Projections
*Self-review drafted 2026-05-06. Address before Nature submission.*

---

## 🔴 MAJOR — Must fix before submission

### 1. Temporal autocorrelation in the coral sub-model

**What the data show:**
- Lag-1 autocorrelation in coral residuals: r = 0.745 (Durbin-Watson ≈ 0.51)
- Rd residuals are fine (lag-1 r = 0.10) — this is coral-specific
- With n = 20 years and lag-1 r ≈ 0.75, effective sample size ≈ 3–5, not 20
- Current SEM has no AR(1) structure anywhere

**Why it matters:**
Posterior credible intervals on all coral-pathway coefficients are almost certainly too narrow. The model is more confident than the data warrant. This affects every downstream indirect-effect calculation and the projection CIs.

**How to fix:**
Add AR(1) to the coral sub-model using the `ar()` helper in the formula (brms >= 2.17):

```r
bf_coraltemp <- bf(logcoral | trunc(ub = 3.39) ~ 1 + temperature + yearresid +
                     ar(time = Year, p = 1))
```

Note: brms handles AR(1) in multivariate models but syntax needs care — test in isolation first before incorporating into the full SEM.

**✅ TESTED 2026-05-06 — Results:**

AR(1) coefficient: φ = 0.916 (95% CI: 0.572–1.181). P(φ ≥ 1) = 28.7% — near unit-root behaviour, consistent with an irreversible ecological decline rather than fluctuation around a stable state.

| Coefficient | No AR(1) [current] | AR(1) added | CI change |
|---|---|---|---|
| Temperature → Coral | −0.545 [−0.808, −0.282] ✓ | −0.465 [−1.100, −0.110] ✓ | +87% wider |
| YearResid → Coral   | −0.664 [−0.932, −0.393] ✓ | −0.730 [−1.919, −0.101] ✓ | +237% wider |

**All six downstream indirect effects remain credible (95% CI excludes zero) under AR(1).** Medians shift ≤10%; CIs widen ~50–100%. Main scientific conclusions are robust — this is an honest-accounting fix, not an overturning one.

The upper bound on Temperature → Coral moves from −0.282 to −0.110, just clearing zero. Report this explicitly in the revision.

**Bonus finding:** φ near 1.0 suggests coral decline behaves close to a random walk with drift — reinforces the "irreversible trajectory" framing and is worth noting in the Discussion.

---

### 2. Path coefficient non-stationarity across the 2010 regime shift

**What the data show:**
```
Coral → Rd slope:     early 2006–2014 = +0.49  |  late 2015–2025 = −0.06
Coral → Fleshy slope: early 2006–2014 = −1.44  |  late 2015–2025 = −0.60
```
The Acanthaster outbreak (2010) and subsequent period dramatically changed how coral mediates ecosystem function. The Coral→Rd relationship appears near-zero (possibly reversed) in the modern post-disturbance regime.

**Why it matters:**
A single pooled coefficient is a weighted average across two qualitatively different states. The projections apply the pooled coefficient to a future that resembles the low-coral post-disturbance state — meaning metabolic responses to further coral loss are likely **overstated**. This is the most serious scientific challenge to the projection conclusions.

**How to fix (options, in order of preference):**
1. **Sensitivity analysis**: re-fit SEM on 2012–2025 only (post-disturbance stable period) and compare projection outputs. If conclusions hold, say so explicitly. If they don't, report both.
2. **Regime indicator**: add a binary disturbance covariate (0 pre-2011, 1 post-2011) to Rd and Fleshy equations to explicitly test for a shift in the coral-mediated pathways.
3. **At minimum**: add a paragraph in Limitations explicitly noting that the pooled coefficients span a regime shift and that the post-2012 period alone gives weaker coral-pathway estimates.

---

## 🟡 MODERATE — Should address

### 3. Temperature extrapolation far outside calibration range under SSP5-8.5

**What the data show:**
```
Observed MCR temperature range (2005–2025): 28.5 – 30.4 °C
SSP5-8.5 projected range (2026–2100):       28.8 – 34.2 °C
Extrapolation beyond observed max:          +3.9 °C
```
The log-linear temperature → coral relationship is applied nearly 4°C beyond its training data by end of century. At 33–34°C, bleaching-driven mass mortality is near-certain and the relationship is fundamentally nonlinear (threshold, not linear). The SSP5-8.5 end-century coral projections (~1%) are likely *optimistic* — the real system would have crossed a tipping point the linear SEM cannot represent.

**How to fix:**
- Low effort: add a sentence to Results (projection paragraph) explicitly flagging that SSP5-8.5 projections beyond ~2060 are applied outside the SEM's training range by 2–4°C and should be interpreted as indicating a fully phase-shifted state rather than quantitative predictions.
- Better: cite bleaching threshold literature (e.g., Hoegh-Guldberg 1999; Hughes et al. 2017) to ground the qualitative interpretation.

---

### 4. `yearresid` is an opaque catch-all predictor

**What the data show:**
- R² of Year ~ Max_temp = 0.187, so yearresid retains **81% of Year's variance**
- Correlation between yearresid and Year itself: r = 0.902
- yearresid ≈ "calendar year with temperature effects removed" — it is nearly a secular time trend

**Why it matters:**
The yearresid coefficient (β = −0.67 on coral) absorbs the Acanthaster outbreak, recovery dynamics, changes in fishing pressure, ocean acidification, any observer effects, and any other secular process — simultaneously. It cannot be interpreted causally in any specific way. Yet Results and Discussion treat it as representing a coherent "secular decline" mechanism distinct from warming.

Additionally, holding yearresid constant at its 2025 value in the projection implicitly assumes all the non-thermal pressures operating during 2005–2025 (including an extraordinary Crown-of-Thorns outbreak) continue indefinitely. This is an unstated and strong assumption.

**How to fix:**
- In Methods, be explicit: "The yearresid coefficient represents an unresolved bundle of secular non-thermal pressures operating 2005–2025 and cannot be attributed to any single mechanism."
- In Limitations, note that freezing yearresid at 2025 assumes the intensity of these non-thermal pressures remains constant at their end-of-record level.
- Consider whether replacing yearresid with a named disturbance indicator (e.g., a Acanthaster index or a post-2010 indicator) would make the causal interpretation cleaner, even if R² drops slightly.

---

### 5. CMIP6 ensemble is a small convenience-selected subset

**What the issue is:**
Five models were selected solely because they were available on the Copernicus CDS for this variable. This is not a performance- or independence-based selection. Some of these models share atmospheric/ocean parameterisation schemes (e.g., ACCESS-CM2 and UKESM1-0-LL both use the UKESM atmospheric component), which inflates apparent inter-model agreement. The full CMIP6 ensemble has 30+ models with `tos` data available via ESGF; the p10–p90 ribbon is narrower than a true ensemble spread.

**How to fix (options):**
1. **Better**: Download 3–5 additional models from ESGF (e.g., GFDL-ESM4, CanESM5, CNRM-CM6-1) to expand to ~8–10 models and widen the ensemble spread.
2. **Acceptable**: Add a sentence in Methods noting this is a CDS-availability subset and that the inter-model spread shown likely **underestimates** the full CMIP6 range. Cite Tokarska et al. (2020, *Sci. Adv.*) or Hausfather et al. (2022) for quantification of CMIP6 spread relative to subsets.

---

## 🟢 MINOR — Nice to have

### 6. Combined uncertainty supplementary figure

The decision to show climate uncertainty (temperature p10/p90) and SEM uncertainty (ecological 95% CI) separately is scientifically defensible and arguably makes both sources more interpretable. However, a reader cannot tell from the ecological panels alone whether the ribbon is narrow because the system is predictable or because climate model spread was dropped.

**Suggested addition:** A supplementary figure showing Coral Cover projections under all 5 models × 3,000 SEM draws (15,000 combinations) for one scenario (e.g., SSP5-8.5) to justify the separation choice and show that the combined uncertainty is dominated by SEM coefficient uncertainty, not inter-model spread.

---

## ✅ What is defensible and should be retained

- Bayesian multivariate SEM with explicit DAG, convergence diagnostics (Rhat < 1.01, ESS > 1,000), and brms missing-data imputation — more rigorous than most time-series ecology papers
- Effect decomposition into direct vs. indirect paths — methodologically sound; the dominance of indirect effects is a genuine publishable contribution
- 20-year time series — genuinely long for this type of integrated analysis
- Delta-anchoring approach for projections — sensible solution to level-bias; correctly motivated in Methods
- Holding yearresid constant rather than tapering to zero — more defensible projection choice
- Removal of 3 non-credible paths (full → reduced model) with stability check — good practice
- Uncertainty source separation in Figure 7 — novel and interpretable
- Bayesian R² reported per response — transparent about which variables are well-explained (coral 0.71) vs. poorly (Pmax 0.19)

---

## Summary priority table

| # | Issue | Priority | Effort |
|---|-------|----------|--------|
| 1 | Temporal autocorrelation (AR1) in coral sub-model | 🔴 Must fix | Medium — re-run SEM |
| 2 | Coefficient non-stationarity / 2010 regime shift | 🔴 Must address | Medium — sensitivity analysis |
| 3 | SSP5-8.5 extrapolation beyond training range | 🟡 Should address | Low — text + citation |
| 4 | `yearresid` interpretability and projection assumption | 🟡 Should address | Low — text clarification |
| 5 | CMIP6 model selection justification | 🟡 Should address | Low–medium |
| 6 | Combined uncertainty supplementary figure | 🟢 Nice to have | Medium |
