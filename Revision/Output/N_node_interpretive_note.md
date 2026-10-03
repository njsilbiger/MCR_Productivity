# Interpretive note: the nitrogen (N) node in the causal DAG

Generated: 2026-10-02

This note summarises the exploratory evidence bearing on interpretation of Turbinaria ornata tissue %N as a nitrogen-availability proxy in the MCR reef causal DAG. All analyses use the raw MCR LTER benthic cover and macroalgal CHN files, with benthic cover normalised to each quadrat's own observed point total and tissue %N restricted to Turbinaria ornata (Sargassum excluded; see `qc_nutrient.md` for rationale).

---

## 1. Within-site coral–N association is confounded by a shared temporal trend

A random-intercept mixed model (`lmer`, Site as random intercept) with a within/between decomposition of tissue %N shows a strong within-site positive association between coral cover and tissue %N across Backreef and Fringing habitats (β_within ≈ 34, 95% CI [21, 47], p < 1e-6). However, adding a linear Year term collapses this to non-significance (β_within ≈ 9, p = 0.29), while Year itself is strongly significant (β_Year ≈ −0.87/year, p < 1e-5). Both coral cover and tissue %N are declining over 2007–2024, and the apparent within-site relationship is largely explained by shared temporal trends rather than a direct link.

Residuals from the Year-free model are substantially autocorrelated within site-habitat series (lag-1 ACF typically 0.4–0.8), meaning all reported p-values and confidence intervals from these models are anti-conservative regardless of the Year issue.

## 2. Measured disturbance covariates do not explain the shared trend

Replacing Year with the available disturbance covariates (satellite DHW_max, COTS density, Cyclone Oli indicator) does not improve model fit over the base model (LRT: χ² = 4.5, df = 3, p = 0.21). Adding disturbance on top of Year contributes nothing further (χ² = 1.9, df = 3, p = 0.59). Year captures the shared coral/N decline far better than any available measured disturbance variable. Whatever secular process drives the co-decline is not well represented by acute heat-stress pulses, annual COTS density, or the single cyclone event.

## 3. Growth dilution partly explains the tissue %N decline

**Turbinaria-specific test:** Turbinaria cover is stable or increasing at most sites while tissue %N declines. The product (Turbinaria cover × %N, a crude standing-N index) is flat or increasing at all site-habitat combinations with enough data. This is consistent with growth dilution: more Turbinaria tissue absorbing the same or increasing total N, lowering per-unit concentration.

**Total-algae test:** The full Algae functional group (turf + fleshy macroalgae) cover is significantly increasing at 10 of 12 site-habitat combinations (slopes +0.3 to +2.1 %/year). However, the total-algae × %N standing-N product is mostly flat or weakly declining, with significant declines at 3 of 12 site-habitats (LTER_4 Backreef and Fringing, LTER_6 Backreef). This means the growth-dilution explanation is sufficient for Turbinaria alone but not clearly sufficient for the whole algal community — a moderate real decline in N availability cannot be ruled out.

## 4. Dissolved N+N at LTER_1 Backreef confirms a supply-side decline

The `WaterColumnN.csv` file contains dissolved N+N (nitrite + nitrate) from LTER_1 Backreef only, 2005–2018. Using only the consistent biannual cruise samples (n = 33 samples across 14 years, 2–4 per year):

- Dissolved N+N declines significantly (slope = −0.033 µmol/L per year, R² = 0.46, p = 0.007).
- The decline correlates with tissue %N (r = 0.72, p = 0.009 on 12 overlapping annual means).
- Dissolved inorganic nitrogen in the water column is not subject to algal biomass dilution. Its decline supports a genuine reduction in N supply or recycling at LTER_1, not purely a growth-dilution artifact of Turbinaria physiology.

## 5. Coral decline precedes the nitrogen decline at LTER_1

At LTER_1 Backreef, coral cover fell from ~31% to ~5% mainly between 2006 and 2015. The sharpest dissolved N+N decline occurred later, 2014–2018. Cross-correlation between annual coral cover and biannual-cruise dissolved N+N (n = 12 years) peaks at lag 0 (r = 0.75) and remains elevated at lag −1 and −2 (coral leading N by 1–2 years; r ≈ 0.55–0.60), but drops off when N leads coral (lag +1: r = 0.48, lag +2: r = 0.20). This asymmetric lead-lag pattern is consistent with — though not proof of — coral loss preceding and potentially driving the nitrogen decline.

## 6. Forereef coral recovery complicates a simple reef-wide recycling story

LTER_1 Forereef (10 m) coral cover crashed to near-zero after Cyclone Oli (2010), then recovered to 60–80% by 2018 before crashing again during the 2019 bleaching. The dissolved N+N decline at the backreef (2014–2018) coincided with forereef coral recovery. If coral-associated N recycling were the dominant driver of backreef water-column N, forereef recovery might have buffered or reversed the decline. Its failure to do so suggests either (a) hydrological separation between habitats (flow is typically offshore-to-lagoon across the reef crest, so forereef coral may not contribute to backreef water-column nutrients), or (b) other drivers contribute to the nitrogen trend alongside coral-associated recycling.

## 7. Lead-lag pattern does not generalise across all sites

Cross-correlation between backreef coral cover and Turbinaria tissue %N (lag ±4 years) was computed at each of the six LTER sites.

| Site | Coral slope (%/yr) | Coral p | %N slope (/yr) | %N p | Peak CCF | At lag |
|---|---|---|---|---|---|---|
| LTER_1 | −1.20 | <0.001 | −0.024 | <0.001 | 0.83 | 0 |
| LTER_2 | −1.87 | <0.001 | −0.014 | 0.005 | 0.63 | −1 |
| LTER_3 | −1.41 | <0.001 | −0.013 | 0.004 | 0.67 | −4 |
| LTER_4 | +0.52 | 0.041 | −0.013 | 0.017 | 0.06 | −4 |
| LTER_5 | −0.16 | 0.194 | −0.008 | 0.128 | 0.32 | −3 |
| LTER_6 | +0.17 | 0.095 | −0.013 | 0.002 | 0.08 | −1 |

The coral–N temporal coupling is strong only at the three sites with major coral declines (LTER_1, 2, 3). At sites where coral is stable or increasing (LTER_4, 5, 6), tissue %N still declines but shows no positive correlation with coral — and at LTER_4 and LTER_6, the lag-0 correlation is weakly to moderately negative. This means:

- The coral → N recycling pathway is plausible *at sites where coral actually declined*, but tissue %N declines at sites with stable coral too, so coral loss alone cannot be the sole driver.
- A reef-wide or oceanic process (declining nutrient supply, warming-driven changes in N cycling) and/or Turbinaria growth dilution likely contribute to the %N decline across all sites, independent of local coral trajectory.
- The LTER_1 lead-lag pattern (Section 5) may partly reflect the coincidence of the strongest coral decline at the site where dissolved N+N data happen to be available, rather than a generalisable causal mechanism.

## 8. Implications for the DAG

1. **Tissue %N is not a clean proxy for reef-scale nitrogen availability.** Its dominant temporal signal includes substantial growth dilution (at least for Turbinaria), and the dissolved N+N validation is limited to one site through 2018. The N node should be labelled as "Turbinaria tissue %N" — an imperfect, species-specific proxy — not "nitrogen availability" or "nitrogen flux."

2. **The Benthos → N direction is more consistent with the temporal evidence than N → Benthos** at LTER_1, based on the lead-lag structure and the timing of coral decline relative to the N decline. However, this is based on n = 12 annual pairs at one site and cannot distinguish the Benthos → N DAG variant from the Rd → N variant (since both predict the same temporal ordering at the metabolism site).

3. **An unmeasured time-varying confounder (the DAG's U node) remains plausible.** No available measured disturbance variable captures the shared trend. Candidates include cumulative or lagged heat stress, long-term warming trends, changes in oceanic nutrient delivery (e.g., internal wave frequency), or other reef-scale processes not represented in the current data.

4. **The spec's warning (line 102) about fitting multiple unmeasured mechanistic links with one observed proxy is directly relevant.** The data cannot separate coral-associated microbial N recycling, fish-mediated N excretion, growth dilution, and oceanic supply changes. Effects estimated through the N node should be reported as conditional on the assumption that tissue %N tracks the intended mechanism, with this limitation stated explicitly.

---

## Summary table

| Finding | Implication for DAG |
|---|---|
| Within-site coral–N association driven by shared Year trend | No direct coral → N link identifiable from this panel alone |
| DHW/COTS/Cyclone don't capture the shared trend | Unmeasured confounder (U) remains active |
| Turbinaria growth dilution explains species-specific %N decline | Tissue %N partly reflects Turbinaria biomass, not just N supply |
| Dissolved N+N declining at LTER_1 (biannual cruise) | A real supply-side decline exists, at least at one site |
| Coral decline precedes N decline (lead-lag) | Consistent with Benthos → N, but n = 12, one site |
| Forereef coral recovery didn't reverse backreef N decline | Habitat-specific story; reef-wide recycling pathway not straightforward |
| Lead-lag pattern doesn't generalise: %N declines at sites with stable coral | Coral loss alone cannot explain the %N decline; reef-wide/oceanic driver or growth dilution likely contributes |
