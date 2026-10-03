# Methods (revised) — sections completed so far

This file accumulates the revised Methods text section by section as the
Revision plan (`Reviewer_Response_Plan.md`) is executed. Sections are
numbered to match the plan, not necessarily the eventual manuscript order.
Each section below is written as manuscript-ready prose, suitable for
pasting into both the revised Methods and the point-by-point response
letter (`Response_to_Reviewers.md`, not yet started).

---

## Section 4. A structural causal model: DAG, adjustment sets, and a data-driven test

**Addresses reviewer comment 1 ("a DAG is not a SEM; parameters estimated
jointly are not causal") and comment 2 ("the DAG is incomplete: fish do
not appear as a cause of respiration or nitrogen, and the coral–algae
feedback is collapsed into a single static arrow").**

### From a simultaneous SEM to one model per causal estimand

The original analysis fit all structural paths jointly in a single
simultaneous-equations model and interpreted every resulting coefficient
as causal. We replaced this with the structural causal model (SCM)
workflow of Pearl (2009) and, for ecological applications specifically,
Arif & MacNeil (2023): a single directed acyclic graph (DAG) encodes our
causal assumptions, the backdoor criterion is applied separately to each
causal question of interest to derive its minimal sufficient adjustment
set, and a **separate statistical model is fit for each estimand**, using
only the covariates that set identifies. Under this workflow, a DAG
records assumptions and testable implications; it does not itself
produce effect estimates, and no coefficient from one estimand's model is
read as the answer to a different estimand's question.

The revised DAG (Fig. 3, replacing the path diagram with β-coefficients
on its arrows) extends the original model in three ways motivated
directly by the reviewer's comments. First, fish community biomass
(partitioned into herbivore, corallivore, and other-fish nodes) is now an
explicit cause of ecosystem respiration and of nitrogen availability, not
an unconnected sink. Second, the coral–algal–herbivore feedback — a cycle
in the original static diagram — is unrolled in time via lagged state
nodes (last year's benthic composition, herbivore biomass, and nitrogen),
which removes the cycle and lets the backdoor criterion be applied to a
genuine DAG. Third, benthic composition is represented as a **single
compositional node** rather than separate coral and algae nodes connected
by an implied (but never justified) causal arrow between them; the
simplex's closure constraint is a mathematical property of the data, not
a causal relationship, and is handled instead by the compositional
likelihoods described in Section 5. One path from the original model was
removed outright: respiration and photosynthetic capacity (Rd and Pmax)
share a direct arrow in the original SEM, but we found no mechanism by
which respiration rate causes photosynthetic capacity; their covariance
instead reflects shared causes (benthic biomass, temperature) that are
already in the model, and is captured as a residual correlation within
the photosynthesis–irradiance model of Section 6 rather than as a causal
edge — a direct illustration of the DAG-versus-SEM distinction the
reviewer raised.

One mechanism could not be settled from first principles: whether
nitrogen availability is better modelled as a downstream consequence of
respiration-linked remineralisation (`Rd → N`) or as driven directly by
standing benthic biomass (`Benthos → N`). Rather than choose on priors
alone, we carried **both DAG variants** in parallel through the
adjustment-set derivation and the data-driven test described below, and
report where the choice does and does not affect downstream estimands.

### Adjustment sets for twelve causal estimands

For each of the twelve causal questions motivating this paper (heat
stress, crown-of-thorns starfish, and herbivory's effects on benthic
composition; benthic composition, fish, and temperature's effects on
respiration and photosynthetic capacity; benthic composition's effects on
fish; and fish's effect on nitrogen), we derived the minimal sufficient
adjustment set via `dagitty`'s implementation of the generalised backdoor
criterion (Textor et al. 2016), distinguishing total and direct effects
where both are of interest (Table S-DAG). Several estimands admit more
than one valid minimal adjustment set — for example, a day's temperature,
flow, and season sit on an interchangeable causal chain, so blocking any
one is sufficient — in which case we selected the most parsimonious set
consistent with the variables available in the disaggregated data tables.
One adjustment set is smaller than earlier drafts of the model assumed:
heat stress's total effect on benthic composition requires conditioning
on a secular-time proxy alone, not also on the previous year's
composition, since the latter is not a common cause of heat stress and
current composition. The two DAG variants (`Rd → N` vs. `Benthos → N`)
return identical adjustment sets for most estimands but diverge for
those concerning photosynthetic capacity, since the `Rd → N` variant
opens an additional causal path through nitrogen that the alternative
does not — a concrete, checkable consequence of the unresolved mechanism
above that the data-driven test below was designed to help adjudicate.

### Testing the DAG against data

Following Shipley (2000) and the model-testing extension of the SCM
framework (Ankan et al. 2021), we tested both DAG variants' implied
conditional independences against the disaggregated data rather than
assuming the graph without checking it. Benthic composition was
represented by its isometric log-ratio (ilr) coordinates (Egozcue et al.
2003) and tested as a full multivariate vector via a canonical-correlation
conditional-independence test (Pillai's trace), avoiding the common
practice of reducing a compositional node to a single scalar proxy when
used inside a causal test. Independences involving ecosystem metabolism
were necessarily tested on the single-site, LTER 1-only subset of the
data using proxy respiration and photosynthesis values (not the
Section 6 model's eventual output, which did not yet exist at this stage
of the analysis), since metabolism is measured at only one of the six
sites; all other independences used the full six-site panel.

The initial test identified 25 of 440 testable missing-edge independences
as violated (Benjamini–Hochberg-corrected p < 0.05), with a clear and
ecologically interpretable pattern: the DAG's lagged state nodes
(previous year's composition, herbivore biomass, and nitrogen) were each
found to be correlated with the secular-time proxy despite having no
edge from it — an omission we recognised as a genuine structural gap
rather than a modelling artefact, since a lagged quantity remains
time-varying. Crown-of-thorns starfish density was likewise found to be
correlated with nitrogen availability despite having no path to it or
connection to secular time in the original graph, consistent with the
outbreak-year clustering documented independently in the disturbance time
series. We added the corresponding edges — secular time as a cause of
each lagged state node, and of crown-of-thorns starfish density — to both
DAG variants and re-tested. Missing-edge violations fell from 25 to 15 of
440, and a global Shipley–Fisher's C test (computed from the graph's
basis set of independences) moved from firmly rejecting both DAG variants
(C = 58.0, df = 14, p < 0.0001) to not rejecting either
(C = 18.1, df = 14, p = 0.20). We treat this global test as corroborating
rather than decisive, since several basis-set tests required conditioning
on seven to ten variables simultaneously at the single-site metabolism
panel's sample size (n = 17–18) — a near-saturated regime in which the
test has little power to distinguish a correctly specified graph from a
misspecified one. Two violations persisted after the edge additions and
are reported as open structural questions rather than patched further
without independent justification: a residual correlation between lagged
composition and lagged herbivore biomass not explained by their shared
time trend, and a residual correlation between crown-of-thorns starfish
density and nitrogen availability for which the nutrient-enrichment
hypothesis of coral-predator outbreaks (Fabricius et al. 2010) offers a
plausible but as-yet-untested mechanism. The persistent correlation
between respiration and photosynthetic capacity, by contrast, is not a
new finding: it is exactly what motivated handling that pair as a
residual correlation rather than a causal edge in the first place (above),
and its persistence after the edge additions confirms that choice rather
than undermining the revised graph.

---

## Section 5. Benthic composition dynamics: a compositional Gompertz model

**Addresses reviewer comments 2 ("the coral–algae feedback is collapsed
into a static arrow") and 4 ("the models treat 20 points as if they were
independent, and autocorrelation is addressed by fitting AR(1) to coral
alone").**

### Model structure

We replaced the original year-level, log-z-scored treatment of coral and
algal cover with a hierarchical, compositional model of benthic dynamics
at **transect resolution** (570 transect-years: 6 backreef sites × 19
years, 2007–2025, × 5 transects per site-year), following the
population-dynamic framing of MacNeil et al. (2019): disturbances drive
losses, lagged state and herbivory drive recovery. Benthic cover is
modelled as a four-part composition (Coral, Algae, CCA, Other) using a
multinomial-logit likelihood on the underlying point counts, which
respects the sum-to-100% constraint exactly and handles the exact zeros
present in some site-years without ad hoc transformation (Aitchison 1986).
Counts were zero-replaced using the multiplicative simple-replacement
method of Palarea-Albaladejo & Martín-Fernández (2015) before computing
additive log-ratio (ALR) coordinates of the previous year's composition
(Other as reference), which enter the current year's model as predictors.
A regression of this year's ALR composition on last year's ALR composition
is a discrete-time multivariate Gompertz/MAR(1) model (Ives et al. 2003):
diagonal coefficients measure within-part density dependence and
off-diagonal coefficients measure between-part competition for space. The
full model is

```
cbind(Other, Coral, Algae, CCA) | trials(n_total) ~
  1 + alr_coral_lag + alr_algae_lag + alr_cca_lag +
  DHW_max + COTS + Cyclone + herb_lag_z + Year_c +
  (1 | Site) + (1 | Site:Year)
```

fit in a Bayesian multinomial-logit framework (`brms`/Stan), with
`(1 | Site)` and `(1 | Site:Year)` absorbing between-site heterogeneity and
site-year process noise respectively. Herbivore biomass (g m⁻²) and heat
stress (DHW, °C-weeks) and crown-of-thorns starfish density enter as the
previous year's lagged, standardised covariates, consistent with the
adjustment sets derived from the causal DAG (Table S-DAG) for the
corresponding estimands (see Section 4). The model converged cleanly
(R̂ ≤ 1.01, bulk/tail effective sample size ≥ ~1,100 for every parameter
across 4,000 post-warmup draws).

### Recovery mediation

To test whether herbivory modifies coral's recovery rate or its
susceptibility to heat stress — rather than simply adding to coral cover —
we fit a second model adding `herb_lag_z : alr_coral_lag` and
`herb_lag_z : DHW_max` interaction terms, and compared the two models by
approximate leave-one-out cross-validation (`loo`, Vehtari et al. 2017).
The interaction model was not favoured (ELPD difference −26.2, SE 13.2,
in favour of the simpler additive model), and neither interaction
coefficient was distinguishable from zero on coral
(`herb_lag_z:alr_coral_lag` = 0.01, 95% CrI [−0.14, 0.15];
`herb_lag_z:DHW_max` = −0.02, 95% CrI [−0.19, 0.14]). We note that the
approximate leave-one-out diagnostic (Pareto-k̂) was poor for both models —
a known limitation of PSIS-LOO for multinomial/binomial likelihoods with
large trial counts (here, 125–250 points per transect-year) — so this
comparison is reported as corroborating rather than decisive evidence; the
individual-parameter result, which does not depend on the LOO
approximation, is the primary basis for concluding that this dataset does
not support a herbivory-modifies-recovery mechanism at this scale.

### Residual temporal dependence

We checked whether the lagged compositional structure above was
sufficient to remove residual temporal autocorrelation (reviewer comment
4), rather than assuming it a priori. We computed the autocorrelation
function (lags 1–3) of both the site-year random effects and Pearson
residuals (binomial-approximation standardisation) for each of the three
non-reference categories at each of the six sites, and tested lag-1
autocorrelation via Ljung–Box tests (36 site × category combinations),
Benjamini–Hochberg corrected for multiple comparisons. Two of 36 raw tests
reached p < 0.05 — at the level expected by chance (≈ 1.8 of 36) — and none
survived correction. We conclude that the lagged state terms adequately
capture the system's temporal dependence, obviating the AR(1)/Gaussian-
process sensitivity checks the a priori plan allowed for in case residual
autocorrelation remained.

### Causal estimands (E1–E3)

Average marginal effects of heat stress (E1), crown-of-thorns starfish
density (E2), and lagged herbivore biomass (E3) on coral cover share were
computed by posterior g-computation: for each predictor in turn, the fitted
model was evaluated at the observed covariates and at a one-unit
counterfactual shift, holding all other covariates at their observed
values, and the resulting change in coral proportion was averaged across
all 570 transect-years and summarised over the posterior. All three
estimands were null, with 95% credible intervals comfortably spanning zero
(DHW: 0.0015 percentage points, 95% CrI [−0.0153, 0.0195]; COTS: 0.0265,
95% CrI [−0.0618, 0.1328]; herbivory: 0.0037, 95% CrI [−0.0162, 0.0242]).
A joint counterfactual removing both heat stress and COTS entirely (set to
zero for every transect-year) likewise showed no detectable change in
coral share (−0.0019, 95% CrI [−0.0103, 0.0068]).

We report these as genuine null results rather than evidence of no
ecological effect. Two structural features of the data limit this model's
power to detect them: DHW and COTS are both heavily zero-inflated (DHW_max
is non-zero in only a minority of site-years, dominated by the single 2019
bleaching event; COTS density is non-zero in a small minority of
site-years — see Section 3), and coral's own lagged-composition term
(1.01, 95% CrI [0.90, 1.12]) already accounts for most of the explainable
year-to-year variance, leaving comparatively little residual variance for
the disturbance covariates to explain even if their true effect is modest
rather than zero. We note separately that coral's near-unit persistence
coefficient is statistically indistinguishable from a random walk at this
temporal resolution (no evidence of density-dependent return to a
long-run mean within the observed 19-year window) — in contrast to algae
(0.45, 95% CrI [0.23, 0.69]) and CCA (0.65, 95% CrI [0.49, 0.82]), both of
which show clear, credible density dependence. This asymmetry, together
with the one unambiguously credible disturbance-adjacent effect in the
model — algal cover's competitive suppression of coral
(`alr_algae_lag` → `muCoral` = −0.35, 95% CrI [−0.59, −0.13]) — is, in our
view, the most defensible population-dynamic signal to draw from this
analysis, and we present it alongside, rather than instead of, the null
E1–E3 results.

---

*(Sections to follow as subsequent plan steps are completed: Section 6,
ecosystem metabolism; Section 7, fish bioenergetics and its contribution
to respiration and nitrogen; the remaining per-estimand models (E4–E10);
a simulation check comparing the original simultaneous-SEM approach
against this one; and robustness/sensitivity diagnostics.)*
