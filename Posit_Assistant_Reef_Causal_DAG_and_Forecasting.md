# Reef causal analysis and ecosystem metabolism forecasting

## Instructions for Posit Assistant

Read this entire file before writing models. Do everything in a new revision document. Treat it as the analysis specification. Work through the numbered steps in order, creating runnable R scripts and a Quarto report. Begin with a data audit, a candidate DAG, and an implementation feasibility assessment. Continue with tasks that do not depend on unresolved scientific choices; flag decisions requiring the investigator's input. Do not fabricate data, fitted results, package functions, environmental covariates, or sampling metadata.

The goal is to evaluate whether coral loss reorganizes reef nutrient cycling, fish communities, and ecosystem metabolism, and then use the fitted generative model to forecast metabolic states. Use `dagitty` for causal structure and identification and `brms` for estimation where the required model is supported. A DAG does not itself estimate parameters or generate forecasts. Bayesian estimation does not remove confounding or guarantee causal identification.

Follow the structural causal model approach of Suchinta Arif and M. Aaron MacNeil, particularly their 2023 Ecological Monographs framework and 2022 coral reef regime shift application. The compositional, measurement, dynamic, and forecasting specifications below are adaptations to this study, not claims that those papers implemented this exact model. References and official software documentation are provided at the end.

## Step 1 Establish the study structure and hypotheses

Use the following study information as given:

- Approximately 20 years at six reef sites, with uneven sampling and missing site-year combinations.
- Benthic coral, macroalgae, CCA, and sand cover sum to 100 percent. Multiple transects provide replicate observations.
- Fish biomass is recorded on multiple transects. Corallivores and herbivores are focal groups, but all observed fish groups are relevant to whole-community fish respiration.
- Macroalgal tissue percent N has approximately five replicate samples per sampled site-year. It is a proxy for nitrogen availability and recycling; it is not a direct measurement of nitrogen flux.
- Ecosystem respiration, ER, and net community production or photosynthesis, NCP, are measured at one site only, on multiple days during one or two seasons depending on logistics.
- The investigator reports that previously estimated fish respiration contributes less than 1 percent of total ER. This is a study-specific finding to verify from existing calculations and their uncertainty, not a universal constraint to impose.

Represent these hypotheses explicitly:

1. Coral loss reduces non-fish ecosystem respiration through changes in coral-associated biomass and microbial processes.
2. Coral loss reduces coral-associated nitrogen recycling, potentially reducing available nitrogen and macroalgal tissue percent N.
3. Coral loss reduces corallivore biomass.
4. Increased macroalgal dominance can increase herbivore biomass through bottom-up pathways. Grazing can affect future macroalgal cover, so consider lagged feedback as an alternative.
5. Benthic composition and nutrient availability affect ecosystem primary production.
6. Fish abundance, sizes, traits, and temperature determine fish respiration, which contributes to total ER.

Do not use stepwise selection, information criteria, or significance to select causal arrows. Predictive validation may compare predeclared forecasting models; it cannot establish causal structure.

Output: study summary and a hypothesis-to-pathway table.

## Step 2 Audit data and create linked observation tables

Inspect available files and create a data dictionary containing units, sampling dates, methods, missing-value codes, identifiers, and spatial support. Determine whether benthic cover comes from point counts, quadrats, or continuous estimates; whether transects are permanent; whether fish and benthos share transects; and whether metabolism estimates have standard errors or joint uncertainty from the original calculation.

Create the following tables without multiplying observations through many-to-many joins:

| Table | Required fields where available |
|---|---|
| State index | site_id, year, optional season, state_id, previous_state_id |
| Benthos | state_id, date, transect_id, four covers or counts, total points, method |
| Fish | state_id, date, transect_id, species, guild, count, size or size class, biomass, surveyed area |
| Tissue N | state_id, date, sample_id, species, percent_N, lab batch, analytical uncertainty |
| Metabolism | state_id, date, season, ER, NCP, uncertainty and covariance if available, light, temperature, flow |
| Environment | site_id, date or year, temperature, disturbance and other actually measured drivers |
| Sampling | site_id, year, variable, sampled indicator, effort, season, reason missing if known |
| Fish traits | species, length-weight relationship, metabolic parameters, source, uncertainty |

Construct the complete six-site annual index spanning the observed record. Empty cells are latent states, not zero-valued observations. Do not add fictitious transects or days. Keep metabolism observations restricted to the actual metabolism site. Forecasts for the other five sites require explicit transportability assumptions and must be labeled accordingly.

Check duplicate keys, composition totals, negative biomass, genuine zeros, detection limits, taxonomic changes, inconsistent units, survey footprint, and season imbalance. Save a missingness plot and effort plot. Explain whether an annual state is defensible or seasonal states are required. Years are crossed with sites; transects, samples, and days are observations linked to their corresponding state.

Output: validated tables, data dictionary, audit report, and state index.

## Step 3 Define the exposure as a compositional intervention

Do not treat the four covers as independent scalar predictors. Closure creates dependence but does not prohibit a causal hypothesis about replacement of coral by particular benthic categories.

Define an intervention on the complete composition. For baseline proportions p = (c, a, k, s), a proportional redistribution after setting coral to c_new is:

```text
a_new = a * (1 - c_new) / (1 - c)
k_new = k * (1 - c_new) / (1 - c)
s_new = s * (1 - c_new) / (1 - c)
```

This is defined when c < 1. Also evaluate explicitly declared alternatives such as coral replaced mainly by macroalgae or mainly by sand. Check feasibility and positivity of each composition. Report effects as effects of these replacement scenarios, not an unspecified effect of coral while all other percentages remain fixed.

For positive proportions, one orthonormal ILR basis is:

```text
z1 = sqrt(3/4) * log(c / (a*k*s)^(1/3))
z2 = sqrt(2/3) * log(a / (k*s)^(1/2))
z3 = sqrt(1/2) * log(k / s)
```

Coordinate changes are not automatically isolated coral interventions. Transform each full replacement composition to ILR coordinates before predicting, and invert the transform for interpretable cover predictions. An ILR basis is a representation, not three distinct ecological processes with independently manipulable causes.

If point counts are available, evaluate a multinomial observation likelihood with a logistic-normal latent composition. If only continuous positive covers exist, evaluate joint log-ratio observation models. Zeros require a justified count, hurdle, detection-limit, or zero-handling model; do not insert an arbitrary pseudocount without sensitivity analysis. If using Dirichlet observations, document restrictions including zero exclusion and covariance assumptions. Do not assume all these models are directly implemented in brms.

Output: chosen compositional representation, inverse transform, zero-handling rationale, and intervention definitions.

## Step 4 Separate ecological states from measurements

Define latent states indexed by site and time: benthic composition B, fish community F, a tissue N state or calibrated nitrogen state, non-fish respiration R_other, and production P or NCP. Distinguish latent process innovations from variation among replicate observations.

For a simple continuous state, the conceptual model is:

```text
X[s,t] = f_X(causal parents, X[s,t-1], environment, site) + eta_X[s,t]
eta_X[s,t] ~ process distribution with scale sigma_process_X
y[s,t,r] ~ observation distribution centered on the appropriate X[s,t]
```

Replicate variation often includes genuine spatial or daily ecological variability, not just instrument error. Separate analytical error, transect heterogeneity, day variability, and state process innovations where data support this. A known measurement standard error can be combined with an additional observation variance, but do not count the same variance twice.

Use recurring transect effects if permanent transects exist. Metabolism days within a campaign may be correlated. Include day-level environmental covariates and campaign effects where justified. Check identifiability of every variance component by simulation; replication alone does not guarantee separation.

Macroalgal percent N alone cannot identify the absolute scale of nitrogen flux or separate coral microbial remineralization from other nitrogen sources. Unless independent calibration exists, fit latent true tissue N on its measured scale and describe microbial recycling as a hypothesized explanation. Account for species, growth dilution, seasonal conditions, and laboratory effects where available. Do not freely fit several unmeasured mechanistic links with one observed proxy.

Output: state definitions and observation equations with a variance-component table.

## Step 5 Build ecological and observation DAGs

Create a biological process DAG, a separate observation diagram, and an expanded temporal DAG. Include measured drivers and plausible unmeasured common causes even when these prevent identification. Site intercepts and temporal autocorrelation do not automatically remove all spatial or time-varying confounding.

Use the following as a starting ecological skeleton, not an accepted identified DAG. B represents the complete composition. The hypothesis about coral and algae is encoded through functions of that composition. FishCommunity contains the guild and size structure necessary for respiration; guild biomass summaries and individual fish observations belong in its measurement layer.

```r
library(dagitty)

g <- dagitty('dag {
  B [exposure]
  ER [outcome]
  U [unobserved]
  Site -> B
  Site -> NState
  Site -> FishCommunity
  Site -> OtherER
  Site -> Production
  Env -> B
  Env -> NState
  Env -> FishCommunity
  Env -> FishResp
  Env -> OtherER
  Env -> Production
  Disturbance -> B
  Disturbance -> NState
  Disturbance -> FishCommunity
  Disturbance -> OtherER
  Disturbance -> Production
  U -> B
  U -> OtherER
  U -> NState
  Bprev -> B
  Bprev -> FishCommunity
  Bprev -> OtherER
  Fprev -> FishCommunity
  Fprev -> B
  Nprev -> NState
  Rprev -> OtherER
  Pprev -> Production
  B -> NState
  B -> FishCommunity
  B -> OtherER
  B -> Production
  NState -> Production
  FishCommunity -> FishResp
  FishResp -> ER
  OtherER -> ER
  Production -> NCP
  ER -> NCP
}')

plot(g)
impliedConditionalIndependencies(g)
adjustmentSets(g, exposure = 'B', outcome = 'ER', effect = 'total')
```

The last two arrows assume Production means gross primary production and NCP = Production - ER on compatible units and integration periods. These are accounting relationships. If the observed variable is instead daytime net photosynthesis or another metabolism metric, revise this structure to its actual definition. Do not estimate a freely varying ER-to-NCP coefficient for an exact accounting identity.

Mark which states are observed only through measurements in the statistical graph; the simplified process graph above suppresses that layer. In particular, an uncalibrated nitrogen state is not an observed confounder or known mediator. The U common cause deliberately leaves important effects potentially unidentified. Never remove it merely to produce an adjustment set. Do not treat sampling indicators as ordinary confounders and adjust for them automatically.

Expand two consecutive time slices explicitly. Each state gets biologically justified persistence and lagged causal arrows. Model grazing as F[t-1] -> B[t], and bottom-up response as B[t] -> F[t] or B[t-1] -> F[t] according to timing. Consider independent common causes affecting earlier and current states. Initial states need their own joint prior or baseline model.

The skeleton omits possible fish nutrient excretion, pelagic food supply, and other pathways. Add them only with ecological justification and available information; distinguish fish respiration from fish-mediated nutrient recycling. Evaluate a small number of scientifically motivated alternative graphs.

Output: executable DAG code, figures, node dictionary, and one row per arrow giving mechanism, timing, measurement status, and source or assumption.

## Step 6 Define causal queries and identification

Create a query table before fitting:

| Query | Target |
|---|---|
| Coral replacement and ER | Total effect of a specified composition replacement on ER at the metabolism site |
| Coral replacement and tissue N | Total effect on measured-scale latent tissue N, with recycling interpretation qualified |
| Coral replacement and corallivores | Total effect on corallivore biomass across sampled sites |
| Algal replacement and herbivores | Total effect of a specified composition change on herbivore biomass |
| Composition and NCP | Total effect on correctly defined NCP at the metabolism site |
| Fish pathway | Contribution of fish respiration to ER and its change under replacement scenarios |
| Dynamic intervention | Effect of a specified intervention trajectory over a declared future horizon |

For each, identify exposure, outcome, target population, time horizon, measured and unmeasured common causes, mediators, colliders, and valid adjustment sets using dagitty. Record explicitly when none exists. Correlated random effects, priors, or a fitted joint likelihood do not restore nonparametric identification under unmeasured confounding.

Distinguish a model's local structural equation conditioned on its parents from a regression targeting a total effect. In a joint structural model, total effects require simulating downstream mediators. Do not hold fish or nitrogen mediators fixed when claiming a total composition effect.

Controlled direct effects, natural direct or indirect effects, and pathway attribution are different quantities. Do not call dagitty's direct-effect option proof that natural mediation is identified. Exposure-induced mediator-outcome confounding and longitudinal feedback may require a sequential g-formula or stronger assumptions. Report assumption-dependent contrasts when necessary.

Output: identification table and predeclared causal estimands.

## Step 7 Add fishflux respiration without double counting

Inspect the installed version, documentation, and source of `fishflux`. Record its version, parameter sources, and exact function arguments. Its metabolic functions use mass, temperature, and other traits; verify the returned units before use. Do not invent a function that accepts total guild biomass as sufficient input.

Calculate fish respiration from species and size or mass distributions, abundances, activity assumptions, temperature, and surveyed area. For draw m:

```text
R_fish[s,t,m] = sum_i(n_i[s,t,m] * r_i[m]) / surveyed_area[s,t]
```

Propagate uncertainty in counts or biomass, length-weight conversion, metabolic parameters, missing traits, activity, and temperature. Preserve shared parameter uncertainty across years and sites; do not redraw a species metabolic calibration independently for every record when it is a shared parameter. Include all surveyed groups. If only aggregate biomass remains, document missing size information and use explicit size-distribution scenarios; biomass alone is insufficient under nonlinear mass scaling.

Convert outputs to the same currency, sign, area, and time basis as ER. If an output is grams C per day and ER is mmol O2 per square meter per day, carbon moles and oxygen demand require a declared respiratory quotient RQ, where mol O2 = mol C / RQ, plus the area normalization. Verify whether the selected output already represents oxygen or energy, and document each conversion. Address fish residence time, movement, nocturnal activity, and mismatch between survey area and the metabolism footprint with sensitivity analyses.

Preferred accounting model, using positive respiration:

```text
R_total[s,t] = R_other[s,t] + kappa[s,t] * R_fish[s,t]
ER_observed[day] ~ observation_model(R_total at that day and season)
```

R_other is non-fish ER. If spatial and temporal supports match, kappa can be fixed at 1. If a scaling factor is necessary, constrain it using independent information and test identifiability; avoid fitting a highly flexible kappa from about 20 annual metabolism states. Apply additive accounting on the original scale even if R_other has a log link.

Do not add fish respiration to a model that already represents total ER and thereby count it twice. Do not treat fishflux-derived predictions as independently observed respiration data or replicate each Monte Carlo draw as a new field observation. Fishflux predictions and fish observations share their underlying information.

Compute posterior R_fish_contribution / R_total and the change in fish respiration relative to the change in total ER over the observed record. A small fraction at individual times does not by itself quantify the contribution to a temporal change. If the denominator of a change ratio is near zero, report absolute changes rather than unstable percentages.

Treat the reported less-than-1-percent result as a check against reproducible existing calculations. Do not truncate the contribution at 1 percent. If prior calculations used these same fish data, do not reuse the result as an independent informative prior. Run a primary model including fish accounting, a model omitting the small component, and plausible high-activity or footprint scenarios. Compare changes in coral contrasts and forecasts. Direct fish respiration can be small even if fish indirectly affect other respiration through ecological pathways.

Output: reproducible fishflux calculation, conversion audit, uncertainty draws, contribution table, and reviewer sensitivity figure.

## Step 8 Specify dynamics and missingness

Use a complete annual process grid when annual states are defensible. Missing observations simply omit their observation likelihood while latent states continue through every year. Do not treat observations five years apart as one AR step. If working directly at irregular times, derive transitions depending on elapsed time, including process variance over the gap; consider a continuous-time process where needed.

A useful starting parameterization is:

```text
mu_X[s,t] = alpha_X + site_effect_X[s] + f_X(parents[s,t], environment[s,t])
X[s,t] = mu_X[s,t] + phi_X * (X[s,t-1] - mu_X[s,t-1]) + eta_X[s,t]
```

This is a candidate residual-persistence model, not a universal ecological equation. Compare it with biologically justified lagged-state dynamics. Avoid adding both a lagged outcome and residual AR structure without identifying their separate roles. No AR term alone guarantees control of shared trends or omitted environmental drivers.

Classify missingness separately for unsampled responses, predictors, sites, years, and seasons. Assess whether storm conditions, reef state, instrument failures, or logistics affect sampling. Use MAR only conditional on a justified measured history; evaluate selection or pattern-mixture sensitivity when MNAR is plausible. Never code missing fish as zero or absent metabolism at five sites as zero.

Missing predictors need a joint model consistent with the causal and temporal structure. Distinguish estimation of past missing states, which can use later data, from honest forecasting, which cannot. Future unknown environment needs scenarios or an explicit forecast distribution.

Output: transition equations, initial-state priors, observation masks, and missingness assumptions and sensitivities.

## Step 9 Check brms feasibility before coding the final model

Create a matrix of each required component, proposed brms syntax, documented support, identifiability, and fallback. Use `mi()` only with a properly modeled uncertain variable and correct indexing. A missing-response model is not automatically a state-space process; `(1 | site_year)` does not by itself create dynamic latent states; correlated multivariate group effects do not automatically share the same latent predictor across likelihoods.

The model needs mixed observation tables, joint composition states, latent predictors reused downstream, missing state transitions, and fish accounting. Verify all of these with a small simulated example using the installed brms version. Inspect generated Stan code with `make_stancode()` and data with `make_standata()` to confirm the intended process and observation likelihoods. Do not claim syntax works merely because a formula looks plausible.

Prefer a joint supported brms model if it genuinely implements the specification. If it cannot, retain dagitty and use one clearly labeled option:

1. A brms-based modular approximation: fit upstream measurement or state models, pass multiple complete posterior state trajectories to downstream models, and combine conditional posterior draws. Preserve each trajectory's temporal and multivariate dependence. This cuts downstream feedback and is not identical to a joint posterior. Compare with simulation and acknowledge the approximation.
2. A simplified brms model with explicitly reduced ambitions, such as annual latent-response models supported by indexed `mi()` terms. Demonstrate that process and observation variance remain distinct.
3. A custom Stan model through cmdstanr if an exact joint state-space model is essential and unsupported. Explain the gap to the investigator before treating this as the final implementation choice. Do not silently replace the requested brms workflow.

Never use posterior mean states as error-free downstream predictors. Never concatenate unrelated model draws and describe them as a coherent joint posterior.

Output: feasibility report and a tested minimal implementation, with every approximation labeled.

## Step 10 Choose priors and verify recoverability

Start parsimoniously: roughly 20 time points at one metabolism site cannot support many weakly identified metabolic paths just because there are many daily measurements. Replicates improve observation precision but do not create independent annual exposure contrasts.

Scale predictors using training data and retain constants for later prediction. Set regularizing coefficient priors and ecologically calibrated positive-scale priors. Example starting priors on standardized Gaussian scales are Normal(0, 0.5) coefficients and half-Normal(0, 1) standard deviations; revise these using prior predictive checks rather than treating them as universal. Specify persistence priors, initial-state priors, and site pooling explicitly. Do not force signs solely because they are the hypothesized result.

Simulate datasets with the actual six-site design, missingness patterns, composition zeros, daily sampling, and small fish contribution. Fit the intended model and assess recovery of process variance, observation variance, persistence, composition contrasts, and fish contribution. Evaluate the degree to which one-site metabolism and the tissue N proxy limit identification. Simplify unsupported mechanisms before fitting the real record.

Output: prior predictive plots and simulation recovery report.

## Step 11 Fit and diagnose

Use at least four chains initially and choose iterations based on effective sample sizes and Monte Carlo precision for target contrasts. Check Rhat near 1, preferably below 1.01, bulk and tail ESS, divergences, treedepth, and chain mixing. Reparameterize or improve priors when needed rather than relying only on higher adapt_delta.

Check posterior predictions at both replicate and state levels: composition totals and zeros, fish distributions, tissue N, daily metabolism, seasonal patterns, temporal autocorrelation, and sites. Examine parameter confounding and prior-to-posterior learning, especially latent N mechanisms and small fish contributions.

Derive testable conditional independencies with dagitty and assess them through models respecting measurement error, repeated sampling, and time. Do not apply ordinary independent-data correlation tests to all transects or interpret absence of rejection as proof of the graph. Explain any ecologically justified graph changes and preserve the original version.

Output: fitted objects, diagnostics, predictive checks, and DAG consistency report.

## Step 12 Calculate causal contrasts by posterior simulation

For each posterior draw, initialize the relevant states, apply the specified composition replacement, and simulate downstream nodes in topological and temporal order. Keep pre-exposure confounders governed by their appropriate distributions; allow descendants such as fish and tissue N to respond. Sum fish and non-fish ER in matched units. If using gross production, derive NCP by subtraction with coherent joint draws.

For a sustained intervention, replace the benthic transition at each declared intervention time. For a one-time intervention, replace the initial composition and allow subsequent ecological dynamics to proceed. These are different questions. A climate or disturbance change is also a different intervention from directly setting benthic composition because it has other downstream pathways.

Report posterior median, 50 and 95 percent intervals, original-scale contrasts, and direction probabilities where helpful. Calculate path contributions within each draw rather than multiplying posterior mean coefficients. A fish-component difference is an accounting attribution; it is not automatically an identified natural indirect effect.

Output: causal contrast tables labeled by identification assumptions and replacement policy.

## Step 13 Build a future metabolic state simulator

Implement and document a function with this interface or an equivalent:

```r
simulate_metabolic_future <- function(
  posterior_draws,
  last_state_draws,
  future_environment,
  horizon,
  benthic_policy = NULL,
  fish_scenario = NULL,
  include_process_noise = TRUE,
  include_observation_noise = FALSE,
  seed = 1
) {
  # Implement after the fitted state equations are verified.
  # This is an interface specification, not a functioning model.
}
```

For each draw and each future year:

1. Draw or select an environmental trajectory, preserving temporal dependence and shared temperature uncertainty.
2. Carry forward all necessary lagged states, not just coral cover.
3. Generate the next valid benthic composition from its transition equation, unless a declared intervention replaces it.
4. Generate nitrogen or tissue N and fish states using the chosen causal timing.
5. Update species or guild size distributions needed for fishflux. If these are not modeled dynamically, use declared trait and size scenarios with uncertainty.
6. Calculate fish respiration using that same draw's fish state, metabolic parameters, and temperature.
7. Generate non-fish ER and production or NCP with new process innovations.
8. Calculate total ER and, when valid, NCP = GPP - ER.
9. Save the entire state to drive the following year.
10. Optionally simulate a future observation with observation noise and a declared sampling design.

Keep three outputs distinct: parameter-uncertain conditional means, future latent ecological states including process uncertainty, and future field measurements including observation uncertainty. Do not present intervals for mean predictions as full future prediction intervals.

Primary forecasts apply to the metabolism site. Other-site outputs are transportability scenarios unless additional metabolism observations justify calibration. Future forecasts require stable mechanisms or explicit scenarios describing how they change. Flag environmental and compositional values outside the historical support. Avoid unsupported extrapolation of a linear temperature response over long climate horizons.

Output: a reusable simulator, documented inputs, and draw-level forecast table.

## Step 14 Declare forecasting scenarios

Use investigator-supplied trajectories or explicit assumptions, without inventing climate projections:

- Continuation under specified future environmental conditions.
- Coral decline with proportional non-coral replacement.
- Coral decline with mainly macroalgal replacement.
- Coral decline with mainly sand replacement.
- Stabilization or recovery under a declared composition policy.
- Disturbance pulse and recovery through the dynamic model.
- Fish respiration sensitivity under plausible metabolic, activity, or footprint assumptions.

If future coral cover is supplied externally, label the result a forecast conditional on that trajectory. To forecast coral endogenously, first fit and validate its environmental and disturbance transition model. Propagate uncertainty from any supplied climate or coral projections rather than replacing them with a single curve.

Save year, site, scenario, draw_id, covers, tissue N state, fish state, R_fish, R_other, ER_total, GPP if modeled, NCP, and extrapolation flags. Report ER and NCP jointly. Only compute P(NCP < 0) when NCP is a comparable whole-day net balance; do not apply that threshold to a daytime metric. Define any low-ER threshold before interpreting its probability.

Output: scenario specifications, trajectories, uncertainty bands, and metabolic state probabilities.

## Step 15 Validate forecasts using historical hindcasts

Use rolling-origin validation: fit only through year t and predict t+1 and longer available horizons. Hold out every replicate from future target site-years. Do not smooth missing training states using data from held-out future years, fit scaling or imputation on the complete record, or use held-out fish and benthos to claim an unconditional metabolism forecast.

Evaluate two distinct tasks if useful: forecasting all future states from past data, and predicting metabolism conditional on newly observed benthos or fish. Label them separately. Compare with persistence and simple environmental benchmarks. Evaluate MAE or RMSE, interval coverage and width, bias, and a distributional score such as CRPS where available. Validation at one site does not establish accuracy across sites.

Predeclare reasonable horizons according to data support and inspect deterioration with lead time. Predictive comparisons can improve forecasting equations within the scientific specification; they do not prove causal mechanisms. Save failed forecasts as well as successes.

Output: hindcast metrics, plots, baseline comparisons, and supported forecast horizons.

## Step 16 Run sensitivity analyses and assemble the project

Assess alternative confounding assumptions, replacement policies, lag choices, annual versus seasonal support, zero handling, priors, missingness, fish size or trait assumptions, and climate extrapolation. Where confounding is unresolved, report model-based scenarios conditional on assumptions rather than identified causal effects.

Create these project outputs:

```text
README.md
renv.lock
data_dictionary.csv
R/01_audit_data.R
R/02_define_composition.R
R/03_build_dags.R
R/04_identify_queries.R
R/05_fishflux_respiration.R
R/06_simulate_and_check_feasibility.R
R/07_fit_models.R
R/08_diagnostics.R
R/09_causal_contrasts.R
R/10_forecast_metabolism.R
R/11_hindcast_validation.R
R/12_sensitivity.R
analysis.qmd
outputs/dags/
outputs/diagnostics/
outputs/causal_contrasts.csv
outputs/fish_respiration_contributions.csv
outputs/forecast_draws.rds
outputs/forecast_summary.csv
outputs/hindcast_metrics.csv
outputs/assumptions_and_limitations.md
```

Use reproducible seeds, relative project paths, explicit dependencies, and saved package versions. Validate joins, composition inversion, unit conversions, additive accounting, gap transitions, and no future-data leakage. Maintain traceable correspondence among arrows, equations, code, and outputs.

The report must clearly distinguish measured evidence, hypothesized mechanisms, identified effects, assumption-dependent effects, and forecast scenarios. Do not claim absolute nitrogen flux from tissue N alone, a six-site metabolic relationship from one site's observations, or a fully joint latent state-space analysis when only modular regressions were fitted.

## Completion criteria

The work is complete when the DAGs parse and have justified arrows; composition interventions preserve the simplex; identification limitations are explicit; the implementation actually separates supported process and observation variation; missing states retain uncertainty; fish accounting is reproducible and does not double count; diagnostics and simulation checks are satisfactory; and the forecasting function generates and validates coherent future state trajectories.

The first response after reading this file should contain the data audit plan, candidate graph critique, unresolved scientific choices, and the brms feasibility strategy. Do not jump immediately to a final `brm()` call.

## References and software documentation

Scientific framework:

- Arif S and MacNeil MA. 2023. Applying the structural causal model framework for observational causal inference in ecology. Ecological Monographs 93, e1554. https://doi.org/10.1002/ecm.1554
- Arif S and MacNeil MA. 2022. Predictive models aren't for causal inference. Ecology Letters. https://doi.org/10.1111/ele.14033
- Arif S, Graham NAJ, Wilson S and MacNeil MA. 2022. Causal drivers of climate-mediated coral reef regime shifts. Ecosphere 13, e3956. https://doi.org/10.1002/ecs2.3956
- Schiettekatte NMD, Barneche DR, Villeger S and colleagues. 2020. Nutrient limitation, bioenergetics and stoichiometry A new model to predict elemental fluxes mediated by fishes. Functional Ecology 34, 1857–1869. https://doi.org/10.1111/1365-2435.13618

Implementation references checked for this specification; recheck the installed versions before coding:

- dagitty adjustment sets: https://search.r-project.org/CRAN/refmans/dagitty/html/adjustmentSets.html
- brms missing values: https://paulbuerkner.com/brms/articles/brms_missings.html
- brms uncertain and missing predictors: https://paulbuerkner.com/brms/reference/mi.html
- brms autoregressive structures: https://paulbuerkner.com/brms/reference/ar.html
- fishflux author repository: https://github.com/nschiett/fishflux
- fishflux author README: https://rdrr.io/github/nschiett/fishflux/f/README.Rmd
- fishflux metabolic rate source and documented inputs: https://rdrr.io/cran/fishflux/src/R/metabolic_rate.R

These sources ground the causal framework and software capabilities. The detailed study model remains a proposed specification requiring the checks above. No real-data analysis has been fitted in this instruction file.
