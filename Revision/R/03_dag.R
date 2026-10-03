# ---------------------------------------------------------------------------
# 03_dag.R
#
# Step 3 of the Revision plan (Revision/Reviewer_Response_Plan.md, Section 4).
# Defines the revised, time-indexed structural causal model (DAG) in
# dagitty, derives minimal adjustment sets for every causal estimand (E1-E12,
# Task 4.2), and produces a coefficient-free DAG figure (Task 4.4).
#
# HARD STOP (Task 4.1 / Section 13): whether N (nitrogen) is caused by Rd
#
# UPDATE 2026-10-01 (post Task 4.3 d-sep test, Revision/R/03b_dsep_test.R):
# added Time -> Benthos_lag, Time -> Herb_lag, Time -> N_lag, Time -> COTS to
# BOTH variants below. The d-sep test found p < 1e-9 marginal correlations
# between each lag node and Time (expected: a lagged quantity is still
# time-varying) and a significant COTS-N/N_lag correlation consistent with
# COTS outbreaks clustering in specific years (Step 2) rather than being
# purely density-dependent on lagged local coral cover. Decided with the PI
# to add all four edges and re-test (see table_s_dsep.csv/.md, re-run after
# this change, for whether the violations resolve).
# (respiration-linked remineralisation) or directly by Benthos. Decided with
# the PI 2026-10-01: draw BOTH, carry both adjustment-set tables through this
# step, and revisit after the Task 4.3 d-separation test against data (not
# yet run -- see the note at the end of this script). Nothing here commits
# to one version; both are written to Output/ for the PI to compare.
#
# This script does NOT yet run Task 4.3 (testing the DAG against data via
# impliedConditionalIndependencies()/localTests()/a Shipley d-sep test).
# That requires assembling several things Step 1/2 did not build yet (ilr
# coordinates for the compositional Benthos node, a nitrogen covariate, and
# Rd/Pmax proxies ahead of the real Section 6 PI-curve fit) and is left for
# a follow-up step so this one stays focused on Tasks 4.1, 4.2 and 4.4.
#
# OUTPUTS:
#   Revision/Output/dag_mcr_rd_to_n.png / .pdf        (DAG figure, Rd -> N variant)
#   Revision/Output/dag_mcr_benthos_to_n.png / .pdf   (DAG figure, Benthos -> N variant)
#   Revision/Output/table_s_dag.csv                   (adjustment sets, both variants)
#   Revision/Output/table_s_dag.md
# ---------------------------------------------------------------------------

source(here::here("Revision", "R", "00_packages.R"))
if (!requireNamespace("dagitty", quietly = TRUE)) install.packages("dagitty")
if (!requireNamespace("ggdag", quietly = TRUE)) install.packages("ggdag")
library(dagitty)
library(ggdag)

# =============================================================================
# A. The two DAG variants
# =============================================================================
# Both share every edge in the plan's draft DAG (Section 4, Task 4.1) except
# the N-mechanism edge. Benthos is a single compositional node (Task 4.1
# note): the closure constraint of the simplex is not drawn as fake causal
# arrows among Coral/Algae/CCA. The coral<->algae feedback is unrolled in
# time via *_lag nodes, which removes the cycle. Rd -> Pmax is deliberately
# absent (Task 4.1): their covariance is modelled as a residual correlation
# inside the Section 6 PI-curve model, not as a causal path -- this is the
# plan's worked example of the DAG-vs-SEM distinction (Comment 1).

dag_str_common <- '
  Time      [pos="0,0"]
  DHW       [exposure, pos="1,0"]
  COTS      [pos="1,1"]
  Cyclone   [pos="1,2"]
  Benthos_lag [pos="0,3"]
  Herb_lag  [pos="0,4"]
  N_lag     [pos="0,5"]
  Benthos   [pos="2,2"]
  Herb      [pos="3,3"]
  Corall    [pos="3,4"]
  OtherFish [pos="3,5"]
  N         [pos="4,5"]
  Season    [pos="3,0"]
  DayTemp   [pos="4,0"]
  Flow      [pos="4,1"]
  Rd        [outcome, pos="5,2"]
  Pmax      [pos="5,3"]

  Time -> DHW
  Time -> Benthos
  Time -> Herb
  Time -> OtherFish
  Time -> Benthos_lag
  Time -> Herb_lag
  Time -> N_lag
  Time -> COTS
  Time -> N
  Benthos_lag -> COTS
  DHW -> Benthos
  COTS -> Benthos
  Cyclone -> Benthos
  Benthos_lag -> Benthos
  Herb_lag -> Benthos
  N_lag -> Benthos
  Benthos_lag -> Herb
  Herb_lag -> Herb
  Benthos -> Herb
  Benthos -> Corall
  Benthos -> OtherFish
  Herb -> N
  OtherFish -> N
  Corall -> N
  DHW -> DayTemp
  Season -> DayTemp
  Season -> Flow
  Benthos -> Rd
  Herb -> Rd
  Corall -> Rd
  OtherFish -> Rd
  DayTemp -> Rd
  Flow -> Rd
  Benthos -> Pmax
  DayTemp -> Pmax
  Flow -> Pmax
  N -> Pmax
'

dag_mcr_rd_to_n <- dagitty(paste0("dag {", dag_str_common, "  Rd -> N\n}"))
dag_mcr_benthos_to_n <- dagitty(paste0("dag {", dag_str_common, "  Benthos -> N\n}"))

dag_variants <- list(
  "Rd -> N"       = dag_mcr_rd_to_n,
  "Benthos -> N"  = dag_mcr_benthos_to_n
)

cat("DAG variants defined:", paste(names(dag_variants), collapse = "; "), "\n")
cat("Nodes:", paste(names(dag_mcr_rd_to_n), collapse = ", "), "\n")

# Sanity check: both should be acyclic (dagitty enforces this at parse time;
# this just confirms parsing succeeded and no edges were mistyped).
stopifnot(all(vapply(dag_variants, isAcyclic, logical(1))))
cat("Both variants parse as valid DAGs (acyclic).\n")

# =============================================================================
# B. Adjustment sets for every estimand (Task 4.2)
# =============================================================================
# Each row is ONE estimand with ITS OWN adjustment set and will get its own
# model in Section 7 (Arif & MacNeil 2023) -- never read one estimand's
# effect off a model built for another. E6/E12 use a vector exposure
# (dagitty's generalised backdoor criterion for a joint exposure set, i.e.
# "the fish groups" treated as a block); E11 is split into two single-
# exposure rows (Benthos -> Herb, Benthos -> Corall) because dagitty
# adjustment sets are defined per exposure-outcome pair, not per a combined
# outcome label as the plan's prose table shorthand suggests.

estimands <- tribble(
  ~id,    ~label,                                          ~exposure,                ~outcome,  ~effect,
  "E1",   "DHW -> Benthos (coral share)",                   "DHW",                    "Benthos", "total",
  "E2",   "COTS -> Benthos",                                "COTS",                   "Benthos", "total",
  "E3",   "Herb_lag -> Benthos (recovery mediation)",       "Herb_lag",                "Benthos", "total",
  "E4",   "Benthos -> Rd (direct, not via fish)",           "Benthos",                 "Rd",      "direct",
  "E5",   "Benthos -> Rd (total)",                          "Benthos",                 "Rd",      "total",
  "E6",   "Fish (Herb+Corall+OtherFish) -> Rd (direct)",    "Herb,Corall,OtherFish",   "Rd",      "direct",
  "E7",   "DayTemp -> Rd (physiological, direct)",          "DayTemp",                 "Rd",      "direct",
  "E8",   "DHW -> Rd (total)",                              "DHW",                     "Rd",      "total",
  "E9",   "Benthos -> Pmax (direct)",                       "Benthos",                 "Pmax",    "direct",
  "E10",  "DayTemp -> Pmax (direct)",                       "DayTemp",                 "Pmax",    "direct",
  "E11a", "Benthos -> Herbivores (total)",                  "Benthos",                 "Herb",    "total",
  "E11b", "Benthos -> Corallivores (total)",                "Benthos",                 "Corall",  "total",
  "E12",  "Fish (Herb+Corall+OtherFish) -> N (total)",      "Herb,Corall,OtherFish",   "N",       "total"
)

format_sets <- function(sets) {
  # dagitty::adjustmentSets() returns a dagitty.sets object (a list of
  # character vectors, one per minimal sufficient adjustment set). Several
  # estimands here have a unique minimal set; a few may have more than one
  # -- list all of them, semicolon-separated, rather than silently picking
  # one, so the PI can see when the choice of covariates is not unique.
  if (length(sets) == 0) return("NONE IDENTIFIABLE (no valid adjustment set)")
  paste(vapply(sets, function(s) {
    if (length(s) == 0) "{}" else paste0("{", paste(s, collapse = ", "), "}")
  }, character(1)), collapse = "; ")
}

get_adjustment_row <- function(dag, id, label, exposure, outcome, effect) {
  exposure_vec <- str_split(exposure, ",")[[1]]
  sets <- tryCatch(
    adjustmentSets(dag, exposure = exposure_vec, outcome = outcome, effect = effect, type = "minimal"),
    error = function(e) structure(list(), error = conditionMessage(e))
  )
  err <- attr(sets, "error")
  tibble(
    id = id, label = label, exposure = exposure, outcome = outcome, effect = effect,
    adjustment_sets = if (!is.null(err)) paste("ERROR:", err) else format_sets(sets),
    n_minimal_sets = if (!is.null(err)) NA_integer_ else length(sets)
  )
}

table_s_dag <- map_dfr(names(dag_variants), function(variant_name) {
  dag <- dag_variants[[variant_name]]
  pmap_dfr(estimands, function(id, label, exposure, outcome, effect) {
    get_adjustment_row(dag, id, label, exposure, outcome, effect)
  }) |>
    mutate(dag_variant = variant_name, .before = 1)
})

write_csv(table_s_dag, here("Revision", "Output", "table_s_dag.csv"))

dag_md <- c(
  "# Table S-DAG: minimal adjustment sets for every causal estimand",
  "",
  "Two DAG variants are carried through this table pending the Task 4.3",
  "d-separation test against data (not yet run -- see note at the end of",
  "`03_dag.R`): `Rd -> N` (respiration-linked remineralisation) vs",
  "`Benthos -> N` (nitrogen driven directly by standing benthic biomass).",
  "Decided with the PI 2026-10-01 to draw and tabulate both rather than",
  "choose one now.",
  "",
  "Each row is one causal estimand (Task 4.2); each gets its own model in",
  "Section 7 using only its own adjustment set (Arif & MacNeil 2023) -- an",
  "estimand's effect is never read off a model fit for a different estimand.",
  "E11 (`Benthos -> Corallivores / Herbivores`) is split into E11a/E11b",
  "because dagitty's adjustment sets are defined per exposure-outcome pair.",
  "",
  "| DAG variant | ID | Estimand | Effect | Adjustment set(s) | # minimal sets |",
  "|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(table_s_dag))) {
  dag_md <- c(dag_md, sprintf(
    "| %s | %s | %s | %s | %s | %s |",
    table_s_dag$dag_variant[i], table_s_dag$id[i], table_s_dag$label[i],
    table_s_dag$effect[i], table_s_dag$adjustment_sets[i],
    ifelse(is.na(table_s_dag$n_minimal_sets[i]), "NA", table_s_dag$n_minimal_sets[i])
  ))
}
writeLines(dag_md, here("Revision", "Output", "table_s_dag.md"))

cat("Done: Section B (adjustment sets). Wrote:\n",
    "  Revision/Output/table_s_dag.csv\n",
    "  Revision/Output/table_s_dag.md\n")
print(table_s_dag, n = 30)

# =============================================================================
# C. Coefficient-free DAG figure (Task 4.4)
# =============================================================================
# Replaces the old Figure 3 ("SEM" with beta labels on arrows). This is a
# pure causal-structure figure -- no coefficients, no fitted values -- per
# the reviewer's DAG-vs-SEM distinction (Comment 1). Effect estimates belong
# in a separate forest plot of estimands (Section 7's new Figure 4, not
# built yet). Node coordinates reuse the pos= hints baked into the dagitty
# string above so the two variants lay out identically except for the one
# differing edge, making them easy to compare side by side.

make_dag_plot <- function(dag, title) {
  tidy_dagitty(dag) |>
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(edge_colour = "grey50") +
    geom_dag_point(colour = "steelblue", size = 8, alpha = 0.6) +
    geom_dag_label_repel(aes(label = name), size = 3, seed = 2451,
                          box.padding = 0.4, label.size = NA,
                          fill = scales::alpha("white", 0.8)) +
    theme_dag() +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5))
}

p_dag_rd_to_n <- make_dag_plot(dag_mcr_rd_to_n, "Hypothesised causal structure (Rd -> N variant)")
p_dag_benthos_to_n <- make_dag_plot(dag_mcr_benthos_to_n, "Hypothesised causal structure (Benthos -> N variant)")

ggsave(here("Revision", "Output", "dag_mcr_rd_to_n.png"), p_dag_rd_to_n, width = 9, height = 7, dpi = 300)
ggsave(here("Revision", "Output", "dag_mcr_rd_to_n.pdf"), p_dag_rd_to_n, width = 9, height = 7)
ggsave(here("Revision", "Output", "dag_mcr_benthos_to_n.png"), p_dag_benthos_to_n, width = 9, height = 7, dpi = 300)
ggsave(here("Revision", "Output", "dag_mcr_benthos_to_n.pdf"), p_dag_benthos_to_n, width = 9, height = 7)

cat("Done: Section C (DAG figures). Wrote:\n",
    "  Revision/Output/dag_mcr_rd_to_n.png / .pdf\n",
    "  Revision/Output/dag_mcr_benthos_to_n.png / .pdf\n")

# =============================================================================
# NEXT STEP (not run here): Task 4.3 -- test the DAG against data
# =============================================================================
# impliedConditionalIndependencies(dag) lists the conditional independences
# implied by each variant; localTests() / a Shipley d-sep (Fisher's C) test
# would check those against the joined site-year data. Running this requires,
# beyond what Steps 1-2 built:
#   - ilr coordinates for the compositional Benthos node (Coral/Algae/CCA/
#     Other proportions -> 3 ilr coordinates), from benthic_site.csv
#   - a nitrogen covariate (N node) -- not yet disaggregated; candidate
#     sources are Data/raw_data/WaterColumnN.csv and
#     Data/raw_data/MCR_LTER_Macroalgal_CHN_2005_to_2024_20250616.csv
#   - Rd/Pmax PROXIES ahead of the real Section 6 PI-curve fit (e.g.
#     R_mean/GP_mean from pp_day.csv), clearly labelled as proxies since the
#     real Rd/Pmax are modelled quantities, not yet estimated
#   - fish trophic-group biomass joined at site-year (fish_transect.csv,
#     summed to Herb/Corall/OtherFish)
# This is left for a dedicated follow-up step rather than folded in here.
