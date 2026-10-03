# ---------------------------------------------------------------------------
# 00_packages.R
# Packages needed for the Revision/ scripts. Run once per session before any
# other Revision/R/*.R script. This file is specific to the revision work and
# does not touch the original analysis environment.
# ---------------------------------------------------------------------------

needed <- c(
  "tidyverse", "here", "lubridate", "janitor", "CCP"
)

# Packages needed by later steps (benthic composition, DAG, PI model, fish
# bioenergetics) are listed here for reference but are not required for
# Step 1 (data disaggregation). Install as later steps are implemented:
#   brms, cmdstanr, posterior, tidybayes, loo, priorsense,
#   dagitty, ggdag, compositions, zCompositions, fishflux, rerddap, sensemakr

missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  install.packages(missing)
}

library(tidyverse)
library(here)
library(lubridate)
library(janitor)

options(dplyr.summarise.inform = FALSE)

# Run brms/Stan chains in parallel rather than sequentially (default rstan
# backend runs chains one at a time otherwise, which is very slow for any
# multi-chain fit). Capped at 4 since that's what every brm() call in this
# project requests.
options(mc.cores = min(4, parallel::detectCores()))
