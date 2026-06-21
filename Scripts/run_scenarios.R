# ============================================================
# Run CSV scenarios for one model: each animal scenario vs its
# matched no-animal baseline (same plants + same climate).
# Run from project root.  Edit `model` and `n_years` as needed.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R"); source("R/compare_functions.R")

model    <- "millennial"                       # "millennial" / "century" / "MIMICS"
n_years  <- 200
scen     <- read_scenarios("Data/scenarios.csv")
scenario_names = names(scen)

# Run each scenario:
out <- lapply(scenario_names, FUN = run_scenario,scen = scen, model = model)

# Get the summary statistics:
lapply(out, function(X) X$eq_compare)
