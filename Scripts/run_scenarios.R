# ============================================================
# Run all scenarios across MULTIPLE models, spin up to equilibrium,
# compare treatment vs baseline, and plot. Run from project root.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/run_models.R");      source("R/plot_eq_compare.R")

scen   <- read_scenarios("Data/scenarios.csv")
models <- c("century", "millennial", "MIMICS")

# ------------------------------------------------------------
# Run every scenario for every model in ONE call.
# Spin-up: BASELINE first, then TREATMENT warm-started from the baseline
# equilibrium (same plants + climate -> fast). method = "runsteady" is the
# robust forward-integration spin-up (replaces stode, which fails on MIMICS).
# To experiment: method = "stode", or spinup_treatment = FALSE, or
# warm_start_treatment = FALSE.
# ------------------------------------------------------------
out <- run_models(models, scen,
                  method               = "runsteady",
                  spinup_treatment     = TRUE,
                  warm_start_treatment = TRUE,
                  stol                 = 1e-8,
                  verbose              = TRUE)

# Convergence at a glance (model x scenario x arm):
print(spin_summary(out))

# Combined tidy comparison (model x scenario x pool):
View(out$combined)

# ------------------------------------------------------------
# Plots: one scenario across models, pools faceted with FREE axes so the
# very different magnitudes all stay readable.
# ------------------------------------------------------------
sc <- "Earthworm"
print(plot_eq_compare(eq_compare_list(out, sc), metric = "percent_change"))
print(plot_eq_compare(eq_compare_list(out, sc), metric = "treatment"))

# All pools for one scenario straight from the combined table:
# print(plot_eq_compare(subset(out$combined, scenario == sc), metric = "difference"))

# Per-model, one-at-a-time still works exactly as before, e.g.:
# res_cent <- lapply(names(scen), run_scenario, scen = scen, model = "century")
