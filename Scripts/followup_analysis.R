# ============================================================
# FAST FOLLOW-UP - reuse saved stable states (from spinup_dynamic.R) for
# short (~100 yr) perturbation experiments:
#   * ADD animals to the spun-up BASELINE   (colonization)
#   * REMOVE animals from the spun-up TREATMENT (extirpation)
# Run from project root.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/dynamic_spinup.R")

scen     <- read_scenarios("Data/scenarios.csv")
model    <- "MIMICS"
scenario <- "Earthworm"

base_saved <- load_spinup(sprintf("Data/spinup/%s_%s_baseline.rds",  model, scenario))
trt_saved  <- load_spinup(sprintf("Data/spinup/%s_%s_treatment.rds", model, scenario))

# Rebuild the setups, then reuse the SAVED parameter lists so any calibration
# done during spin-up is preserved.
treatment_setup <- setup_scenario(model, scen, scenario, animals = TRUE)
treatment_setup$parms <- trt_saved$parms
baseline_setup  <- setup_scenario(model, scen, scenario, animals = FALSE)
baseline_setup$parms  <- base_saved$parms

# ---- ADD animals to the baseline limit cycle (seed = input biomass) ----
add <- followup_add_animals(base_saved, treatment_setup, n_years = 100, by = 30)

# ---- REMOVE animals from the treatment limit cycle ----
rem <- followup_remove_animals(trt_saved, baseline_setup, n_years = 100, by = 30)

# ---- inspect ----
cat("\nADD-animals end state:\n");    print(final_state(add$out))
cat("\nREMOVE-animals end state:\n"); print(final_state(rem$out))
# plot_ode_output(add$out)
# plot_ode_output(rem$out)

# Effect of the manipulation = end state vs the saved pre-manipulation state:
cat("\nAdd-animals effect (end - baseline limit cycle), shared pools:\n")
print(compare_vectors(final_state(add$out), base_saved$state))
cat("\nRemove-animals effect (end - treatment limit cycle), shared pools:\n")
print(compare_vectors(final_state(rem$out), trt_saved$state))
