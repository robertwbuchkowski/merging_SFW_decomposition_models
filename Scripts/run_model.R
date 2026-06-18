# ============================================================
# RUN ONE MODEL — manual, linear, transparent
# ------------------------------------------------------------
# Pick a model in the SELECT block below, then step through the script.
# Switch plants/animals OFF by setting their initial values to 0 (a pool
# with zero biomass produces zero flux, so it stays off). No config DSL.
#
# Run from the project ROOT (so the relative paths work).
# ============================================================

library(pacman)
p_load(deSolve, rootSolve, tidyverse, yaml)

# ---- shared utilities ----
source("R/climate_forcing.R")
source("R/plot_ode_output.R")
source("R/spinup.R")
source("R/setup.R")


# # ---- switch MODEL by name, switch SCENARIO via `off` ----
run <- setup_model("millennial", off = c("Earthworm","Detritivore","DetPredator","RootHerb"))
print(run$working_state)


# ============================================================
# QUICK MASS-BALANCE CHECK at t = 0 (should be ~0)
# ============================================================
run$parms$climate_forcing <- make_climate_forcing(run$parms)
mb0 <- run$wrapped_model(0, run$working_state, run$parms)[[2]]
cat(sprintf("mass_balance_check at t=0: %.3e\n", mb0))


# ============================================================
# Fast warm-start to constant-forcing steady state.
# ============================================================

cfss <- stode(y = run$working_state,
              func = run$wrapped_model,
              parms = run$parms)

if(attr(cfss, "steady")){
  cat("Result is steady")
  run$init_state_spin = cfss[[1]]
  
}else{
  cat("Result did not reach stability. Using input state.")
  run$init_state_spin = run$working_state
}

# ============================================================
# Comparison at steady state with animal scenario
# ============================================================

run_animal <- setup_model("millennial", off = c("Earthworm","DetPredator","RootHerb"))
print(run_animal$working_state)

# QUICK MASS-BALANCE CHECK at t = 0 (should be ~0)
run_animal$parms$climate_forcing <- make_climate_forcing(run_animal$parms)
mb0 <- run_animal$wrapped_model(0, run_animal$working_state, run_animal$parms)[[2]]
cat(sprintf("mass_balance_check at t=0: %.3e\n", mb0))


# Fast warm-start to constant-forcing steady state.
cfss_animal <- stode(y = run_animal$working_state,
              func = run_animal$wrapped_model,
              parms = run_animal$parms)

if(attr(cfss_animal, "steady")){
  cat("Result is steady")
  run_animal$init_state_spin = cfss[[1]]
  
}else{
  cat("Result did not reach stability. Using input state.")
  rrun_animalun$init_state_spin = run_animal$working_state
}

# Compare with and without animal at equilibrium:

source("R/compare_functions.R")

compare_vectors(
  cfss$y,
  cfss_animal$y
)

# ============================================================
# SPIN-UP  — ONE semi-continuous run. Increase n_years until stable.
#   by = 1   -> daily output (see within-year dynamics)
#   by = 365 -> yearly snapshots (lighter output for very long runs)
# ============================================================

stable_dynamic = spinup_until_stable(run$init_state, run$parms)

deSolve::ode(y = run$init_state, times = 1:2, func = run$wrapped_model, parms = run$parms)
