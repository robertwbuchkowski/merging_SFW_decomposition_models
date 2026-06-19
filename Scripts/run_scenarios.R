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

results <- vector('list', length(names(scen)))
names(results) = names(scen)

for (nm in names(scen)) {
  message("\n==================  ", nm, "  ==================")
  pair  <- setup_scenario_pair(model, scen, nm)

  # ============================================================
  # QUICK MASS-BALANCE CHECK at t = 0 (should be ~0)
  # ============================================================
  mb0 <- pair$treatment$wrapped_model(0, pair$treatment$working_state, pair$treatment$parms)[[2]]
  cat(sprintf("mass_balance_check at t=0: %.3e\n", mb0))
  
  mb0 <- pair$baseline$wrapped_model(0, pair$baseline$working_state, pair$baseline$parms)[[2]]
  cat(sprintf("mass_balance_check at t=0: %.3e\n", mb0))
  
  
  # ============================================================
  # Fast warm-start to constant-forcing steady state.
  # ============================================================
  
  pair$treatment$parms$climate_forcing = make_climate_forcing_equilibrium(pair$treatment$parms)
  cfss <- stode(y = pair$treatment$working_state,
                func = pair$treatment$wrapped_model,
                parms = pair$treatment$parms)
  pair$treatment$parms$climate_forcing = make_climate_forcing(pair$treatment$parms)
  
  if(attr(cfss, "steady")){
    cat("Result is steady")
    pair$treatment$init_state_spin = cfss[[1]]
    
  }else{
    cat("Result did not reach stability. Using input state.")
    pair$treatment$init_state_spin = pair$treatment$working_state
  }
  
  # ============================================================
  # QUICK MASS-BALANCE CHECK at t = 0 (should be ~0)
  # ============================================================
  mb0 <- pair$baseline$wrapped_model(0, pair$baseline$working_state, pair$baseline$parms)[[2]]
  cat(sprintf("mass_balance_check at t=0: %.3e\n", mb0))
  
  mb0 <- pair$baseline$wrapped_model(0, pair$baseline$working_state, pair$baseline$parms)[[2]]
  cat(sprintf("mass_balance_check at t=0: %.3e\n", mb0))
  
  
  # ============================================================
  # Fast warm-start to constant-forcing steady state.
  # ============================================================
  
  pair$baseline$parms$climate_forcing = make_climate_forcing_equilibrium(pair$baseline$parms)
  cfss <- stode(y = pair$baseline$working_state,
                func = pair$baseline$wrapped_model,
                parms = pair$baseline$parms)
  pair$baseline$parms$climate_forcing = make_climate_forcing(pair$baseline$parms)
  
  if(attr(cfss, "steady")){
    cat("Result is steady")
    pair$baseline$init_state_spin = cfss[[1]]
    
  }else{
    cat("Result did not reach stability. Using input state.")
    pair$baseline$init_state_spin = pair$baseline$working_state
  }
  
  results[[nm]] = list(
    pair = pair,
    eq_compare = compare_vectors(pair$treatment$init_state_spin, pair$baseline$init_state_spin)
  )
}

lapply(results, function(X) X["eq_compare"])
