# ============================================================
# SLOW STEP - run once per scenario, then reuse the saved states.
# Equilibrium spin-up -> (optional calibration) -> SEASONAL dynamic spin-up
# -> save the stable limit-cycle state to Data/spinup/.
# Run from project root. This can take a while; that is why it is separate.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")

scen   <- read_scenarios("Data/scenarios.xlsx")
models <- c("century", "millennial", "MIMICS")
do_fit   <- F                # calibrate treatment animal params first?
do_treatment <- F            # also spin up the treatment arm now?

for (model in models) {
  for (scenario in names(scen)) {
    
    pair <- setup_scenario_pair(model, scen, scenario)
    
    # ------------------------------------------------------------
    # (part 1) CALIBRATION - fit the treatment animal to a reasonable biomass
    # (defaults to its input value) and a user-defined effect on a pool.
    # ------------------------------------------------------------
    pair$baseline <- spinup_equilibrium(pair$baseline)        # needed as the effect reference

    # For response curves along a parameter gradient (manual tuning across
    # stable/unstable regions), see Scripts/fit_all_animals.R (scan_animal_param).

    if (do_fit) {
      pair$treatment <- fit_animal_params(
        pair$treatment, pair$baseline,
        animal            = "Detritivore",       # the animal in this scenario
        target_biomass    = NULL,         # NULL = match its input (starting) value
        effect_pool       = NULL,         # user-defined effect target:
        target_effect_pct = NULL              #   target percent change
      )
      cat("\nCalibration history:\n"); print(pair$treatment$fit$history)
    }

    # ------------------------------------------------------------
    # QUICK CHECK - equilibrium animal biomass and the animal's effect on ALL
    # state variables (treatment vs baseline), before the long seasonal spin-up.
    # ------------------------------------------------------------
    if (is.null(pair$treatment$init_state_spin))
      pair$treatment <- spinup_equilibrium(pair$treatment,
                                           warm_start = pair$baseline$init_state_spin)
    eq_t    <- pair$treatment$init_state_spin
    animals <- intersect(c("Earthworm", "Detritivore", "DetPredator", "RootHerb"),
                         names(eq_t))
    cat("\nEquilibrium animal biomass (", model, "/", scenario, "):\n", sep = "")
    print(round(eq_t[animals], 4))
    cat("\nAnimal effect on all state variables (treatment vs baseline):\n")
    print(compare_vectors(eq_t, pair$baseline$init_state_spin), digits = 4)
    
    # ------------------------------------------------------------
    # (part 2) BASELINE FIRST: equilibrium -> seasonal dynamic spin-up -> save
    # ------------------------------------------------------------
    cat("\n--- Baseline dynamic spin-up ---\n")
    dyn_b <- dynamic_spinup(pair$baseline, n_years = 300, by = 1, tol = 1e-4)
    cat("baseline converged:", dyn_b$converged, "\n")
    save_spinup(pair$baseline, dyn_b$final_state, scenario, "baseline")
    
    # ------------------------------------------------------------
    # TREATMENT (flexible - can be run later in a separate session)
    # ------------------------------------------------------------
    if (do_treatment) {
      cat("\n--- Treatment dynamic spin-up ---\n")
      if (is.null(pair$treatment$init_state_spin))
        pair$treatment <- spinup_equilibrium(pair$treatment,
                                             warm_start = pair$baseline$init_state_spin)
      dyn_t <- dynamic_spinup(pair$treatment, n_years = 300, by = 1, tol = 1e-4)
      cat("treatment converged:", dyn_t$converged, "\n")
      save_spinup(pair$treatment, dyn_t$final_state, scenario, "treatment")
    }
    
    cat("\nDone. Stable states saved under Data/spinup/. Use Scripts/followup_analysis.R next.\n")
    
  }
}
