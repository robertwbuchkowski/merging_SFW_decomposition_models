# ============================================================
# SLOW STEP - run once per scenario, then reuse the saved states.
# Equilibrium spin-up -> (optional calibration) -> SEASONAL dynamic spin-up
# -> save the stable limit-cycle state to Data/spinup/.
# Run from project root. This can take a while; that is why it is separate.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")

scen     <- read_scenarios("Data/scenarios.csv")
model    <- "MIMICS"            # which model
scenario <- "Mite"         # which scenario
do_fit   <- T                   # calibrate treatment animal params first?
do_treatment <- TRUE            # also spin up the treatment arm now?

pair <- setup_scenario_pair(model, scen, scenario)

# ------------------------------------------------------------
# (part 1) CALIBRATION - fit the treatment animal to a reasonable biomass
# (defaults to its input value) and a user-defined effect on a pool.
# ------------------------------------------------------------
pair$baseline <- spinup_equilibrium(pair$baseline)        # needed as the effect reference

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
  dyn_t <- dynamic_spinup(pair$treatment, n_years = 300, by = 30, tol = 1e-4)
  cat("treatment converged:", dyn_t$converged, "\n")
  save_spinup(pair$treatment, dyn_t$final_state, scenario, "treatment")
}

cat("\nDone. Stable states saved under Data/spinup/. Use Scripts/followup_analysis.R next.\n")
