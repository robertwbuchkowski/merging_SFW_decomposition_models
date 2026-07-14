# ============================================================
# SLOW STEP, IN PARALLEL - same job as Scripts/spinup_dynamic.R, but the
# scenarios are spread over multiple cores.
#
# Scenarios are independent (each is its own ODE solve writing its own .rds),
# so this is close to a linear speed-up: 5 scenarios on 5 cores ~= the time of
# the single slowest scenario, instead of the sum of all of them.
#
#   ncores <- 1              sequential (identical to spinup_dynamic.R; debug)
#   ncores <- detect_cores() all physical cores but one (default)
#
# Works on macOS/Linux (fork) and Windows (PSOCK) -- see R/parallel_utils.R.
# Run from the project root.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")
source("R/parallel_utils.R")

scen   <- read_scenarios("Data/scenarios.xlsx")
models <- c("millennial")

use_fitted_params <- TRUE     # apply saved fitted params (from fit_all_animals.R)?
do_spinup         <- TRUE     # run the long seasonal dynamic spin-ups?
do_treatment      <- TRUE     # also spin up the treatment arm?
ncores            <- detect_cores()   # e.g. set to 4 to cap it; 1 = sequential

scen$MitePredator <- NULL

fitted_params <- if (use_fitted_params)
  load_fitted_params("Results/fitted_animal_params.csv") else NULL

cat("scenarios:", paste(names(scen), collapse = ", "), "\n")
cat("cores    :", ncores, "of", parallel::detectCores(logical = FALSE), "physical\n\n")

# ------------------------------------------------------------
# THE PARALLEL BIT - one scenario per core. Each worker runs the baseline
# equilibrium, applies the fitted params, computes the total vs direct-only
# animal effect, does the seasonal dynamic spin-ups and saves its own
# Data/spinup/*.rds. Failures are reported, not fatal.
# ------------------------------------------------------------
results <- list()
for (model in models) {
  res <- spinup_scenarios_parallel(
    model, scen,
    scenarios     = names(scen),
    ncores        = ncores,
    fitted_params = fitted_params,
    do_spinup     = do_spinup,
    do_treatment  = do_treatment,
    n_years = 600, by = 1, tol = 1e-4)
  results <- c(results, par_ok(res))
}

# ---- convergence report -------------------------------------------------
conv <- do.call(rbind, lapply(results, function(r)
  data.frame(model = r$model, scenario = r$scenario,
             baseline = r$converged[["baseline"]],
             treatment = r$converged[["treatment"]])))
cat("\n--- dynamic spin-up convergence ---\n"); print(conv)

# ---- equilibrium animal biomass ----------------------------------------
for (r in results) {
  a <- intersect(c("Earthworm", "Detritivore", "DetPredator", "RootHerb"),
                 names(r$eq_treatment))
  cat("\nEquilibrium animal biomass (", r$model, "/", r$scenario, "):\n", sep = "")
  print(round(r$eq_treatment[a], 4))
}

# ---- animal effects (total vs direct-only) ------------------------------
animal_eq_effect <- do.call(rbind, lapply(results, `[[`, "eq_effect"))
dir.create("Results", showWarnings = FALSE)
write_csv(animal_eq_effect, "Results/animal_eq_effect.csv")

# sanity check: any state variable with no baseline counterpart?
animal_eq_effect %>% filter(is.na(baseline))

cat("\nDone. Stable states saved under Data/spinup/.",
    "Use Scripts/followup_analysis.R next.\n")
