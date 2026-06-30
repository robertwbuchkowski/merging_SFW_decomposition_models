# ============================================================
# Fit animal parameters for ALL models x ALL scenarios and summarise how the
# fitted parameters differ across models (with the baseline/default value for
# reference). Run from project root. This calibrates on the (fast) constant-
# forcing equilibrium only -- no seasonal spin-up -- but still does many
# runsteady solves, so it takes a few minutes.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")

scen   <- read_scenarios("Data/scenarios.xlsx")
models <- c("century", "millennial", "MIMICS")

# ------------------------------------------------------------
# OPTIONAL effect targets (user-defined). Empty list = biomass-only fitting,
# which is robust and always well-posed. To also fit an effect, add an entry
# per animal with an ACHIEVABLE, correctly-signed target, e.g.:
#   effect_spec <- list(Detritivore = list(pool = "SOM_1", pct = +10))
# ------------------------------------------------------------
effect_spec <- list()

animals_order <- c("Earthworm", "Detritivore", "DetPredator", "RootHerb")  # prey before predator
tol_biomass   <- 0.02

scen$MitePredator = NULL

models <- c("century", "millennial")

rows <- list()
for (model in models) {
  for (scenario in names(scen)) {

    pair <- tryCatch(setup_scenario_pair(model, scen, scenario),
                     error = function(e) { message("setup failed ", model, "/", scenario,
                                                    ": ", conditionMessage(e)); NULL })
    if (is.null(pair)) next

    pair$baseline <- spinup_equilibrium(pair$baseline, verbose = FALSE)  # effect reference

    animals <- animals_order[animals_order %in% pair$treatment$active]
    if (!length(animals)) next

    # capture DEFAULT (baseline) parameter values up front, before any fitting
    defaults <- lapply(animals, function(a) {
      bp <- animal_fit_defaults[[a]]$biomass_param
      ep <- animal_fit_defaults[[a]]$effect_param
      list(bp = bp, bp_default = if (!is.na(bp)) pair$treatment$parms[[bp]] else NA,
           ep = ep, ep_default = if (!is.na(ep)) pair$treatment$parms[[ep]] else NA)
    })
    names(defaults) <- animals

    for (a in animals) {
      es <- effect_spec[[a]]
      fit <- tryCatch(
        fit_animal_params(pair$treatment, pair$baseline, animal = a,
                          effect_pool = es$pool, target_effect_pct = es$pct,
                          verbose = FALSE),
        error = function(e) { message("fit failed ", model, "/", scenario, "/", a,
                                      ": ", conditionMessage(e)); NULL })
      if (is.null(fit)) next
      pair$treatment <- fit                      # carry forward for later animals
      f <- fit$fit; d <- defaults[[a]]

      # biomass (feeding-rate) parameter row
      rows[[length(rows) + 1]] <- data.frame(
        model = model, scenario = scenario, animal = a,
        param = d$bp, role = "biomass (feeding rate)",
        baseline = d$bp_default, fitted = f$fitted_biomass_param,
        ratio = f$fitted_biomass_param / d$bp_default,
        target = f$target_biomass, achieved = f$achieved_biomass,
        converged = abs(f$achieved_biomass - f$target_biomass) /
                    max(abs(f$target_biomass), 1e-8) < tol_biomass,
        stringsAsFactors = FALSE)

      # effect parameter row (only if an effect target was given for this animal)
      if (!is.null(es) && !is.na(d$ep)) {
        rows[[length(rows) + 1]] <- data.frame(
          model = model, scenario = scenario, animal = a,
          param = d$ep, role = "effect on pool",
          baseline = d$ep_default, fitted = f$fitted_effect_param,
          ratio = f$fitted_effect_param / d$ep_default,
          target = es$pct, achieved = f$achieved_effect_pct,
          converged = abs(f$achieved_effect_pct - es$pct) < 1,
          stringsAsFactors = FALSE)
      }
    }
    cat("Done", scenario, "for", model, "\n")
  }
}

if (!length(rows)) stop("No fits succeeded - check the model/scenario setup.")
summary_long <- do.call(rbind, rows)
dir.create("Results", showWarnings = FALSE)
write.csv(summary_long, "Results/animal_fit_summary_long.csv", row.names = FALSE)

cat("\n================= FITTED ANIMAL PARAMETERS (long) =================\n")
print(summary_long, digits = 4)

# ------------------------------------------------------------
# Wide view: fitted feeding rate by MODEL (baseline shown once), so you can
# read across a row to see how different the models are.
# ------------------------------------------------------------
bm <- summary_long[summary_long$role == "biomass (feeding rate)", ]
wide <- reshape(
  bm[, c("scenario", "animal", "param", "baseline", "model", "fitted")],
  idvar = c("scenario", "animal", "param", "baseline"),
  timevar = "model", direction = "wide")
names(wide) <- sub("^fitted\\.", "fitted_", names(wide))
write.csv(wide, "Results/animal_fit_feeding_rate_by_model.csv", row.names = FALSE)

cat("\n========= FITTED FEEDING RATE BY MODEL (baseline = default) =========\n")
print(wide, digits = 4)

cat("\nSaved Results/animal_fit_summary_long.csv and Results/animal_fit_feeding_rate_by_model.csv\n")
