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
models <- c("millennial")


# ------------------------------------------------------------
# EFFECT fitting -- OPT-IN and EXPLICIT. Nothing happens in the background:
# an effect is fit ONLY for animals that `effect_spec` gives a pool AND a
# target size for. Pick ONE of:
#
#   effect_spec <- list()                    # biomass-only fitting (fast, robust)
#   effect_spec <- load_effect_targets()     # ALSO fit effects, reading the saved
#                                            #   sizes from config/effect_targets.csv
#
# The CSV has one row per model x scenario x animal (model, scenario, animal,
# pool, pct, param) -- edit it to change which pool and how big an effect to
# fit. Blank/NA pool or pct = no effect fit for that animal. To (re)generate
# the file from the built-in defaults, run once:  save_effect_targets()
# ------------------------------------------------------------
effect_spec <- list()
# effect_spec <- load_effect_targets("Data/effect_targets.csv")

animals_order <- c("Earthworm", "Detritivore", "DetPredator", "RootHerb")  # prey before predator
tol_biomass   <- 0.02

scen$MitePredator = NULL

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

    # capture DEFAULT (baseline) parameter values up front, before any fitting.
    # effect_param is resolved per model x scenario (animal_fit_spec).
    defaults <- lapply(animals, function(a) {
      sp <- animal_fit_spec(a, model, scenario)
      bp <- sp$biomass_param
      ep <- sp$effect_param
      list(bp = bp, bp_default = if (!is.na(bp)) pair$treatment$parms[[bp]] else NA,
           ep = ep, ep_default = if (!is.na(ep) && ep %in% names(pair$treatment$parms))
                                  pair$treatment$parms[[ep]] else NA)
    })
    names(defaults) <- animals

    for (a in animals) {
      # Effect target for THIS model x scenario x animal, straight from
      # effect_spec. NULL (e.g. effect_spec <- list()) => biomass-only fit.
      es  <- get_effect_target(effect_spec, model, scenario, a)
      fit <- tryCatch(
        fit_animal_params(pair$treatment, pair$baseline, animal = a, scenario = scenario,
                          effect_pool       = es$pool,
                          target_effect_pct = es$pct,
                          effect_param      = es$param,
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

      # effect parameter row -- written whenever an effect was actually fit
      # (pool + target resolved for this model x scenario x animal).
      did_effect <- !is.na(f$effect_param) && !is.null(f$effect_pool) &&
                    !is.na(f$target_effect_pct)
      if (did_effect) {
        rows[[length(rows) + 1]] <- data.frame(
          model = model, scenario = scenario, animal = a,
          param = f$effect_param, role = paste0("effect on ", f$effect_pool),
          baseline = d$ep_default, fitted = f$fitted_effect_param,
          ratio = if (is.na(d$ep_default) || d$ep_default == 0) NA
                  else f$fitted_effect_param / d$ep_default,
          target = f$target_effect_pct, achieved = f$achieved_effect_pct,
          converged = is.finite(f$achieved_effect_pct) &&
                      abs(f$achieved_effect_pct - f$target_effect_pct) < 1,
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

# Save the fitted parameters (keyed by MODEL x scenario x param) for reuse in
# Scripts/spinup_dynamic.R -- so the spin-up reads them instead of re-fitting.
save_fitted_params(summary_long, "Results/fitted_animal_params.csv")

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

# ------------------------------------------------------------
# Interactive single scan + plot (run by hand while tuning). scan_animal_param
# returns an object plot_animal_scan() understands directly; the batch loop
# above just flattens many such scans into one CSV.
# ------------------------------------------------------------
if (FALSE) {
  model <- "millennial"; scenario <- "RootHerbivore"; a <- "RootHerb"
  pair <- setup_scenario_pair(model, scen, scenario)
  pair$baseline <- spinup_equilibrium(pair$baseline, verbose = FALSE)
  spec <- animal_fit_spec(a, model, scenario)
  grid <- fit_param_grid(pair$treatment$parms[[spec$biomass_param]], buffer = 2, n = 15)
  sc <- scan_animal_param(pair$treatment, param = spec$biomass_param, values = grid,
                          animal = a, baseline = pair$baseline,
                          effect_pool = spec$effect_pool)
  plot_animal_scan(sc, target_biomass = pair$treatment$working_state[a], log_x = TRUE)
  
  # Add in the fitted parameter:
  pair$treatment <- apply_fitted_params(pair$treatment, load_fitted_params("Results/fitted_animal_params.csv"), model, scenario)
  
  grid <- fit_param_grid(pair$treatment$parms[[spec$effect_param]], buffer = 2, n = 15, scale = "linear")
  
  grid <- seq(0, 1.3e-05, length = 10)
  sc <- scan_animal_param(pair$treatment, param = spec$effect_param, values = grid,
                          animal = a, baseline = pair$baseline,
                          effect_pool = spec$effect_pool)
  plot_animal_scan(sc, target_biomass = pair$treatment$working_state[a], log_x = TRUE)
  
  sc <- scan_animal_param(pair$treatment, param = spec$effect_param, values = grid,
                          animal = a, baseline = pair$baseline,
                          effect_pool = "C_root_herb")
  plot_animal_scan(sc, target_biomass = pair$treatment$working_state[a], log_x = TRUE)
}
