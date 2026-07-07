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

manual_tune_add = T

# ------------------------------------------------------------
# EFFECT fitting. Effect targets now come from animal_fit_spec(), which
# resolves a pool + target size PER model x scenario x animal from
# effect_pool_overrides / animal_fit_defaults in R/fit_animals.R -- so you can
# have a different effect pool for every model x scenario. Edit them there.
#
#   fit_effects <- TRUE   also fit each animal's effect (pool + size) using
#                         those per model x scenario defaults.
#   fit_effects <- FALSE  biomass-only fitting (robust, always well-posed).
#
# `effect_spec` below still lets you FORCE an explicit pool/size for a given
# animal, overriding the per-scenario default for every scenario (rarely
# needed now; leave empty to use the per model x scenario table).
# ------------------------------------------------------------
fit_effects <- TRUE
effect_spec <- list()

animals_order <- c("Earthworm", "Detritivore", "DetPredator", "RootHerb")  # prey before predator
tol_biomass   <- 0.02

scen$MitePredator = NULL

rows <- list()
for (model in models) {
  for (scenario in names(scen)) {

    if(manual_tune_add){
      if(model == "MIMICS" & scenario == "Mite") next # Skip ones I did manually
    }
    
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
      # explicit effect_spec (if any) overrides the per model x scenario default;
      # otherwise fit_animal_params() pulls pool + size from animal_fit_spec().
      es  <- effect_spec[[a]]
      fit <- tryCatch(
        fit_animal_params(pair$treatment, pair$baseline, animal = a, scenario = scenario,
                          effect_pool = if (fit_effects) es$pool else NULL,
                          target_effect_pct = if (fit_effects) es$pct else NULL,
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

if(manual_tune_add){
  rows[[length(rows) + 1]] <- data.frame(
    model = "MIMICS", scenario = "Mite", animal = "Detritivore",
    param = "adj_detritivores", role = "biomass (feeding rate)",
    baseline = 1, fitted = 0.01529,
    ratio = 0.01529,
    target = 0.1, achieved = 0.1,
    converged = TRUE,
    stringsAsFactors = FALSE)
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

# ============================================================
# PARAMETER-GRADIENT SCANS (manual-tuning aid)
# ------------------------------------------------------------
# For each model x scenario x animal, sweep the feeding rate over a grid that
# BUFFERS around its default value (fit_param_grid) and record equilibrium
# biomass + the effect on the per-model effect_pool (from animal_fit_defaults
# via animal_fit_spec). Non-converged / unstable points are flagged, not fatal.
# Saved to Results/animal_scan_long.csv so you can see biomass and effect along
# the gradient and tune by hand where the automatic optimum lands in an
# unstable region.
# ============================================================
do_scan      <- TRUE
scan_buffer  <- 2       # grid spans default x 10^(-buffer) .. x 10^(+buffer)
scan_n       <- 13      # points in the grid
track_effect <- TRUE    # also record the effect on the per-model effect_pool

if (do_scan) {
  scan_rows <- list()
  for (model in models) {
    for (scenario in names(scen)) {

      pair <- tryCatch(setup_scenario_pair(model, scen, scenario),
                       error = function(e) NULL)
      if (is.null(pair)) next
      pair$baseline <- spinup_equilibrium(pair$baseline, verbose = FALSE)

      animals <- animals_order[animals_order %in% pair$treatment$active]
      for (a in animals) {
        spec <- animal_fit_spec(a, model, scenario)   # per model x scenario effect_pool
        bp   <- spec$biomass_param
        if (is.na(bp)) next
        # gradient buffered around the DEFAULT feeding rate
        grid  <- fit_param_grid(pair$treatment$parms[[bp]], buffer = scan_buffer, n = scan_n)
        epool <- if (track_effect && !is.na(spec$effect_pool)) spec$effect_pool else NULL

        sc <- tryCatch(
          scan_animal_param(pair$treatment, param = bp, values = grid, animal = a,
                            baseline = pair$baseline, effect_pool = epool,
                            verbose = FALSE),
          error = function(e) { message("scan failed ", model, "/", scenario, "/", a,
                                        ": ", conditionMessage(e)); NULL })
        if (is.null(sc)) next

        scan_rows[[length(scan_rows) + 1]] <- data.frame(
          model = model, scenario = scenario, animal = a, param = bp,
          effect_pool = if (is.null(epool)) NA_character_ else epool,
          value = sc[[bp]], biomass = sc$biomass, effect_pct = sc$effect_pct,
          converged = sc$converged, max_deriv = sc$max_deriv,
          stringsAsFactors = FALSE)
      }
      cat("Scanned", scenario, "for", model, "\n")
    }
  }

  if (length(scan_rows)) {
    scan_long <- do.call(rbind, scan_rows)
    write.csv(scan_long, "Results/animal_scan_long.csv", row.names = FALSE)
    cat("\nSaved Results/animal_scan_long.csv (", nrow(scan_long), " rows ).\n", sep = "")
  }
}

# ------------------------------------------------------------
# Interactive single scan + plot (run by hand while tuning). scan_animal_param
# returns an object plot_animal_scan() understands directly; the batch loop
# above just flattens many such scans into one CSV.
# ------------------------------------------------------------
if (FALSE) {
  model <- "MIMICS"; scenario <- "MitePredator"; a <- "DetPredator"
  pair <- setup_scenario_pair(model, scen, scenario)
  pair$baseline <- spinup_equilibrium(pair$baseline, verbose = FALSE)
  spec <- animal_fit_spec(a, model, scenario)
  grid <- fit_param_grid(pair$treatment$parms[[spec$biomass_param]], buffer = 2, n = 15)
  grid <- seq(0.001, 0.01, length.out = 20)
  sc <- scan_animal_param(pair$treatment, param = spec$biomass_param, values = grid,
                          animal = a, baseline = pair$baseline,
                          effect_pool = spec$effect_pool)
  plot_animal_scan(sc, target_biomass = pair$treatment$working_state[a], log_x = TRUE)
  
  
  grid <- fit_param_grid(pair$treatment$parms[[spec$effect_param]], buffer = 2, n = 15)
  sc <- scan_animal_param(pair$treatment, param = spec$effect_param, values = grid,
                          animal = a, baseline = pair$baseline,
                          effect_pool = spec$effect_pool)
  plot_animal_scan(sc, target_biomass = pair$treatment$working_state[a], log_x = TRUE)
}
