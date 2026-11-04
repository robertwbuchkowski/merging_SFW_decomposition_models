#------------------------------------------------------------------------#
# Seasonal spin-up and animal follow-up (invasion) simulations ----
#------------------------------------------------------------------------#
# The only part of the project that simulates the SEASONAL dynamics.
#   SPIN-UP  each scenario's no-animal soil is brought to its seasonal limit
#           cycle (Newton shooting, warm-started from the constant-forcing
#           equilibrium) and saved to Data/spinup/. Optional: the treatment
#           arm too (do_treatment). Set do_spinup <- FALSE to reuse saved states.
# Then FOLLOW-UP experiments from the saved stable states (~50 yr), looped
# over ALL models x scenarios:
#   LOOP 1  ADD animals to the spun-up BASELINE (colonization), plus a
#           time-matched CONTINUED-BASELINE control (same setup, no animals
#           added) so the two can be plotted against each other.
#   LOOP 2  REMOVE animals from the spun-up TREATMENT (extirpation) -- kept
#           separate so it can be run on its own, independent of Loop 1.
# Each run is saved to Data/followup/ keyed by model_scenario_kind.rds, so
# either loop can be (re)run independently and plotting can happen later,
# in a fresh session, without re-running any simulation.
# Run from project root. Requires Results/fitted_animal_params.csv from
# Scripts/1_fit_all_animals.R when use_fitted_params = TRUE.
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")
source("R/millennial_model.R"); source("R/derive_millennial_parms.R")
source("R/init_millennial_state.R")
source("R/run_utils.R")
debugonce_functions()   # only active when options(sfw.debugonce = TRUE); see Scripts/0_run_all.R

scen   <- read_scenarios("Data/scenarios.xlsx")
models <- c("millennial")
use_fitted_params <- TRUE    # apply saved fitted params (from fit_all_animals.R)?

n_years <- 50                    # length of every follow-up run below
by      <- 1

# Fitted animal parameters saved by Scripts/fit_all_animals.R, keyed BY MODEL
# (x scenario x param). Read once here; applied per model/scenario in the loop
# instead of re-fitting.
fitted_params <- if (use_fitted_params)
  load_fitted_params("Results/fitted_animal_params.csv") else NULL

#------------------------------------------------------------------------#
# SPIN-UP: seasonal limit cycle, saved to Data/spinup/ ----
#------------------------------------------------------------------------#
do_spinup    <- TRUE    # run the seasonal spin-up (FALSE = reuse saved Data/spinup/ states)
do_treatment <- FALSE   # also spin up the treatment arm (needed only by LOOP 2)
use_newton   <- TRUE    # TRUE = Newton shooting (fast, exact limit cycle);
                        # FALSE = forward-integration spin-up

if (do_spinup) {
  for (model in models) {
    for (scenario in names(scen)) {
      cat("\n==== SPIN-UP:", model, "/", scenario, "====\n")
      pair <- setup_scenario_pair(model, scen, scenario)

      # constant-forcing baseline equilibrium: warm start for the seasonal spin-up
      pair$baseline <- spinup_equilibrium(pair$baseline)

      # baseline: seasonal spin-up to the limit cycle, then save
      cat("\n--- Baseline dynamic spin-up ---\n")
      dyn_b <- if (use_newton) dynamic_spinup_newton(pair$baseline)
               else dynamic_spinup(pair$baseline, n_years = 600, by = 1, tol = 1e-4)
      cat("baseline converged:", dyn_b$converged, "\n")
      save_spinup(pair$baseline, dyn_b$final_state, scenario, "baseline")

      # treatment: optional (fitted animal parameters applied first)
      if (do_treatment) {
        cat("\n--- Treatment dynamic spin-up ---\n")
        if (use_fitted_params)
          pair$treatment <- apply_fitted_params(pair$treatment, fitted_params, model, scenario)
        pair$treatment <- spinup_equilibrium(pair$treatment,
                                             warm_start = pair$baseline$init_state_spin)
        dyn_t <- if (use_newton) dynamic_spinup_newton(pair$treatment)
                 else dynamic_spinup(pair$treatment, n_years = 600, by = 1, tol = 1e-4)
        cat("treatment converged:", dyn_t$converged, "\n")
        save_spinup(pair$treatment, dyn_t$final_state, scenario, "treatment")
      }
      rm(pair, dyn_b)
    }
  }
  cat("\nStable states saved under Data/spinup/.\n")
}

#------------------------------------------------------------------------#
# LOOP 1: ADD animals + CONTINUED-BASELINE control ----
#------------------------------------------------------------------------#
# For each model x scenario with saved spin-ups, reuses the SAVED parameter
# lists (so any calibration from fit_all_animals.R / the spin-up above is
# preserved), then runs:
#   add               animals introduced into the baseline limit cycle
#   continue_baseline the same baseline limit cycle continued with no animals
# Both are saved to Data/followup/ for plotting (see plot_followup_add()).
add_results <- list()

for (model in models) {
  for (scenario in names(scen)) {

    base_file <- sprintf("Data/spinup/%s_%s_baseline.rds",  model, scenario)
    if (!file.exists(base_file)) {
      message("skip (no saved spin-up): ", model, " / ", scenario)
      next
    }
    cat("\n==== ADD + continue-baseline:", model, "/", scenario, "====\n")

    base_saved <- load_spinup(base_file)

    treatment_setup <- setup_scenario(model, scen, scenario, animals = TRUE)
    baseline_setup  <- setup_scenario(model, scen, scenario, animals = FALSE)
    baseline_setup$parms  <- base_saved$parms
    treatment_setup$parms <- base_saved$parms
    
    if (use_fitted_params) {
      baseline_setup <- apply_fitted_params(baseline_setup, fitted_params, model, scenario)
      treatment_setup <- apply_fitted_params(treatment_setup, fitted_params, model, scenario)
    }
    
    add     <- followup_add_animals(base_saved, treatment_setup, n_years = n_years, by = by)
    control <- followup_continue_baseline(base_saved, baseline_setup, n_years = n_years, by = by)

    save_followup(model, scenario, "add", add)
    save_followup(model, scenario, "continue_baseline", control)

    key <- paste(model, scenario, sep = ".")
    add_results[[key]] <- list(add = add, control = control)

    cat("Add-animals effect (end - baseline limit cycle), shared pools:\n")
    print(compare_vectors(final_state(add$out), base_saved$state))
  }
}

#------------------------------------------------------------------------#
# LOOP 2: REMOVE animals ----
#------------------------------------------------------------------------#
# LOOP 2: REMOVE animals (fully independent of Loop 1 -- reads only the
# saved spin-ups from disk, so it can be run on its own / in a fresh session).
remove_results <- list()

for (model in models) {
  for (scenario in names(scen)) {

    base_file <- sprintf("Data/spinup/%s_%s_baseline.rds",  model, scenario)
    trt_file  <- sprintf("Data/spinup/%s_%s_treatment.rds", model, scenario)
    if (!file.exists(base_file) || !file.exists(trt_file)) {
      message("skip (no saved spin-up): ", model, " / ", scenario)
      next
    }
    cat("\n==== REMOVE:", model, "/", scenario, "====\n")

    base_saved <- load_spinup(base_file)
    trt_saved  <- load_spinup(trt_file)

    baseline_setup <- setup_scenario(model, scen, scenario, animals = FALSE)
    baseline_setup$parms <- base_saved$parms

    rem <- followup_remove_animals(trt_saved, baseline_setup, n_years = n_years, by = by)
    save_followup(model, scenario, "remove", rem)

    key <- paste(model, scenario, sep = ".")
    remove_results[[key]] <- rem

    cat("Remove-animals effect (end - treatment limit cycle), shared pools:\n")
    print(compare_vectors(final_state(rem$out), trt_saved$state))
  }
}

#------------------------------------------------------------------------#
# PLOTTING: continued baseline vs. animals added ----
#------------------------------------------------------------------------#
# PLOTTING: continued baseline (no animals) vs. animals-added, per
# model/scenario. Works directly from the saved Data/followup/*.rds files, so
# this can be run later without re-running either loop above.
# Example, one model/scenario:
plot_followup_add("millennial", "RootHerbivore", by = NULL)

dir.create("Results/figures", showWarnings = FALSE, recursive = TRUE)
pdf("Results/figures/followup_trajectories.pdf", width = 8, height = 8)
# All model x scenario combos that have saved add + continue_baseline runs:
for (scenario in names(scen)) {
  for (model in models) {
    if (file.exists(sprintf("Data/followup/%s_%s_add.rds", model, scenario)) &&
        file.exists(sprintf("Data/followup/%s_%s_continue_baseline.rds", model, scenario))) {
      print(plot_followup_add(model, scenario, by = 365))
    }
  }
}
dev.off()

#------------------------------------------------------------------------#
# COMBINED animal-effect graphic ----
#------------------------------------------------------------------------#
# COMBINED animal-effect graphic: the stacked C-change trajectory over the
# follow-up PLUS the equilibrium effect (from Scripts/2_animal_equilibrium_effect.R ->
# Results/animal_eq_effect.csv) as a stacked bar past a broken x-axis, with
# matching pool colours and a dashed net-change line. One facet per scenario.
# This single figure replaces the separate equilibrium-effect plot.
p_combined <- plot_followup_with_eq("millennial", by = 365,
                                    eq_csv = "Results/animal_eq_effect.csv",
                                    eq_type = "total")
ggsave("Results/figures/followup_added_with_equilibrium.png", p_combined,
       width = 11, height = 7, dpi = 300, bg = "white")
print(p_combined)
