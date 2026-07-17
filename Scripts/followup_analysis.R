# ============================================================
# FAST FOLLOW-UP - reuse saved stable states (from spinup_dynamic.R) for
# short (~100 yr) perturbation experiments, looped over ALL models x scenarios:
#   LOOP 1  ADD animals to the spun-up BASELINE (colonization), plus a
#           time-matched CONTINUED-BASELINE control (same setup, no animals
#           added) so the two can be plotted against each other.
#   LOOP 2  REMOVE animals from the spun-up TREATMENT (extirpation) -- kept
#           separate so it can be run on its own, independent of Loop 1.
# Each run is saved to Data/followup/ keyed by model_scenario_kind.rds, so
# either loop can be (re)run independently and plotting can happen later,
# in a fresh session, without re-running any simulation.
# Run from project root; requires Data/spinup/*.rds from spinup_dynamic.R.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")

scen   <- read_scenarios("Data/scenarios.xlsx")
models <- c("millennial")
use_fitted_params <- TRUE    # apply saved fitted params (from fit_all_animals.R)?

scen$MitePredator = NULL          # match the scenarios actually spun up

n_years <- 25                    # length of every follow-up run below
by      <- 1

# Fitted animal parameters saved by Scripts/fit_all_animals.R, keyed BY MODEL
# (x scenario x param). Read once here; applied per model/scenario in the loop
# instead of re-fitting.
fitted_params <- if (use_fitted_params)
  load_fitted_params("Results/fitted_animal_params.csv") else NULL

# ------------------------------------------------------------
# LOOP 1: ADD animals + CONTINUED-BASELINE control
# For each model x scenario with saved spin-ups, reuses the SAVED parameter
# lists (so any calibration from fit_all_animals.R / spinup_dynamic.R is
# preserved), then runs:
#   add               animals introduced into the baseline limit cycle
#   continue_baseline the same baseline limit cycle continued with no animals
# Both are saved to Data/followup/ for plotting (see plot_followup_add()).
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# LOOP 2: REMOVE animals (fully independent of Loop 1 -- reads only the
# saved spin-ups from disk, so it can be run on its own / in a fresh session).
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# PLOTTING: continued baseline (no animals) vs. animals-added, per
# model/scenario. Works directly from the saved Data/followup/*.rds files, so
# this can be run later without re-running either loop above.
# ------------------------------------------------------------
# Example, one model/scenario:
plot_followup_add("millennial", "RootHerbivore", by = NULL)

pdf("Plots/output.pdf", width = 8, height = 8)
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

# All models for one scenario combos that have saved add + continue_baseline runs:
for (model in models) {
  scenario = "Mite"
  if (file.exists(sprintf("Data/followup/%s_%s_add.rds", model, scenario)) &&
      file.exists(sprintf("Data/followup/%s_%s_continue_baseline.rds", model, scenario))) {
    print(plot_followup_add(model, scenario))
  }
}

# ------------------------------------------------------------
# STACKED, GROUPED animal-effect graphic (added animals): net C change over the
# follow-up, pools combined into groups and relabelled, one facet per scenario,
# annual time steps (by = 365). Auto-detects scenarios with saved add +
# continue_baseline runs.
# ------------------------------------------------------------
dir.create("Plots", showWarnings = FALSE)
p_stacked <- plot_followup_stacked("millennial", by = 365)
ggsave("Plots/followup_added_stacked.png", p_stacked,
       width = 10, height = 7, dpi = 300)
print(p_stacked)
