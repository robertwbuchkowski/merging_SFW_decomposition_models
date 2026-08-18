# ============================================================
# RUN ALL  -  regenerate the whole project end to end.
#
# Sources Scripts/1..4 in order, so a single call reproduces every saved
# output (fitted params -> spin-up states -> follow-up runs -> sensitivity).
# Run from the project root:  source("Scripts/0_run_all.R")
#
# HOW TO USE
#   * Toggle which steps run below. Later steps depend on earlier ones:
#       1 fit_all_animals        -> Results/fitted_animal_params.csv
#       2 spinup_dynamic         -> Data/spinup/*.rds, Results/animal_eq_effect.csv
#       3 followup_analysis      -> Data/followup/*.rds, Plots/*
#       4 sensitivity_animal_effects -> Results/animal_effect_sensitivity*.csv, figures
#       5 animal_param_uncertainty    -> Results/animal_param_uncertainty_effect.csv, figure
#   * Each script sources the same R/ helpers and is self-contained; this
#     wrapper just runs them in a clean environment and times each one.
#   * The slow step is 2 (spin-up). Set run_step2 <- FALSE to reuse saved
#     spin-ups when you only need to refresh downstream steps.
# ============================================================

stopifnot(file.exists("R/setup.R"))        # guard: must be run from project root

run_step1 <- TRUE    # fit animal parameters
run_step2 <- TRUE    # equilibrium + seasonal spin-up (slow)
run_step3 <- TRUE    # follow-up add/remove experiments
run_step4 <- TRUE    # parameter sensitivity of the animal effect
run_step5 <- TRUE    # range of animal effects across animal-parameter uncertainty

steps <- c(
  "1" = "Scripts/1_fit_all_animals.R",
  "2" = "Scripts/2_spinup_dynamic.R",
  "3" = "Scripts/3_followup_analysis.R",
  "4" = "Scripts/4_sensitivity_animal_effects.R",
  "5" = "Scripts/5_animal_param_uncertainty.R")
run <- c(run_step1, run_step2, run_step3, run_step4, run_step5)

run_one <- function(path) {
  message("\n", strrep("=", 60), "\n== RUN: ", path, "\n", strrep("=", 60))
  t0 <- Sys.time()
  # each script runs in its own environment so their globals don't collide
  sys.source(path, envir = new.env(parent = globalenv()))
  message(sprintf("== DONE: %s  (%.1f min)", path,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

t_all <- Sys.time()
for (i in seq_along(steps)) if (run[i]) run_one(steps[[i]])
message(sprintf("\nAll requested steps finished in %.1f min.",
                as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
