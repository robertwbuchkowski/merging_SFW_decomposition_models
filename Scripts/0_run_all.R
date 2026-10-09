#------------------------------------------------------------------------#
# Run all: regenerate the whole project ----
#------------------------------------------------------------------------#
# RUN ALL  -  regenerate the whole project end to end.
#
# Sources Scripts/1..6 in order, so a single call reproduces every saved
# output (fitted params -> spin-up states -> follow-up runs -> sensitivity).
# Run from the project root:  source("Scripts/0_run_all.R")
#
# HOW TO USE
#   * Toggle which steps run below. Later steps depend on earlier ones:
#       1 fit_all_animals        -> Results/fitted_animal_params.csv
#       2 spinup_dynamic         -> Data/spinup/*.rds, Results/animal_eq_effect.csv
#       3 followup_analysis      -> Data/followup/*.rds, Plots/*
#       4 sensitivity_animal_effects -> Results/animal_effect_morris*.csv, figures
#                                      (Morris: main = SD rule; supplemental =
#                                       standardized: CV = cv_default for all)
#       5 uncertainty_animal_effect  -> Results/animal_effect_uncertainty_*.csv, figure
#                                      (Latin hypercube on the headline effect)
#       6 multistability_test        -> Results/multistability_*.csv, figure (diagnostic)
#   * Each script sources the R/ helpers it needs and is self-contained. This
#     wrapper runs each one in the GLOBAL environment (as if you ran it by
#     hand), optionally clearing it first so every step starts fresh, and
#     times each step. All custom functions live in R/ and are sourced at the
#     top of each script.
#   * Step-through debugging: set debug_functions <- TRUE. Every function the
#     script has just sourced is flagged with debugonce(), so running the steps
#     top to bottom opens the browser the first time each function is called.
#     Model evaluations then run serially (no parallel workers), so the
#     debugger can reach them. Set back to FALSE for normal runs.
#   * The slow step is 2 (spin-up). Set run_step2 <- FALSE to reuse saved
#     spin-ups when you only need to refresh downstream steps.

stopifnot(file.exists("R/setup.R"))        # guard: must be run from project root
source("R/run_utils.R")                     # run_one(), debugonce_functions()

debug_functions         <- FALSE  # TRUE = debugonce() on every function loaded by each script
clean_env_between_steps <- TRUE   # TRUE = clear the global environment before each step

run_step1 <- TRUE    # fit animal parameters
run_step2 <- TRUE    # equilibrium + seasonal spin-up (slow)
run_step3 <- TRUE    # follow-up add/remove experiments
run_step4 <- TRUE    # Morris sensitivity of the animal effect
run_step5 <- FALSE    # uncertainty of the headline animal effect (Latin hypercube)
run_step6 <- FALSE   # multiple-stable-states test (diagnostic; off by default)

steps <- c(
  "1" = "Scripts/1_fit_all_animals.R",
  "2" = "Scripts/2_spinup_dynamic.R",
  "3" = "Scripts/3_followup_analysis.R",
  "4" = "Scripts/4_sensitivity_animal_effects.R",
  "5" = "Scripts/5_uncertainty_animal_effect.R",
  "6" = "Scripts/6_multistability_test.R")
run <- c(run_step1, run_step2, run_step3, run_step4, run_step5, run_step6)


options(sfw.debugonce = debug_functions)   # read by debugonce_functions() in each script
keep_objs <- c("run_one", "debugonce_functions", "steps", "run", "i", "t_all",
               "debug_functions", "clean_env_between_steps", "keep_objs",
               paste0("run_step", 1:6))

t_all <- Sys.time()
for (i in seq_along(steps)) if (run[i]) run_one(steps[[i]], clean = clean_env_between_steps,
                                                keep = keep_objs)
message(sprintf("\nAll requested steps finished in %.1f min.",
                as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
