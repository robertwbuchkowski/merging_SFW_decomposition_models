#------------------------------------------------------------------------#
# Run all: regenerate the whole project ----
#------------------------------------------------------------------------#
# RUN ALL  -  regenerate the whole project end to end.
#
# Sources Scripts/1..8 in order, so a single call reproduces every saved
# output (fitted params -> spin-up states -> invasion runs -> sensitivity).
# Run from the project root:  source("Scripts/0_run_all.R")
#
# HOW TO USE
#   * Toggle which steps run below. Later steps depend on earlier ones:
#       1 fit_all_animals        -> Results/fitted_animal_params.csv
#                                      (feeding rates fitted at the constant-forcing
#                                       equilibrium)
#       2 spinup_dynamic         -> Data/spinup/*.rds, Results/animal_eq_effect.csv,
#                                      Plots/flux_network_*.png
#                                      (equilibrium effects, then the SEASONAL
#                                       spin-up of the no-animal soil to its
#                                       annual limit cycle)
#       3 followup_analysis      -> Data/followup/*.rds, Plots/*
#                                      (animal INVASION: animals added to the
#                                       seasonal no-animal soil, ~50-yr trajectories)
#       4 sensitivity_animal_effects -> Results/animal_effect_morris*.csv, figures
#                                      (Morris: main = SD rule; supplemental =
#                                       standardized: CV = cv_default for all)
#       5 uncertainty_animal_effect  -> Results/animal_effect_uncertainty_*.csv, figure
#                                      (Latin hypercube on the headline effect)
#       6 multistability_test        -> Results/multistability_*.csv, figure (diagnostic)
#       7 parameter_identifiability  -> Results/identifiability_*.csv, figures
#                                      (FME: which parameters artificial steady-state
#                                       and animal-introduction data can identify)
#       8 sobol_variance_partition   -> Results/sobol_*.csv, figures/sobol_*.png
#                                      (Morris screening of step 4 -> Sobol' indices;
#                                       needs steps 1 and 4 first)
#   * FORCING: steps 1 and 4-8 use the EQUILIBRIUM under constant mean climate
#     forcing. Only steps 2 and 3 simulate the seasonal dynamics (the seasonal
#     spin-up and the animal-invasion trajectories), with the feeding rates
#     fitted in step 1.
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
#   * Run time: steps 1-3 take ~10 min in total; the sensitivity steps 4, 5
#     and 8 are the long ones (see their headers).

stopifnot(file.exists("R/setup.R"))        # guard: must be run from project root
source("R/run_utils.R")                     # run_one(), debugonce_functions()

debug_functions         <- FALSE  # TRUE = debugonce() on every function loaded by each script
clean_env_between_steps <- TRUE   # TRUE = clear the global environment before each step

run_step1 <- TRUE    # fit animal parameters (constant-forcing equilibrium)
run_step2 <- TRUE    # equilibrium effects + seasonal spin-up of the no-animal soil
run_step3 <- TRUE    # animal invasion trajectories (seasonal)
run_step4 <- TRUE    # Morris sensitivity of the animal effect
run_step5 <- TRUE    # uncertainty of the headline animal effect (Latin hypercube)
run_step6 <- FALSE   # multiple-stable-states test (diagnostic; off by default)
run_step7 <- TRUE    # parameter identifiability with artificial data (FME)
run_step8 <- TRUE    # Sobol' variance partition after Morris screening

steps <- c(
  "1" = "Scripts/1_fit_all_animals.R",
  "2" = "Scripts/2_spinup_dynamic.R",
  "3" = "Scripts/3_followup_analysis.R",
  "4" = "Scripts/4_sensitivity_animal_effects.R",
  "5" = "Scripts/5_uncertainty_animal_effect.R",
  "6" = "Scripts/6_multistability_test.R",
  "7" = "Scripts/7_parameter_identifiability.R",
  "8" = "Scripts/8_sobol_variance_partition.R")
run <- c(run_step1, run_step2, run_step3, run_step4, run_step5, run_step6, run_step7,
         run_step8)

options(sfw.debugonce = debug_functions)   # read by debugonce_functions() in each script
keep_objs <- c("run_one", "debugonce_functions", "steps", "run", "i", "t_all",
               "debug_functions", "clean_env_between_steps", "keep_objs",
               paste0("run_step", 1:8))

t_all <- Sys.time()
for (i in seq_along(steps)) if (run[i]) run_one(steps[[i]], clean = clean_env_between_steps,
                                                keep = keep_objs)
message(sprintf("\nAll requested steps finished in %.1f min.",
                as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
