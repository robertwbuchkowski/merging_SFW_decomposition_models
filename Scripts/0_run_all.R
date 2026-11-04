#------------------------------------------------------------------------#
# Run all: regenerate the whole project ----
#------------------------------------------------------------------------#
# RUN ALL  -  regenerate the whole project end to end.
#
# Sources Scripts/1..7 in order, so a single call reproduces every saved
# output (fitted params -> spin-up states -> invasion runs -> uncertainty ->
# identifiability -> Sobol' -> multistability).
# Run from the project root:  source("Scripts/0_run_all.R")
#
# HOW TO USE
#   * Toggle which steps run below. Later steps depend on earlier ones:
#       1 fit_all_animals        -> Results/fitted_animal_params.csv
#                                      (feeding rates fitted at the constant-forcing
#                                       equilibrium)
#       2 animal_equilibrium_effect -> Results/animal_eq_effect.csv,
#                                      Results/animal_effect_summary.csv,
#                                      Results/figures/flux_network_*.png
#                                      (equilibrium animal effects at default
#                                       parameters; flux networks per scenario +
#                                       all-scenario change figure)
#       3 followup_analysis      -> Data/spinup/*.rds, Data/followup/*.rds,
#                                      Results/figures/followup_*
#                                      (SEASONAL spin-up of the no-animal soil to
#                                       its annual limit cycle, then animal
#                                       INVASION, ~50-yr trajectories)
#       4 uncertainty_animal_effect  -> Results/animal_effect_uncertainty_*.csv,
#                                      Results/animal_effect_table.csv (paper table), figure
#                                      (Latin hypercube on the headline effect)
#       5 parameter_identifiability  -> Results/identifiability_*.csv, figures
#                                      (FME: which parameters artificial steady-state
#                                       and animal-introduction data can identify)
#       6 sobol_variance_partition   -> Results/sobol_*.csv, figures/sobol_*.png
#                                      (Sobol' indices of EVERY parameter, no
#                                       screening; needs step 1 first)
#       7 multistability_test        -> Results/multistability_*.csv, figure (diagnostic)
#   * FIGURES: every figure goes to Results/figures/.
#   * PARAMETER RANGES (steps 4 and 6): SD reported in scenarios.xlsx, else
#     CV = cv_default (0.25) in R/sensitivity_settings.R.
#   * FORCING: steps 1, 2 and 4-7 use the EQUILIBRIUM under constant mean
#     climate forcing. Only step 3 simulates the seasonal dynamics (the
#     seasonal spin-up and the animal-invasion trajectories), with the feeding
#     rates fitted in step 1.
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
#   * Run time: steps 1-3 take ~10 min in total; steps 4 and 6 are the long
#     ones (step 6, Sobol' on every parameter, ~5-6 CPU-hours per scenario
#     at N = 1000; see its header).

stopifnot(file.exists("R/setup.R"))        # guard: must be run from project root
source("R/run_utils.R")                     # run_one(), debugonce_functions()

debug_functions         <- FALSE  # TRUE = debugonce() on every function loaded by each script
clean_env_between_steps <- TRUE   # TRUE = clear the global environment before each step

run_step1 <- TRUE    # fit animal parameters (constant-forcing equilibrium)
run_step2 <- TRUE    # equilibrium animal effects + flux networks (constant forcing)
run_step3 <- TRUE    # seasonal spin-up + animal invasion trajectories
run_step4 <- TRUE    # uncertainty of the headline animal effect (Latin hypercube)
run_step5 <- TRUE    # parameter identifiability with artificial data (FME)
run_step6 <- TRUE    # Sobol' variance partition (all parameters)
run_step7 <- FALSE   # multiple-stable-states test (diagnostic; off by default)

steps <- c(
  "1" = "Scripts/1_fit_all_animals.R",
  "2" = "Scripts/2_animal_equilibrium_effect.R",
  "3" = "Scripts/3_followup_analysis.R",
  "4" = "Scripts/4_uncertainty_animal_effect.R",
  "5" = "Scripts/5_parameter_identifiability.R",
  "6" = "Scripts/6_sobol_variance_partition.R",
  "7" = "Scripts/7_multistability_test.R")
run <- c(run_step1, run_step2, run_step3, run_step4, run_step5, run_step6, run_step7)

options(sfw.debugonce = debug_functions)   # read by debugonce_functions() in each script
keep_objs <- c("run_one", "debugonce_functions", "steps", "run", "i", "t_all",
               "debug_functions", "clean_env_between_steps", "keep_objs",
               paste0("run_step", 1:7))

t_all <- Sys.time()
for (i in seq_along(steps)) if (run[i]) run_one(steps[[i]], clean = clean_env_between_steps,
                                                keep = keep_objs)
message(sprintf("\nAll requested steps finished in %.1f min.",
                as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
