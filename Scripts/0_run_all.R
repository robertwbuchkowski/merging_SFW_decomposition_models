#------------------------------------------------------------------------#
# Run all: regenerate the whole project ----
#------------------------------------------------------------------------#
# RUN ALL  -  regenerate the whole project end to end.
#
# Sources Scripts/1..9 in order, so a single call reproduces every saved
# output (fitted params -> spin-up states -> follow-up runs -> sensitivity).
# Run from the project root:  source("Scripts/0_run_all.R")
#
# HOW TO USE
#   * Toggle which steps run below. Later steps depend on earlier ones:
#       1 fit_all_animals        -> Results/fitted_animal_params.csv
#                                      (fit at the seasonal steady state by default;
#                                       fit_forcing in the script)
#       2 spinup_dynamic         -> Data/spinup/*.rds, Results/animal_eq_effect.csv
#       3 followup_analysis      -> Data/followup/*.rds, Plots/*
#       4 sensitivity_animal_effects -> Results/animal_effect_morris*.csv, figures
#                                      (Morris: main = SD rule; supplemental =
#                                       standardized: CV = cv_default for all)
#       5 uncertainty_animal_effect  -> Results/animal_effect_uncertainty_*.csv, figure
#                                      (Latin hypercube on the headline effect)
#       6 multistability_test        -> Results/multistability_*.csv, figure (diagnostic)
#       7 parameter_identifiability  -> Results/identifiability_*.csv, figures
#                                      (FME: which parameters artificial steady-state
#                                       and animal-introduction data can identify)
#       8 compare_forcing            -> Results/forcing_comparison/ (seasonal vs
#                                      constant: tables, figures, one PDF)
#       9 sobol_variance_partition   -> <results>/sobol_*.csv, figures/sobol_*.png
#                                      (Morris screening of step 4 -> Sobol' indices
#                                       of the shortlist + all animal parameters;
#                                       constant forcing; needs steps 1 and 4 first)
#     Steps 1, 4, 5, 7 and 9 evaluate the model at the steady state chosen by
#     forcing_run (below): "constant" mean forcing (current default) or the
#     "seasonal" annual limit cycle. Seasonal is ~50-100x slower: step 4
#     needs ~400 CPU-hours at morris_r = 100 and step 5 ~80 CPU-hours at
#     up_n = 2000 (divide by n_cores); step 7 takes ~10 min per scenario.
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

# FORCING of the steady state in steps 1, 4, 5, 7 and 9 (and where they save):
#   "seasonal" -> Results/                    (annual limit cycle)
#   "constant" -> Results/constant_forcing/   (constant mean forcing; current
#                                              default, used for the Sobol' step 9)
# To compare the two: run once with "seasonal", once with "constant" (only
# steps 1, 4, 5, 7 use it; steps 2, 3, 6 are skipped under "constant" so the
# seasonal spin-ups and follow-ups in Data/ are never overwritten), then step 8.
forcing_run <- "constant"

debug_functions         <- FALSE  # TRUE = debugonce() on every function loaded by each script
clean_env_between_steps <- TRUE   # TRUE = clear the global environment before each step

run_step1 <- TRUE    # fit animal parameters
run_step2 <- TRUE    # equilibrium + seasonal spin-up (slow)
run_step3 <- TRUE    # follow-up add/remove experiments
run_step4 <- TRUE    # Morris sensitivity of the animal effect
run_step5 <- TRUE    # uncertainty of the headline animal effect (Latin hypercube)
run_step6 <- FALSE   # multiple-stable-states test (diagnostic; off by default)
run_step7 <- TRUE    # parameter identifiability with artificial data (FME)
run_step8 <- FALSE   # compare seasonal vs constant forcing results (needs both runs)
run_step9 <- TRUE    # Sobol' variance partition after Morris screening (constant forcing)

steps <- c(
  "1" = "Scripts/1_fit_all_animals.R",
  "2" = "Scripts/2_spinup_dynamic.R",
  "3" = "Scripts/3_followup_analysis.R",
  "4" = "Scripts/4_sensitivity_animal_effects.R",
  "5" = "Scripts/5_uncertainty_animal_effect.R",
  "6" = "Scripts/6_multistability_test.R",
  "7" = "Scripts/7_parameter_identifiability.R",
  "8" = "Scripts/8_compare_forcing.R",
  "9" = "Scripts/9_sobol_variance_partition.R")
run <- c(run_step1, run_step2, run_step3, run_step4, run_step5, run_step6, run_step7,
         run_step8, run_step9)
stopifnot(forcing_run %in% c("seasonal", "constant"))
if (forcing_run == "constant" && any(run[c(2, 3, 6)])) {
  message("forcing_run = \"constant\": skipping steps 2, 3 and 6 (seasonal spin-up, ",
          "follow-up and multistability are not forcing-specific and would overwrite Data/).")
  run[c(2, 3, 6)] <- FALSE
}
options(sfw.forcing = forcing_run)          # read by R/sensitivity_settings.R (steps 1, 4, 5, 7, 9)


options(sfw.debugonce = debug_functions)   # read by debugonce_functions() in each script
keep_objs <- c("run_one", "debugonce_functions", "steps", "run", "i", "t_all",
               "debug_functions", "clean_env_between_steps", "keep_objs",
               "forcing_run", paste0("run_step", 1:9))

t_all <- Sys.time()
for (i in seq_along(steps)) if (run[i]) run_one(steps[[i]], clean = clean_env_between_steps,
                                                keep = keep_objs)
message(sprintf("\nAll requested steps finished in %.1f min.",
                as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
