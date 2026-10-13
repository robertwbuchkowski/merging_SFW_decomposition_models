#------------------------------------------------------------------------#
# PARAMETER IDENTIFIABILITY WITH DIFFERENT TYPES OF DATA (FME) ----
#------------------------------------------------------------------------#
# Which parameters could be estimated from which kind of field data?
# Local identifiability analysis (Brun, Reichert & Kuensch 2001) with the FME
# package (sensFun + collin), on ARTIFICIAL data produced by the model itself
# at the current parameter values (fitted animal parameters applied).
#
# ARTIFICIAL DATA TYPES (every model pool is "measured" unless obs_pools is set):
#   ss_noanimal  steady-state pools without the animal (constant mean forcing)
#   ss_animal    steady-state pools with the animal, incl. animal biomass
#   ts           animal INTRODUCTION: the animal is added to the no-animal soil
#                (at ts_seed_mod x its input biomass) and all pools incl. the
#                animal are measured once per year for ts_years (= 3) years.
#                ts_forcing = "seasonal" starts from the no-animal seasonal
#                limit cycle and samples on the same calendar day each year
#                (realistic, slower); "constant" uses constant mean forcing.
# DATA DESIGNS (combinations of data types; edit `designs` below).
# PARAMETER SETS:
#   all     every parameter of the Morris screening (general, scenario-
#           specific and animal) is to be estimated from the data
#   animal  only the animal parameters are estimated; all other parameters
#           are assumed known (fixed at their values)
#
# METHOD: FME::sensFun gives the normalized sensitivities
#   S_ij = (dy_i/dp_j) (p_j / y_i) of every observation to every parameter (forward
#   differences, relative step fd_rel). A parameter carries information in a
#   design when its RMS sensitivity (FME's L2) is >= min_sens. FME::collin gives
#   the collinearity index gamma of a parameter set (gamma > ident_threshold:
#   the parameters can compensate each other and cannot be estimated jointly).
#   The identifiable subset is built greedily: parameters that carry
#   information are added in order of decreasing sensitivity while gamma stays
#   <= ident_threshold. For small sets (<= exhaustive_max parameters) every
#   subset is also tested and the largest identifiable subset size reported.
#   This is LOCAL identifiability at the current parameter values; relative
#   sensitivities imply equal relative measurement error for all observations.
#
# Run from the project root (needs Results/fitted_animal_params.csv from
# Scripts/1). Outputs (Results/):
#   identifiability_summary.csv        one row per scenario x parameter set x design
#   identifiability_parameters.csv     status of every parameter (identifiable /
#                                      collinear / no information), L1, L2
#   identifiability_collinear_pairs.csv  most collinear parameter pairs
#   identifiability_sensitivities.rds  the sensFun objects (to redo the
#                                      collinearity step without re-solving)
#   figures/identifiability_summary.png, identifiability_animal_params.png,
#   figures/identifiability_all_params.png
#
# HELPER FUNCTIONS live in R/identifiability.R (sourced below).
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl, FME)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/derive_millennial_parms.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")
source("R/scenario_uncertainty.R"); source("R/morris_sensitivity.R")   # morris_cluster()
source("R/sensitivity_settings.R")    # shared: parameter list, SD rule, bounds, solver
source("R/sensitivity_helpers.R")     # parameter sets, make_base_pair(), infeasible_soil()
source("R/identifiability.R")         # artificial data, sensitivities, collinearity
source("R/millennial_model.R"); source("R/init_millennial_state.R")
source("R/run_utils.R")
debugonce_functions()   # only active when options(sfw.debugonce = TRUE); see Scripts/0_run_all.R

model     <- "millennial"
scen      <- read_scenarios("Data/scenarios.xlsx")
scenarios <- names(scen)

use_fitted_params <- TRUE
fitted_params <- if (use_fitted_params && file.exists("Results/fitted_animal_params.csv"))
  load_fitted_params("Results/fitted_animal_params.csv") else NULL

fig_dir <- "Results/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
res_dir <- "Results";         dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

derive_fn <- match.fun(model_table[[model]]$derive)

#------------------------------------------------------------------------#
# Artificial data and identifiability settings ----
#------------------------------------------------------------------------#
ident_settings <- list(
  ts_years    = 3,           # introduction time series: measured once per year for 3 years
  ts_forcing  = "seasonal",  # "seasonal" (from the no-animal limit cycle) or "constant"
  ts_seed_mod = 0.01,        # animals introduced at 1% of their input biomass
                             # (same as followup_add_animals() in Scripts/3)
  obs_pools   = NULL,        # NULL = every pool is measured; or e.g.
                             # c("P", "M", "MIC", "B", "Earthworm", "Detritivore", "RootHerb")
  obs_floor   = 1e-9,        # observations below this (numerically zero) are dropped
  fd_rel      = 1e-3,        # relative parameter step for the finite differences
  eq_max_time = eq_max_time, # equilibrium solver horizon (R/sensitivity_settings.R)
  eq_stol     = 1e-10,       # equilibrium derivative tolerance (tight: FD needs precision)
  lc_tol      = 1e-7         # Newton limit-cycle tolerance (seasonal time series; error
                             # in the sensitivities ~ lc_tol / fd_rel)
)

# data designs: which data types are used together
designs <- list(
  "Steady state, no animal"         = "ss_noanimal",
  "Steady state, with animal"       = "ss_animal",
  "Steady state, both"              = c("ss_noanimal", "ss_animal"),
  "Introduction, years 1-3"         = "ts",
  "No-animal steady state + introduction" = c("ss_noanimal", "ts"),
  "All data"                        = c("ss_noanimal", "ss_animal", "ts"))

param_sets      <- c("all", "animal")  # (1) all parameters, (2) animal parameters only
ident_threshold <- 20      # collinearity index above which a set is not identifiable
                           # (Brun et al. 2001: 10-15 critical; 20 common)
min_sens        <- 1e-3    # minimum RMS normalized sensitivity to carry information
exhaustive_max  <- 12      # test all subsets when at most this many parameters

#------------------------------------------------------------------------#
# Parallel execution ----
#------------------------------------------------------------------------#
# Scenarios are processed in parallel (one worker per scenario).
n_cores <- min(length(scenarios), max(1, parallel::detectCores(logical = FALSE) - 1))
if (isTRUE(getOption("sfw.debugonce", FALSE))) n_cores <- 1   # debugger cannot reach workers

#------------------------------------------------------------------------#
# Parameter sets (same as the Morris screening) ----
#------------------------------------------------------------------------#
unc_all <- tryCatch(read_param_uncertainty("Data/scenarios.xlsx", warn_cv = FALSE),
                    error = function(e) { message("uncertainty read failed: ",
                                                  conditionMessage(e)); NULL })
animal_unc <- if (!is.null(unc_all))
  unc_all[unc_all$category %in% c("Animal", "Effect"), , drop = FALSE] else NULL
scenspec_unc <- if (!is.null(unc_all))
  unc_all[unc_all$category %in% c("Environment", "Productivity"), , drop = FALSE] else NULL

#------------------------------------------------------------------------#
# Sensitivities of the artificial data (FME::sensFun) ----
#------------------------------------------------------------------------#
# One nominal run plus one run per parameter, per scenario. With seasonal
# introduction data each run includes a Newton limit-cycle solve.
cat("Sensitivities:", length(scenarios), "scenarios on", n_cores, "worker(s); time series:",
    ident_settings$ts_forcing, "forcing\n")
one_scenario <- function(s, settings) tryCatch(ident_sensitivity(s, settings),
  error = function(e) { message("identifiability failed for ", s, ": ", conditionMessage(e)); NULL })
cl <- morris_cluster(
  n_cores,
  source_files = c("R/climate_forcing.R", "R/spinup.R", "R/derive_millennial_parms.R",
                   "R/millennial_model.R", "R/init_millennial_state.R",
                   "R/setup.R", "R/compare_functions.R", "R/fit_animals.R",
                   "R/sensitivity_settings.R", "R/sensitivity_helpers.R",
                   "R/identifiability.R"),
  packages = c("deSolve", "rootSolve", "yaml", "FME"),
  export = c("derive_fn", "model", "scen", "fitted_params", "unc_all", "animal_unc",
             "scenspec_unc"))
sens_list <- if (is.null(cl)) lapply(scenarios, one_scenario, settings = ident_settings) else
  parallel::parLapplyLB(cl, scenarios, one_scenario, settings = ident_settings)
if (!is.null(cl)) parallel::stopCluster(cl)
names(sens_list) <- scenarios
sens_list <- Filter(Negate(is.null), sens_list)
if (!length(sens_list)) stop("Sensitivities failed for every scenario (see messages above).")
for (sr in sens_list)
  cat(sprintf("  %-14s %d parameters (%d animal), %d observations, %.1f min%s\n",
              sr$scenario, length(sr$ptype), length(sr$animal), nrow(sr$obs), sr$time_s / 60,
              if (length(sr$skipped_zero)) paste0("; skipped (value 0): ",
                                                  paste(sr$skipped_zero, collapse = ", ")) else ""))
saveRDS(list(sens = sens_list, settings = ident_settings),
        file.path(res_dir, "identifiability_sensitivities.rds"))

#------------------------------------------------------------------------#
# Collinearity and identifiable subsets (FME::collin) ----
#------------------------------------------------------------------------#
idf <- run_identifiability(sens_list, designs, param_sets, threshold = ident_threshold,
                           min_sens = min_sens, exhaustive_max = exhaustive_max)
write_csv(idf$summary, file.path(res_dir, "identifiability_summary.csv"))
write_csv(idf$params,  file.path(res_dir, "identifiability_parameters.csv"))
write_csv(idf$pairs,   file.path(res_dir, "identifiability_collinear_pairs.csv"))

cat("\nIdentifiability (collinearity index threshold ", ident_threshold, "):\n", sep = "")
print(as.data.frame(idf$summary %>%
  transmute(scenario, param_set, design, n_params, n_informative,
            gamma_all = round(collinearity_all, 1), all_identifiable,
            n_ident = n_identifiable_greedy, n_ident_exh = n_identifiable_exhaustive)),
  row.names = FALSE)

cat("\nAnimal parameters, identifiable subset per design:\n")
print(as.data.frame(idf$summary %>% filter(param_set == "animal") %>%
  select(scenario, design, identifiable_greedy, no_information)), row.names = FALSE)

#------------------------------------------------------------------------#
# Figures ----
#------------------------------------------------------------------------#
p_sum <- plot_ident_summary(idf$summary, names(designs))
ggsave(file.path(fig_dir, "identifiability_summary.png"), p_sum, width = 13, height = 7, dpi = 150, bg = "white")
p_an <- plot_ident_params(idf$params, names(designs), set = "animal")
ggsave(file.path(fig_dir, "identifiability_animal_params.png"), p_an, width = 11, height = 7, dpi = 150, bg = "white")
p_all <- plot_ident_params(idf$params, names(designs), set = "all")
ggsave(file.path(fig_dir, "identifiability_all_params.png"), p_all, width = 13, height = 14, dpi = 150, bg = "white")
if (interactive()) print(p_sum)
cat("\nWrote Results/identifiability_*.csv and figures/identifiability_*.png\n")
