#------------------------------------------------------------------------#
# MORRIS SENSITIVITY OF THE ANIMAL EFFECT TO ALL PARAMETERS ----
#------------------------------------------------------------------------#
# MORRIS SENSITIVITY OF THE ANIMAL EFFECT TO ALL PARAMETERS (steady state)
# One combined Morris elementary-effects screening over the GENERAL model
# parameters, the SCENARIO-SPECIFIC parameters (Environment / Productivity rows
# of scenarios.xlsx) and the ANIMAL parameters (Animal / Effect rows, only in the
# scenario where that animal is active). Parameters are varied in COMBINATION
# along random Morris trajectories; per parameter it reports:
#   mu_star  mean |elementary effect|  -> TOTAL sensitivity (rank on this)
#   sigma    sd of elementary effect   -> interactions / non-linearity
# The scalar output is the steady-state animal effect on total soil + root C
# (treatment - baseline), as a PERCENTAGE of the no-animal value by default
# (effect_output = "pct" in R/sensitivity_settings.R; "abs" = g C m-2). With forcing_mode = "seasonal" (R/sensitivity_settings.R)
# the steady state is the annual limit cycle under the seasonal climate forcing
# and the effect is the difference of annual means; "constant" uses the
# equilibrium under constant mean forcing (~50-100x faster, but the model is
# nonlinear, so it can differ strongly from the seasonal result).
# RUN TIME (seasonal, measured): ~10-50 s per evaluation (both arms; mean
#   ~20-35 s by scenario; ~0.4 s with constant forcing), r x (k + 1) evaluations per scenario
#   and range mode (k ~ 50 parameters). r = 100 is ~50 CPU-hours per scenario
#   and mode, i.e. ~400 CPU-hours for 4 scenarios x 2 modes, divided by n_cores.
#   Points that are impossible soils are skipped at no cost.
#
# TWO RANGE DEFINITIONS (toggle each below):
#   MAIN (run_main_morris) - every parameter varies over its model value
#     +/- 2 SD, where SD = the SD reported in scenarios.xlsx, else a CV of
#     cv_default (0.5: SD = half the value), reduced when a reported Min/Max implies less spread
#     (labelled SD / CV-MinMax / CV-default). Same rule as the uncertainty
#     propagation in Scripts/5_uncertainty_animal_effect.R.
#   SUPPLEMENTAL (run_supp_morris) - STANDARDIZED, apples-to-apples:
#     CV = cv_default for EVERY parameter (SD = cv_default x |value|),
#     ignoring reported SDs and Min/Max entirely; the range (+/- 2 SD) is cut
#     only at the physical bounds. Comparing it with MAIN shows which rankings
#     come from the reported SDs and Min/Max limits.
#   Both are cut to the physical bounds (proportions in [0, 1] including p_a
#   and p_b, sign kept, pH, clay+silt). Every cut is reported in the
#   bound_hit column and printed per scenario.
#   Both start from the same seed per scenario, so the trajectories are the
#   same unless one mode has to redraw trajectories (below).
# IMPOSSIBLE SOILS are excluded BEFORE any simulation:
#   * the seasonal moisture amplitude is varied RELATIVE to the room available,
#     theta_amp_rel = theta_amp / min(MAtheta, porosity - MAtheta) in [0, 1]
#     (moisture_amp_relative in R/sensitivity_settings.R), so moisture never
#     leaves [0, porosity] because of the amplitude;
#   * every trajectory is checked point by point (mean moisture < porosity,
#     porosity > 0; only BD, MAtheta and the amplitude matter) and redrawn if
#     any point is impossible. The count is printed and kept in n_traj_redrawn.
#   Shared settings (parameter list, SD rule, bounds, solver) are in
#   R/sensitivity_settings.R.
#
# Run from the project root. Outputs (Results/ with seasonal forcing,
# Results/constant_forcing/ with constant forcing):
#   animal_effect_morris.csv                + figures/animal_effect_morris*.png   (main)
#   animal_effect_morris_standardized.csv   + figures/animal_effect_morris_standardized*.png
#   animal_effect_morris_main_vs_standardized.csv  (rank comparison, if both run)
#   *_elementary_effects.csv, figures/*_bootstrap.png  (bootstrap of the ranking)
#   *_convergence.csv  (convergence status per scenario)
#   figures/animal_effect_morris_paired.png  (current knowledge vs standardized)
#
# MORRIS BOOTSTRAP AND CONVERGENCE: the r trajectories are resampled with
#   replacement (boot_B times) to give 95% CIs on mu* and sigma, a rank interval
#   and P(top-k). Convergence is judged per scenario as in Sarrazin, Pianosi &
#   Wagener (2016), on all and on unflagged ("clean") elementary effects:
#     index convergence     every top-k parameter has CI width <= conv_ci_tol x mu*
#     screening convergence the top-k parameters are SECURE (CI above every
#                           parameter outside the top-k); parameters whose CIs
#                           overlap across the top-k boundary are reported as
#                           TIED (near-equal mu*), not as a failure.
#   Status: converged (both) / tied (index converged, ties at the boundary) /
#   screening only (top-k set secure, CIs still wide) / unconverged (neither:
#   raise morris_r). Written to animal_effect_morris*_convergence.csv.
#
# REGIME FLAGS: every evaluation is tagged if the treatment steady state has
#   (a) microbial collapse (B < collapse_rel x baseline), or
#   (b) aggregate breakdown at the model's 0.001 floor (kb_floor, earthworms).
#   Flagged points are kept; Morris also reports mu*/sigma/rank from unflagged
#   steps only (*_clean columns).
#
# PARALLEL: model evaluations run on n_cores worker processes (PSOCK, works on
#   Windows). Set n_cores <- 1 to run serially.
#
# HELPER FUNCTIONS live in R/sensitivity_helpers.R and R/sensitivity_morris_run.R
#   (sourced below); each lists where it is used and which settings it reads.
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl, ggrepel)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/derive_millennial_parms.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")
source("R/scenario_uncertainty.R"); source("R/morris_sensitivity.R")
source("R/sensitivity_settings.R")    # shared: parameter list, SD rule, bounds, solver
source("R/seasonal_steady_state.R")  # steady state under seasonal / constant forcing
source("R/sensitivity_helpers.R")     # parameter sets, ranges, model evaluation
source("R/sensitivity_morris_run.R")  # run_morris(), outputs, figures
source("R/millennial_model.R"); source("R/init_millennial_state.R")
source("R/run_utils.R")
debugonce_functions()   # only active when options(sfw.debugonce = TRUE); see Scripts/0_run_all.R

model   <- "millennial"
scen    <- read_scenarios("Data/scenarios.xlsx")
scenarios <- names(scen)

use_fitted_params <- TRUE
# results go to a folder per forcing (R/sensitivity_settings.R): seasonal ->
# Results/, constant -> Results/constant_forcing/; the fitted animal parameters
# are read from the same folder (Scripts/1 run under the same forcing)
res_dir <- res_dir_forcing;           dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(res_dir, "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
fitted_file <- file.path(res_dir, "fitted_animal_params.csv")
fitted_params <- if (use_fitted_params && file.exists(fitted_file))
  load_fitted_params(fitted_file, fallback = file.path(res_dir, "animal_fit_summary_long.csv")) else NULL
if (use_fitted_params && is.null(fitted_params))
  stop("No ", fitted_file, ": run Scripts/1 with the same forcing first ",
       "(or set use_fitted_params <- FALSE to use the uncalibrated animal parameters).")
cat("Forcing:", forcing_mode, "| results ->", res_dir, "\n")


derive_fn    <- match.fun(model_table[[model]]$derive)

#------------------------------------------------------------------------#
# Morris + range settings ----
#------------------------------------------------------------------------#
morris_r      <- 200     # trajectories per scenario (raise for stable mu*/sigma)
morris_levels <- 4L     # grid levels p in the Morris design
morris_seed   <- 27082026      # base seed; each scenario uses morris_seed + its index

run_main_morris <- TRUE   # main-text analysis (existing uncertainty)
run_supp_morris <- TRUE   # supplemental analysis (standardized: CV = cv_default for every parameter)

#------------------------------------------------------------------------#
# Morris bootstrap (convergence of the ranking) ----
#------------------------------------------------------------------------#
boot_B     <- 1000        # bootstrap replicates over trajectories
boot_top_k <- 5           # size of the "most influential" set checked for convergence
conv_ci_tol <- 0.2        # index convergence: 95% CI width of mu* <= 0.2 x mu* for
                          # every top-k parameter (Sarrazin et al. 2016)

paired_n_top <- 8          # parameters per scenario in the paired Morris figure

#------------------------------------------------------------------------#
# Parallel execution ----
#------------------------------------------------------------------------#
# Model evaluations (Morris trajectories) run on n_cores
# worker processes; 1 = serial. The bootstrap itself is vectorised and takes
# seconds, so it stays on the main process.
n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
if (isTRUE(getOption("sfw.debugonce", FALSE))) n_cores <- 1   # debugger cannot reach workers

#------------------------------------------------------------------------#
# Reported parameter uncertainty ----
#------------------------------------------------------------------------#
# Reported uncertainty from scenarios.xlsx. Parameter ranges are built from it by
# range_for() in R/sensitivity_helpers.R (main and standardized modes).
unc_all <- tryCatch(read_param_uncertainty("Data/scenarios.xlsx", warn_cv = FALSE),
                    error = function(e) { message("uncertainty read failed: ",
                                                  conditionMessage(e)); NULL })


# Parameter classes from scenarios.xlsx:
#   * ANIMAL            - Category Animal/Effect (only where that animal is active)
#   * SCENARIO-SPECIFIC - Category Environment/Productivity (non-animal, but with
#                         per-scenario values/uncertainty in the sheet)
# Everything else in sweep_params is a GENERAL model parameter.
animal_unc <- if (!is.null(unc_all))
  unc_all[unc_all$category %in% c("Animal", "Effect"), , drop = FALSE] else NULL
scenspec_unc <- if (!is.null(unc_all))
  unc_all[unc_all$category %in% c("Environment", "Productivity"), , drop = FALSE] else NULL


#------------------------------------------------------------------------#
# Figure colours and labels ----
#------------------------------------------------------------------------#
# Colour / legend labels by range source, used by write_morris_outputs().
main_cols <- c(SD = "#1b7837", `CV-MinMax` = "#2166ac", `CV-default` = "#b2182b")
main_labs <- c(SD = "Reported SD", `CV-MinMax` = "CV reduced to fit Min/Max",
               `CV-default` = paste0("CV = ", cv_default, " (no SD or Min/Max)"))
supp_cols <- main_cols["CV-default"]
supp_labs <- c(`CV-default` = paste0("CV = ", cv_default, " (all parameters)"))


#------------------------------------------------------------------------#
# RUN ----
#------------------------------------------------------------------------#
# start the worker pool (NULL = serial). Workers source the R/ helpers (so they
# have every function) and receive the script-level settings the model
# evaluations read.
cl <- morris_cluster(
  n_cores,
  source_files = c("R/climate_forcing.R", "R/spinup.R", "R/derive_millennial_parms.R",
                   "R/millennial_model.R", "R/init_millennial_state.R",
                   "R/setup.R", "R/compare_functions.R", "R/fit_animals.R",
                   "R/sensitivity_settings.R", "R/seasonal_steady_state.R",
                   "R/sensitivity_helpers.R"),
  export = c("derive_fn"))
if (!is.null(cl)) cat("Running model evaluations on", length(cl), "worker processes.\n")

morris_main <- morris_supp <- NULL
if (run_main_morris) {
  morris_main <- run_morris("main")
  write_morris_outputs(morris_main, "", main_cols, main_labs)
}
if (run_supp_morris) {
  morris_supp <- run_morris("standardized")
  write_morris_outputs(morris_supp, "_standardized", supp_cols, supp_labs)
}

#------------------------------------------------------------------------#
# MAIN vs STANDARDIZED comparison ----
#------------------------------------------------------------------------#
# MAIN vs STANDARDIZED comparison: per-parameter mu* and rank side by side, and
# a per-scenario Spearman rank correlation of mu*. A parameter whose rank
# drops a lot under the standardized range was important mainly because its
# reported uncertainty is wide (and vice versa).
if (!is.null(morris_main) && !is.null(morris_supp)) {
  cmp <- full_join(
    morris_main %>% select(scenario, parameter, parameter_label, param_type,
                           source_main = range_source, mu_star_main = mu_star,
                           sigma_main = sigma, rank_main = rank),
    morris_supp %>% select(scenario, parameter,
                           source_std = range_source, mu_star_std = mu_star,
                           sigma_std = sigma, rank_std = rank),
    by = c("scenario", "parameter")) %>%
    mutate(rank_change = rank_std - rank_main) %>%
    arrange(scenario, rank_main)
  write_csv(cmp, file.path(res_dir, "animal_effect_morris_main_vs_standardized.csv"))

  rho <- cmp %>% filter(is.finite(mu_star_main), is.finite(mu_star_std)) %>%
    group_by(scenario) %>%
    summarise(spearman_mu_star = cor(mu_star_main, mu_star_std, method = "spearman"),
              n_params = n(), .groups = "drop")
  cat("\nRank agreement of mu* (main vs standardized):\n")
  print(as.data.frame(rho), row.names = FALSE)
  cat("Wrote:", file.path(res_dir, "animal_effect_morris_main_vs_standardized.csv"), "\n")
}


#------------------------------------------------------------------------#
# Paired Morris ranking: current knowledge vs standardized ----
#------------------------------------------------------------------------#
# Compares the two rankings with bootstrap intervals (plot_morris_paired() in
# R/sensitivity_morris_run.R). Uses this session's runs, or the saved CSVs when
# the Morris runs were switched off.
paired_main <- if (!is.null(morris_main)) morris_main else {
  f <- file.path(res_dir, "animal_effect_morris.csv"); if (file.exists(f)) read_csv(f, show_col_types = FALSE) }
paired_std <- if (!is.null(morris_supp)) morris_supp else {
  f <- file.path(res_dir, "animal_effect_morris_standardized.csv"); if (file.exists(f)) read_csv(f, show_col_types = FALSE) }
if (!is.null(paired_main) && !is.null(paired_std)) {
  p_paired <- plot_morris_paired(paired_main, paired_std, n_top = paired_n_top)
  ggsave(file.path(fig_dir, "animal_effect_morris_paired.png"), p_paired,
         width = 11, height = 8, dpi = 300)
  cat("Wrote:", file.path(fig_dir, "animal_effect_morris_paired.png"), "\n")
}


# shut down the worker pool
if (!is.null(cl)) parallel::stopCluster(cl)
