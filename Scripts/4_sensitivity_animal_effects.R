#------------------------------------------------------------------------#
# MORRIS SENSITIVITY OF THE ANIMAL EFFECT TO ALL PARAMETERS ----
#------------------------------------------------------------------------#
# MORRIS SENSITIVITY OF THE ANIMAL EFFECT TO ALL PARAMETERS (equilibrium)
# One combined Morris elementary-effects screening over the GENERAL model
# parameters, the SCENARIO-SPECIFIC parameters (Environment / Productivity rows
# of scenarios.xlsx) and the ANIMAL parameters (Animal / Effect rows, only in the
# scenario where that animal is active). Parameters are varied in COMBINATION
# along random Morris trajectories; per parameter it reports:
#   mu_star  mean |elementary effect|  -> TOTAL sensitivity (rank on this)
#   sigma    sd of elementary effect   -> interactions / non-linearity
# The scalar output is the equilibrium animal effect on total soil + root C
# (treatment - baseline), constant forcing, no dynamic spin-up.
#
# TWO RANGE DEFINITIONS (toggle each below):
#   MAIN (run_main_morris) - every parameter varies over its model value
#     +/- 2 SD, where SD = the SD reported in scenarios.xlsx, else a CV of
#     cv_default (0.5: SD = half the value), reduced when a reported Min/Max implies less spread
#     (labelled SD / CV-MinMax / CV-default). Same rule as the uncertainty
#     propagation in Scripts/5_uncertainty_animal_effect.R.
#   SUPPLEMENTAL (run_supp_morris) - STANDARDIZED, apples-to-apples:
#     every parameter gets [value x 0.5, value x 2], cut to a narrower
#     reported Min/Max.
#   Both are cut to the physical bounds (proportions in [0, 1] including p_a
#   and p_b, sign kept, pH, clay+silt). Every cut is reported in the
#   bound_hit column and printed per scenario.
#   Both use the SAME Morris trajectories (same seed per scenario), so any rank
#   change between them is caused by the ranges alone.
#   Shared settings (parameter list, SD rule, bounds, solver) are in
#   R/sensitivity_settings.R.
#
# Run from the project root. Outputs (Results/):
#   animal_effect_morris.csv                + figures/animal_effect_morris*.png   (main)
#   animal_effect_morris_standardized.csv   + figures/animal_effect_morris_standardized*.png
#   animal_effect_morris_main_vs_standardized.csv  (rank comparison, if both run)
#   *_elementary_effects.csv, figures/*_bootstrap.png  (bootstrap of the ranking)
#   figures/animal_effect_morris_paired.png  (current knowledge vs standardized)
#
# MORRIS BOOTSTRAP: the r trajectories are resampled with replacement (boot_B
#   times) to give 95% CIs on mu* and sigma, a rank interval, and P(top-k). A
#   scenario's ranking is called converged when every top-k parameter stays in
#   the top-k in >= 80% of replicates (otherwise a warning: raise morris_r).
#
# REGIME FLAGS: every evaluation is tagged if the treatment equilibrium has
#   (a) microbial collapse (MIC or B < collapse_rel x baseline), or
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
source("R/sensitivity_helpers.R")     # parameter sets, ranges, model evaluation
source("R/sensitivity_morris_run.R")  # run_morris(), outputs, figures
source("R/millennial_model.R"); source("R/init_millennial_state.R")
source("R/run_utils.R")
debugonce_functions()   # only active when options(sfw.debugonce = TRUE); see Scripts/0_run_all.R

model   <- "millennial"
scen    <- read_scenarios("Data/scenarios.xlsx")
scenarios <- names(scen)

use_fitted_params <- TRUE
fitted_params <- if (use_fitted_params && file.exists("Results/fitted_animal_params.csv"))
  load_fitted_params("Results/fitted_animal_params.csv") else NULL

fig_dir <- "Results/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
res_dir <- "Results";         dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

derive_fn    <- match.fun(model_table[[model]]$derive)

#------------------------------------------------------------------------#
# Morris + range settings ----
#------------------------------------------------------------------------#
morris_r      <- 200     # trajectories per scenario (raise for stable mu*/sigma)
morris_levels <- 4L     # grid levels p in the Morris design
supp_frac     <- c(0.5, 2)     # SUPPLEMENTAL standardized range = 50-200% of value (half..double)
morris_seed   <- 27082026      # base seed; each scenario uses morris_seed + its index

run_main_morris <- TRUE   # main-text analysis (existing uncertainty)
run_supp_morris <- TRUE   # supplemental analysis (standardized half-double, Min/Max-bound)

#------------------------------------------------------------------------#
# Morris bootstrap (convergence of the ranking) ----
#------------------------------------------------------------------------#
boot_B     <- 1000        # bootstrap replicates over trajectories
boot_top_k <- 5           # report P(parameter is in the top-k by mu*)

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
supp_cols <- c(`50-200%` = "#b2182b", `MinMax-bound` = "#2166ac",
               `guard-rail-bound` = "#762a83", `MinMax-excludes-default` = "#e08214")
supp_labs <- c(`50-200%` = "50% to 200% (half to double)",
               `MinMax-bound` = "50-200% cut to Min/Max",
               `guard-rail-bound` = "50-200% cut to physical bound",
               `MinMax-excludes-default` = "50-200% (Min/Max excludes default)")


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
                   "R/sensitivity_settings.R", "R/sensitivity_helpers.R"),
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
