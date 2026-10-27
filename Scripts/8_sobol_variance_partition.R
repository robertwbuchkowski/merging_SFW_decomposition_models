#------------------------------------------------------------------------#
# VARIANCE PARTITION OF THE ANIMAL EFFECT (Morris screening -> Sobol') ----
#------------------------------------------------------------------------#
# Two-stage global sensitivity analysis of the steady-state animal effect on
# total soil + root C (treatment - baseline), per scenario:
#   1. SCREENING (sobol_screen = TRUE) with the Morris results of Scripts/4 (main AND standardized
#      ranges): keep every parameter whose mu* is >= screen_rel_cut x the
#      largest mu* in either analysis, or whose bootstrap CI reaches that cut
#      (tied), at most max_screened of them, plus ALL animal parameters.
#   OUTPUT: the animal effect as a PERCENTAGE of the no-animal soil + root C
#      ("pct", primary) and in g C m-2 ("abs"), both from the same runs. The
#      absolute effect scales with the baseline soil C, which spans orders of
#      magnitude across the parameter distributions (runaway carbon at
#      implausible combinations), so its variance is dominated by a handful of
#      extreme runs (output_kurtosis in sobol_summary.csv); the relative effect
#      is much better behaved. "rank" = Sobol' indices of the RANKS of the %
#      effect, robust to the remaining heavy tails (influence on the order of
#      the outcomes). The Morris screening (Scripts/4) must use the same effect
#      (effect_output in R/sensitivity_settings.R, default "pct"); the
#      fixed-group index checks the screening for each output.
#   2. SOBOL' INDICES of the shortlisted parameters (first order S: share of
#      the output variance explained by the parameter alone; total ST: incl.
#      all its interactions), with the parameters drawn from the MAIN
#      (knowledge-based) distributions of Scripts/5. All other parameters vary
#      together as ONE group; its total index ST_rest is the variance lost by
#      fixing them, i.e. a direct test of the screening (screening_ok when the
#      upper 95% bound of ST_rest <= rest_st_max).
#   3. PAWN (distribution-based) indices for EVERY parameter from the same runs,
#      robust to the heavy tails caused by regime shifts; a cross-check.
#   4. REALISM FILTER reported separately (option A): Sobol' indices assume
#      independent inputs, so they are computed on the unfiltered design; the
#      share of samples with a realistic baseline (targets of Scripts/5) and the
#      effect among them are reported next to them.
# Methods: Saltelli et al. (2010) first-order and Jansen (1999) total-effect
# estimators, N (k + 2) runs (k = shortlisted parameters + 1 group); 95% CIs by
# bootstrapping the N base rows; rows with an unconverged run are dropped
# from all indices (complete cases, count reported). Impossible soils are
# excluded BEFORE any simulation (design rows redrawn).
#
# FORCING: equilibria under constant mean climate forcing. The Morris
#   screening (Scripts/4) and the fitted animal parameters (Scripts/1) are read
#   from Results/, so run Scripts/1 and 4 first.
# RUN TIME: N x (k + 2) runs at ~0.4 s, k ~ 10-20: N = 1000 is
#   ~1.5-2.5 CPU-hours per scenario, divided by n_cores.
#
# Outputs (Results/):
#   sobol_shortlist.csv   screening: relative mu* per parameter and the reason
#   sobol_indices.csv     S, ST with 95% CIs per factor (incl. the fixed group),
#                         per output (pct / abs)
#   sobol_summary.csv     per scenario x output: runs, failures, kurtosis, sum S,
#                         interaction share, ST_rest, screening_ok, realism filter
#   sobol_pawn.csv        PAWN for every parameter, per output
#   sobol_runs.rds        design and model output (to recompute indices)
#   figures/sobol_indices.png (primary output), sobol_indices_<other outputs>.png,
#   figures/sobol_pawn.png, sobol_realism.png
#
# HELPER FUNCTIONS live in R/sobol_variance.R (sourced below).
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/derive_millennial_parms.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")
source("R/scenario_uncertainty.R"); source("R/morris_sensitivity.R")   # morris_cluster()
source("R/sensitivity_settings.R")    # shared: parameter list, SD rule, bounds, realism filter
source("R/sensitivity_helpers.R")     # parameter sets, ranges, model evaluation
source("R/sensitivity_uncertainty.R") # sampling distributions, realism targets
source("R/sobol_variance.R")          # screening, Sobol' design and indices, PAWN
source("R/millennial_model.R"); source("R/init_millennial_state.R")
source("R/run_utils.R")
debugonce_functions()   # only active when options(sfw.debugonce = TRUE); see Scripts/0_run_all.R

model     <- "millennial"
scen      <- read_scenarios("Data/scenarios.xlsx")
scenarios <- names(scen)

#------------------------------------------------------------------------#
# Settings ----
#------------------------------------------------------------------------#
# screening
sobol_screen       <- TRUE   # FALSE = skip the screening and analyse ALL parameters
                             # (N x (k + 2) runs with k ~ 50: ~6 CPU-hours per scenario
                             # at N = 1000 under constant forcing). Use it when the
                             # screening check fails (screening_ok = FALSE).
screen_rel_cut     <- 0.10   # keep mu* >= 10% of the largest mu* (main or standardized)
max_screened       <- 12     # at most this many non-animal parameters
include_all_animal <- TRUE   # always analyse every animal parameter
# Sobol'
sobol_N      <- 1000         # base rows; N (k + 2) model runs per scenario
sobol_seed   <- 20261020
sobol_boot_B <- 1000         # bootstrap replicates for the 95% CIs
rest_st_max  <- 0.05         # screening accepted if the fixed group's ST (upper CI) <= 5%
pawn_n_cond  <- 10           # PAWN conditioning classes
sobol_primary <- "pct"       # output for the console checks and main figures:
                             # "pct" (% of no-animal soil + root C), "rank" (ranks
                             # of the % effect; robust) or "abs" (g C m-2)
res_dir <- res_dir_default;               dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(res_dir, "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

use_fitted_params <- TRUE
fitted_file <- file.path(res_dir, "fitted_animal_params.csv")
fitted_params <- if (use_fitted_params && file.exists(fitted_file))
  load_fitted_params(fitted_file, fallback = file.path(res_dir, "animal_fit_summary_long.csv")) else NULL
if (use_fitted_params && is.null(fitted_params))
  stop("No ", fitted_file, ": run Scripts/1 first ",
       "(or set use_fitted_params <- FALSE to use the uncalibrated animal parameters).")

derive_fn <- match.fun(model_table[[model]]$derive)

#------------------------------------------------------------------------#
# Morris screening results (Scripts/4) ----
#------------------------------------------------------------------------#
read_morris <- function(tag) {
  f <- file.path(res_dir, paste0("animal_effect_morris", tag, ".csv"))
  if (file.exists(f)) read_csv(f, show_col_types = FALSE) else NULL
}
morris_main <- read_morris("")
morris_std  <- read_morris("_standardized")
if (is.null(morris_main))
  stop("No Morris results in ", res_dir, ": run Scripts/4 first.")
if (is.null(morris_std))
  warning("No standardized Morris results: the screening uses the main analysis only.")
screened_on <- if ("effect_output" %in% names(morris_main)) unique(morris_main$effect_output) else "abs"
if (!identical(screened_on, if (sobol_primary == "abs") "abs" else "pct"))
  warning("The Morris screening used the '", screened_on, "' effect but the Sobol' primary output is '",
          sobol_primary, "': rerun Scripts/4 with effect_output = \"",
          if (sobol_primary == "abs") "abs" else "pct", "\" (R/sensitivity_settings.R).")

#------------------------------------------------------------------------#
# Reported parameter uncertainty and parameter classes ----
#------------------------------------------------------------------------#
unc_all <- tryCatch(read_param_uncertainty("Data/scenarios.xlsx", warn_cv = FALSE),
                    error = function(e) { message("uncertainty read failed: ",
                                                  conditionMessage(e)); NULL })
animal_unc <- if (!is.null(unc_all))
  unc_all[unc_all$category %in% c("Animal", "Effect"), , drop = FALSE] else NULL
scenspec_unc <- if (!is.null(unc_all))
  unc_all[unc_all$category %in% c("Environment", "Productivity"), , drop = FALSE] else NULL

#------------------------------------------------------------------------#
# Parallel execution ----
#------------------------------------------------------------------------#
n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
if (isTRUE(getOption("sfw.debugonce", FALSE))) n_cores <- 1   # debugger cannot reach workers
cl <- morris_cluster(
  n_cores,
  source_files = c("R/climate_forcing.R", "R/spinup.R", "R/derive_millennial_parms.R",
                   "R/millennial_model.R", "R/init_millennial_state.R",
                   "R/setup.R", "R/compare_functions.R", "R/fit_animals.R",
                   "R/sensitivity_settings.R",
                   "R/sensitivity_helpers.R", "R/sensitivity_uncertainty.R",
                   "R/sobol_variance.R"),
  export = c("derive_fn"))
if (!is.null(cl)) cat("Running model evaluations on", length(cl), "worker processes.\n")

#------------------------------------------------------------------------#
# Sobol' analysis per scenario ----
#------------------------------------------------------------------------#
sob <- list()
for (si in seq_along(scenarios))
  sob[[scenarios[si]]] <- run_sobol_scenario(scenarios[si], si)
if (!is.null(cl)) parallel::stopCluster(cl)

shortlist <- bind_rows(lapply(sob, `[[`, "shortlist")) %>%
  mutate(parameter_label = pretty_param(parameter), .after = parameter)
indices   <- bind_rows(lapply(sob, `[[`, "indices")) %>%
  select(scenario, output, factor, label, type, is_animal, n_params, S, S_lo, S_hi,
         ST, ST_lo, ST_hi, n_used, variance)
summ      <- bind_rows(lapply(sob, `[[`, "summary"))
pawn      <- bind_rows(lapply(sob, `[[`, "pawn")) %>%
  select(scenario, output, parameter, label, shortlisted, pawn_median, pawn_max, ks_crit)

write_csv(shortlist, file.path(res_dir, "sobol_shortlist.csv"))
write_csv(indices,   file.path(res_dir, "sobol_indices.csv"))
write_csv(summ,      file.path(res_dir, "sobol_summary.csv"))
write_csv(pawn,      file.path(res_dir, "sobol_pawn.csv"))
saveRDS(lapply(sob, function(x) x[c("design", "samples", "targets")]),
        file.path(res_dir, "sobol_runs.rds"))

#------------------------------------------------------------------------#
# Checks printed to the console ----
#------------------------------------------------------------------------#
cat("\nSobol' summary (pct = % of no-animal soil C, abs = g C m-2):\n")
print(as.data.frame(summ %>% transmute(scenario, output, n_shortlisted, n_runs, n_rows_used,
                                       kurtosis = round(output_kurtosis, 1),
                                       sum_S = round(sum_S, 2),
                                       interaction_share = round(interaction_share, 2),
                                       ST_rest = round(ST_rest, 3), screening_ok,
                                       realistic_rate = round(realistic_rate, 2))),
      row.names = FALSE)
bad <- summ$scenario[summ$output == sobol_primary & summ$screening_ok %in% FALSE]
if (length(bad))
  warning("Screening NOT validated (", sobol_primary, " output) for: ", paste(bad, collapse = ", "),
          " -- the fixed parameters explain more than ", 100 * rest_st_max,
          "% of the variance; lower screen_rel_cut, raise max_screened, or set sobol_screen <- FALSE ",
          "(all parameters).")
miss <- pawn %>% filter(output == sobol_primary, !shortlisted, pawn_median > ks_crit)
if (nrow(miss)) {
  cat("\nParameters fixed by the screening but with a PAWN index above the KS critical value:\n")
  print(as.data.frame(miss %>% select(scenario, parameter, pawn_median, ks_crit)), row.names = FALSE)
}
fails <- summ %>% filter(output == sobol_primary, n_rows_used < 0.9 * N)
if (nrow(fails))
  warning("More than 10% of the Sobol' rows had an unconverged run for: ",
          paste(fails$scenario, collapse = ", "), " (indices use complete rows only).")

#------------------------------------------------------------------------#
# Figures ----
#------------------------------------------------------------------------#
ggsave(file.path(fig_dir, "sobol_indices.png"), plot_sobol_indices(indices, summ, sobol_primary),
       width = 13, height = 9, dpi = 150, bg = "white")
for (other in setdiff(c("pct", "rank", "abs"), sobol_primary))
  ggsave(file.path(fig_dir, paste0("sobol_indices_", other, ".png")),
         plot_sobol_indices(indices, summ, other), width = 13, height = 9, dpi = 150, bg = "white")
ggsave(file.path(fig_dir, "sobol_pawn.png"), plot_pawn(pawn, output_use = sobol_primary),
       width = 13, height = 9, dpi = 150, bg = "white")
ggsave(file.path(fig_dir, "sobol_realism.png"),
       plot_sobol_realism(bind_rows(lapply(sob, `[[`, "samples")), summ, sobol_primary),
       width = 10, height = 7, dpi = 150, bg = "white")
cat("\nWrote sobol_*.csv, sobol_runs.rds and figures/sobol_*.png to", res_dir, "\n")
