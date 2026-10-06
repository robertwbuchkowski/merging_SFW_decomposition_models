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
#   MAIN (run_main_morris) - reflects EXISTING uncertainty:
#     * reported range from scenarios.xlsx: value +/- 2 SD, else [Min, Max];
#     * otherwise the buffer range [value x 0.5, value x 2] (half to double).
#   SUPPLEMENTAL (run_supp_morris) - STANDARDIZED, apples-to-apples:
#     * every parameter gets [value x 0.5, value x 2] (50-200%, half to double);
#     * if scenarios.xlsx gives a Min/Max that is narrower, the range is cut
#       to that Min/Max ("Min/Max binds");
#   Both are clamped to the physical guard rails (efficiencies [0,1], etc.).
#   Both use the SAME Morris trajectories (same seed per scenario), so any rank
#   change between them is caused by the ranges alone.
#
# Run from the project root. Outputs (Results/):
#   animal_effect_morris.csv                + figures/animal_effect_morris*.png   (main)
#   animal_effect_morris_standardized.csv   + figures/animal_effect_morris_standardized*.png
#   animal_effect_morris_main_vs_standardized.csv  (rank comparison, if both run)
#   *_elementary_effects.csv, figures/*_bootstrap.png  (bootstrap of the ranking)
#
# MORRIS BOOTSTRAP: the r trajectories are resampled with replacement (boot_B
#   times) to give 95% CIs on mu* and sigma, a rank interval, and P(top-k). A
#   scenario's ranking is called converged when every top-k parameter stays in
#   the top-k in >= 80% of replicates (otherwise a warning: raise morris_r).
#
# UNCERTAINTY PROPAGATION (run_uncertainty_propagation): all screened
#   parameters sampled jointly by Latin hypercube (truncated normal for SD,
#   triangular for Min/Max, log-normal for unmeasured parameters); the total
#   and direct equilibrium animal effects are recomputed per sample, and only
#   samples whose baseline stays near realistic stocks are ACCEPTED (baseline
#   realism filter; targets editable in Data/baseline_targets.csv).
#   Outputs: animal_effect_uncertainty_*.csv, baseline_targets_used.csv,
#            figures/animal_effect_uncertainty.png
#
# PAIRED MORRIS FIGURE: figures/animal_effect_morris_paired.png compares the
#   current-knowledge and standardized rankings with bootstrap intervals.
#
# REGIME FLAGS: every evaluation is tagged if the treatment equilibrium has
#   (a) microbial collapse (MIC or B < collapse_rel x baseline), or
#   (b) aggregate breakdown at the model's 0.001 floor (kb_floor, earthworms).
#   Flagged points are kept; Morris also reports mu*/sigma/rank from unflagged
#   steps only (*_clean columns) and the uncertainty summary gives intervals
#   with and without flagged samples.
#
# PARALLEL: model evaluations run on n_cores worker processes (PSOCK, works on
#   Windows). Set n_cores <- 1 to run serially.
#
# HELPER FUNCTIONS live in R/sensitivity_helpers.R, R/sensitivity_morris_run.R
#   and R/sensitivity_uncertainty.R (sourced below); each lists where it is
#   used and which settings from this script it reads.
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl, ggrepel)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/derive_millennial_parms.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")
source("R/scenario_uncertainty.R"); source("R/morris_sensitivity.R")
source("R/sensitivity_helpers.R")     # parameter sets, ranges, model evaluation
source("R/sensitivity_morris_run.R")  # run_morris(), outputs, figures
source("R/sensitivity_uncertainty.R") # Latin-hypercube uncertainty propagation
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

animal_pools <- c("Earthworm", "Detritivore", "RootHerb")
derive_fn    <- match.fun(model_table[[model]]$derive)

#------------------------------------------------------------------------#
# Morris + range settings ----
#------------------------------------------------------------------------#
morris_r      <- 200     # trajectories per scenario (raise for stable mu*/sigma)
morris_levels <- 4L     # grid levels p in the Morris design
buffer        <- 0.5    # MAIN no-uncertainty range = [value*buffer, value/buffer]
                        #                            = [0.5*value, 2*value] (half..double)
supp_frac     <- c(0.5, 2)     # SUPPLEMENTAL standardized range = 50-200% of value (half..double)
morris_seed   <- 27082026      # base seed; each scenario uses morris_seed + its index

run_main_morris <- TRUE   # main-text analysis (existing uncertainty)
run_supp_morris <- TRUE   # supplemental analysis (standardized half-double, Min/Max-bound)

#------------------------------------------------------------------------#
# Morris bootstrap (convergence of the ranking) ----
#------------------------------------------------------------------------#
boot_B     <- 1000        # bootstrap replicates over trajectories
boot_top_k <- 5           # report P(parameter is in the top-k by mu*)

#------------------------------------------------------------------------#
# Uncertainty propagation to the headline animal effect ----
#------------------------------------------------------------------------#
run_uncertainty_propagation <- TRUE
up_n        <- 2000       # Latin-hypercube samples per scenario (3 equilibria each);
                          # ~10-30% pass the realism filter, so keep this large
up_seed     <- 20260827
up_gsd          <- 1.5           # geometric SD of the log-normal used for parameters
                                 # with no reported uncertainty (1.5 = ~68% within x/1.5)
up_minmax_dist  <- "triangular"  # Min/Max parameters: "triangular" (mode = model value) or "uniform"

# Baseline realism filter (see the section below). Defaults: each scenario's
# all-default baseline equilibrium, accepted range = reference / factor .. x factor.
up_filter_pools  <- c("TotalC", "P", "M", "MIC", "B")   # TotalC = all non-animal pools
up_filter_factor <- 2                 # default accepted range: reference / 2 .. x 2
up_filter_factor_by_pool <- c(MIC = 10) # per-pool overrides: organic-horizon microbial
                                        # biomass is tiny and volatile, so allow / 10 .. x 10
up_targets_file  <- "Data/baseline_targets.csv"   # optional; overrides the defaults

paired_n_top <- 8          # parameters per scenario in the paired Morris figure

#------------------------------------------------------------------------#
# Parallel execution ----
#------------------------------------------------------------------------#
# Model evaluations (Morris trajectories, uncertainty samples) run on n_cores
# worker processes; 1 = serial. The bootstrap itself is vectorised and takes
# seconds, so it stays on the main process.
n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
if (isTRUE(getOption("sfw.debugonce", FALSE))) n_cores <- 1   # debugger cannot reach workers

#------------------------------------------------------------------------#
# Regime flags (kept in the results, reported with and without) ----
#------------------------------------------------------------------------#
#   microbe_collapse : MIC or B in the treatment below collapse_rel x baseline
#   kb_floor         : aggregate breakdown rate k_b*(1 + Earthworm*k_b_slope_pint)
#                      at or below the model's 0.001 floor (breakdown switched off)
collapse_rel <- 1e-6
kb_floor     <- 0.001

#------------------------------------------------------------------------#
# Equilibrium solver horizon (days) ----
#------------------------------------------------------------------------#
# Parameter sets near microbial collapse approach equilibrium very slowly
# (> 1e7 d). The solver stops as soon as the steady-state test is met, so a
# long horizon only costs time where it is needed.
eq_max_time <- 1e9

#------------------------------------------------------------------------#
# GENERAL model parameters to include ----
#------------------------------------------------------------------------#
# GENERAL model parameters to include (animal parameters are added per scenario
# from scenarios.xlsx below). Only those present in a scenario's parms are used.
sweep_params <- c("k_frag_litter", "k_frag_organic", "k_l_o", "k_l", "k_b", "k_pa",
                  "k_ma", "k_litterfall_ann", "k_litterfall_herb_ann", "k_MICd",
                  "k_bd", "root_to_organic", "a_root_herb", "k_exudate_tree",
                  "k_exudate_herb", "k_frag_CWD", "K_ol", "alpha_ol", "alpha_ob",
                  "K_ob", "p_c", "psi_matric", "lambda_mat", "k_a_min",
                  "alpha_pl", "alpha_lb", "K_pl", "K_lb", "p1", "p2", "K_ld",
                  "CUE_T", "p_a", "p_b", "k_mort_root_tree", "k_mort_root_herb")
# rho_p (particle density) is deliberately NOT varied: it is close to constant
# for mineral soils, and with porosity = 1 - BD / rho_p a half-to-double range
# produces impossible soils (porosity <= 0 or below the soil moisture).
# Porosity uncertainty enters through bulk density (BD) instead.

#------------------------------------------------------------------------#
# Pretty display names ----
#------------------------------------------------------------------------#
# Pretty display names (code -> label), general + animal parameters combined.
# pretty_param() falls back to the raw code for anything not listed.
param_labels <- c(
  k_frag_litter = "Litter fragmentation", k_frag_organic = "Organic fragmentation",
  k_frag_CWD = "CWD fragmentation", k_l_o = "DOM -> POC turnover",
  k_l = "LMWC leaching", k_b = "Aggregate breakdown", k_pa = "POC -> aggregate",
  k_ma = "MAOC -> aggregate", k_MICd = "Microbial turnover (organic)",
  k_bd = "Microbial turnover (mineral)", k_a_min = "Min. aeration factor",
  k_exudate_tree = "Tree root exudation", k_exudate_herb = "Herb root exudation",
  pct_claysilt = "Clay + silt (%)", root_to_organic = "Root death -> Organic",
  a_root_herb = "Herb root allocation", MAT = "Mean annual temperature",
  MAtheta = "Mean annual soil moisture", NPP_herb = "Herbaceous NPP", NPP_tree = "Tree NPP",
  BD = "Bulk density", rho_p = "Particle density", pH = "Soil pH",
  psi_matric = "Matric potential", lambda_mat = "Matric sensitivity",
  CUE_T = "CUE temperature slope", K_ol = "Organic->DOM half-saturation",
  K_ob = "DOM->microbe half-saturation", K_pl = "POC->LMWC half-saturation",
  K_lb = "LMWC->microbe half-saturation", K_ld = "LMWC->MAOC max sorption",
  alpha_ol = "Organic->DOM pre-exponential", alpha_ob = "DOM->microbe pre-exponential",
  alpha_pl = "POC->LMWC pre-exponential", alpha_lb = "LMWC->microbe pre-exponential",
  p1 = "pH sorption coef. 1", p2 = "pH sorption coef. 2",
  p_a = "Aggregate->POC partition", p_b = "Necromass->MAOC partition",
  p_c = "Clay+silt protection coef.", k_litterfall_ann = "Tree annual litterfall prop.",
  k_litterfall_herb_ann = "Herb. annual litterfall prop.",
  k_mort_root_tree = "Tree root mortality", k_mort_root_herb = "Herb. root mortality",
  T_amp = "Temperature amplitude", theta_amp = "Moisture amplitude",
  fCLAY = "Clay fraction", LigFrac = "Lignin fraction",
  # animal parameters
  c_earthworm_om = "Earthworm organic matter consumption",
  k_b_slope_pint = "Earthworm effect on aggregate stability (proportion)",
  c_earthworm_soil = "Earthworm soil consumption", p_earthworm = "Earthworm production efficiency",
  d_earthworm = "Earthworm mortality", prop_feaces_earthworm_LMWC = "Earthworm feces to LMWC fraction",
  a_earthworm_soil = "Earthworm soil assimilation", c_earthworm_litter = "Earthworm litter consumption",
  a_earthworm = "Earthworm litter assimilation efficiency",
  d_detritivores = "Detritivore mortality",
  slope_pint_det_k_frag_organic = "Detritivore effect of organic fragmentation (proportion)",
  c_detritivores = "Detritivore consumption", p_detritivores = "Detritivore production efficiency",
  a_detritivores = "Detritivore assimilation efficiency",
  slope_pint_det_k_frag_litter = "Detritivore litter fragmentation response",
  a_rootherb = "Root herbivore assimilation efficiency",
  p_rootherb = "Root herbivore production efficiency", c_rootherb = "Root herbivore consumption",
  k_exudate_intercept = "Root exudation intercept",
  k_exudate_slope = "Root herbivore effect on exudation", d_rootherb = "Root herbivore mortality")


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
main_cols <- c(SD = "#1b7837", MinMax = "#2166ac", buffer = "#b2182b")
main_labs <- c(SD = "+/- 2 standard deviations", MinMax = "Min/Max",
               buffer = "50% to 200% (half to double)")
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
                   "R/sensitivity_helpers.R", "R/sensitivity_uncertainty.R"),
  export = c("derive_fn", "animal_pools", "collapse_rel", "kb_floor", "eq_max_time"))
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


#------------------------------------------------------------------------#
# UNCERTAINTY PROPAGATION TO THE HEADLINE ANIMAL EFFECT ----
#------------------------------------------------------------------------#
# Morris ranks parameters; it does not say how uncertain the animal effect
# itself is. Here ALL screened parameters are sampled jointly (Latin hypercube,
# up_n samples per scenario) and for every sample the equilibrium TOTAL and
# DIRECT animal effects on total soil + root C are recomputed.
#
# Sampling distributions (concentrated on plausible values, not uniform):
#   SD reported      normal(sheet value, SD), truncated at +/- 2 SD
#   Min/Max reported triangular(Min, mode = model value, Max)   [up_minmax_dist]
#   nothing reported log-normal centred on the model value with geometric SD
#                    up_gsd, truncated to half..double
#   All truncated to the physical guard rails. Parameters are independent.
#
# BASELINE REALISM FILTER: a sample is ACCEPTED only if its no-animal
#   equilibrium stays within the target range for every listed pool. Default
#   targets = the all-default baseline equilibrium of each scenario, with
#   lower/upper = reference / factor and reference * factor, where factor is
#   up_filter_factor unless the pool is listed in up_filter_factor_by_pool.
#   To use your own targets, copy Results/baseline_targets_used.csv to
#   Data/baseline_targets.csv and edit lower/upper (or add rows; pool = any
#   state variable or TotalC). Rows in that file override the defaults.
#
# Outputs: animal_effect_uncertainty_samples.csv      every sample
#          animal_effect_uncertainty_summary.csv      per scenario x subset
#                                                     (all / unflagged / accepted)
#          animal_effect_uncertainty_diagnostics.csv  counts, acceptance, rejections
#          animal_effect_uncertainty_inputs.csv       distribution of every parameter
#          baseline_targets_used.csv                  realism targets actually used
#          figures/animal_effect_uncertainty.png


if (run_uncertainty_propagation) {
  up <- run_uncertainty_propagation_fn()
  write_csv(up$samples,     file.path(res_dir, "animal_effect_uncertainty_samples.csv"))
  write_csv(up$summary,     file.path(res_dir, "animal_effect_uncertainty_summary.csv"))
  write_csv(up$diagnostics, file.path(res_dir, "animal_effect_uncertainty_diagnostics.csv"))
  write_csv(up$inputs,      file.path(res_dir, "animal_effect_uncertainty_inputs.csv"))
  write_csv(up$targets,     file.path(res_dir, "baseline_targets_used.csv"))

  cat("\nHeadline animal effect with propagated parameter uncertainty (g C m-2),",
      "accepted samples:\n")
  print(up$summary %>% filter(subset == "accepted") %>%
          transmute(scenario, n,
                    total  = sprintf("%.1f (%.1f to %.1f)", total_median, total_q025, total_q975),
                    direct = sprintf("%.1f (%.1f to %.1f)", direct_median, direct_q025, direct_q975),
                    p_sign_as_default = round(p_sign_as_default, 3)),
        row.names = FALSE)

  # FIGURE: accepted samples (filled) vs all converged samples (dashed outline);
  # solid line = all-default effect; shaded band = central 95% of accepted.
  long <- up$samples %>% filter(converged == 1) %>%
    select(scenario, accepted, total, direct) %>%
    pivot_longer(c(total, direct), names_to = "effect", values_to = "value")
  long_acc <- long %>% filter(accepted)
  bands <- up$summary %>% filter(subset == "accepted") %>%
    transmute(scenario,
              total_lo = total_q025, total_hi = total_q975, total_def = total_default,
              direct_lo = direct_q025, direct_hi = direct_q975, direct_def = direct_default) %>%
    pivot_longer(-scenario, names_to = c("effect", ".value"), names_sep = "_")
  acc_txt <- up$diagnostics %>%
    transmute(scenario, lab = sprintf("accepted %d of %d", n_accepted, n_converged)) %>%
    tidyr::crossing(effect = c("total", "direct"))
  # total effect first, direct effect second
  eff_lv   <- c("total", "direct")
  long     <- mutate(long,     effect = factor(effect, levels = eff_lv))
  long_acc <- mutate(long_acc, effect = factor(effect, levels = eff_lv))
  bands    <- mutate(bands,    effect = factor(effect, levels = eff_lv))
  acc_txt  <- mutate(acc_txt,  effect = factor(effect, levels = eff_lv))
  p_up <- ggplot(long_acc, aes(value, fill = effect, colour = effect)) +
    geom_rect(data = bands, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf, fill = effect),
              inherit.aes = FALSE, alpha = 0.12) +
    geom_density(alpha = 0.35) +
    geom_density(data = long, fill = NA, linetype = "dashed", linewidth = 0.4) +
    geom_vline(data = bands, aes(xintercept = def, colour = effect), linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
    geom_text(data = acc_txt, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE,
              hjust = 1.05, vjust = 1.5, size = 2.8, colour = "grey30") +
    facet_wrap(~ scenario + effect, scales = "free", ncol = 2, labeller = scenario_labeller) +
    scale_fill_manual(values = c(total = "#1b7837", direct = "#762a83"), name = "Effect") +
    scale_colour_manual(values = c(total = "#1b7837", direct = "#762a83"), name = "Effect") +
    labs(x = expression("Animal effect on total soil + root C (g C m"^-2*")"),
         y = "Density",
         caption = paste("Filled = samples passing the baseline realism filter; dashed = all converged samples.",
                         "Solid line = all-default effect; shaded band = central 95% of accepted; dotted = no effect.")) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "animal_effect_uncertainty.png"), p_up,
         width = 10, height = 2.6 * length(unique(long$scenario)) + 1.2, dpi = 150)
  cat("Wrote uncertainty outputs to", res_dir, "and", file.path(fig_dir, "animal_effect_uncertainty.png"), "\n")
}

# shut down the worker pool
if (!is.null(cl)) parallel::stopCluster(cl)
