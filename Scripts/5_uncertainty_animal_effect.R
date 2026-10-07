#------------------------------------------------------------------------#
# UNCERTAINTY OF THE HEADLINE ANIMAL EFFECT (Latin hypercube) ----
#------------------------------------------------------------------------#
# Propagates parameter uncertainty to the headline result: the equilibrium
# TOTAL and DIRECT animal effect on total soil + root C, per scenario.
# All parameters screened by the Morris analysis (Scripts/4) are sampled
# JOINTLY by Latin hypercube; every sample is solved to equilibrium (baseline,
# treatment, treatment with indirect effects off).
#
# SAMPLING: each parameter is drawn from a normal distribution centred on its
#   model value, with SD = the SD reported in scenarios.xlsx, else a CV of
#   cv_default (0.5: SD = half the value), reduced when a reported Min/Max implies less spread
#   (labelled SD / CV-MinMax / CV-default). It is truncated at +/- 2 SD, at the
#   reported Min/Max and at the physical bounds (proportions in [0, 1]
#   including p_a and p_b, sign kept, pH, clay+silt). Parameters reaching a bound are
#   listed in animal_effect_uncertainty_inputs.csv (bound_hit). Same rule as
#   the Morris main analysis; shared settings are in R/sensitivity_settings.R.
#
# BASELINE REALISM FILTER: a sample is ACCEPTED only if its no-animal
#   equilibrium (a) keeps EVERY pool positive and non-zero (> up_pool_floor) and
#   (b) stays within target ranges for the listed pools (defaults: each
#   scenario's all-default baseline, reference / factor .. x factor; editable
#   in Data/baseline_targets.csv).
#
# Run from the project root (needs Results/fitted_animal_params.csv from
# Scripts/1). Helper functions: R/sensitivity_helpers.R, R/sensitivity_uncertainty.R.
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/derive_millennial_parms.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")
source("R/scenario_uncertainty.R"); source("R/morris_sensitivity.R")   # morris_cluster()
source("R/sensitivity_settings.R")    # shared: parameter list, SD rule, bounds, solver
source("R/sensitivity_helpers.R")     # parameter sets, SD rule, model evaluation
source("R/sensitivity_uncertainty.R") # sampling, headline effects, realism filter
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
# Sampling and filter settings ----
#------------------------------------------------------------------------#
up_n        <- 2000       # Latin-hypercube samples per scenario (3 equilibria each);
                          # ~10-30% pass the realism filter, so keep this large
up_seed     <- 20260827

# Baseline realism filter. Defaults: each scenario's all-default baseline
# equilibrium, accepted range = reference / factor .. x factor.
up_filter_pools  <- c("TotalC", "P", "M", "MIC", "B")   # TotalC = all non-animal pools
up_filter_factor <- 2                 # default accepted range: reference / 2 .. x 2
up_filter_factor_by_pool <- c(MIC = 10) # per-pool overrides: organic-horizon microbial
                                        # biomass is tiny and volatile, so allow / 10 .. x 10
up_targets_file  <- "Data/baseline_targets.csv"   # optional; overrides the defaults
up_pool_floor    <- 1e-9  # every baseline pool must also be above this (g C m-2):
                          # positive and non-zero (1e-9 = numerically zero)

#------------------------------------------------------------------------#
# Parallel execution ----
#------------------------------------------------------------------------#
n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
if (isTRUE(getOption("sfw.debugonce", FALSE))) n_cores <- 1   # debugger cannot reach workers

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
# Start the worker pool ----
#------------------------------------------------------------------------#
cl <- morris_cluster(
  n_cores,
  source_files = c("R/climate_forcing.R", "R/spinup.R", "R/derive_millennial_parms.R",
                   "R/millennial_model.R", "R/init_millennial_state.R",
                   "R/setup.R", "R/compare_functions.R", "R/fit_animals.R",
                   "R/sensitivity_settings.R", "R/sensitivity_helpers.R",
                   "R/sensitivity_uncertainty.R"),
  export = c("derive_fn"))
if (!is.null(cl)) cat("Running model evaluations on", length(cl), "worker processes.\n")

#------------------------------------------------------------------------#
# Run the uncertainty propagation ----
#------------------------------------------------------------------------#
up <- run_uncertainty_propagation_fn()
write_csv(up$samples,     file.path(res_dir, "animal_effect_uncertainty_samples.csv"))
write_csv(up$summary,     file.path(res_dir, "animal_effect_uncertainty_summary.csv"))
write_csv(up$diagnostics, file.path(res_dir, "animal_effect_uncertainty_diagnostics.csv"))
write_csv(up$inputs,      file.path(res_dir, "animal_effect_uncertainty_inputs.csv"))
write_csv(up$targets,     file.path(res_dir, "baseline_targets_used.csv"))

# parameters whose sampling range is cut by a bound (proportions, sign, Min/Max)
hits <- up$inputs %>% filter(nzchar(bound_hit))
if (nrow(hits)) {
  cat("\nParameters whose sampling range reaches a bound:\n")
  print(hits %>% transmute(scenario, parameter, sd_source, cv = round(cv, 3),
                           lo = signif(lo, 4), hi = signif(hi, 4), bound_hit),
        row.names = FALSE)
}

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

# shut down the worker pool
if (!is.null(cl)) parallel::stopCluster(cl)
