# ============================================================
# RANGE OF ANIMAL EFFECTS ACROSS ANIMAL-PARAMETER UNCERTAINTY (equilibrium)
# ------------------------------------------------------------
# For each scenario, take the ANIMAL parameters listed in scenarios.xlsx
# (Category "Animal" or "Effect") and, one at a time, re-evaluate the
# equilibrium animal effect at the parameter's reported low and high bound. The
# range of the resulting effect shows how uncertain each animal parameter makes
# the predicted animal effect.
#
# Uncertainty bounds come from R/scenario_uncertainty.R with the precedence:
#   SD (value +/- SD)  ->  else Min/Max  ->  else CV = 2 (value +/- 2|value|)
# CV = 2 rows are placeholders and raise a warning.
#
# Output: Results/animal_param_uncertainty_effect.csv and
#         Results/figures/animal_param_uncertainty_range.png
# Run from the project root (or via Scripts/0_run_all.R).
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")
source("R/scenario_uncertainty.R")
source("R/derive_millennial_parms.R")

model <- "millennial"
scen  <- read_scenarios("Data/scenarios.xlsx")

fitted_params <- if (file.exists("Results/fitted_animal_params.csv"))
  load_fitted_params("Results/fitted_animal_params.csv") else NULL

fig_dir <- "Results/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
res_dir <- "Results";         dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

animal_pools <- c("Earthworm", "Detritivore", "DetPredator", "RootHerb")
derive_fn    <- match.fun(model_table[[model]]$derive)
`%||%` <- function(a, b) if (is.null(a)) b else a

set_param <- function(obj, param, value) {
  obj$parms[[param]]        <- value
  obj$parms                 <- derive_fn(obj$parms)
  obj$parms$climate_forcing <- make_climate_forcing(obj$parms)
  obj
}

# equilibrium TOTAL animal effect on total soil+root C (a single scalar summary)
animal_effect_totalC <- function(pair) {
  pair$baseline  <- spinup_equilibrium(pair$baseline, verbose = FALSE)
  pair$treatment <- spinup_equilibrium(pair$treatment,
                                       warm_start = pair$baseline$init_state_spin,
                                       verbose = FALSE)
  eq_b <- pair$baseline$init_state_spin
  eq_t <- pair$treatment$init_state_spin
  soil <- setdiff(intersect(names(eq_b), names(eq_t)), animal_pools)
  sum(eq_t[soil]) - sum(eq_b[soil])          # g C m-2 (treatment - baseline)
}

# evaluate the effect for one scenario with ONE animal parameter set to `value`
effect_at <- function(scenario, param = NULL, value = NULL) {
  pair <- setup_scenario_pair(model, scen, scenario)
  if (!is.null(fitted_params))
    pair$treatment <- apply_fitted_params(pair$treatment, fitted_params, model, scenario,
                                          verbose = FALSE)
  if (!is.null(param)) {
    pair$treatment <- set_param(pair$treatment, param, value)
    pair$baseline  <- set_param(pair$baseline,  param, value)   # keep arms consistent
  }
  tryCatch(animal_effect_totalC(pair), error = function(e) {
    message("  ", scenario, " ", param %||% "default", "=", signif(value, 4),
            " failed: ", conditionMessage(e)); NA_real_ })
}

# ------------------------------------------------------------
# For each scenario, evaluate the effect at each animal parameter's lo / default
# / hi bound. Only parameters that actually exist in the model's parms are used.
# ------------------------------------------------------------
unc <- read_param_uncertainty("Data/scenarios.xlsx")

rows <- list()
for (scenario in unique(unc$scenario)) {
  if (!scenario %in% names(scen)) next
  base_parms <- setup_scenario_pair(model, scen, scenario)$treatment$parms
  eff_default <- effect_at(scenario)                       # all params at default
  sub <- unc[unc$scenario == scenario, , drop = FALSE]
  for (i in seq_len(nrow(sub))) {
    p <- sub$parameter[i]
    if (!p %in% names(base_parms)) next                    # not a model parameter
    e_lo <- effect_at(scenario, p, sub$lo[i])
    e_hi <- effect_at(scenario, p, sub$hi[i])
    rows[[length(rows) + 1]] <- data.frame(
      scenario = scenario, parameter = p, unc_source = sub$unc_source[i],
      value = sub$value[i], lo = sub$lo[i], hi = sub$hi[i],
      effect_default = eff_default, effect_lo = e_lo, effect_hi = e_hi,
      effect_min = min(e_lo, e_hi, na.rm = TRUE),
      effect_max = max(e_lo, e_hi, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }
}
res <- bind_rows(rows) %>%
  mutate(effect_span = effect_max - effect_min) %>%
  arrange(scenario, desc(effect_span))
write_csv(res, file.path(res_dir, "animal_param_uncertainty_effect.csv"))


# Renaming the parameter names:
param_labels <- c(
  c_earthworm_om             = "Earthworm organic matter consumption",
  k_b_slope_pint             = "Earthworm effect on aggregate stability (proportion)",
  MAT                        = "Mean annual temperature",
  c_earthworm_soil           = "Earthworm soil consumption",
  p_earthworm                = "Earthworm production efficiency",
  d_earthworm                = "Earthworm mortality",
  prop_feaces_earthworm_LMWC = "Earthworm feces to LMWC fraction",
  a_earthworm_soil           = "Earthworm soil assimilation",
  c_earthworm_litter         = "Earthworm litter consumption",
  BD                         = "Bulk density",
  a_earthworm                = "Earthworm litter assimilation efficiency",
  NPP_tree                   = "Tree NPP",
  NPP_herb                   = "Herbaceous NPP",
  pct_claysilt               = "Clay + silt (%)",
  pH                         = "Soil pH",
  T_amp                      = "Temperature amplitude",
  MAtheta                    = "Mean annual soil moisture",
  theta_amp                  = "Moisture amplitude",
  fCLAY                      = "Clay fraction",
  LigFrac                    = "Lignin fraction",
  d_detritivores             = "Detritivore mortality",
  slope_pint_det_k_frag_organic = "Detritivore effect of organic fragmentation (proportion)",
  c_detritivores             = "Detritivore consumption",
  p_detritivores             = "Detritivore production efficiency",
  a_detritivores             = "Detritivore assimilation efficiency",
  slope_pint_det_k_frag_litter = "Detritivore litter fragmentation response",
  a_rootherb                 = "Root herbivore assimilation efficiency",
  p_rootherb                 = "Root herbivore production efficiency",
  c_rootherb                 = "Root herbivore consumption",
  k_exudate_intercept        = "Root exudation intercept",
  k_exudate_slope            = "Root herbivore effect on exudation",
  d_rootherb                 = "Root herbivore mortality"
)

pretty_param <- function(x) {
  ifelse(x %in% names(param_labels), param_labels[x], x)
}

# ------------------------------------------------------------
# FIGURE: for each scenario, a horizontal range (lo-hi effect) per animal
# parameter, ordered by span, with the all-default effect marked. Line-type /
# colour flags where the bound came from (SD vs Min-Max vs CV2 placeholder).
# ------------------------------------------------------------
plot_df <- res %>%
  group_by(scenario) %>%
  mutate(parameter = pretty_param(parameter)) %>%
  mutate(parameter = fct_reorder(parameter, effect_span)) %>%
  ungroup() %>%
  mutate(unc_source = ifelse(unc_source == "CV2", "CV = 2", unc_source))

p_range <- ggplot(plot_df, aes(y = parameter)) +
  geom_vline(aes(xintercept = effect_default), linetype = "dotted", colour = "grey50") +
  geom_segment(aes(x = effect_min, xend = effect_max,
                   yend = parameter, colour = unc_source), linewidth = 1.1) +
  geom_point(aes(x = effect_min, colour = unc_source), size = 1.6) +
  geom_point(aes(x = effect_max, colour = unc_source), size = 1.6) +
  facet_wrap(~scenario, scales = "free") +
  scale_colour_manual(values = c(SD = "#1b7837", MinMax = "#2166ac", `CV = 2` = "#b2182b"),
                      name = "Uncertainty source") +
  labs(
       x = expression("Animal effect on total soil + root C (g C m"^-2*")"),
       y = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "animal_param_uncertainty_range.png"),
       p_range, width = 12, height = 8, dpi = 150)
print(p_range)

cat("\nWrote:\n  ", file.path(res_dir, "animal_param_uncertainty_effect.csv"),
    "\n  ", file.path(fig_dir, "animal_param_uncertainty_range.png"), "\n")
