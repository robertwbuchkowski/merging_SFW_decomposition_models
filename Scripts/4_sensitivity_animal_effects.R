# ============================================================
# SENSITIVITY OF THE ANIMAL EFFECT TO MODEL + ANIMAL PARAMETERS (equilibrium)
# ------------------------------------------------------------
# Merged sensitivity + parameter-uncertainty analysis. For each scenario and
# each parameter it re-computes the equilibrium animal effect on total soil +
# root C (treatment - baseline) at chosen parameter values, with no dynamic
# spin-up.
#
# Two ranges per scenario x parameter:
#   BASELINE RANGE (always): effect at 50% and 150% of the default value (i.e.
#       value x 0.5 and x 1.5), clamped to physical guard rails. If a guard rail
#       cuts the point off (e.g. 1.5 x root_to_organic = 1.05 > 1), the range is
#       reported up to the nearest allowed value (1.0). This is the ranking
#       metric: parameters are ranked per scenario by |effect(150%)-effect(50%)|.
#   REPORTED RANGE (when available): effect at the parameter's SD or Min/Max
#       bounds from scenarios.xlsx, drawn ON TOP of the baseline range.
#
# Parameters: the soil/model list (sweep_params, all scenarios) PLUS the animal
# parameters from scenarios.xlsx (Category Animal/Effect), which are added only
# to the scenario they belong to (an animal parameter is irrelevant to a
# scenario without that animal).
#
# Output: Results/animal_effect_sensitivity.csv (long, all points),
#         Results/animal_effect_sensitivity_summary.csv (ranked ranges),
#         Results/figures/animal_effect_sensitivity_range.png
# Run from the project root.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/derive_millennial_parms.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")
source("R/scenario_uncertainty.R")

model   <- "millennial"
scen    <- read_scenarios("Data/scenarios.xlsx")
scen$MitePredator <- NULL
scenarios <- names(scen)

use_fitted_params <- TRUE
fitted_params <- if (use_fitted_params && file.exists("Results/fitted_animal_params.csv"))
  load_fitted_params("Results/fitted_animal_params.csv") else NULL

fig_dir <- "Results/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
res_dir <- "Results";         dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

animal_pools <- c("Earthworm", "Detritivore", "DetPredator", "RootHerb")
derive_fn    <- match.fun(model_table[[model]]$derive)
`%||%` <- function(a, b) if (is.null(a)) b else a

baseline_frac <- c(0.9, 1.1)     # always-evaluated baseline range (90% / 110%)

# ------------------------------------------------------------
# PARAMETERS. Soil/model parameters (all scenarios) + a pretty-name lookup.
# ------------------------------------------------------------
sweep_params <- c("k_frag_litter", "k_frag_organic", "k_l_o", "k_l",
                  "k_b", "k_pa", "k_ma", "pct_claysilt",
                  "k_MICd", "k_bd", "root_to_organic", "a_root_herb",
                  "MAT", "MAtheta", "NPP_herb", "NPP_tree", "k_exudate_tree",
                  "k_exudate_herb", "k_frag_CWD", "K_ol", "alpha_ol", "alpha_ob",
                  "K_ob", "BD", "pH", "p_c", "rho_p", "psi_matric", "lambda_mat",
                  "k_a_min", "alpha_pl", "alpha_lb", "K_pl", "K_lb", "p1", "p2",
                  "K_ld", "CUE_T", "p_a", "p_b")

param_labels <- c(
  k_frag_litter = "Litter fragmentation", k_frag_organic = "Organic fragmentation",
  k_frag_CWD = "CWD fragmentation", k_l_o = "DOM -> POC turnover",
  k_l = "LMWC leaching", k_b = "Aggregate breakdown", k_pa = "POC -> aggregate",
  k_ma = "MAOC -> aggregate", k_MICd = "Microbial turnover (organic)",
  k_bd = "Microbial turnover (mineral)", k_a_min = "Min. aeration factor",
  k_exudate_tree = "Tree root exudation", k_exudate_herb = "Herb root exudation",
  pct_claysilt = "Clay + silt (%)", root_to_organic = "Root death -> Organic",
  a_root_herb = "Herb root allocation", MAT = "Mean annual temperature",
  MAtheta = "Mean annual moisture", NPP_herb = "Herbaceous NPP", NPP_tree = "Tree NPP",
  BD = "Bulk density", rho_p = "Particle density", pH = "Soil pH",
  psi_matric = "Matric potential", lambda_mat = "Matric sensitivity",
  CUE_T = "CUE temperature slope", K_ol = "Organic->DOM half-saturation",
  K_ob = "DOM->microbe half-saturation", K_pl = "POC->LMWC half-saturation",
  K_lb = "LMWC->microbe half-saturation", K_ld = "LMWC->MAOC max sorption",
  alpha_ol = "Organic->DOM pre-exponential", alpha_ob = "DOM->microbe pre-exponential",
  alpha_pl = "POC->LMWC pre-exponential", alpha_lb = "LMWC->microbe pre-exponential",
  p1 = "pH sorption coef. 1", p2 = "pH sorption coef. 2",
  p_a = "Aggregate->POC partition", p_b = "Necromass->MAOC partition",
  p_c = "Clay+silt protection coef.",
  # animal parameters (labels for common ones; unlisted fall back to the code)
  a_earthworm = "Earthworm assimilation eff.", p_earthworm = "Earthworm production eff.",
  d_earthworm = "Earthworm death rate", E_earthworm = "Earthworm activation energy",
  c_earthworm_litter = "Earthworm litter feeding", c_earthworm_soil = "Earthworm soil feeding",
  c_earthworm_om = "Earthworm OM feeding", prop_feaces_earthworm_LMWC = "Earthworm faeces -> LMWC",
  a_detritivores = "Detritivore assimilation eff.", p_detritivores = "Detritivore production eff.",
  d_detritivores = "Detritivore death rate", E_detritivores = "Detritivore activation energy",
  c_detritivores = "Detritivore feeding rate",
  a_rootherb = "Root herbivore assimilation eff.", p_rootherb = "Root herbivore production eff.",
  d_rootherb = "Root herbivore death rate", E_root_herb = "Root herbivore activation energy",
  c_rootherb = "Root herbivore feeding rate",
  slope_pint_det_k_frag_litter = "Detritivore -> litter frag.",
  slope_pint_det_k_frag_organic = "Detritivore -> organic frag.",
  k_b_slope_pint = "Earthworm -> aggregate breakdown",
  k_exudate_slope = "Root herbivore -> exudation",
  k_exudate_intercept = "Root exudation intercept")

pretty_param <- function(x) ifelse(x %in% names(param_labels), param_labels[x], x)

# ------------------------------------------------------------
# PHYSICAL GUARD RAILS. a_*/p_* efficiencies in [0,1] (except the soil partition
# coefficients p_a, p_b [0,1] and p_c unbounded); named fractions in [0,1];
# pct_claysilt in [0,100]. clamp_to_bounds() pulls a value onto the nearest
# allowed edge; in_bounds() tells whether it was inside.
# ------------------------------------------------------------
param_bounds <- function(param) {
  if (grepl("^(a_|p_)", param) && !param %in% c("p_a", "p_b", "p_c"))
    return(c(0, 1))
  switch(param,
         pct_claysilt = c(0, 100),
         p_a = , p_b = c(0, 1),
         LigFrac = , a_root_herb = , root_to_organic = ,
         prop_feaces_earthworm_LMWC = c(0, 1),
         c(-Inf, Inf))
}
clamp_to_bounds <- function(v, param) { b <- param_bounds(param); pmin(pmax(v, b[1]), b[2]) }

# ------------------------------------------------------------
# set_param(): change ONE parameter and rebuild derived params + forcing.
# ------------------------------------------------------------
set_param <- function(obj, param, value) {
  obj$parms[[param]]        <- value
  obj$parms                 <- derive_fn(obj$parms)
  obj$parms$climate_forcing <- make_climate_forcing(obj$parms)
  obj
}

# equilibrium TOTAL animal effect on total soil + root C (scalar, g C m-2)
animal_effect_totalC <- function(pair) {
  pair$baseline  <- spinup_equilibrium(pair$baseline, verbose = FALSE)
  pair$treatment <- spinup_equilibrium(pair$treatment,
                                       warm_start = pair$baseline$init_state_spin,
                                       verbose = FALSE)
  eq_b <- pair$baseline$init_state_spin
  eq_t <- pair$treatment$init_state_spin
  soil <- setdiff(intersect(names(eq_b), names(eq_t)), animal_pools)
  sum(eq_t[soil]) - sum(eq_b[soil])
}

# build the scenario pair once per scenario (with fitted params applied)
make_pair <- function(scenario) {
  p <- setup_scenario_pair(model, scen, scenario)
  if (!is.null(fitted_params))
    p$treatment <- apply_fitted_params(p$treatment, fitted_params, model, scenario, verbose = FALSE)
  p
}

# evaluate the effect for a prebuilt pair with ONE parameter set to `value`
effect_at <- function(base_pair, scenario, param = NULL, value = NULL) {
  pair <- base_pair
  if (!is.null(param)) {
    pair$treatment <- set_param(pair$treatment, param, value)
    pair$baseline  <- set_param(pair$baseline,  param, value)
  }
  tryCatch(animal_effect_totalC(pair), error = function(e) {
    message("  ", scenario, " ", param %||% "default", "=", signif(value, 4),
            " failed: ", conditionMessage(e)); NA_real_ })
}

# ------------------------------------------------------------
# Assemble the parameter set per scenario: soil params (all scenarios) + the
# animal parameters listed for THAT scenario in scenarios.xlsx. Attach the
# reported uncertainty (SD / Min-Max / CV2) so it can be drawn on top.
# ------------------------------------------------------------
unc_all <- tryCatch(read_param_uncertainty("Data/scenarios.xlsx"),
                    error = function(e) { message("uncertainty read failed: ",
                                                  conditionMessage(e)); NULL })
animal_unc <- if (!is.null(unc_all))
  unc_all[unc_all$category %in% c("Animal", "Effect"), , drop = FALSE] else
  unc_all[0, , drop = FALSE]

# ------------------------------------------------------------
# MAIN LOOP: for each scenario x parameter, evaluate the baseline (10/110%)
# range and, when present, the reported SD/Min-Max range.
# ------------------------------------------------------------
records <- list()
for (scenario in scenarios) {
  base_pair  <- make_pair(scenario)
  parms_here <- base_pair$treatment$parms

  # which animal params belong to this scenario (present in its parms)
  a_here <- if (nrow(animal_unc)) {
    a <- animal_unc[animal_unc$scenario == scenario, , drop = FALSE]
    a[a$parameter %in% names(parms_here), , drop = FALSE]
  } else animal_unc[0, , drop = FALSE]

  params_here <- unique(c(intersect(sweep_params, names(parms_here)), a_here$parameter))
  eff_default <- effect_at(base_pair, scenario)          # all params at default

  for (p in params_here) {
    d0 <- parms_here[[p]]
    if (!is.finite(d0)) next
    is_animal <- p %in% a_here$parameter

    # --- BASELINE range: 50% / 150%, clamped to guard rails ---
    raw   <- d0 * baseline_frac
    lo_b  <- clamp_to_bounds(min(raw), p)
    hi_b  <- clamp_to_bounds(max(raw), p)
    cut   <- any(abs(c(lo_b, hi_b) - range(raw)) > 1e-12)   # guard rail truncated a point
    e_lo_b <- effect_at(base_pair, scenario, p, lo_b)
    e_hi_b <- effect_at(base_pair, scenario, p, hi_b)

    # --- REPORTED range (SD / Min-Max / CV2), if this param has one here ---
    src <- NA_character_; e_lo_r <- NA_real_; e_hi_r <- NA_real_; lo_r <- NA_real_; hi_r <- NA_real_
    ur <- if (!is.null(unc_all))
      unc_all[unc_all$scenario == scenario & unc_all$parameter == p, , drop = FALSE] else NULL
    if (!is.null(ur) && nrow(ur)) {
      src  <- ur$unc_source[1]
      lo_r <- clamp_to_bounds(ur$lo[1], p); hi_r <- clamp_to_bounds(ur$hi[1], p)
      e_lo_r <- effect_at(base_pair, scenario, p, lo_r)
      e_hi_r <- effect_at(base_pair, scenario, p, hi_r)
    }

    records[[length(records) + 1]] <- tibble(
      scenario = scenario, parameter = p, is_animal = is_animal, default = d0,
      base_lo_val = lo_b, base_hi_val = hi_b, base_clamped = cut,
      effect_default = eff_default, effect_base_lo = e_lo_b, effect_base_hi = e_hi_b,
      rep_source = src, rep_lo_val = lo_r, rep_hi_val = hi_r,
      effect_rep_lo = e_lo_r, effect_rep_hi = e_hi_r)
  }
}
res <- bind_rows(records) %>%
  mutate(base_min = pmin(effect_base_lo, effect_base_hi),
         base_max = pmax(effect_base_lo, effect_base_hi),
         base_span = base_max - base_min,
         rep_min  = pmin(effect_rep_lo, effect_rep_hi),
         rep_max  = pmax(effect_rep_lo, effect_rep_hi))
write_csv(res, file.path(res_dir, "animal_effect_sensitivity.csv"))

# ------------------------------------------------------------
# RANKING: within each scenario, rank parameters by the BASELINE (10/110%) span.
# A parameter "does not vary by scenario" if it is NOT in scenarios.xlsx (its
# value is the invariant model default) -> flagged with a leading dot, as in the
# old heat map.
# ------------------------------------------------------------
in_sheet <- if (!is.null(unc_all)) unique(unc_all$parameter) else character(0)
res <- res %>%
  mutate(varies   = parameter %in% in_sheet,
         lab_core = pretty_param(parameter),
         lab      = ifelse(varies, lab_core, paste0("\u00b7 ", lab_core)))

summary_tbl <- res %>%
  group_by(scenario) %>% arrange(scenario, base_span) %>%
  mutate(rank = row_number()) %>% ungroup() %>%
  select(scenario, parameter, lab, is_animal, varies, rep_source,
         default, base_lo_val, base_hi_val, base_clamped,
         effect_default, base_min, base_max, base_span,
         rep_lo_val, rep_hi_val, rep_min, rep_max, rank)
write_csv(summary_tbl, file.path(res_dir, "animal_effect_sensitivity_summary.csv"))

# ------------------------------------------------------------
# FIGURE: one facet per scenario. Grey baseline (10/110%) range bar per
# parameter (ranked, most sensitive at top), with the reported SD/Min-Max range
# drawn ON TOP in colour when available. Dotted line = all-default effect.
# A leading dot on a parameter label = not in scenarios.xlsx (does not vary by
# scenario). An open tick marks a baseline point truncated by a guard rail.
# ------------------------------------------------------------
plot_df <- res %>%
  group_by(scenario) %>%
  mutate(lab = fct_reorder(lab, base_span)) %>%
  ungroup()

rep_df <- plot_df %>% filter(!is.na(rep_source), is.finite(rep_min), is.finite(rep_max), rep_source != "CV2")
clamp_df <- plot_df %>% filter(base_clamped)

p_range <- ggplot(plot_df, aes(y = lab)) +
  geom_vline(aes(xintercept = effect_default), linetype = "dotted", colour = "grey50") +
  # baseline 10/110% range (grey, always)
  geom_segment(aes(x = base_min, xend = base_max, yend = lab),
               colour = "grey70", linewidth = 2.4, lineend = "round") +
  # reported SD / Min-Max range on top (coloured)
  geom_segment(data = rep_df,
               aes(x = rep_min, xend = rep_max, yend = lab, colour = rep_source),
               linewidth = 1.0, lineend = "round") +
  geom_point(data = rep_df, aes(x = rep_min, colour = rep_source), size = 1.4) +
  geom_point(data = rep_df, aes(x = rep_max, colour = rep_source), size = 1.4) +
  # mark guard-rail-truncated baseline points
  geom_point(data = clamp_df, aes(x = base_min), shape = 1, size = 2, colour = "grey30") +
  geom_point(data = clamp_df, aes(x = base_max), shape = 1, size = 2, colour = "grey30") +
  facet_wrap(~scenario, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c(SD = "#1b7837", MinMax = "#2166ac", CV2 = "#b2182b"),
                      name = "Reported uncertainty", na.translate = FALSE) +
  labs(
    # title = "Sensitivity of the animal effect to model + animal parameters",
    #    subtitle = paste0("Grey = effect range across 10-110% of default (ranking metric); ",
    #                      "colour = reported SD / Min-Max.\n",
    #                      "Leading dot = parameter not in scenarios.xlsx (does not vary by scenario); ",
    #                      "open circle = 10/110% point hit a guard rail."),
       x = expression("Animal effect on total soil + root C (g C m"^-2*")"), y = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "animal_effect_sensitivity_range.png"),
       p_range, width = 13, height = 10, dpi = 600)
print(p_range)

cat("\nWrote:\n  ", file.path(res_dir, "animal_effect_sensitivity.csv"),
    "\n  ", file.path(res_dir, "animal_effect_sensitivity_summary.csv"),
    "\n  ", file.path(fig_dir, "animal_effect_sensitivity_range.png"), "\n")
