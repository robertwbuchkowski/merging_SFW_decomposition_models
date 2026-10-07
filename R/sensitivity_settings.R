#------------------------------------------------------------------------#
# Shared sensitivity / uncertainty settings ----
#------------------------------------------------------------------------#
# Settings shared by Scripts/4_sensitivity_animal_effects.R (Morris) and
# Scripts/5_uncertainty_animal_effect.R (Latin hypercube), so both analyses use
# the same parameters, the same standard-deviation rule and the same bounds.
# Edit here; both scripts source this file.

#------------------------------------------------------------------------#
# Standard-deviation rule ----
#------------------------------------------------------------------------#
# Every parameter gets a standard deviation (see param_sd_info() in
# R/sensitivity_helpers.R):
#   reported SD in scenarios.xlsx                     -> that SD      ("SD")
#   otherwise CV = cv_default (SD = cv_default*|value|)          ("CV-default")
#   ... reduced when a reported Min/Max implies less: SD = (Max-Min)/4 ("CV-MinMax")
# The analysed range is value +/- 2 SD, cut to the reported Min/Max and to the
# physical bounds below. Every cut is recorded (bound_hit) in the outputs.
cv_default <- 0.5   # CV = 0.5: the SD is half the model value

#------------------------------------------------------------------------#
# Physical bounds ----
#------------------------------------------------------------------------#
# Proportions are bounded to [0, 1]: all a_* and p_* efficiencies, plus the
# names below. The soil coefficients p_a, p_b, p_c are not efficiencies, so the
# prefix rule skips them; p_a and p_b are partition fractions (aggregate
# breakdown to POC vs MAOC; necromass to MAOC vs LMWC) and are bounded to
# [0, 1] in the table, while p_c (a Qmax scaling factor) is not bounded.
bound_exempt_prefix <- c("p_a", "p_b", "p_c")
param_bound_table <- list(
  p_a = c(0, 1), p_b = c(0, 1),
  root_to_organic = c(0, 1), a_root_herb = c(0, 1), LigFrac = c(0, 1), fCLAY = c(0, 1),
  prop_feaces_earthworm_LMWC = c(0, 1), winter_root_act_prop = c(0, 1),
  k_litterfall_ann = c(0, 1), k_litterfall_herb_ann = c(0, 1),
  MAtheta = c(0, 1),
  pct_claysilt = c(0, 100),
  pH = c(1, 14))
# All other parameters keep their sign: a positive value cannot go below
# positive_floor_frac x value (a rate of exactly 0 has no equilibrium), a
# negative value cannot exceed 0. Parameters listed in sign_free_params may
# change sign.
positive_floor_frac <- 0.01
sign_free_params    <- c("MAT")

#------------------------------------------------------------------------#
# GENERAL model parameters to include ----
#------------------------------------------------------------------------#
# GENERAL model parameters (scenario-specific and animal parameters are added
# per scenario from scenarios.xlsx). Only those present in a scenario's parms
# are used.
sweep_params <- c("k_frag_litter", "k_frag_organic", "k_l_o", "k_l", "k_b", "k_pa",
                  "k_ma", "k_litterfall_ann", "k_litterfall_herb_ann", "k_MICd",
                  "k_bd", "root_to_organic", "a_root_herb", "k_exudate_tree",
                  "k_exudate_herb", "k_frag_CWD", "K_ol", "alpha_ol", "alpha_ob",
                  "K_ob", "p_c", "psi_matric", "lambda_mat", "k_a_min",
                  "alpha_pl", "alpha_lb", "K_pl", "K_lb", "p1", "p2", "K_ld",
                  "CUE_T", "p_a", "p_b", "k_mort_root_tree", "k_mort_root_herb")
# rho_p (particle density) is deliberately NOT varied: it is close to constant
# for mineral soils, and with porosity = 1 - BD / rho_p a wide range produces
# impossible soils. Porosity uncertainty enters through bulk density (BD).

#------------------------------------------------------------------------#
# Pretty display names ----
#------------------------------------------------------------------------#
# Pretty display names (code -> label); pretty_param() falls back to the code.
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
# Model evaluation: pools, regime flags, solver horizon ----
#------------------------------------------------------------------------#
animal_pools <- c("Earthworm", "Detritivore", "RootHerb")
#   microbe_collapse : MIC or B in the treatment below collapse_rel x baseline
#   kb_floor         : aggregate breakdown k_b*(1 + Earthworm*k_b_slope_pint) at
#                      or below the model's 0.001 floor (breakdown switched off)
collapse_rel <- 1e-6
kb_floor     <- 0.001
# Parameter sets near microbial collapse approach equilibrium very slowly
# (> 1e7 d); the solver stops as soon as it is steady, so a long horizon only
# costs time where it is needed.
eq_max_time <- 1e9
