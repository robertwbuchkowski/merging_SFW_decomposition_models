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
  a_root_herb = c(0, 1),
  prop_feaces_earthworm_LMWC = c(0, 1), winter_root_act_prop = c(0, 1),
  MAtheta = c(0, 1), theta_amp_rel = c(0, 1),
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
sweep_params <- c("k_frag_litter", "k_frag_CWD", "k_l", "k_b", "k_pa", "k_ma",
                  "k_bd", "a_root_herb", "k_exudate_tree",
                  "p_c", "psi_matric", "lambda_mat", "k_a_min",
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
  k_frag_litter = "Litter fragmentation", k_frag_CWD = "CWD fragmentation",
  k_l = "LMWC leaching", k_b = "Aggregate breakdown", k_pa = "POC -> aggregate",
  k_ma = "MAOC -> aggregate",
  k_bd = "Microbial turnover", k_a_min = "Min. aeration factor",
  k_exudate_tree = "Tree root exudation",
  pct_claysilt = "Clay + silt (%)",
  a_root_herb = "Herb root allocation", MAT = "Mean annual temperature",
  MAtheta = "Mean annual soil moisture", NPP_herb = "Herbaceous NPP", NPP_tree = "Tree NPP",
  BD = "Bulk density", rho_p = "Particle density", pH = "Soil pH",
  psi_matric = "Matric potential", lambda_mat = "Matric sensitivity",
  CUE_T = "CUE temperature slope",
  K_pl = "Depolymerization half-saturation (POC, litter, CWD)",
  K_lb = "LMWC->microbe half-saturation", K_ld = "LMWC->MAOC max sorption",
  alpha_pl = "Depolymerization pre-exponential (POC, litter, CWD)",
  alpha_lb = "LMWC->microbe pre-exponential",
  p1 = "pH sorption coef. 1", p2 = "pH sorption coef. 2",
  p_a = "Aggregate->POC partition", p_b = "Necromass->MAOC partition",
  p_c = "Clay+silt protection coef.",
  k_mort_root_tree = "Tree root mortality", k_mort_root_herb = "Herb. root mortality",
  T_amp = "Temperature amplitude", theta_amp = "Moisture amplitude",
  theta_amp_rel = "Moisture amplitude (fraction of room)",
  # animal parameters
  k_b_slope_pint = "Earthworm effect on aggregate stability (proportion)",
  c_earthworm_soil = "Earthworm soil consumption", p_earthworm = "Earthworm production efficiency",
  d_earthworm = "Earthworm mortality", prop_feaces_earthworm_LMWC = "Earthworm feces to LMWC fraction",
  a_earthworm_soil = "Earthworm soil assimilation", c_earthworm_litter = "Earthworm litter consumption",
  a_earthworm = "Earthworm litter assimilation efficiency",
  d_detritivores = "Detritivore mortality",
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
#   microbe_collapse : B in the treatment below collapse_rel x baseline
#   kb_floor         : aggregate breakdown k_b*(1 + Earthworm*k_b_slope_pint) at
#                      or below the model's 0.001 floor (breakdown switched off)
collapse_rel <- 1e-6
kb_floor     <- 0.001
# Parameter sets near microbial collapse approach equilibrium very slowly
# (> 1e7 d); the solver stops as soon as it is steady, so a long horizon only
# costs time where it is needed.
eq_max_time <- 1e9

#------------------------------------------------------------------------#
# Forcing used for the steady state (Scripts/4, 5, 7) ----
#------------------------------------------------------------------------#
# "seasonal": the steady state is the annual LIMIT CYCLE under the seasonal
#             climate forcing (R/seasonal_steady_state.R); effects and the
#             realism filter use its ANNUAL MEAN. Realistic; ~10-30x slower.
# "constant": equilibrium under constant mean forcing (runsteady), as before.
# Set it for a whole run from Scripts/0_run_all.R (forcing_run), or for one
# script with options(sfw.forcing = "constant") before sourcing it.
forcing_mode <- getOption("sfw.forcing", "seasonal")
if (!forcing_mode %in% c("seasonal", "constant"))
  stop("option sfw.forcing must be \"seasonal\" or \"constant\"")

# Where Scripts/1, 4, 5 and 7 read and write their results, by forcing, so a
# run under one forcing never overwrites the other. Seasonal results stay in
# Results/ (as before); Scripts/8_compare_forcing.R compares the two.
results_dir_by_forcing <- c(seasonal = "Results", constant = "Results/constant_forcing")
res_dir_forcing <- results_dir_by_forcing[[forcing_mode]]

#------------------------------------------------------------------------#
# Seasonal soil moisture amplitude (Scripts/4 and 5) ----
#------------------------------------------------------------------------#
# Soil moisture swings MAtheta +/- theta_amp over the year and must stay
# within [0, porosity]. Varying MAtheta, theta_amp and BD (porosity =
# 1 - BD / rho_p) independently creates impossible soils. With
# moisture_amp_relative = TRUE the Morris screening and the uncertainty
# propagation vary the amplitude RELATIVE to the room available:
#   theta_amp_rel = theta_amp / min(MAtheta, porosity - MAtheta),  in [0, 1]
# so every sampled amplitude is physically possible; theta_amp is recomputed
# from theta_amp_rel whenever MAtheta or BD change. Its SD is the SD of
# theta_amp (reported or CV rule) divided by the default room, i.e. the same
# CV. The remaining impossible combination (mean moisture >= porosity) is
# removed before any simulation: Morris trajectories containing such a point
# are redrawn (see run_morris()). The identifiability analysis (Scripts/7)
# keeps the absolute amplitude.
moisture_amp_relative <- TRUE
# parameters that decide whether a soil is physically possible
soil_feasibility_params <- c("BD", "rho_p", "MAtheta", "theta_amp", "theta_amp_rel")
# default room per scenario, filled by make_base_pair()
moisture_room_default <- new.env()
lc_tol       <- 1e-6   # limit cycle: max relative |Phi(y) - y| over one year
lc_ode_tol   <- 1e-7   # ODE rtol/atol of the one-year solves (must be < lc_tol)
lc_reuse_jacobian <- FALSE  # FALSE: fresh Jacobian per evaluation (Morris / LHS points
                            # are far from nominal); identifiability uses TRUE
# (the identifiability analysis, Scripts/7, tightens both: it differentiates
#  the solution numerically)

#------------------------------------------------------------------------#
# Baseline realism filter (Scripts/5 and 9) ----
#------------------------------------------------------------------------#
# A parameter sample is "realistic" when its no-animal steady state (annual
# mean if seasonal) keeps every pool positive (> up_pool_floor) and the pools
# listed here inside a target range. Defaults: each scenario's all-default
# baseline, accepted range = reference / factor .. x factor; rows in
# up_targets_file (columns scenario, pool, lower, upper) override them.
up_filter_pools  <- c("TotalC", "P", "M", "B")   # TotalC = all non-animal pools
up_filter_factor <- 2                 # default accepted range: reference / 2 .. x 2
up_filter_factor_by_pool <- c()       # per-pool overrides, e.g. c(B = 5) for / 5 .. x 5
up_targets_file  <- "Data/baseline_targets.csv"   # optional; overrides the defaults
up_pool_floor    <- 1e-9  # every baseline pool must also be above this (g C m-2):
                          # positive and non-zero (1e-9 = numerically zero)

#------------------------------------------------------------------------#
# Output of the Morris screening and Sobol' analysis (Scripts/4 and 9) ----
#------------------------------------------------------------------------#
# "pct": animal effect on total soil + root C as a PERCENTAGE of the
#        no-animal value (default). Across the parameter distributions the
#        baseline soil C spans orders of magnitude (runaway carbon at
#        implausible parameter combinations); the absolute effect scales with
#        it, so its elementary effects and variance are dominated by a few
#        extreme runs. The relative effect is far better behaved.
# "abs": the effect in g C m-2 (the earlier behaviour of Scripts/4).
# The Sobol' step (Scripts/9) should screen on the same output it analyses.
effect_output <- "pct"
