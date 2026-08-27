# ============================================================
# MORRIS SENSITIVITY OF THE ANIMAL EFFECT TO ALL PARAMETERS (equilibrium)
# ------------------------------------------------------------
# One combined Morris elementary-effects screening over BOTH the general model
# parameters AND the animal parameters (the latter only in the scenarios where
# that animal is active). Parameters are varied in COMBINATION along random
# Morris trajectories; per parameter it reports:
#   mu_star  mean |elementary effect|  -> TOTAL sensitivity (rank on this)
#   sigma    sd of elementary effect   -> interactions / non-linearity
#
# ACCEPTABLE RANGE per parameter (one rule for everything):
#   * If scenarios.xlsx reports SD or Min/Max for it -> use that range.
#   * Otherwise (no reported uncertainty -- whether an animal parameter without
#     SD/Min-Max, or a general model parameter absent from the sheet) -> use the
#     BUFFER range [value x 0.5, value x 2] (half to double the current value).
#   Both are clamped to the physical guard rails (efficiencies [0,1], etc.).
#
# The scalar output is the equilibrium animal effect on total soil + root C
# (treatment - baseline). No dynamic spin-up. Run from the project root.
# Output: Results/animal_effect_morris.csv, Results/figures/animal_effect_morris.png
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl, ggrepel)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/derive_millennial_parms.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")
source("R/scenario_uncertainty.R"); source("R/morris_sensitivity.R")

model   <- "millennial"
scen    <- read_scenarios("Data/scenarios.xlsx")
scenarios <- names(scen)

use_fitted_params <- TRUE
fitted_params <- if (use_fitted_params && file.exists("Results/fitted_animal_params.csv"))
  load_fitted_params("Results/fitted_animal_params.csv") else NULL

fig_dir <- "Results/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
res_dir <- "Results";         dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

animal_pools <- c("Earthworm", "Detritivore", "DetPredator", "RootHerb")
derive_fn    <- match.fun(model_table[[model]]$derive)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Morris + range settings ---------------------------------------------
morris_r      <- 50     # trajectories per scenario (raise for stable mu*/sigma)
morris_levels <- 4L     # grid levels p in the Morris design
buffer        <- 0.5    # no-uncertainty range = [value*buffer, value/buffer]
                        #                       = [0.5*value, 2*value] (half..double)

# ------------------------------------------------------------
# GENERAL model parameters to include (animal parameters are added per scenario
# from scenarios.xlsx below). Only those present in a scenario's parms are used.
# ------------------------------------------------------------
sweep_params <- c("k_frag_litter", "k_frag_organic", "k_l_o", "k_l", "k_b", "k_pa",
                  "k_ma", "k_litterfall_ann", "k_litterfall_herb_ann", "k_MICd",
                  "k_bd", "root_to_organic", "a_root_herb", "k_exudate_tree",
                  "k_exudate_herb", "k_frag_CWD", "K_ol", "alpha_ol", "alpha_ob",
                  "K_ob", "p_c", "rho_p", "psi_matric", "lambda_mat", "k_a_min",
                  "alpha_pl", "alpha_lb", "K_pl", "K_lb", "p1", "p2", "K_ld",
                  "CUE_T", "p_a", "p_b", "k_mort_root_tree", "k_mort_root_herb")

# ------------------------------------------------------------
# Pretty display names (code -> label), general + animal parameters combined.
# pretty_param() falls back to the raw code for anything not listed.
# ------------------------------------------------------------
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

pretty_param <- function(x) ifelse(x %in% names(param_labels), param_labels[x], x)

# ------------------------------------------------------------
# PHYSICAL GUARD RAILS. a_*/p_* efficiencies in [0,1] (except soil partition
# coefficients p_a, p_b [0,1] and p_c unbounded); named fractions in [0,1];
# pct_claysilt in [0,100]. clamp_to_bounds() pulls a value onto the nearest edge.
# ------------------------------------------------------------
param_bounds <- function(param) {
  if (grepl("^(a_|p_)", param) && !param %in% c("p_a", "p_b", "p_c"))
    return(c(0, 1))
  switch(param,
         winter_root_act_prop = , root_to_organic = ,
         k_litterfall_ann = , k_litterfall_herb_ann = c(0, 1),
         pct_claysilt = c(0, 100),
         p_a = , p_b = c(0, 1),
         LigFrac = , a_root_herb = , prop_feaces_earthworm_LMWC = c(0, 1),
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
# equilibrium TOTAL animal effect on total soil + root C (scalar, g C m-2).
# Attaches `converged` (TRUE only if BOTH arms reached a stable steady state) so
# the Morris driver can flag evaluations that did not settle.
animal_effect_totalC <- function(pair) {
  pair$baseline  <- spinup_equilibrium(pair$baseline, verbose = FALSE)
  pair$treatment <- spinup_equilibrium(pair$treatment,
                                       warm_start = pair$baseline$init_state_spin,
                                       verbose = FALSE)
  eq_b <- pair$baseline$init_state_spin
  eq_t <- pair$treatment$init_state_spin
  conv <- isTRUE(pair$baseline$spin_info$converged) &&
          isTRUE(pair$treatment$spin_info$converged)
  soil <- setdiff(intersect(names(eq_b), names(eq_t)), animal_pools)
  val  <- sum(eq_t[soil]) - sum(eq_b[soil])
  attr(val, "converged") <- conv
  val
}

# per-Morris-run counter of equilibrium evaluations that did NOT reach a stable
# state. Reset before each morris_run() call; effect_scalar() increments it.
.morris_nonconv <- new.env(parent = emptyenv())
.morris_nonconv$n <- 0L; .morris_nonconv$total <- 0L

# scalar output for a named parameter vector (set all, both arms, then evaluate).
# Returns NA (dropped from the elementary effects) when an arm fails to reach a
# stable equilibrium, and records that in .morris_nonconv.
effect_scalar <- function(base_pair, pv) {
  pair <- base_pair
  for (nm in names(pv)) {
    pair$treatment <- set_param(pair$treatment, nm, pv[[nm]])
    pair$baseline  <- set_param(pair$baseline,  nm, pv[[nm]])
  }
  val <- animal_effect_totalC(pair)
  .morris_nonconv$total <- .morris_nonconv$total + 1L
  if (!isTRUE(attr(val, "converged"))) {
    .morris_nonconv$n <- .morris_nonconv$n + 1L
    return(NA_real_)
  }
  as.numeric(val)
}

# ------------------------------------------------------------
# ACCEPTABLE RANGE per (scenario, parameter): reported SD/Min-Max if available,
# else the buffer range [value*buffer, value/buffer]. Returns lo, hi (clamped),
# and the source label used ("SD", "MinMax", or "buffer").
# ------------------------------------------------------------
unc_all <- tryCatch(read_param_uncertainty("Data/scenarios.xlsx", warn_cv = FALSE),
                    error = function(e) { message("uncertainty read failed: ",
                                                  conditionMessage(e)); NULL })

range_for <- function(scenario, param, default) {
  src <- "buffer"; lo <- default * buffer; hi <- default / buffer
  if (!is.null(unc_all)) {
    r <- unc_all[unc_all$scenario == scenario & unc_all$parameter == param, , drop = FALSE]
    if (nrow(r) && r$unc_source[1] %in% c("SD", "MinMax")) {
      src <- r$unc_source[1]; lo <- r$lo[1]; hi <- r$hi[1]        # reported range
    }
  }
  rng <- sort(c(clamp_to_bounds(lo, param), clamp_to_bounds(hi, param)))
  list(lo = rng[1], hi = rng[2], source = src)
}

# Parameter classes from scenarios.xlsx:
#   * ANIMAL           - Category Animal/Effect (only where that animal is active)
#   * SCENARIO-SPECIFIC- Category Environment/Productivity (non-animal, but still
#                        scenario-specific because the sheet gives per-scenario
#                        values/uncertainty). These were previously dropped.
# Everything else in sweep_params is a GENERAL model parameter.
animal_unc <- if (!is.null(unc_all))
  unc_all[unc_all$category %in% c("Animal", "Effect"), , drop = FALSE] else NULL
scenspec_unc <- if (!is.null(unc_all))
  unc_all[unc_all$category %in% c("Environment", "Productivity"), , drop = FALSE] else NULL

# ------------------------------------------------------------
# MORRIS per scenario over the combined parameter set.
# ------------------------------------------------------------
set.seed(27082026)
morris_rows <- list()
for (scenario in scenarios) {
  base_pair <- setup_scenario_pair(model, scen, scenario)
  if (!is.null(fitted_params))
    base_pair$treatment <- apply_fitted_params(base_pair$treatment, fitted_params,
                                               model, scenario, verbose = FALSE)
  parms_here <- base_pair$treatment$parms

  # animal parameters for this scenario (present in its parms)
  a_here <- if (!is.null(animal_unc)) {
    a <- animal_unc[animal_unc$scenario == scenario, , drop = FALSE]
    a$parameter[a$parameter %in% names(parms_here)]
  } else character(0)
  # scenario-specific non-animal parameters (Environment / Productivity)
  s_here <- if (!is.null(scenspec_unc)) {
    s <- scenspec_unc[scenspec_unc$scenario == scenario, , drop = FALSE]
    s$parameter[s$parameter %in% names(parms_here)]
  } else character(0)
  # general model parameters = sweep_params that are NOT scenario-specific here
  g_here <- setdiff(intersect(sweep_params, names(parms_here)), c(a_here, s_here))
  pnames <- unique(c(g_here, s_here, a_here))

  ptype <- setNames(rep("general", length(pnames)), pnames)
  ptype[pnames %in% s_here] <- "scenario-specific"
  ptype[pnames %in% a_here] <- "animal"

  lo <- hi <- setNames(numeric(length(pnames)), pnames)
  src <- setNames(character(length(pnames)), pnames)
  for (p in pnames) {
    d0 <- parms_here[[p]]
    if (!is.finite(d0) || d0 == 0) { lo[p] <- hi[p] <- NA; next }
    rg <- range_for(scenario, p, d0)
    lo[p] <- rg$lo; hi[p] <- rg$hi; src[p] <- rg$source
  }
  keep <- names(lo)[is.finite(lo) & is.finite(hi) & (hi > lo)]
  lo <- lo[keep]; hi <- hi[keep]; src <- src[keep]

  cat("Morris:", scenario, "-", length(keep), "parameters,", morris_r, "trajectories\n")
  .morris_nonconv$n <- 0L; .morris_nonconv$total <- 0L        # reset before this run
  trajs <- morris_trajectories(length(keep), morris_r, levels = morris_levels)
  ee    <- morris_run(trajs, lo, hi,
                      eval_fun = function(pv) effect_scalar(base_pair, pv),
                      verbose = FALSE)
  n_bad <- .morris_nonconv$n; n_tot <- .morris_nonconv$total
  if (n_bad > 0)
    warning(sprintf("Morris %s: %d of %d equilibrium evaluations did NOT reach a stable state (%.1f%%).",
                    scenario, n_bad, n_tot, 100 * n_bad / max(n_tot, 1)))
  cat(sprintf("  stable-state check: %d/%d evaluations converged%s\n",
              n_tot - n_bad, n_tot, if (n_bad > 0) "  <-- SEE WARNING" else ""))

  sm    <- morris_summary(ee)
  sm$scenario     <- scenario
  sm$range_source <- src[sm$parameter]
  sm$param_type   <- ptype[sm$parameter]
  sm$is_animal    <- sm$parameter %in% a_here
  sm$n_nonconverged <- n_bad
  sm$n_evaluations  <- n_tot
  morris_rows[[scenario]] <- sm
}

morris_tbl <- bind_rows(morris_rows) %>%
  mutate(parameter_label = pretty_param(parameter)) %>%
  group_by(scenario) %>% arrange(scenario, desc(mu_star)) %>%
  mutate(rank = row_number()) %>% ungroup() %>%
  select(scenario, parameter, parameter_label, param_type, is_animal,
         range_source, mu_star, sigma, n_ee, rank, n_nonconverged, n_evaluations)
write_csv(morris_tbl, file.path(res_dir, "animal_effect_morris.csv"))

# ------------------------------------------------------------
# FIGURE: traditional Morris mu* vs sigma map, one facet per scenario. Colour =
# how the range was set (reported SD / Min-Max, or the half-double buffer);
# shape = general vs animal parameter.
# ------------------------------------------------------------
p_morris <- ggplot(morris_tbl, aes(mu_star, sigma)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(aes(colour = range_source, shape = param_type), alpha = 0.85, size = 2) +
  ggrepel::geom_text_repel(aes(label = parameter_label), size = 2.2, max.overlaps = 14) +
  facet_wrap(~scenario, scales = "free") +
  scale_colour_manual(values = c(SD = "#1b7837", MinMax = "#2166ac", buffer = "#b2182b"),
                      labels = c(SD = "Two standard deviations", MinMax = "Min/Max", buffer = "50% to 150%"),
                      name = "Range source") +
  scale_shape_manual(values = c(general = 16, `scenario-specific` = 15, animal = 17),
                     name = "Parameter type") +
  labs(
    # title = "Morris screening of the animal effect (all parameters, varied in combination)",
    #    subtitle = "mu* = total sensitivity (mean |EE|); sigma = interactions / non-linearity.  buffer range = half to double.",
       x = expression(mu*"* (mean |elementary effect|, g C m"^-2*")"),
       y = expression(sigma*" (sd of elementary effect)")) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "animal_effect_morris.png"), p_morris,
       width = 13, height = 9, dpi = 150)
print(p_morris)

p_morris_animal <- ggplot(morris_tbl %>% filter(is_animal), aes(mu_star, sigma)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(aes(colour = range_source, shape = param_type), alpha = 0.85, size = 2) +
  ggrepel::geom_text_repel(aes(label = parameter_label), size = 2.2, max.overlaps = 14) +
  facet_wrap(~scenario, scales = "free") +
  scale_colour_manual(values = c(SD = "#1b7837", MinMax = "#2166ac", buffer = "#b2182b"),
                      labels = c(SD = "Standard deviation", MinMax = "Min/Max", buffer = "50% to 150%"),
                      name = "Range source") +
  scale_shape_manual(values = c(general = 16, `scenario-specific` = 15, animal = 17),
                     name = "Parameter type") +
  labs(
    # title = "Morris screening of the animal effect (all parameters, varied in combination)",
    #    subtitle = "mu* = total sensitivity (mean |EE|); sigma = interactions / non-linearity.  buffer range = half to double.",
    x = expression(mu*"* (mean |elementary effect|, g C m"^-2*")"),
    y = expression(sigma*" (sd of elementary effect)")) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "animal_effect_morris_animal.png"), p_morris_animal,
       width = 13, height = 9, dpi = 150)

cat("\nWrote:\n  ", file.path(res_dir, "animal_effect_morris.csv"),
    "\n  ", file.path(fig_dir, "animal_effect_morris.png"),
    "\n  ", file.path(fig_dir, "animal_effect_morris_animal.png"), "\n")