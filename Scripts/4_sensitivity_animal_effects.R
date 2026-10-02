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
#   parameters sampled jointly by Latin hypercube over the MAIN ranges; the
#   total and direct equilibrium animal effects are recomputed per sample.
#   Outputs: animal_effect_uncertainty_samples.csv / _summary.csv,
#            figures/animal_effect_uncertainty.png
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

animal_pools <- c("Earthworm", "Detritivore", "RootHerb")
derive_fn    <- match.fun(model_table[[model]]$derive)
`%||%` <- function(a, b) if (is.null(a)) b else a

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
up_n        <- 500        # Latin-hypercube samples per scenario (3 equilibria each)
up_seed     <- 20260827

#------------------------------------------------------------------------#
# Parallel execution ----
#------------------------------------------------------------------------#
# Model evaluations (Morris trajectories, uncertainty samples) run on n_cores
# worker processes; 1 = serial. The bootstrap itself is vectorised and takes
# seconds, so it stays on the main process.
n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)

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

pretty_param <- function(x) ifelse(x %in% names(param_labels), param_labels[x], x)

#------------------------------------------------------------------------#
# PHYSICAL GUARD RAILS ----
#------------------------------------------------------------------------#
# PHYSICAL GUARD RAILS. a_*/p_* efficiencies in [0,1] (except soil partition
# coefficients p_a, p_b [0,1] and p_c unbounded); named fractions in [0,1];
# pct_claysilt in [0,100]. clamp_to_bounds() pulls a value onto the nearest edge.
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

#------------------------------------------------------------------------#
# set_param() ----
#------------------------------------------------------------------------#
# set_param(): change ONE parameter and rebuild derived params + forcing.
set_param <- function(obj, param, value) {
  obj$parms[[param]]        <- value
  obj$parms                 <- derive_fn(obj$parms)
  obj$parms$climate_forcing <- make_climate_forcing(obj$parms)
  obj
}

# infeasible_soil(): parameter sets that describe an impossible soil. Porosity
# (1 - BD / rho_p) must be positive and exceed the mean soil moisture, otherwise
# the moisture functions take the square root of a negative number. Returns ""
# or a short reason; such evaluations are skipped and reported, not solved.
infeasible_soil <- function(parms) {
  phi <- parms$phi_por
  if (!is.finite(phi) || phi <= 0) return("infeasible soil: porosity <= 0 (BD >= rho_p)")
  if (phi <= parms$MAtheta)        return("infeasible soil: porosity <= mean soil moisture")
  ""
}

# regime_flags(): which special regime (if any) an equilibrium pair is in.
# Returns "" (normal), "microbe_collapse", "kb_floor" or both joined by "+".
regime_flags <- function(eq_b, eq_t, parms) {
  f   <- character(0)
  mic <- intersect(c("MIC", "B"), intersect(names(eq_b), names(eq_t)))
  if (length(mic) && any(is.finite(eq_t[mic]) &
                         eq_t[mic] < collapse_rel * pmax(eq_b[mic], 1e-12)))
    f <- c(f, "microbe_collapse")
  if ("Earthworm" %in% names(eq_t)) {
    kb <- parms$k_b + eq_t[["Earthworm"]] * parms$k_b_slope_pint * parms$k_b
    if (is.finite(kb) && kb <= kb_floor) f <- c(f, "kb_floor")
  }
  paste(f, collapse = "+")
}

# equilibrium TOTAL animal effect on total soil + root C (scalar, g C m-2).
# Attaches `converged` (TRUE only if BOTH arms reached a stable steady state)
# and `regime` (see regime_flags) for the Morris driver.
animal_effect_totalC <- function(pair) {
  pair$baseline  <- spinup_equilibrium(pair$baseline, max_time = eq_max_time,
                                       verbose = FALSE)
  pair$treatment <- spinup_equilibrium(pair$treatment,
                                       warm_start = pair$baseline$init_state_spin,
                                       max_time = eq_max_time, verbose = FALSE)
  eq_b <- pair$baseline$init_state_spin
  eq_t <- pair$treatment$init_state_spin
  conv <- isTRUE(pair$baseline$spin_info$converged) &&
          isTRUE(pair$treatment$spin_info$converged)
  soil <- setdiff(intersect(names(eq_b), names(eq_t)), animal_pools)
  val  <- sum(eq_t[soil]) - sum(eq_b[soil])
  attr(val, "converged") <- conv
  attr(val, "regime")    <- if (conv) regime_flags(eq_b, eq_t, pair$treatment$parms) else ""
  val
}

# make_effect_fun(): the Morris evaluator for one scenario. Built by a
# top-level factory so its environment carries only base_pair to the workers.
make_effect_fun <- function(base_pair) {
  force(base_pair)            # evaluate now so workers receive the object itself
  function(pv) effect_scalar(base_pair, pv)
}

# scalar output for a named parameter vector (set all, both arms, then
# evaluate). The converged / regime attributes travel with the value;
# morris_run() drops non-converged points and tags flagged ones.
effect_scalar <- function(base_pair, pv) {
  pair <- base_pair
  for (nm in names(pv)) {
    pair$treatment <- set_param(pair$treatment, nm, pv[[nm]])
    pair$baseline  <- set_param(pair$baseline,  nm, pv[[nm]])
  }
  why <- infeasible_soil(pair$baseline$parms)
  if (nzchar(why)) return(structure(NA_real_, converged = FALSE, error = why))
  animal_effect_totalC(pair)
}

#------------------------------------------------------------------------#
# PARAMETER RANGES ----
#------------------------------------------------------------------------#
# PARAMETER RANGES. range_for() returns lo, hi (clamped to guard rails) and a
# label for how the range was set, under either mode:
#   mode = "main"         reported +/-2 SD or [Min, Max] if available, else the
#                         half-to-double buffer.            labels: SD / MinMax / buffer
#   mode = "standardized" 50-200% for every parameter, cut to the reported
#                         Min/Max wherever that is narrower.
#                         labels: 50-200% / MinMax-bound / guard-rail-bound /
#                                 MinMax-excludes-default (Min/Max does not
#                                 overlap 50-200%; the unbound 50-200% is kept)
unc_all <- tryCatch(read_param_uncertainty("Data/scenarios.xlsx", warn_cv = FALSE),
                    error = function(e) { message("uncertainty read failed: ",
                                                  conditionMessage(e)); NULL })

unc_row <- function(scenario, param) {
  if (is.null(unc_all)) return(NULL)
  r <- unc_all[unc_all$scenario == scenario & unc_all$parameter == param, , drop = FALSE]
  if (nrow(r)) r[1, , drop = FALSE] else NULL
}

range_for <- function(scenario, param, default, mode = c("main", "standardized")) {
  mode <- match.arg(mode)
  r    <- unc_row(scenario, param)

  if (mode == "main") {
    src <- "buffer"; lo <- default * buffer; hi <- default / buffer
    if (!is.null(r) && r$unc_source %in% c("SD", "MinMax")) {
      src <- r$unc_source; lo <- r$lo; hi <- r$hi                 # reported range
    }
    rng <- sort(c(clamp_to_bounds(lo, param), clamp_to_bounds(hi, param)))
    return(list(lo = rng[1], hi = rng[2], source = src))
  }

  # standardized: 50-200% (sorted so negative defaults stay ordered)
  rng <- sort(default * supp_frac); lo <- rng[1]; hi <- rng[2]
  src <- "50-200%"
  if (!is.null(r) && is.finite(r$min) && is.finite(r$max) && r$max > r$min) {
    lo2 <- max(lo, r$min); hi2 <- min(hi, r$max)
    if (hi2 > lo2) {
      if (lo2 > lo || hi2 < hi) src <- "MinMax-bound"
      lo <- lo2; hi <- hi2
    } else {
      src <- "MinMax-excludes-default"                            # keep 50-200%
    }
  }
  clo <- clamp_to_bounds(lo, param); chi <- clamp_to_bounds(hi, param)
  if ((clo != lo || chi != hi) && src == "50-200%") src <- "guard-rail-bound"
  rng <- sort(c(clo, chi))
  list(lo = rng[1], hi = rng[2], source = src)
}

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
# param_set_for() ----
#------------------------------------------------------------------------#
# param_set_for(): the parameters screened in one scenario and their class
# (general / scenario-specific / animal). Shared by Morris and the uncertainty
# propagation so both use exactly the same parameter set.
param_set_for <- function(scenario, parms_here) {
  a_here <- if (!is.null(animal_unc)) {
    a <- animal_unc[animal_unc$scenario == scenario, , drop = FALSE]
    a$parameter[a$parameter %in% names(parms_here)]
  } else character(0)
  s_here <- if (!is.null(scenspec_unc)) {
    s <- scenspec_unc[scenspec_unc$scenario == scenario, , drop = FALSE]
    s$parameter[s$parameter %in% names(parms_here)]
  } else character(0)
  g_here <- setdiff(intersect(sweep_params, names(parms_here)), c(a_here, s_here))
  pnames <- unique(c(g_here, s_here, a_here))
  ptype  <- setNames(rep("general", length(pnames)), pnames)
  ptype[pnames %in% s_here] <- "scenario-specific"
  ptype[pnames %in% a_here] <- "animal"
  list(pnames = pnames, ptype = ptype, animal = a_here)
}

# base scenario pair with the fitted animal parameters applied
make_base_pair <- function(scenario) {
  bp <- setup_scenario_pair(model, scen, scenario)
  if (!is.null(fitted_params))
    bp$treatment <- apply_fitted_params(bp$treatment, fitted_params,
                                        model, scenario, verbose = FALSE)
  bp
}

#------------------------------------------------------------------------#
# run_morris() ----
#------------------------------------------------------------------------#
# run_morris(): Morris screening for every scenario under one range mode.
# Each scenario is seeded with morris_seed + its index, so the main and
# standardized runs share the same unit-hypercube trajectories.
run_morris <- function(mode = c("main", "standardized")) {
  mode <- match.arg(mode)
  morris_rows <- list()
  ee_all      <- list()
  for (si in seq_along(scenarios)) {
    scenario  <- scenarios[si]
    base_pair  <- make_base_pair(scenario)
    parms_here <- base_pair$treatment$parms

    ps     <- param_set_for(scenario, parms_here)
    pnames <- ps$pnames; ptype <- ps$ptype; a_here <- ps$animal

    lo <- hi <- setNames(numeric(length(pnames)), pnames)
    src <- setNames(character(length(pnames)), pnames)
    for (p in pnames) {
      d0 <- parms_here[[p]]
      if (!is.finite(d0) || d0 == 0) { lo[p] <- hi[p] <- NA; next }
      rg <- range_for(scenario, p, d0, mode)
      lo[p] <- rg$lo; hi[p] <- rg$hi; src[p] <- rg$source
    }
    keep <- names(lo)[is.finite(lo) & is.finite(hi) & (hi > lo)]
    lo <- lo[keep]; hi <- hi[keep]; src <- src[keep]

    cat("Morris [", mode, "]: ", scenario, " - ", length(keep), " parameters, ",
        morris_r, " trajectories\n", sep = "")
    set.seed(morris_seed + si)
    trajs <- morris_trajectories(length(keep), morris_r, levels = morris_levels)
    ee    <- morris_run(trajs, lo, hi,
                        eval_fun = make_effect_fun(base_pair),
                        verbose = FALSE, cl = cl)
    n_bad <- attr(ee, "n_nonconverged"); n_tot <- attr(ee, "n_eval")
    # impossible soils are skipped before solving: count them separately
    sr    <- attr(ee, "skip_reasons")
    n_inf <- if (length(sr)) sum(sr[grepl("^infeasible soil", names(sr))]) else 0L
    n_bad <- n_bad - n_inf
    if (n_inf > 0)
      message(sprintf("Morris [%s] %s: %d of %d points are unrealistic soils (porosity <= 0 or <= soil moisture) and were removed.",
                      mode, scenario, n_inf, n_tot))
    rc    <- attr(ee, "regime_counts")
    n_mic <- if ("microbe_collapse" %in% names(rc)) rc[["microbe_collapse"]] else 0L
    n_kb  <- if ("kb_floor" %in% names(rc)) rc[["kb_floor"]] else 0L
    if (n_mic + n_kb > 0)
      message(sprintf("Morris [%s] %s: flagged evaluations - microbe_collapse %d, kb_floor %d (of %d).",
                      mode, scenario, n_mic, n_kb, n_tot))
    if (n_bad > 0)
      warning(sprintf("Morris [%s] %s: %d of %d equilibrium evaluations did NOT reach a stable state (%.1f%%).",
                      mode, scenario, n_bad, n_tot, 100 * n_bad / max(n_tot, 1)))
    out_rng <- keep[!(vapply(keep, function(p) parms_here[[p]], numeric(1)) >= lo &
                      vapply(keep, function(p) parms_here[[p]], numeric(1)) <= hi)]
    if (length(out_rng))
      warning(sprintf("Morris [%s] %s: model value lies OUTSIDE its range for: %s",
                      mode, scenario, paste(out_rng, collapse = ", ")))
    cat(sprintf("  stable-state check: %d/%d realistic evaluations converged%s; %d unrealistic soils removed\n",
                n_tot - n_inf - n_bad, n_tot - n_inf, if (n_bad > 0) "  <-- SEE WARNING" else "", n_inf))

    sm <- morris_summary(ee)
    bt <- morris_bootstrap(ee, B = boot_B, top_k = boot_top_k,
                           seed = morris_seed + 1000 + si)
    if (!is.null(bt)) sm <- merge(sm, bt, by = "parameter", all.x = TRUE, sort = FALSE)
    # same metrics using only steps where neither end is in a flagged regime
    ee_clean <- ee[!nzchar(ee$regime), , drop = FALSE]
    sm_clean <- if (nrow(ee_clean)) morris_summary(ee_clean) else
      data.frame(parameter = character(), mu_star = numeric(), sigma = numeric(), n_ee = integer())
    names(sm_clean)[names(sm_clean) != "parameter"] <-
      paste0(names(sm_clean)[names(sm_clean) != "parameter"], "_clean")
    sm <- merge(sm, sm_clean, by = "parameter", all.x = TRUE, sort = FALSE)
    sm$n_ee_flagged <- vapply(sm$parameter, function(p)
      sum(ee$parameter == p & nzchar(ee$regime)), integer(1))
    ee$scenario <- scenario
    ee_all[[scenario]] <- ee
    sm$scenario       <- scenario
    sm$range_source   <- src[sm$parameter]
    sm$range_lo       <- lo[sm$parameter]
    sm$range_hi       <- hi[sm$parameter]
    sm$default        <- vapply(sm$parameter, function(p) parms_here[[p]], numeric(1))
    # flag ranges that do not contain the value actually used in the model (e.g. a
    # fitted animal parameter whose reported SD/Min-Max was built around the
    # sheet value, not the fitted value)
    sm$default_in_range <- sm$default >= sm$range_lo & sm$default <= sm$range_hi
    sm$param_type     <- ptype[sm$parameter]
    sm$is_animal      <- sm$parameter %in% a_here
    # elementary effects lost because an end point was removed or unsolved
    sm$n_ee_removed <- vapply(sm$parameter, function(p)
      sum(ee$parameter == p & !is.finite(ee$ee)), integer(1))
    sm$n_eval_infeasible <- n_inf
    sm$n_nonconverged <- n_bad
    sm$n_evaluations  <- n_tot
    sm$n_eval_microbe_collapse <- n_mic
    sm$n_eval_kb_floor         <- n_kb
    morris_rows[[scenario]] <- sm
  }

  tbl <- bind_rows(morris_rows) %>%
    mutate(parameter_label = pretty_param(parameter)) %>%
    group_by(scenario) %>% arrange(scenario, desc(mu_star)) %>%
    mutate(rank = row_number()) %>%
    mutate(rank_clean = rank(-mu_star_clean, ties.method = "first", na.last = "keep")) %>%
    ungroup() %>%
    select(scenario, parameter, parameter_label, param_type, is_animal,
           range_source, default, range_lo, range_hi, default_in_range,
           mu_star, mu_star_lo, mu_star_hi, sigma, sigma_lo, sigma_hi,
           n_ee, rank, rank_median, rank_lo, rank_hi, p_top_k,
           mu_star_clean, sigma_clean, n_ee_clean, rank_clean, n_ee_flagged,
           n_ee_removed, n_eval_infeasible, n_nonconverged, n_evaluations,
           n_eval_microbe_collapse, n_eval_kb_floor)
  attr(tbl, "ee") <- bind_rows(ee_all)          # raw elementary effects
  tbl
}

#------------------------------------------------------------------------#
# FIGURE helper: Morris mu* vs sigma map ----
#------------------------------------------------------------------------#
# FIGURE helper: traditional Morris mu* vs sigma map, one facet per scenario.
# Colour = how the range was set; shape = parameter type.
plot_morris <- function(tbl, colours, labels, legend_name = "Range source") {
  ggplot(tbl, aes(mu_star, sigma)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
    geom_point(aes(colour = range_source, shape = param_type), alpha = 0.85, size = 2) +
    ggrepel::geom_text_repel(aes(label = parameter_label), size = 2.2, max.overlaps = 14) +
    facet_wrap(~scenario, scales = "free", labeller = scenario_labeller) +
    scale_colour_manual(values = colours, labels = labels, name = legend_name, drop = FALSE) +
    scale_shape_manual(values = c(general = 16, `scenario-specific` = 15, animal = 17),
                       name = "Parameter type") +
    labs(x = expression(mu*"* (mean |elementary effect|, g C m"^-2*")"),
         y = expression(sigma*" (sd of elementary effect)")) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
}

# order a label within each facet (suffix trick, stripped again on the axis)
tidytext_reorder <- function(x, by, within) {
  key <- paste(x, within, sep = "___")
  factor(key, levels = unique(key[order(within, by)]))
}

main_cols <- c(SD = "#1b7837", MinMax = "#2166ac", buffer = "#b2182b")
main_labs <- c(SD = "+/- 2 standard deviations", MinMax = "Min/Max",
               buffer = "50% to 200% (half to double)")
supp_cols <- c(`50-200%` = "#b2182b", `MinMax-bound` = "#2166ac",
               `guard-rail-bound` = "#762a83", `MinMax-excludes-default` = "#e08214")
supp_labs <- c(`50-200%` = "50% to 200% (half to double)",
               `MinMax-bound` = "50-200% cut to Min/Max",
               `guard-rail-bound` = "50-200% cut to physical bound",
               `MinMax-excludes-default` = "50-200% (Min/Max excludes default)")

write_morris_outputs <- function(tbl, tag, colours, labels) {
  csv <- file.path(res_dir, paste0("animal_effect_morris", tag, ".csv"))
  write_csv(tbl, csv)
  p_all <- plot_morris(tbl, colours, labels)
  ggsave(file.path(fig_dir, paste0("animal_effect_morris", tag, ".png")), p_all,
         width = 13, height = 9, dpi = 150)
  p_an <- plot_morris(tbl %>% filter(is_animal), colours, labels)
  ggsave(file.path(fig_dir, paste0("animal_effect_morris", tag, "_animal.png")), p_an,
         width = 13, height = 9, dpi = 150)
  # raw elementary effects (lets the bootstrap / summaries be recomputed later)
  ee <- attr(tbl, "ee")
  if (!is.null(ee))
    write_csv(ee, file.path(res_dir, paste0("animal_effect_morris", tag, "_elementary_effects.csv")))

  # bootstrap figure: mu* with 95% CI for the top 15 parameters per scenario
  top <- tbl %>% group_by(scenario) %>% slice_min(rank, n = 15) %>% ungroup() %>%
    mutate(lab = tidytext_reorder(parameter_label, mu_star, scenario))
  p_ci <- ggplot(top, aes(mu_star, lab, colour = range_source, shape = param_type)) +
    geom_errorbar(aes(xmin = mu_star_lo, xmax = mu_star_hi), width = 0.25, orientation = "y") +
    geom_point(size = 2) +
    facet_wrap(~scenario, scales = "free", labeller = scenario_labeller) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_colour_manual(values = colours, labels = labels, name = "Range source", drop = FALSE) +
    scale_shape_manual(values = c(general = 16, `scenario-specific` = 15, animal = 17),
                       name = "Parameter type") +
    labs(x = expression(mu*"* with "*95*"% bootstrap CI (g C m"^-2*")"), y = NULL) +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, paste0("animal_effect_morris", tag, "_bootstrap.png")), p_ci,
         width = 13, height = 9, dpi = 150)

  # convergence verdict per scenario: is the top-k set stable across replicates?
  conv <- tbl %>% group_by(scenario) %>% slice_min(rank, n = boot_top_k) %>%
    summarise(min_p_top_k    = min(p_top_k, na.rm = TRUE),
              max_rel_ci_mu  = max((mu_star_hi - mu_star_lo) / mu_star, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(top_k_stable = min_p_top_k >= 0.8)
  cat("\nMorris bootstrap convergence (", tag, "): top-", boot_top_k,
      " stable if every member is in the top-", boot_top_k, " in >= 80% of replicates\n", sep = "")
  print(as.data.frame(conv), row.names = FALSE)
  if (any(!conv$top_k_stable))
    warning("Morris ranking not converged for: ",
            paste(conv$scenario[!conv$top_k_stable], collapse = ", "),
            " -- increase morris_r.")

  if (interactive()) print(tbl %>% group_by(scenario) %>% slice_max(mu_star, n = 5))
  cat("Wrote:", csv, "and figures animal_effect_morris", tag, "[ _animal].png\n", sep = "")
}

#------------------------------------------------------------------------#
# RUN ----
#------------------------------------------------------------------------#
# start the worker pool (NULL = serial). Workers source the R/ helpers and get
# the script-level objects that the model evaluations need.
cl <- morris_cluster(
  n_cores,
  source_files = c("R/climate_forcing.R", "R/spinup.R", "R/derive_millennial_parms.R",
                   "R/millennial_model.R", "R/init_millennial_state.R",
                   "R/setup.R", "R/compare_functions.R", "R/fit_animals.R"),
  export = c("set_param", "derive_fn", "animal_pools", "regime_flags",
             "animal_effect_totalC", "effect_scalar", "collapse_rel", "kb_floor",
             "eq_max_time", "infeasible_soil"))
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
# UNCERTAINTY PROPAGATION TO THE HEADLINE ANIMAL EFFECT ----
#------------------------------------------------------------------------#
# Morris ranks parameters; it does not say how uncertain the animal effect
# itself is. Here ALL screened parameters are sampled jointly (Latin hypercube,
# up_n samples per scenario) over the MAIN-text ranges, and for every sample the
# equilibrium TOTAL and DIRECT animal effects on total soil + root C are
# recomputed. Sampling distributions by range source:
#   SD      normal(sheet value, SD) truncated to value +/- 2 SD (and guard rails)
#   MinMax  uniform on [Min, Max]
#   buffer  log-uniform on [value/2, 2*value] (symmetric in ratio around value)
# Parameters are sampled independently (no correlation structure).
# Outputs: Results/animal_effect_uncertainty_samples.csv  (every sample)
#          Results/animal_effect_uncertainty_summary.csv  (per scenario)
#          Results/figures/animal_effect_uncertainty.png

# inverse-CDF sampler for one parameter given u in (0,1)
up_quantile <- function(u, rg, r) {
  lo <- rg$lo; hi <- rg$hi
  if (rg$source == "SD" && !is.null(r) && is.finite(r$sd) && r$sd > 0) {
    a <- pnorm(lo, r$value, r$sd); b <- pnorm(hi, r$value, r$sd)
    return(qnorm(a + u * (b - a), r$value, r$sd))
  }
  if (rg$source == "buffer" && lo * hi > 0) {               # same sign: log-uniform
    s <- sign(lo); l <- log(abs(lo)); h <- log(abs(hi))
    v <- s * exp(l + u * (h - l))
    return(v)
  }
  lo + u * (hi - lo)                                         # uniform
}

# Latin hypercube in (0,1)^k: one stratified draw per row in each column
lhs_unit <- function(n, k) {
  vapply(seq_len(k), function(j) (sample.int(n) - runif(n)) / n, numeric(n))
}

# total and direct equilibrium effects for one parameter vector
headline_effects <- function(base_pair, pv) {
  pair <- base_pair
  for (nm in names(pv)) {
    pair$treatment <- set_param(pair$treatment, nm, pv[[nm]])
    pair$baseline  <- set_param(pair$baseline,  nm, pv[[nm]])
  }
  out <- c(total = NA_real_, direct = NA_real_, converged = 0,
           microbe_collapse = 0, kb_floor = 0, infeasible = 0)
  if (nzchar(infeasible_soil(pair$baseline$parms))) { out[["infeasible"]] <- 1; return(out) }
  res <- tryCatch({
    b  <- spinup_equilibrium(pair$baseline, max_time = eq_max_time, verbose = FALSE)
    eb <- b$init_state_spin
    tr <- spinup_equilibrium(pair$treatment, warm_start = eb, max_time = eq_max_time,
                             verbose = FALSE)
    td <- spinup_equilibrium(zero_indirect_effects(pair$treatment, verbose = FALSE),
                             warm_start = eb, max_time = eq_max_time, verbose = FALSE)
    soil <- setdiff(intersect(names(eb), names(tr$init_state_spin)), animal_pools)
    conv <- isTRUE(b$spin_info$converged) && isTRUE(tr$spin_info$converged) &&
            isTRUE(td$spin_info$converged)
    # regime flags: total-effect state (both flags) and direct-effect state
    # (microbes only; its indirect slopes are zero, so no k_b floor)
    rg_t <- regime_flags(eb, tr$init_state_spin, pair$treatment$parms)
    rg_d <- regime_flags(eb, td$init_state_spin, pair$treatment$parms)
    c(total  = sum(tr$init_state_spin[soil]) - sum(eb[soil]),
      direct = sum(td$init_state_spin[soil]) - sum(eb[soil]),
      converged = as.numeric(conv),
      microbe_collapse = as.numeric(grepl("microbe_collapse", paste(rg_t, rg_d))),
      kb_floor         = as.numeric(grepl("kb_floor", rg_t)),
      infeasible       = 0)
  }, error = function(e) out)
  res
}

# factory so a parallel worker receives only base_pair with the evaluator
make_headline_fun <- function(base_pair) {
  force(base_pair)
  function(pv) headline_effects(base_pair, pv)
}

run_uncertainty_propagation_fn <- function() {
  samp_rows <- list(); summ_rows <- list()
  for (si in seq_along(scenarios)) {
    scenario   <- scenarios[si]
    base_pair  <- make_base_pair(scenario)
    parms_here <- base_pair$treatment$parms
    ps         <- param_set_for(scenario, parms_here)

    rgs <- list(); rows <- list()
    for (p in ps$pnames) {
      d0 <- parms_here[[p]]
      if (!is.finite(d0) || d0 == 0) next
      rg <- range_for(scenario, p, d0, "main")
      if (is.finite(rg$lo) && is.finite(rg$hi) && rg$hi > rg$lo) {
        rgs[[p]] <- rg; rows[[p]] <- unc_row(scenario, p)
      }
    }
    pn <- names(rgs); k <- length(pn)

    cat("Uncertainty propagation: ", scenario, " - ", k, " parameters, ",
        up_n, " LHS samples\n", sep = "")
    set.seed(up_seed + si)
    U <- lhs_unit(up_n, k)
    X <- vapply(seq_len(k), function(j) up_quantile(U[, j], rgs[[pn[j]]], rows[[pn[j]]]),
                numeric(up_n))
    if (up_n == 1) X <- matrix(X, nrow = 1)
    colnames(X) <- pn

    default <- headline_effects(base_pair, setNames(numeric(0), character(0)))
    rows_list <- lapply(seq_len(up_n), function(i) setNames(X[i, ], pn))
    Y <- if (is.null(cl)) lapply(rows_list, function(pv) headline_effects(base_pair, pv)) else
      parallel::parLapply(cl, rows_list, make_headline_fun(base_pair))
    Y <- do.call(rbind, Y)

    n_inf <- sum(Y[, "infeasible"] == 1)
    n_bad <- sum(Y[, "converged"] < 1) - n_inf
    if (n_inf > 0)
      warning(sprintf("Uncertainty propagation %s: %d of %d samples describe an impossible soil (porosity <= moisture); excluded.",
                      scenario, n_inf, up_n))
    if (n_bad > 0)
      warning(sprintf("Uncertainty propagation %s: %d of %d samples did NOT reach a stable state; excluded.",
                      scenario, n_bad, up_n))
    cat(sprintf("  stable-state check: %d/%d realistic samples converged%s; %d unrealistic soils removed\n",
                up_n - n_inf - n_bad, up_n - n_inf, if (n_bad > 0) "  <-- SEE WARNING" else "", n_inf))

    s_df <- data.frame(scenario = scenario, sample = seq_len(up_n), Y, X,
                       check.names = FALSE)
    samp_rows[[scenario]] <- s_df

    ok   <- Y[, "converged"] == 1 & is.finite(Y[, "total"]) & is.finite(Y[, "direct"])
    flag <- Y[, "microbe_collapse"] == 1 | Y[, "kb_floor"] == 1
    n_mic <- sum(ok & Y[, "microbe_collapse"] == 1); n_kb <- sum(ok & Y[, "kb_floor"] == 1)
    if (n_mic + n_kb > 0)
      message(sprintf("Uncertainty propagation %s: flagged samples - microbe_collapse %d, kb_floor %d (of %d converged).",
                      scenario, n_mic, n_kb, sum(ok)))
    tot <- Y[ok, "total"]; dir <- Y[ok, "direct"]
    tot_c <- Y[ok & !flag, "total"]; dir_c <- Y[ok & !flag, "direct"]
    pct <- 100 * dir / tot
    q   <- function(v, p) if (length(v)) unname(stats::quantile(v, p, na.rm = TRUE)) else NA_real_
    summ_rows[[scenario]] <- data.frame(
      scenario = scenario, n_samples = up_n, n_converged = sum(ok),
      n_infeasible = n_inf, n_microbe_collapse = n_mic, n_kb_floor = n_kb,
      n_clean = sum(ok & !flag),
      total_default  = default[["total"]],
      total_median   = q(tot, 0.5), total_mean = mean(tot),
      total_q025     = q(tot, 0.025), total_q975 = q(tot, 0.975),
      p_total_positive = mean(tot > 0),
      p_sign_as_default = mean(sign(tot) == sign(default[["total"]])),
      direct_default = default[["direct"]],
      direct_median  = q(dir, 0.5),
      direct_q025    = q(dir, 0.025), direct_q975 = q(dir, 0.975),
      pct_direct_default = 100 * default[["direct"]] / default[["total"]],
      pct_direct_median  = q(pct[is.finite(pct)], 0.5),
      pct_direct_q025    = q(pct[is.finite(pct)], 0.025),
      pct_direct_q975    = q(pct[is.finite(pct)], 0.975),
      total_median_clean = q(tot_c, 0.5),
      total_q025_clean   = q(tot_c, 0.025), total_q975_clean = q(tot_c, 0.975),
      p_sign_as_default_clean = if (length(tot_c)) mean(sign(tot_c) == sign(default[["total"]])) else NA_real_,
      direct_median_clean = q(dir_c, 0.5),
      direct_q025_clean   = q(dir_c, 0.025), direct_q975_clean = q(dir_c, 0.975))
  }
  list(samples = bind_rows(samp_rows), summary = bind_rows(summ_rows))
}

if (run_uncertainty_propagation) {
  if (!is.null(cl)) parallel::clusterExport(cl, "headline_effects")   # defined after cluster start
  up <- run_uncertainty_propagation_fn()
  write_csv(up$samples, file.path(res_dir, "animal_effect_uncertainty_samples.csv"))
  write_csv(up$summary, file.path(res_dir, "animal_effect_uncertainty_summary.csv"))
  cat("\nHeadline animal effect with propagated parameter uncertainty (g C m-2):\n")
  print(up$summary %>%
          transmute(scenario, n_converged,
                    total = sprintf("%.1f (%.1f to %.1f)", total_median, total_q025, total_q975),
                    direct = sprintf("%.1f (%.1f to %.1f)", direct_median, direct_q025, direct_q975),
                    p_sign_as_default = round(p_sign_as_default, 3),
                    flagged = n_microbe_collapse + n_kb_floor,
                    total_unflagged = sprintf("%.1f (%.1f to %.1f)", total_median_clean,
                                              total_q025_clean, total_q975_clean)),
        row.names = FALSE)

  # FIGURE: distribution of total and direct effects per scenario; solid line =
  # all-default effect, shaded band = central 95% of the samples.
  long <- up$samples %>% filter(converged == 1) %>%
    select(scenario, total, direct, microbe_collapse, kb_floor) %>%
    pivot_longer(c(total, direct), names_to = "effect", values_to = "value")
  long_clean <- long %>% filter(microbe_collapse == 0, kb_floor == 0)
  flag_txt <- up$summary %>%
    transmute(scenario, lab = sprintf("flagged: %d collapse, %d k_b floor",
                                      n_microbe_collapse, n_kb_floor))
  bands <- up$summary %>%
    transmute(scenario,
              total_lo = total_q025, total_hi = total_q975, total_def = total_default,
              direct_lo = direct_q025, direct_hi = direct_q975, direct_def = direct_default) %>%
    pivot_longer(-scenario, names_to = c("effect", ".value"), names_sep = "_")
  p_up <- ggplot(long, aes(value, fill = effect, colour = effect)) +
    geom_rect(data = bands, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf, fill = effect),
              inherit.aes = FALSE, alpha = 0.12) +
    geom_density(alpha = 0.35) +
    geom_density(data = long_clean, fill = NA, linetype = "dashed", linewidth = 0.5) +
    geom_text(data = flag_txt, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE,
              hjust = 1.05, vjust = 1.5, size = 2.8, colour = "grey30") +
    geom_vline(data = bands, aes(xintercept = def, colour = effect), linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
    facet_wrap(~scenario, scales = "free", labeller = scenario_labeller) +
    scale_fill_manual(values = c(total = "#1b7837", direct = "#762a83"), name = "Effect") +
    scale_colour_manual(values = c(total = "#1b7837", direct = "#762a83"), name = "Effect") +
    labs(x = expression("Animal effect on total soil + root C (g C m"^-2*")"),
         y = "Density",
         caption = paste("Filled = all converged samples; dashed outline = excluding flagged regimes",
                         "(microbial collapse, k_b at floor). Solid line = all-default effect;",
                         "shaded band = central 95%; dotted = no effect.")) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "animal_effect_uncertainty.png"), p_up,
         width = 12, height = 8, dpi = 150)
  cat("Wrote:", file.path(res_dir, "animal_effect_uncertainty_samples.csv"),
      file.path(res_dir, "animal_effect_uncertainty_summary.csv"),
      file.path(fig_dir, "animal_effect_uncertainty.png"), "\n")
}

# shut down the worker pool
if (!is.null(cl)) parallel::stopCluster(cl)
