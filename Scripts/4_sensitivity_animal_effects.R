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
# Paired Morris ranking: current knowledge vs standardized ----
#------------------------------------------------------------------------#
# Per scenario, the top parameters with bootstrap 95% intervals for BOTH range
# definitions. mu* is scaled to the most influential parameter within each
# scenario x analysis (raw mu* depends on range width, so it is not comparable
# between analyses). "*" after a label = the bootstrap rank intervals of the two
# analyses do not overlap (rank shift beyond sampling noise). Uses this
# session's runs, or the saved CSVs when the Morris runs were switched off.
plot_morris_paired <- function(main_tbl, std_tbl, n_top = 8) {
  prep <- function(d, tag) d %>%
    transmute(scenario, parameter, parameter_label, analysis = tag,
              mu_star, mu_star_lo, mu_star_hi, rank, rank_lo, rank_hi)
  mm <- bind_rows(prep(main_tbl, "main"), prep(std_tbl, "std")) %>%
    group_by(scenario, analysis) %>%
    mutate(s_max = max(mu_star, na.rm = TRUE), rel = mu_star / s_max,
           rel_lo = mu_star_lo / s_max, rel_hi = mu_star_hi / s_max) %>%
    ungroup()
  keep <- mm %>% filter(rank <= n_top) %>% distinct(scenario, parameter)
  pd   <- mm %>% semi_join(keep, by = c("scenario", "parameter"))
  shift <- pd %>% select(scenario, parameter, analysis, rank_lo, rank_hi) %>%
    pivot_wider(names_from = analysis, values_from = c(rank_lo, rank_hi)) %>%
    mutate(star = rank_hi_main < rank_lo_std | rank_hi_std < rank_lo_main) %>%
    select(scenario, parameter, star)
  ord <- pd %>% filter(analysis == "main") %>% select(scenario, parameter, rel_main = rel)
  pd <- pd %>% left_join(ord, by = c("scenario", "parameter")) %>%
    left_join(shift, by = c("scenario", "parameter")) %>%
    mutate(label = paste0(parameter_label, ifelse(!is.na(star) & star, " *", "")),
           key   = paste(label, scenario, sep = "___"))
  pd$key <- factor(pd$key, levels = unique(pd$key[order(pd$scenario, pd$rel_main)]))

  dodge <- c(main = 0.18, std = -0.18)       # current knowledge above, standardized below
  ggplot(pd, aes(y = key, colour = analysis)) +
    lapply(names(dodge), function(a) list(
      geom_errorbar(data = filter(pd, analysis == a),
                    aes(xmin = rel_lo, xmax = rel_hi), width = 0.25,
                    orientation = "y", position = position_nudge(y = dodge[[a]])),
      geom_point(data = filter(pd, analysis == a), aes(x = rel), size = 2,
                 position = position_nudge(y = dodge[[a]])))) +
    facet_wrap(~scenario, scales = "free_y", labeller = scenario_labeller) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_colour_manual(values = c(main = "#2166ac", std = "#b2182b"),
                        labels = c(main = "Current knowledge", std = "Standardized (50-200%)"),
                        name = NULL) +
    labs(x = expression(mu*"* relative to the most influential parameter (95% bootstrap CI)"),
         y = NULL,
         caption = paste0("Top ", n_top, " parameters from either analysis per scenario, ordered by ",
                          "current-knowledge mu*. * = rank shift beyond bootstrap noise.")) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

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

# inverse-CDF samplers (u in (0,1))
q_trunc_norm <- function(u, mean, sd, lo, hi) {
  a <- pnorm(lo, mean, sd); b <- pnorm(hi, mean, sd)
  qnorm(a + u * (b - a), mean, sd)
}
q_trunc_lnorm <- function(u, centre, gsd, lo, hi) {        # centre, lo, hi > 0
  ml <- log(centre); sl <- log(gsd)
  a <- pnorm(log(lo), ml, sl); b <- pnorm(log(hi), ml, sl)
  exp(qnorm(a + u * (b - a), ml, sl))
}
q_triangular <- function(u, lo, mode, hi) {
  mode <- min(max(mode, lo), hi); fc <- (mode - lo) / (hi - lo)
  ifelse(u < fc, lo + sqrt(u * (hi - lo) * (mode - lo)),
                 hi - sqrt((1 - u) * (hi - lo) * (hi - mode)))
}

# up_distribution(): sampling distribution for one scenario x parameter
up_distribution <- function(scenario, p, d0) {
  rg <- range_for(scenario, p, d0, "main"); r <- unc_row(scenario, p)
  if (rg$source == "SD" && !is.null(r) && is.finite(r$sd) && r$sd > 0)
    return(list(dist = "truncated normal", lo = rg$lo, hi = rg$hi,
                centre = r$value, spread = r$sd))
  if (rg$source == "MinMax")
    return(list(dist = if (up_minmax_dist == "uniform") "uniform" else "triangular",
                lo = rg$lo, hi = rg$hi, centre = min(max(d0, rg$lo), rg$hi), spread = NA_real_))
  if (rg$lo * rg$hi > 0)                                    # same sign: log-normal on magnitude
    return(list(dist = "truncated log-normal", lo = rg$lo, hi = rg$hi,
                centre = d0, spread = up_gsd))
  list(dist = "uniform", lo = rg$lo, hi = rg$hi, centre = d0, spread = NA_real_)
}

up_quantile <- function(u, ds) {
  switch(ds$dist,
    "truncated normal"     = q_trunc_norm(u, ds$centre, ds$spread, ds$lo, ds$hi),
    "triangular"           = q_triangular(u, ds$lo, ds$centre, ds$hi),
    "truncated log-normal" = {
      s  <- sign(ds$centre); mg <- sort(abs(c(ds$lo, ds$hi)))
      s * q_trunc_lnorm(u, abs(ds$centre), ds$spread, mg[1], mg[2])
    },
    ds$lo + u * (ds$hi - ds$lo))                            # uniform
}

# Latin hypercube in (0,1)^k: one stratified draw per row in each column
lhs_unit <- function(n, k) {
  vapply(seq_len(k), function(j) (sample.int(n) - runif(n)) / n, numeric(n))
}

# baseline values used by the realism filter ("TotalC" = all non-animal pools)
baseline_values <- function(eb, pools) {
  soil <- setdiff(names(eb), animal_pools)
  v <- c(TotalC = sum(eb[soil]), eb[soil])
  setNames(vapply(pools, function(p) if (p %in% names(v)) unname(v[[p]]) else NA_real_,
                  numeric(1)), paste0("base_", pools))
}

# total and direct equilibrium effects + baseline values for one parameter vector
headline_effects <- function(base_pair, pv, pools) {
  pair <- base_pair
  for (nm in names(pv)) {
    pair$treatment <- set_param(pair$treatment, nm, pv[[nm]])
    pair$baseline  <- set_param(pair$baseline,  nm, pv[[nm]])
  }
  out <- c(total = NA_real_, direct = NA_real_, converged = 0,
           microbe_collapse = 0, kb_floor = 0, infeasible = 0,
           setNames(rep(NA_real_, length(pools)), paste0("base_", pools)))
  if (nzchar(infeasible_soil(pair$baseline$parms))) { out[["infeasible"]] <- 1; return(out) }
  tryCatch({
    b  <- spinup_equilibrium(pair$baseline, max_time = eq_max_time, verbose = FALSE)
    eb <- b$init_state_spin
    tr <- spinup_equilibrium(pair$treatment, warm_start = eb, max_time = eq_max_time,
                             verbose = FALSE)
    td <- spinup_equilibrium(zero_indirect_effects(pair$treatment, verbose = FALSE),
                             warm_start = eb, max_time = eq_max_time, verbose = FALSE)
    soil <- setdiff(intersect(names(eb), names(tr$init_state_spin)), animal_pools)
    conv <- isTRUE(b$spin_info$converged) && isTRUE(tr$spin_info$converged) &&
            isTRUE(td$spin_info$converged)
    rg_t <- regime_flags(eb, tr$init_state_spin, pair$treatment$parms)
    rg_d <- regime_flags(eb, td$init_state_spin, pair$treatment$parms)
    c(total  = sum(tr$init_state_spin[soil]) - sum(eb[soil]),
      direct = sum(td$init_state_spin[soil]) - sum(eb[soil]),
      converged = as.numeric(conv),
      microbe_collapse = as.numeric(grepl("microbe_collapse", paste(rg_t, rg_d))),
      kb_floor         = as.numeric(grepl("kb_floor", rg_t)),
      infeasible       = 0,
      baseline_values(eb, pools))
  }, error = function(e) out)
}

# factory so a parallel worker receives only what the evaluator needs
make_headline_fun <- function(base_pair, pools) {
  force(base_pair); force(pools)
  function(pv) headline_effects(base_pair, pv, pools)
}

# realism targets for one scenario: defaults from the all-default baseline,
# overridden / extended by rows in up_targets_file
targets_for <- function(scenario, base_pair, user_targets) {
  eb  <- spinup_equilibrium(base_pair$baseline, max_time = eq_max_time,
                            verbose = FALSE)$init_state_spin
  ref <- baseline_values(eb, up_filter_pools)
  fac <- ifelse(up_filter_pools %in% names(up_filter_factor_by_pool),
                up_filter_factor_by_pool[up_filter_pools], up_filter_factor)
  tg  <- data.frame(scenario = scenario, pool = up_filter_pools,
                    reference = unname(ref), lower = unname(ref) / fac,
                    upper = unname(ref) * fac, factor = unname(fac), source = "default",
                    stringsAsFactors = FALSE)
  tg <- tg[is.finite(tg$reference), , drop = FALSE]
  if (!is.null(user_targets)) {
    u <- user_targets[user_targets$scenario == scenario, , drop = FALSE]
    for (i in seq_len(nrow(u))) {
      j <- which(tg$pool == u$pool[i])
      if (length(j)) { tg$lower[j] <- u$lower[i]; tg$upper[j] <- u$upper[i]; tg$source[j] <- "user" }
      else tg <- rbind(tg, data.frame(scenario = scenario, pool = u$pool[i],
                                      reference = if ("reference" %in% names(u)) u$reference[i] else NA_real_,
                                      lower = u$lower[i], upper = u$upper[i], factor = NA_real_,
                                      source = "user"))
    }
  }
  tg
}

run_uncertainty_propagation_fn <- function() {
  user_targets <- if (file.exists(up_targets_file)) {
    u <- read.csv(up_targets_file, stringsAsFactors = FALSE)
    if (!all(c("scenario", "pool", "lower", "upper") %in% names(u)))
      stop(up_targets_file, " needs columns scenario, pool, lower, upper")
    u$scenario <- rename_scenarios(u$scenario)
    cat("Baseline realism targets: using", up_targets_file, "\n"); u
  } else NULL

  samp_rows <- summ_rows <- diag_rows <- input_rows <- target_rows <- list()
  q <- function(v, p) if (length(v)) unname(stats::quantile(v, p, na.rm = TRUE)) else NA_real_

  for (si in seq_along(scenarios)) {
    scenario   <- scenarios[si]
    base_pair  <- make_base_pair(scenario)
    parms_here <- base_pair$treatment$parms
    ps         <- param_set_for(scenario, parms_here)

    # sampling distributions
    dists <- list()
    for (p in ps$pnames) {
      d0 <- parms_here[[p]]
      if (!is.finite(d0) || d0 == 0) next
      ds <- up_distribution(scenario, p, d0)
      if (is.finite(ds$lo) && is.finite(ds$hi) && ds$hi > ds$lo) dists[[p]] <- ds
    }
    pn <- names(dists); k <- length(pn)
    input_rows[[scenario]] <- data.frame(
      scenario = scenario, parameter = pn, parameter_label = pretty_param(pn),
      param_type = ps$ptype[pn], default = vapply(pn, function(p) parms_here[[p]], 0),
      distribution = vapply(dists, `[[`, "", "dist"),
      lo = vapply(dists, `[[`, 0, "lo"), hi = vapply(dists, `[[`, 0, "hi"),
      centre = vapply(dists, `[[`, 0, "centre"), spread = vapply(dists, `[[`, 0, "spread"),
      row.names = NULL)

    # realism targets
    tg <- targets_for(scenario, base_pair, user_targets)
    target_rows[[scenario]] <- tg
    pools <- tg$pool

    cat("Uncertainty propagation: ", scenario, " - ", k, " parameters, ",
        up_n, " LHS samples\n", sep = "")
    set.seed(up_seed + si)
    U <- lhs_unit(up_n, k)
    X <- vapply(seq_len(k), function(j) up_quantile(U[, j], dists[[pn[j]]]), numeric(up_n))
    if (up_n == 1) X <- matrix(X, nrow = 1)
    colnames(X) <- pn

    default <- headline_effects(base_pair, setNames(numeric(0), character(0)), pools)
    rows_list <- lapply(seq_len(up_n), function(i) setNames(X[i, ], pn))
    hf <- make_headline_fun(base_pair, pools)
    Y  <- if (is.null(cl)) lapply(rows_list, hf) else parallel::parLapply(cl, rows_list, hf)
    Y  <- do.call(rbind, Y)

    # classification
    inf  <- Y[, "infeasible"] == 1
    ok   <- Y[, "converged"] == 1 & is.finite(Y[, "total"]) & is.finite(Y[, "direct"])
    flag <- Y[, "microbe_collapse"] == 1 | Y[, "kb_floor"] == 1
    pass_pool <- vapply(seq_len(nrow(tg)), function(j) {
      v <- Y[, paste0("base_", tg$pool[j])]
      is.finite(v) & v >= tg$lower[j] & v <= tg$upper[j]
    }, logical(up_n))
    if (up_n == 1) pass_pool <- matrix(pass_pool, nrow = 1)
    accepted <- ok & apply(pass_pool, 1, all)

    n_inf <- sum(inf); n_bad <- sum(!ok) - n_inf
    cat(sprintf("  %d/%d realistic samples converged; %d unrealistic soils removed\n",
                sum(ok), up_n - n_inf, n_inf))
    cat(sprintf("  baseline realism filter: %d of %d converged samples accepted (%.0f%%)\n",
                sum(accepted), sum(ok), 100 * sum(accepted) / max(sum(ok), 1)))
    if (sum(ok) > 0 && sum(accepted) / sum(ok) < 0.1)
      warning(sprintf("Uncertainty propagation %s: only %d of %d samples pass the baseline realism filter -- the parameter distributions may be wider than plausible, or the targets too strict.",
                      scenario, sum(accepted), sum(ok)))
    if (n_bad > 0)
      warning(sprintf("Uncertainty propagation %s: %d samples did NOT reach a stable state; excluded.",
                      scenario, n_bad))

    samp_rows[[scenario]] <- data.frame(scenario = scenario, sample = seq_len(up_n),
                                        accepted = accepted, Y, X, check.names = FALSE)

    rej <- setNames(colSums(!pass_pool[ok, , drop = FALSE]), paste0("n_reject_", tg$pool))
    diag_rows[[scenario]] <- data.frame(
      scenario = scenario, n_samples = up_n, n_infeasible = n_inf,
      n_nonconverged = n_bad, n_converged = sum(ok),
      n_microbe_collapse = sum(ok & Y[, "microbe_collapse"] == 1),
      n_kb_floor = sum(ok & Y[, "kb_floor"] == 1),
      n_accepted = sum(accepted), acceptance_rate = sum(accepted) / max(sum(ok), 1),
      t(rej), check.names = FALSE)

    summ <- function(idx, label) {
      tot <- Y[idx, "total"]; dir <- Y[idx, "direct"]
      pct <- 100 * dir / tot; pct <- pct[is.finite(pct)]
      data.frame(
        scenario = scenario, subset = label, n = sum(idx),
        total_default = default[["total"]],
        total_median = q(tot, 0.5), total_mean = if (length(tot)) mean(tot) else NA_real_,
        total_q025 = q(tot, 0.025), total_q975 = q(tot, 0.975),
        p_total_positive = if (length(tot)) mean(tot > 0) else NA_real_,
        p_sign_as_default = if (length(tot)) mean(sign(tot) == sign(default[["total"]])) else NA_real_,
        direct_default = default[["direct"]],
        direct_median = q(dir, 0.5), direct_q025 = q(dir, 0.025), direct_q975 = q(dir, 0.975),
        pct_direct_default = 100 * default[["direct"]] / default[["total"]],
        pct_direct_median = q(pct, 0.5), pct_direct_q025 = q(pct, 0.025),
        pct_direct_q975 = q(pct, 0.975))
    }
    summ_rows[[scenario]] <- rbind(summ(ok, "all converged"),
                                   summ(ok & !flag, "unflagged"),
                                   summ(accepted, "accepted"))
  }
  list(samples = bind_rows(samp_rows), summary = bind_rows(summ_rows),
       diagnostics = bind_rows(diag_rows), inputs = bind_rows(input_rows),
       targets = bind_rows(target_rows))
}

if (run_uncertainty_propagation) {
  if (!is.null(cl)) parallel::clusterExport(cl, c("headline_effects", "baseline_values"))
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
  fx <- function(d) mutate(d, effect = factor(effect, levels = c("total", "direct")))
  long <- fx(long); long_acc <- fx(long_acc); bands <- fx(bands); acc_txt <- fx(acc_txt)
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
