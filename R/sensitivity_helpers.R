#------------------------------------------------------------------------#
# Sensitivity analysis: parameter sets, ranges and model evaluation ----
#------------------------------------------------------------------------#
# Helper functions for Scripts/4_sensitivity_animal_effects.R. They were moved here from the script so
# that every custom function is defined BEFORE the script's analysis code
# runs (sourced at the top of the script), which also lets
# debugonce_functions() flag them all for step-through debugging.
#
# These functions read settings and data that the script creates at top level
# (listed per function under 'Needs'). They are looked up when the function is
# CALLED, so source this file first and run the script's settings before use.

#------------------------------------------------------------------------#
# pretty_param() ----
#------------------------------------------------------------------------#
# Used in: called by run_morris(), run_uncertainty_propagation_fn()
# Needs (set in 4_sensitivity_animal_effects.R): param_labels
pretty_param <- function(x) ifelse(x %in% names(param_labels), param_labels[x], x)

#------------------------------------------------------------------------#
# param_bounds() ----
#------------------------------------------------------------------------#
# Used in: called by clamp_to_bounds()
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

#------------------------------------------------------------------------#
# clamp_to_bounds() ----
#------------------------------------------------------------------------#
# Used in: called by range_for()
clamp_to_bounds <- function(v, param) { b <- param_bounds(param); pmin(pmax(v, b[1]), b[2]) }

#------------------------------------------------------------------------#
# set_param() ----
#------------------------------------------------------------------------#
# Used in: called by effect_scalar(), headline_effects()
# Needs (set in 4_sensitivity_animal_effects.R): derive_fn
# set_param(): change ONE parameter and rebuild derived params + forcing.
set_param <- function(obj, param, value) {
  obj$parms[[param]]        <- value
  obj$parms                 <- derive_fn(obj$parms)
  obj$parms$climate_forcing <- make_climate_forcing(obj$parms)
  obj
}

#------------------------------------------------------------------------#
# infeasible_soil() ----
#------------------------------------------------------------------------#
# Used in: called by effect_scalar(), headline_effects()
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

#------------------------------------------------------------------------#
# regime_flags() ----
#------------------------------------------------------------------------#
# Used in: called by animal_effect_totalC(), headline_effects()
# Needs (set in 4_sensitivity_animal_effects.R): collapse_rel, kb_floor
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

#------------------------------------------------------------------------#
# animal_effect_totalC() ----
#------------------------------------------------------------------------#
# Used in: called by effect_scalar()
# Needs (set in 4_sensitivity_animal_effects.R): animal_pools, eq_max_time
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

#------------------------------------------------------------------------#
# make_effect_fun() ----
#------------------------------------------------------------------------#
# Used in: called by run_morris()
# make_effect_fun(): the Morris evaluator for one scenario. Built by a
# top-level factory so its environment carries only base_pair to the workers.
make_effect_fun <- function(base_pair) {
  force(base_pair)            # evaluate now so workers receive the object itself
  function(pv) effect_scalar(base_pair, pv)
}

#------------------------------------------------------------------------#
# effect_scalar() ----
#------------------------------------------------------------------------#
# Used in: called by make_effect_fun()
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
# unc_row() ----
#------------------------------------------------------------------------#
# Used in: called by range_for(), up_distribution()
# Needs (set in 4_sensitivity_animal_effects.R): unc_all
unc_row <- function(scenario, param) {
  if (is.null(unc_all)) return(NULL)
  r <- unc_all[unc_all$scenario == scenario & unc_all$parameter == param, , drop = FALSE]
  if (nrow(r)) r[1, , drop = FALSE] else NULL
}

#------------------------------------------------------------------------#
# range_for() ----
#------------------------------------------------------------------------#
# Used in: called by run_morris(), up_distribution()
# Needs (set in 4_sensitivity_animal_effects.R): buffer, supp_frac
# range_for(): returns lo, hi (clamped to the physical guard rails) and a label
# for how the range was set, under either mode:
#   mode = "main"         reported +/-2 SD or [Min, Max] if available, else the
#                         half-to-double buffer.            labels: SD / MinMax / buffer
#   mode = "standardized" 50-200% for every parameter, cut to the reported
#                         Min/Max wherever that is narrower.
#                         labels: 50-200% / MinMax-bound / guard-rail-bound /
#                                 MinMax-excludes-default (Min/Max does not
#                                 overlap 50-200%; the unbound 50-200% is kept)
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

#------------------------------------------------------------------------#
# param_set_for() ----
#------------------------------------------------------------------------#
# Used in: called by run_morris(), run_uncertainty_propagation_fn()
# Needs (set in 4_sensitivity_animal_effects.R): animal_unc, scenspec_unc, sweep_params
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

#------------------------------------------------------------------------#
# make_base_pair() ----
#------------------------------------------------------------------------#
# Used in: called by run_morris(), run_uncertainty_propagation_fn()
# Needs (set in 4_sensitivity_animal_effects.R): fitted_params, model, scen
# base scenario pair with the fitted animal parameters applied
make_base_pair <- function(scenario) {
  bp <- setup_scenario_pair(model, scen, scenario)
  if (!is.null(fitted_params))
    bp$treatment <- apply_fitted_params(bp$treatment, fitted_params,
                                        model, scenario, verbose = FALSE)
  bp
}
