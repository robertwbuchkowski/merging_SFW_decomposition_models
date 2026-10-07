#------------------------------------------------------------------------#
# Sensitivity analysis: parameter sets, ranges and model evaluation ----
#------------------------------------------------------------------------#
# Helper functions for the sensitivity scripts. They were moved here from the scripts so
# that every custom function is defined BEFORE the script's analysis code
# runs (sourced at the top of the script), which also lets
# debugonce_functions() flag them all for step-through debugging.
#
# Shared settings (parameter list, SD rule, bounds) are in R/sensitivity_settings.R.
# These functions read settings and data that the script creates at top level
# (listed per function under 'Needs'). They are looked up when the function is
# CALLED, so source this file first and run the script's settings before use.

#------------------------------------------------------------------------#
# pretty_param() ----
#------------------------------------------------------------------------#
# Used in: called by run_morris(), run_uncertainty_propagation_fn()
# Needs (set in R/sensitivity_settings.R): param_labels
pretty_param <- function(x) ifelse(x %in% names(param_labels), param_labels[x], x)

#------------------------------------------------------------------------#
# param_bounds() ----
#------------------------------------------------------------------------#
# Used in: called by clamp_to_bounds(), param_sd_info(), range_for()
# Needs (set in R/sensitivity_settings.R): param_bound_table, bound_exempt_prefix,
#   positive_floor_frac, sign_free_params
# Physical bounds for one parameter, as c(lower, upper) with attribute "why"
# naming the reason for each side (used to label bound hits):
#   proportion  a_* / p_* efficiencies in [0, 1], except p_a, p_b, p_c
#   table       named entries in param_bound_table (e.g. root_to_organic,
#               pct_claysilt in [0, 100], pH in [1, 14])
#   positive    a positive value stays >= positive_floor_frac x value
#   sign        a negative value stays <= 0
# `default` (the parameter's model value) is needed for the sign rules.
param_bounds <- function(param, default = NA_real_) {
  lo <- -Inf; hi <- Inf; why <- c("", "")
  if (grepl("^(a_|p_)", param) && !param %in% bound_exempt_prefix) {
    lo <- 0; hi <- 1; why <- c("proportion", "proportion")
  }
  if (param %in% names(param_bound_table)) {
    lo <- param_bound_table[[param]][1]; hi <- param_bound_table[[param]][2]
    why <- c("table", "table")
  }
  if (is.finite(default) && !param %in% sign_free_params) {
    if (default > 0 && lo < positive_floor_frac * default) {
      lo <- positive_floor_frac * default; why[1] <- "positive"
    }
    if (default < 0 && hi > 0) { hi <- 0; why[2] <- "sign" }
  }
  structure(c(lo, hi), why = why)
}

#------------------------------------------------------------------------#
# clamp_to_bounds() ----
#------------------------------------------------------------------------#
# Used in: called by range_for()
clamp_to_bounds <- function(v, param, default = NA_real_) {
  b <- param_bounds(param, default); pmin(pmax(v, b[1]), b[2])
}

#------------------------------------------------------------------------#
# set_param() ----
#------------------------------------------------------------------------#
# Used in: called by effect_scalar(), headline_effects()
# Needs (set in Scripts/4 and Scripts/5): derive_fn
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
# Needs (set in R/sensitivity_settings.R): collapse_rel, kb_floor
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
# Needs (set in R/sensitivity_settings.R): animal_pools, eq_max_time
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
# Needs (set in Scripts/4 and Scripts/5): unc_all
unc_row <- function(scenario, param) {
  if (is.null(unc_all)) return(NULL)
  r <- unc_all[unc_all$scenario == scenario & unc_all$parameter == param, , drop = FALSE]
  if (nrow(r)) r[1, , drop = FALSE] else NULL
}

#------------------------------------------------------------------------#
# param_sd_info() ----
#------------------------------------------------------------------------#
# Used in: called by range_for() (Morris, Scripts/4) and up_distribution()
#   (uncertainty propagation, Scripts/5)
# Needs (set in R/sensitivity_settings.R): cv_default
# The standard deviation and analysed range for one scenario x parameter,
# centred on the model value `default`:
#   sd     reported SD ("SD"); else cv_default x |default| ("CV-default"), reduced to
#          (Max - Min) / 4 when a reported Min/Max implies less ("CV-MinMax")
#   lo/hi  default +/- 2 sd, cut to the reported Min/Max, then to the physical
#          bounds of param_bounds()
#   bound_hit  which side(s) were cut and why, e.g. "lower:positive",
#          "upper:Max", "lower:proportion; upper:proportion" ("" = none)
# use_reported_sd = FALSE ignores reported SDs and use_minmax = FALSE ignores the
# reported Min/Max entirely (no SD reduction, no cut), so every parameter gets
# SD = cv_default x |default| cut only at the physical bounds -- the
# standardized Morris analysis.
param_sd_info <- function(scenario, param, default, use_reported_sd = TRUE,
                          use_minmax = TRUE) {
  r <- unc_row(scenario, param)
  has_mm <- use_minmax && !is.null(r) && is.finite(r$min) && is.finite(r$max) && r$max > r$min
  if (use_reported_sd && !is.null(r) && is.finite(r$sd_reported)) {
    sd <- r$sd_reported; src <- "SD"
  } else {
    sd <- cv_default * abs(default); src <- "CV-default"
    if (has_mm && (r$max - r$min) / 4 < sd) { sd <- (r$max - r$min) / 4; src <- "CV-MinMax" }
  }
  lo <- default - 2 * sd; hi <- default + 2 * sd; hit <- character(0)
  if (has_mm) {
    if (lo < r$min) { lo <- r$min; hit <- c(hit, "lower:Min") }
    if (hi > r$max) { hi <- r$max; hit <- c(hit, "upper:Max") }
  }
  b <- param_bounds(param, default); why <- attr(b, "why")
  if (lo < b[1]) { lo <- b[1]; hit <- c(hit, paste0("lower:", why[1])) }
  if (hi > b[2]) { hi <- b[2]; hit <- c(hit, paste0("upper:", why[2])) }
  if (lo > hi) {                      # model value outside the reported Min/Max
    lo <- hi <- default; hit <- c(hit, "default outside Min/Max")
  }
  list(sd = sd, cv = if (default != 0) sd / abs(default) else NA_real_,
       source = src, lo = lo, hi = hi, centre = default,
       bound_hit = paste(hit, collapse = "; "))
}

#------------------------------------------------------------------------#
# range_for() ----
#------------------------------------------------------------------------#
# Used in: called by run_morris()
# range_for(): lo, hi, the range-source label and the bound hits for one
# parameter, under either Morris mode (both via param_sd_info()):
#   mode = "main"         reported SD where available, else the CV rule:
#                         default +/- 2 SD.        labels: SD / CV-MinMax / CV-default
#   mode = "standardized" SD = cv_default x |default| for EVERY parameter,
#                         ignoring reported SDs and Min/Max entirely; the range
#                         (+/- 2 SD) is cut only at the physical bounds.
#                                                  label: CV-default
range_for <- function(scenario, param, default, mode = c("main", "standardized")) {
  mode <- match.arg(mode)
  s <- param_sd_info(scenario, param, default, use_reported_sd = (mode == "main"),
                     use_minmax = (mode == "main"))
  list(lo = s$lo, hi = s$hi, source = s$source, sd = s$sd, cv = s$cv,
       bound_hit = s$bound_hit)
}

#------------------------------------------------------------------------#
# param_set_for() ----
#------------------------------------------------------------------------#
# Used in: called by run_morris(), run_uncertainty_propagation_fn()
# Needs (set in Scripts/4 and Scripts/5): animal_unc, scenspec_unc
# Needs (set in R/sensitivity_settings.R): sweep_params
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
# Needs (set in Scripts/4 and Scripts/5): fitted_params, model, scen
# base scenario pair with the fitted animal parameters applied
make_base_pair <- function(scenario) {
  bp <- setup_scenario_pair(model, scen, scenario)
  if (!is.null(fitted_params))
    bp$treatment <- apply_fitted_params(bp$treatment, fitted_params,
                                        model, scenario, verbose = FALSE)
  bp
}
