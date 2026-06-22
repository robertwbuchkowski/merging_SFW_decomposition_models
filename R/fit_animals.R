# ============================================================
# fit_animals.R - calibrate animal parameters for a TREATMENT scenario
# ------------------------------------------------------------
# Two goals, fit one animal at a time:
#   (1) BIOMASS: tune a feeding-rate parameter so the animal's EQUILIBRIUM
#       biomass matches a target (default = its starting/input value).
#   (2) EFFECT (optional, user-defined): tune an "effect" parameter so the
#       animal changes a chosen soil pool by a target % relative to the
#       no-animal baseline equilibrium.
#
# Both responses are monotonic (verified), so each goal is a 1-D root-find
# (stats::uniroot) with an automatically expanded bracket. The two are fit by
# coordinate descent because they interact weakly. Equilibria are obtained
# with the robust runsteady spin-up (spinup_equilibrium), warm-started across
# iterations for speed.
#
# Returns the calibrated treatment object plus a full iteration $history so
# you can see exactly what happened. Requires setup.R, spinup.R sourced.
# ============================================================

# default parameter levers per animal (override via args if your model differs)
animal_fit_defaults <- list(
  Earthworm   = list(biomass_param = "c_earthworm_soil",  effect_param = "k_b_slope_pint"),
  Detritivore = list(biomass_param = "c_detritivores",    effect_param = "slope_pint_det_k_frag_litter"),
  DetPredator = list(biomass_param = "c_detpredator",     effect_param = NA),
  RootHerb    = list(biomass_param = "c_rootherb",        effect_param = "k_exudate_slope")
)

# expand around x0 until f changes sign; returns c(lo, hi) or NULL (monotone,
# no root in range -> target not achievable, e.g. a saturating effect).
.bracket_root <- function(f, x0, step, max_expand = 40) {
  f0 <- f(x0); if (!is.finite(f0)) return(NULL)
  if (f0 == 0) return(c(x0, x0))
  lo <- x0; hi <- x0; flo <- f0; fhi <- f0
  for (k in seq_len(max_expand)) {
    lo <- lo - step; flo <- f(lo)
    if (is.finite(flo) && flo * f0 <= 0) return(c(lo, hi))
    hi <- hi + step; fhi <- f(hi)
    if (is.finite(fhi) && fhi * f0 <= 0) return(c(lo, hi))
  }
  NULL
}

fit_animal_params <- function(treatment, baseline,
                              animal            = NULL,
                              target_biomass    = NULL,   # default: the animal's input value
                              biomass_param     = NULL,
                              effect_pool       = NULL,   # set these two to ALSO fit an effect
                              target_effect_pct = NULL,
                              effect_param      = NULL,
                              tol_biomass = 0.02,   # relative
                              tol_effect  = 1.0,    # percentage points
                              max_outer   = 6,
                              max_time = 1e7, stol = 1e-8,
                              verbose = TRUE) {

  animals_all    <- c("Earthworm", "Detritivore", "DetPredator", "RootHerb")
  active_animals <- intersect(animals_all, treatment$active)
  if (is.null(animal)) {
    if (length(active_animals) != 1)
      stop("Specify `animal`; treatment has active animals: ",
           paste(active_animals, collapse = ", "))
    animal <- active_animals
  }
  if (is.null(biomass_param)) biomass_param <- animal_fit_defaults[[animal]]$biomass_param
  if (is.null(target_biomass)) target_biomass <- unname(treatment$working_state[animal])
  if (is.null(biomass_param) || is.na(biomass_param))
    stop("No biomass_param available for ", animal, " - pass biomass_param=.")

  do_effect <- !is.null(effect_pool) && !is.null(target_effect_pct)
  if (do_effect && is.null(effect_param)) effect_param <- animal_fit_defaults[[animal]]$effect_param
  if (do_effect && (is.null(effect_param) || is.na(effect_param)))
    stop("Effect fitting requested but no effect_param for ", animal, " - pass effect_param=.")

  # baseline equilibrium (animal params do not affect it) -- compute once
  if (is.null(baseline$init_state_spin))
    baseline <- spinup_equilibrium(baseline, max_time = max_time, stol = stol, verbose = FALSE)
  base_eq <- baseline$init_state_spin

  # warm-started equilibrium evaluator for the treatment
  last_eq <- treatment$working_state
  eq_now <- function() {
    sp <- spinup_equilibrium(treatment, warm_start = last_eq,
                             max_time = max_time, stol = stol, verbose = FALSE)
    last_eq <<- sp$init_state_spin
    sp$init_state_spin
  }
  effect_of <- function(eq) 100 * (eq[effect_pool] - base_eq[effect_pool]) /
                            max(abs(base_eq[effect_pool]), 1e-8)

  history <- list()
  for (outer in seq_len(max_outer)) {

    # ---- (1) biomass via feeding rate (log10 scale, monotone increasing) ----
    fb <- function(logc) {
      treatment$parms[[biomass_param]] <<- 10^logc
      unname(eq_now()[animal]) - target_biomass
    }
    x0 <- log10(treatment$parms[[biomass_param]])
    br <- .bracket_root(fb, x0, step = 0.4)
    if (is.null(br)) {
      warning("Could not bracket a feeding rate for target biomass ", target_biomass,
              " (target may be outside the achievable range).")
    } else {
      root <- stats::uniroot(fb, lower = br[1], upper = br[2], tol = 1e-4)$root
      treatment$parms[[biomass_param]] <- 10^root
    }

    # ---- (2) effect via effect parameter (linear scale, monotone) ----
    if (do_effect) {
      fe <- function(val) {
        treatment$parms[[effect_param]] <<- val
        unname(effect_of(eq_now())) - target_effect_pct
      }
      e0 <- treatment$parms[[effect_param]]
      step_e <- max(abs(e0), 1) * 0.5
      bre <- .bracket_root(fe, e0, step = step_e)
      if (is.null(bre)) {
        warning("Target effect ", target_effect_pct, "% on ", effect_pool,
                " not achievable (the parameter saturates or the target is too large).")
      } else {
        roote <- stats::uniroot(fe, lower = bre[1], upper = bre[2], tol = 1e-6)$root
        treatment$parms[[effect_param]] <- roote
      }
    }

    eq  <- eq_now()
    B   <- unname(eq[animal])
    eff <- if (do_effect) unname(effect_of(eq)) else NA_real_
    history[[outer]] <- data.frame(
      outer = outer,
      biomass_param = treatment$parms[[biomass_param]],
      biomass = B,
      effect_param = if (do_effect) treatment$parms[[effect_param]] else NA_real_,
      effect_pct = eff
    )
    if (verbose) {
      msg <- sprintf("outer %d: %s = %.4g -> biomass = %.4g (target %.4g)",
                     outer, biomass_param, treatment$parms[[biomass_param]], B, target_biomass)
      if (do_effect)
        msg <- paste0(msg, sprintf(" | %s = %.4g -> %s %+.1f%% (target %+.1f%%)",
                                   effect_param, treatment$parms[[effect_param]],
                                   effect_pool, eff, target_effect_pct))
      message(msg)
    }
    ok_b <- abs(B - target_biomass) / max(abs(target_biomass), 1e-8) < tol_biomass
    ok_e <- !do_effect || abs(eff - target_effect_pct) < tol_effect
    if (ok_b && ok_e) break
  }

  treatment$init_state_spin <- eq_now()
  treatment$fit <- list(
    animal = animal,
    biomass_param = biomass_param, fitted_biomass_param = treatment$parms[[biomass_param]],
    target_biomass = target_biomass, achieved_biomass = unname(treatment$init_state_spin[animal]),
    effect_param = if (do_effect) effect_param else NA,
    fitted_effect_param = if (do_effect) treatment$parms[[effect_param]] else NA,
    effect_pool = effect_pool, target_effect_pct = target_effect_pct,
    achieved_effect_pct = if (do_effect) unname(effect_of(treatment$init_state_spin)) else NA,
    history = do.call(rbind, history)
  )
  treatment
}
