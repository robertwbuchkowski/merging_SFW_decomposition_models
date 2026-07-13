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

# default parameter levers per animal (override via args if your model differs).
#   biomass_param  feeding-rate ADJUSTMENT FACTOR (adj_*) tuned to hit a target
#                  biomass. adj_* scales ALL of an animal's individual feeding
#                  coefficients together (see the models), so multi-food animals
#                  (e.g. earthworms, detritivores) keep their relative food
#                  preferences while overall feeding scales up/down.
#   effect_param   parameter tuned to hit a target pool effect (same name in
#                  all three models here)
#   effect_pool    the pool whose response defines the "effect" -- specified
#                  the DEFAULT pool (Millennial pools: Litter/CWD/Organic/
#                  DOM/MIC/P/L/A/M/B); override it per scenario in
#                  effect_pool_overrides below.
#   effect_pct     default target effect SIZE (percent change vs the no-animal
#                  baseline equilibrium) for this animal, used when you fit an
#                  effect but do not pass an explicit target. Sign matters:
#                  detritivores/earthworms typically DEPLETE the litter/pool
#                  they feed on (negative), root herbivores boost exudate-fed
#                  pools (positive). These are starting guesses -- confirm with
#                  a scan_animal_param() sweep and adjust to your data.
animal_fit_defaults <- list(
  Earthworm = list(
    biomass_param = "adj_earthworm",     # scales all earthworm feeding rates
    effect_param  = "k_b_slope_pint",                 # acts on desorb / binding / f_PASSIVE
    effect_pool   = list(millennial = "M"),
    effect_pct    = -10),
  Detritivore = list(
    biomass_param = "adj_detritivores",  # scales all detritivore feeding rates
    effect_param  = "slope_pint_det_k_frag_litter",   # acts on FI / f_MetLitter / fragmentation
    effect_pool   = list(millennial = "Litter"),
    effect_pct    = -15),
  DetPredator = list(
    biomass_param = "adj_detpredator",   # scales predator feeding rate
    effect_param  = NA,                               # predator: no direct pool effect
    effect_pool   = list(),
    effect_pct    = NA),
  RootHerb = list(
    biomass_param = "adj_rootherb",      # scales root-herbivore feeding rate
    effect_param  = "k_exudate_slope",                # acts on root exudation
    effect_pool   = list(millennial = "DOM"),
    effect_pct    = +10)
)

# ------------------------------------------------------------
# PER model x scenario x animal OVERRIDES of the effect target (pool, param,
# and/or size). This is the key knob you asked for: the effect pool no longer
# has to be the same for a given animal across scenarios -- set a different
# pool (and, if you like, a different target size or even a different effect
# parameter) for any specific model x scenario x animal here. Anything left
# out falls back to the per-MODEL default in animal_fit_defaults above.
#
# Keyed as effect_pool_overrides[[model]][[scenario]][[animal]] = list(
#     pool = "<pool name in THIS model>",   # which pool the effect is measured on
#     pct  = <numeric percent change>,      # optional; else animal's effect_pct
#     param = "<effect parameter>"          # optional; else animal's effect_param
#   ). Scenario names are matched case/punctuation-insensitively (normalize()),
# so "EarthwormTemperate" and "Earthworm - Temperate" both work.
#
# SUGGESTED DEFAULTS (starting points -- verify with scan_animal_param and edit):
#   Detritivores shred/consume litter, so they DEPLETE the metabolic-litter /
#   structural-litter pool they feed on (target the litter pool, negative pct).
#   Earthworms rework and stabilise organo-mineral C, so target a protected /
#   passive soil pool (SOM_1 / PASSIVE / M) -- sign depends on whether your
#   hypothesis is net protection (+) or accelerated turnover (-).
#   Root herbivores raise root exudation, feeding the fast/active DOM pool
#   (SOM_1 / ACTIVE / DOM, positive).
effect_pool_overrides <- list(
  millennial = list(
    Isopod        = list(Detritivore = list(pool = "Litter", pct = -10)),
    Mite          = list(Detritivore = list(pool = "Litter", pct = -10)),
    RootHerbivore = list(RootHerb    = list(pool = "DOM",    pct = +10)),
    Earthworm     = list(Earthworm   = list(pool = "M",      pct = +10))
  )
)

# resolve effect_pool_overrides[[model]][[scenario]][[animal]] with scenario
# names normalized (so punctuation/case differences don't matter).
.norm_scen <- function(x) tolower(gsub("[^[:alnum:]]+", "", trimws(as.character(x))))
.lookup_effect_override <- function(model, scenario, animal) {
  mo <- effect_pool_overrides[[model]]
  if (is.null(mo) || is.null(scenario)) return(NULL)
  key <- names(mo)[.norm_scen(names(mo)) == .norm_scen(scenario)]
  if (!length(key)) return(NULL)
  mo[[key[1]]][[animal]]
}

# ------------------------------------------------------------
# animal_fit_spec(): resolve biomass_param, effect_param, effect_pool, and the
# default effect size (effect_pct) for an animal. Resolution order for the
# effect pool/param/size is:
#   (1) a per model x scenario override in effect_pool_overrides, then
#   (2) the per-MODEL default in animal_fit_defaults, then
#   (3) NA (nothing to fit).
# Pass `scenario` to enable per model x scenario resolution; omit it to get the
# per-MODEL default only (backwards compatible).
# ------------------------------------------------------------
animal_fit_spec <- function(animal, model = NULL, scenario = NULL) {
  d <- animal_fit_defaults[[animal]]
  if (is.null(d)) stop("No animal_fit_defaults entry for '", animal, "'.")

  ep_default <- if (!is.null(model)) d$effect_pool[[model]] else NULL
  ov <- if (!is.null(model)) .lookup_effect_override(model, scenario, animal) else NULL

  effect_pool  <- if (!is.null(ov$pool))  ov$pool  else ep_default
  effect_param <- if (!is.null(ov$param)) ov$param else d$effect_param
  effect_pct   <- if (!is.null(ov$pct))   ov$pct   else d$effect_pct

  list(biomass_param = if (is.null(d$biomass_param)) NA else d$biomass_param,
       effect_param  = if (is.null(effect_param)) NA else effect_param,
       effect_pool   = if (is.null(effect_pool))  NA else effect_pool,
       effect_pct    = if (is.null(effect_pct))   NA else effect_pct)
}

# ------------------------------------------------------------
# fit_param_grid(): build a gradient that BUFFERS around a parameter's default
# (or any center) value, for scan_animal_param().
#   center  the value to buffer around (e.g. parms[[biomass_param]])
#   buffer  half-width: decades each side on a log scale (default), or a
#           fraction of |center| on a linear scale
#   scale   "log" for positive rate parameters; "linear" for params that can
#           be <= 0 (e.g. k_b_slope_pint)
# ------------------------------------------------------------
fit_param_grid <- function(center, buffer = 1, n = 11, scale = c("log", "linear")) {
  scale <- match.arg(scale)
  if (scale == "log") {
    if (!is.finite(center) || center <= 0)
      stop("log grid needs a positive center; use scale='linear' for <= 0 params.")
    center * 10^seq(-buffer, buffer, length.out = n)
  } else {
    half <- if (center != 0) abs(center) * buffer else buffer
    seq(center - half, center + half, length.out = n)
  }
}

# expand around x0 until f changes sign; returns c(lo, hi) or NULL (monotone,
# no root in range -> target not achievable, e.g. a saturating effect).
.bracket_root <- function(f, x0, step, max_expand = 40) {
  f0 <- f(x0); if (!is.finite(f0)) return(NULL)
  if (f0 == 0) return(c(x0, x0))
  lo <- x0; hi <- x0
  for (k in seq_len(max_expand)) {
    lo <- lo - step; flo <- f(lo)
    # tight bracket between the sign-change point and x0 (BOTH finite) -- never
    # return a far endpoint that may sit in a blown-up / non-converged region.
    if (is.finite(flo) && flo * f0 <= 0) return(c(lo, x0))
    hi <- hi + step; fhi <- f(hi)
    if (is.finite(fhi) && fhi * f0 <= 0) return(c(x0, hi))
  }
  NULL
}

fit_animal_params <- function(treatment, baseline,
                              animal            = NULL,
                              scenario          = NULL,   # enables per model x scenario effect pool
                              target_biomass    = NULL,   # default: the animal's input value
                              biomass_param     = NULL,
                              effect_pool       = NULL,   # set to ALSO fit an effect (else uses
                              target_effect_pct = NULL,   #   the per model x scenario default)
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
  if (is.null(scenario)) scenario <- treatment$scenario   # set by setup_scenario()
  spec <- animal_fit_spec(animal, treatment$model, scenario)
  if (is.null(biomass_param)) biomass_param <- spec$biomass_param
  if (is.null(target_biomass)) target_biomass <- unname(treatment$working_state[animal])
  if (is.null(biomass_param) || is.na(biomass_param))
    stop("No biomass_param available for ", animal, " - pass biomass_param=.")

  # EFFECT FITTING IS EXPLICIT: an effect is fit ONLY when the caller supplies
  # BOTH an effect_pool AND a target_effect_pct (e.g. from an effect_spec read
  # in by Scripts/fit_all_animals.R). Nothing is pulled from the defaults table
  # in the background -- so leaving them NULL gives clean biomass-only fitting.
  # (effect_param, the LEVER, still defaults from the spec once you ask for an
  # effect; pass effect_param= to override it.)
  do_effect <- !is.null(effect_pool) && !is.null(target_effect_pct) &&
               !all(is.na(effect_pool)) && !all(is.na(target_effect_pct))
  if (do_effect && is.null(effect_param)) effect_param <- spec$effect_param
  if (do_effect && (is.null(effect_param) || is.na(effect_param)))
    stop("Effect fitting requested but no effect_param for ", animal, " - pass effect_param=.")

  # baseline equilibrium (animal params do not affect it) -- compute once
  if (is.null(baseline$init_state_spin))
    baseline <- spinup_equilibrium(baseline, max_time = max_time, stol = stol, verbose = FALSE)
  base_eq <- baseline$init_state_spin

  # warm-started equilibrium evaluator for the treatment (error-safe)
  last_eq   <- treatment$working_state
  last_conv <- TRUE
  eq_now <- function() {
    sp <- spinup_equilibrium(treatment, warm_start = last_eq,
                             max_time = max_time, stol = stol, verbose = FALSE)
    if (all(is.finite(sp$init_state_spin))) last_eq <<- sp$init_state_spin
    last_conv <<- isTRUE(sp$spin_info$converged)
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
      root <- tryCatch(stats::uniroot(fb, lower = br[1], upper = br[2], tol = 1e-4)$root,
                       error = function(e) { warning("biomass uniroot failed: ",
                                                     conditionMessage(e)); NA_real_ })
      if (is.finite(root)) treatment$parms[[biomass_param]] <- 10^root
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
        roote <- tryCatch(stats::uniroot(fe, lower = bre[1], upper = bre[2], tol = 1e-6)$root,
                          error = function(e) { warning("effect uniroot failed: ",
                                                        conditionMessage(e)); NA_real_ })
        if (is.finite(roote)) treatment$parms[[effect_param]] <- roote
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
    solver_converged = last_conv,
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

# ============================================================
# scan_animal_param() - response curves along a USER-DEFINED parameter gradient
# ------------------------------------------------------------
# fit_animal_params() returns a single optimum; this returns the WHOLE curve so
# you can see how equilibrium biomass and the pool effect respond to a
# parameter, spot non-convergent / unstable regions (flagged, not fatal), and
# tune by hand. Each point uses the error-safe runsteady spin-up and is
# warm-started from the previous point for speed and continuity.
#
#   treatment    a treatment setup object (animals on)
#   param        name of the parameter to vary (e.g. "c_detritivores",
#                "slope_pint_det_k_frag_litter")
#   values       numeric vector to sweep (the gradient you define)
#   animal       which animal's biomass to record (auto if exactly one active)
#   baseline     baseline setup/equilibrium; needed only if effect_pool is set
#   effect_pool  optional pool whose % change vs the baseline equilibrium to record
#   warm_start   start each point from the previous equilibrium (default TRUE)
#
# Returns a data.frame with columns: <param>, biomass, effect_pct, converged,
# max_deriv. attr(, "param"/"animal"/"effect_pool") describe the scan.
# ============================================================
scan_animal_param <- function(treatment, param, values,
                              animal = NULL, baseline = NULL, effect_pool = NULL,
                              warm_start = TRUE,
                              max_time = 1e7, stol = 1e-8, verbose = TRUE) {

  if (!param %in% names(treatment$parms))
    warning("Parameter '", param, "' is not currently in parms; it will be added.")

  animals_all    <- c("Earthworm", "Detritivore", "DetPredator", "RootHerb")
  active_animals <- intersect(animals_all, treatment$active)
  if (is.null(animal)) {
    if (length(active_animals) != 1)
      stop("Specify `animal`; treatment has active animals: ",
           paste(active_animals, collapse = ", "))
    animal <- active_animals
  }

  base_eq <- NULL
  if (!is.null(effect_pool)) {
    if (is.null(baseline)) stop("effect_pool set but no `baseline` supplied.")
    if (is.null(baseline$init_state_spin))
      baseline <- spinup_equilibrium(baseline, max_time = max_time, stol = stol, verbose = FALSE)
    base_eq <- baseline$init_state_spin
  }

  last_eq <- treatment$working_state
  rows <- vector("list", length(values))
  for (i in seq_along(values)) {
    treatment$parms[[param]] <- values[i]
    sp <- spinup_equilibrium(treatment, warm_start = if (warm_start) last_eq else NULL,
                             max_time = max_time, stol = stol, verbose = FALSE)
    eq <- sp$init_state_spin
    if (all(is.finite(eq))) last_eq <- eq
    B   <- unname(eq[animal])
    eff <- if (!is.null(effect_pool))
             100 * (unname(eq[effect_pool]) - unname(base_eq[effect_pool])) /
               max(abs(unname(base_eq[effect_pool])), 1e-8)
           else NA_real_
    rows[[i]] <- data.frame(value = values[i], biomass = B, effect_pct = eff,
                            converged = isTRUE(sp$spin_info$converged),
                            max_deriv = sp$spin_info$max_deriv)
    if (verbose)
      cat(sprintf("  %s = %.4g -> %s = %.4g%s  [%s]\n",
                  param, values[i], animal, B,
                  if (!is.null(effect_pool)) sprintf(", %s %+.1f%%", effect_pool, eff) else "",
                  if (isTRUE(sp$spin_info$converged)) "ok" else "NOT converged"))
  }
  res <- do.call(rbind, rows)
  names(res)[1] <- param
  attr(res, "param") <- param
  attr(res, "animal") <- animal
  attr(res, "effect_pool") <- effect_pool
  res
}

# ------------------------------------------------------------
# plot_animal_scan(): quick base-R view of a scan_animal_param() result.
# Biomass (and effect, if present) vs the parameter; non-converged points are
# drawn in red so unstable regions are obvious. Optional reference lines for a
# target biomass and a fitted optimum.
# ------------------------------------------------------------
plot_animal_scan <- function(scan, target_biomass = NULL, fitted_value = NULL,
                             log_x = FALSE) {
  param <- attr(scan, "param"); animal <- attr(scan, "animal")
  ep    <- attr(scan, "effect_pool")
  x     <- scan[[param]]
  col   <- ifelse(scan$converged, "black", "red")
  has_eff <- !all(is.na(scan$effect_pct))
  lx <- if (log_x) "x" else ""

  op <- graphics::par(mfrow = c(1, if (has_eff) 2 else 1), mar = c(4.2, 4.2, 2.4, 1))
  on.exit(graphics::par(op))

  graphics::plot(x, scan$biomass, type = "b", log = lx, pch = 19, col = col,
                 xlab = param, ylab = paste(animal, "equilibrium biomass"),
                 main = "Biomass vs parameter")
  if (!is.null(target_biomass)) graphics::abline(h = target_biomass, lty = 2, col = "blue")
  if (!is.null(fitted_value))   graphics::abline(v = fitted_value,  lty = 3, col = "darkgreen")
  graphics::legend("topleft", bty = "n", pch = c(19, 19), col = c("black", "red"),
                   legend = c("converged", "not converged"), cex = 0.8)

  if (has_eff) {
    graphics::plot(x, scan$effect_pct, type = "b", log = lx, pch = 19, col = col,
                   xlab = param, ylab = paste0("% change in ", ep, " vs baseline"),
                   main = "Effect vs parameter")
    graphics::abline(h = 0, col = "grey60")
    if (!is.null(fitted_value)) graphics::abline(v = fitted_value, lty = 3, col = "darkgreen")
  }
  invisible(NULL)
}

# ============================================================
# SAVE / LOAD / APPLY fitted animal parameters
# ------------------------------------------------------------
# Fit once (Scripts/fit_all_animals.R), save to a tidy file keyed by MODEL, and
# reuse elsewhere (e.g. Scripts/spinup_dynamic.R) instead of re-fitting.
# ============================================================

# save_fitted_params(): write one row per model x scenario x fitted parameter.
# `summary_long` is the data frame built by fit_all_animals.R (needs columns
# model, scenario, animal, param, fitted, role).
save_fitted_params <- function(summary_long,
                               file = "Results/fitted_animal_params.csv") {
  need <- c("model", "scenario", "animal", "param", "fitted", "role")
  if (!all(need %in% names(summary_long)))
    stop("save_fitted_params(): summary_long needs columns ", paste(need, collapse = ", "))
  keep <- summary_long[, need]
  names(keep)[names(keep) == "fitted"] <- "value"
  keep <- keep[!is.na(keep$param) & is.finite(keep$value), , drop = FALSE]
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(keep, file, row.names = FALSE)
  message("saved fitted parameters -> ", file, " (", nrow(keep), " rows)")
  invisible(file)
}

# load_fitted_params(): read the fitted-parameter table back. If the dedicated
# file is missing, fall back to the full summary (animal_fit_summary_long.csv),
# so the spin-up still works from an earlier fit run.
load_fitted_params <- function(file = "Results/fitted_animal_params.csv",
                               fallback = "Results/animal_fit_summary_long.csv") {
  if (!file.exists(file)) {
    if (file.exists(fallback)) {
      message("'", file, "' not found; using ", fallback)
      d <- utils::read.csv(fallback, stringsAsFactors = FALSE)
      if ("fitted" %in% names(d)) names(d)[names(d) == "fitted"] <- "value"
      return(d[, intersect(c("model", "scenario", "animal", "param", "value", "role"),
                           names(d))])
    }
    stop("No fitted-parameter file found (looked for ", file, " and ", fallback,
         ")  -- run Scripts/fit_all_animals.R first.")
  }
  utils::read.csv(file, stringsAsFactors = FALSE)
}

# apply_fitted_params(): overwrite a setup object's parms with the fitted values
# for a given MODEL and scenario. Matching is BY MODEL first, then scenario, so
# each model gets its own fitted parameters. Returns the object unchanged (with
# a message) if no saved fit matches.
apply_fitted_params <- function(obj, fitted, model, scenario, verbose = TRUE) {
  if (is.null(fitted) || !nrow(fitted)) {
    if (verbose) message("  no fitted-parameter table supplied -- parms unchanged")
    return(obj)
  }
  sub <- fitted[fitted$model == model & fitted$scenario == scenario, , drop = FALSE]
  sub <- sub[!is.na(sub$param) & is.finite(sub$value), , drop = FALSE]
  if (!nrow(sub)) {
    if (verbose) message("  no saved fit for ", model, " / ", scenario, " -- parms unchanged")
    return(obj)
  }
  for (i in seq_len(nrow(sub))) obj$parms[[ sub$param[i] ]] <- sub$value[i]
  if (verbose)
    cat(sprintf("  applied %d fitted param(s) for %s / %s: %s\n",
                nrow(sub), model, scenario,
                paste(sub$param, signif(sub$value, 4), sep = "=", collapse = ", ")))
  obj
}

# ============================================================
# EFFECT TARGETS - saved effect sizes, read in explicitly for fitting
# ------------------------------------------------------------
# Effect fitting is now OPT-IN and EXPLICIT. In Scripts/fit_all_animals.R:
#
#   effect_spec <- list()                          # biomass-only fitting
#   effect_spec <- load_effect_targets()           # ALSO fit effects, using the
#                                                  #   saved sizes on disk
#
# The saved file (default config/effect_targets.csv) has one row per
# model x scenario x animal:
#     model, scenario, animal, pool, pct, param
#   pool   the pool the effect is measured on (a pool name in THIS model)
#   pct    target % change vs the no-animal baseline equilibrium
#   param  (optional, may be blank) the effect parameter to tune; blank falls
#          back to the animal's effect_param in animal_fit_defaults
# Edit that CSV to change which pool / how big an effect you fit, per
# model x scenario. Rows with a blank or NA pool/pct are skipped (no effect fit
# for that animal), so you can turn individual animals off without deleting them.
# ============================================================

# write the built-in animal_fit_defaults / effect_pool_overrides table out to a
# CSV you can then edit by hand. Run once to (re)generate the starting file.
save_effect_targets <- function(file = "config/effect_targets.csv",
                                models = names(effect_pool_overrides)) {
  rows <- list()
  for (m in models) {
    for (sc in names(effect_pool_overrides[[m]])) {
      for (a in names(effect_pool_overrides[[m]][[sc]])) {
        sp <- animal_fit_spec(a, m, sc)
        rows[[length(rows) + 1]] <- data.frame(
          model = m, scenario = sc, animal = a,
          pool  = if (is.na(sp$effect_pool)) NA_character_ else sp$effect_pool,
          pct   = if (is.na(sp$effect_pct))  NA_real_      else sp$effect_pct,
          param = if (is.na(sp$effect_param)) NA_character_ else sp$effect_param,
          stringsAsFactors = FALSE)
      }
    }
  }
  out <- do.call(rbind, rows)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(out, file, row.names = FALSE)
  message("wrote effect targets -> ", file, " (", nrow(out), " rows)")
  invisible(out)
}

# read the saved effect targets. Returns a data.frame that get_effect_target()
# understands. Errors if the file is missing (so a typo can't silently turn
# effect fitting off).
load_effect_targets <- function(file = "Data/effect_targets.csv") {
  if (!file.exists(file))
    stop("Effect-target file not found: ", file,
         "  -- run save_effect_targets() once to create it, or use ",
         "effect_spec <- list() for biomass-only fitting.")
  d <- utils::read.csv(file, stringsAsFactors = FALSE)
  need <- c("model", "scenario", "animal", "pool", "pct")
  if (!all(need %in% names(d)))
    stop("effect-target file needs columns: ", paste(need, collapse = ", "))
  message("effect targets loaded from ", file, " (", nrow(d), " rows)")
  d
}

# look up one effect target. Works with:
#   spec = list()            -> NULL (no effect fitting)  <- the biomass-only case
#   spec = <data.frame>      -> the matching row, or NULL if absent/blank
# Returns list(pool=, pct=, param=) or NULL.
get_effect_target <- function(spec, model, scenario, animal) {
  if (is.null(spec) || !length(spec)) return(NULL)          # effect_spec <- list()
  if (!is.data.frame(spec)) return(spec[[animal]])          # hand-written list, as before
  sub <- spec[.norm_scen(spec$model)    == .norm_scen(model) &
              .norm_scen(spec$scenario) == .norm_scen(scenario) &
              .norm_scen(spec$animal)   == .norm_scen(animal), , drop = FALSE]
  if (!nrow(sub)) return(NULL)
  r <- sub[1, ]
  if (is.na(r$pool) || !nzchar(as.character(r$pool)) || is.na(r$pct)) return(NULL)
  list(pool  = as.character(r$pool),
       pct   = as.numeric(r$pct),
       param = if (!is.null(r$param) && !is.na(r$param) && nzchar(as.character(r$param)))
                 as.character(r$param) else NULL)
}
