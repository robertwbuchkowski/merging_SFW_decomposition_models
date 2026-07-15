# ============================================================
# SPIN-UP & STABILITY HELPERS
# ------------------------------------------------------------
#   spinup_equilibrium()  FAST warm-start to steady state under CONSTANT forcing
#                         (rootSolve::runsteady forward integration). Now
#                         ERROR-SAFE: if the solver blows up (e.g. an extreme
#                         feeding rate during fitting overgrazes the microbes),
#                         it returns NA + converged=FALSE instead of throwing.
#   spinup_until_stable() SEASONAL integration to the annual limit cycle, using
#                         a phase-robust stability test (annual means) that
#                         ignores negligible pools. Writes PNG trajectory plots
#                         sized for many panels.
# ============================================================

# ------------------------------------------------------------
# runsteady_spinup(): integrate y forward to steady state under the forcing in
# `parms`. Error-safe: a solver failure or non-finite result -> converged=FALSE
# and an NA state, so callers (e.g. the fitter's bracket search) can skip it.
# ------------------------------------------------------------
runsteady_spinup <- function(y, model_fn, parms,
                             max_time = 1e7, stol = 1e-8, verbose = TRUE) {
  ss <- tryCatch(
    rootSolve::runsteady(y = y, time = c(0, max_time),
                         func = model_fn, parms = parms, stol = stol),
    error = function(e) { if (verbose) message("  runsteady error: ", conditionMessage(e)); NULL }
  )

  if (is.null(ss) || any(!is.finite(ss$y))) {
    return(list(final_state = setNames(rep(NA_real_, length(y)), names(y)),
                converged = FALSE, time_to_steady = NA_real_, max_deriv = Inf, raw = ss))
  }

  y_eq <- ss$y; names(y_eq) <- names(y)
  converged <- isTRUE(attr(ss, "steady"))
  t_steady  <- attr(ss, "time")
  d_end   <- tryCatch(model_fn(t_steady, y_eq, parms)[[1]], error = function(e) Inf)
  max_drv <- max(abs(d_end))

  if (verbose)
    cat(sprintf("  runsteady: %s | t = %.3g d | max|dy/dt| = %.2e\n",
                if (converged) "CONVERGED" else "NOT converged (raise max_time / relax stol)",
                t_steady, max_drv))

  list(final_state = y_eq, converged = converged,
       time_to_steady = t_steady, max_deriv = max_drv, raw = ss)
}

# ------------------------------------------------------------
# spinup_equilibrium(): runsteady on a setup object. Stores obj$init_state_spin
# and obj$spin_info. warm_start: shared pools override the object's init.
# ------------------------------------------------------------
spinup_equilibrium <- function(obj, warm_start = NULL,
                               max_time = 1e7, stol = 1e-8, verbose = TRUE) {
  y0 <- obj$working_state
  if (!is.null(warm_start)) {
    shared <- intersect(names(y0), names(warm_start))
    ok <- shared[is.finite(warm_start[shared])]          # never warm-start from NA
    y0[ok] <- warm_start[ok]
    if (verbose && length(ok))
      cat(sprintf("  warm-started %d shared pool(s) from a previous equilibrium\n", length(ok)))
  }

  obj$parms$climate_forcing <- make_climate_forcing_equilibrium(obj$parms)
  res <- runsteady_spinup(y0, obj$wrapped_model, obj$parms,
                          max_time = max_time, stol = stol, verbose = verbose)
  obj$parms$climate_forcing <- make_climate_forcing(obj$parms)   # restore seasonal

  obj$init_state_spin <- res$final_state
  obj$spin_info <- list(method = "runsteady", converged = res$converged,
                        time_to_steady = res$time_to_steady, max_deriv = res$max_deriv)
  obj
}

# ------------------------------------------------------------
# save_trajectory_png(): write ALL pools of a deSolve run to a single PNG,
# sized for the number of panels (default 4 columns). Replaces the old PDF.
# ------------------------------------------------------------
save_trajectory_png <- function(out, file, ncol = 4, panel_px = 320, res = 110) {
  df    <- as.data.frame(out)
  pools <- setdiff(names(df), "time")
  np    <- length(pools)
  ncol  <- min(ncol, np)
  nrow  <- ceiling(np / ncol)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  grDevices::png(file, width = ncol * panel_px, height = nrow * panel_px, res = res)
  on.exit(grDevices::dev.off())
  op <- graphics::par(mfrow = c(nrow, ncol), mar = c(3.2, 3.6, 2, 1), mgp = c(2, 0.6, 0))
  on.exit(graphics::par(op), add = TRUE)
  yr <- df$time / 365
  for (p in pools)
    graphics::plot(yr, df[[p]], type = "l", xlab = "year", ylab = p, main = p)
  invisible(file)
}

# ------------------------------------------------------------
# check_stability(): is the seasonal limit cycle stationary? Compares the
# ANNUAL MEAN of each pool over the last full period vs the previous one
# (means are phase-independent, so an uneven `by` no longer creates spurious
# drift). Pools whose mean is below `abs_floor` are flagged negligible and
# excluded from the convergence decision (so a near-zero pool can't dominate).
# Returns a table sorted by relative drift.
# ------------------------------------------------------------
# ------------------------------------------------------------
# check_by(): the output step MUST divide the year, or the seasonal cycle is
# sampled at a drifting phase and the annual-mean stability test never settles.
# 365 = 5 x 73, so the only clean steps are 1, 5, 73 and 365. (check_stability()
# now interpolates to a daily grid, so a bad `by` no longer breaks the test --
# but it still costs seasonal resolution, hence the warning.)
# ------------------------------------------------------------
valid_by <- c(1, 5, 73, 365)
check_by <- function(by) {
  if (!isTRUE(365 %% by == 0))
    warning("by = ", by, " does not divide 365 (365 = 5 x 73). The seasonal cycle ",
            "will be sampled at a drifting phase and sharp pools (e.g. C_leaf_herb, ",
            "whose litterfall pulse is ~30 d wide) will be poorly resolved. ",
            "Use one of: ", paste(valid_by, collapse = ", "), ".", call. = FALSE)
  invisible(by)
}

check_stability <- function(out, period = 365, abs_floor = 1e-3, dt = 1, use_old = FALSE) {
  df   <- as.data.frame(out)
  tt   <- df[[1]]
  cols <- setdiff(names(df), c("time", "mass_balance_check"))
  tmax <- max(tt)

  if(!use_old){
    abs_drift <- (df[tmax-period,] - df[(tmax-2*period),])[cols]
    rel_drift <- (df[tmax-period,] - df[(tmax-2*period),])[cols] / pmax(df[tmax-period,][cols], abs_floor)
    res <- data.frame(pool = cols, abs_drift = as.numeric(abs_drift), rel_drift = as.numeric(rel_drift),
                      above_floor = abs(as.numeric(abs_drift)) >= abs_floor)
  }else{
    # ------------------------------------------------------------
    # GRID-INDEPENDENT annual means. The raw output grid (seq(0, 365*n_years,
    # by = by)) only lines up with the year when `by` divides 365 (i.e. by is
    # 1, 5, 73 or 365). With, say, by = 30 (365/30 = 12.167) consecutive 365-day
    # windows sample DIFFERENT phases of the seasonal cycle, so the annual mean
    # wobbles forever and a perfectly converged limit cycle is reported as
    # drifting -- worst for pools with a big, sharp seasonal swing (C_leaf_herb).
    # Fix: interpolate every pool onto a uniform `dt`-day grid FIRST, then take
    # the trapezoidal means. The metric is then correct for any `by`.
    # ------------------------------------------------------------
    mean_over <- function(col, t0, t1) {
      t0 <- max(t0, min(tt))
      if (t1 - t0 <= 0 || sum(tt >= t0 - 1e-9 & tt <= t1 + 1e-9) < 2) return(NA_real_)
      g <- seq(t0, t1, by = dt)                       # uniform, year-aligned grid
      y <- stats::approx(tt, df[[col]], xout = g, rule = 2)$y
      if (length(g) < 2 || anyNA(y)) return(NA_real_)
      sum(diff(g) * (utils::head(y, -1) + utils::tail(y, -1)) / 2) / (max(g) - min(g))
    }
    
    if (tmax - min(tt) < 2 * period) {            # fallback: not two full periods yet
      i2 <- which.min(abs(tt - tmax)); i1 <- which.min(abs(tt - (tmax - period)))
      m_last <- as.numeric(df[i2, cols]); m_prev <- as.numeric(df[i1, cols])
    } else {
      m_last <- vapply(cols, mean_over, numeric(1), t0 = tmax - period,     t1 = tmax)
      m_prev <- vapply(cols, mean_over, numeric(1), t0 = tmax - 2 * period, t1 = tmax - period)
    }
    
    abs_drift <- abs(m_last - m_prev)
    rel_drift <- abs_drift / pmax(abs(m_prev), abs_floor)
    res <- data.frame(pool = cols, mean_prev = m_prev, mean_last = m_last,
                      abs_drift = abs_drift, rel_drift = rel_drift,
                      above_floor = abs(m_prev) >= abs_floor)
  }

  res[order(-res$rel_drift), ]
}

# ------------------------------------------------------------
# spinup_until_stable(): SEASONAL integration, lengthening until the annual
# means stop drifting (only meaningful pools count). PNG plots each iteration.
#   tol        max allowed year-over-year relative drift of annual means
#   abs_floor  pools with mean below this (g C m-2) are ignored in the test
# ------------------------------------------------------------
spinup_until_stable <- function(init_state, parms, model_fn,
                                n_years = 100, by = 5,
                                max_iter = 10, tol = 1e-3, abs_floor = 1e-3,
                                verbose = TRUE, plot_spinup = TRUE,
                                plot_dir = "Data/spinup") {
  if (missing(model_fn) || is.null(model_fn))
    stop("spinup_until_stable(): supply model_fn (e.g. obj$wrapped_model).")
  check_by(by)
  state <- init_state
  parms$climate_forcing <- make_climate_forcing(parms)

  for (i in seq_len(max_iter)) {
    if (verbose) cat("\n--- Seasonal spin-up iteration", i, "(", round(n_years), "yr ) ---\n")
    times <- seq(0, 365 * n_years, by = by)
    out   <- deSolve::ode(y = state, times = times, func = model_fn, parms = parms)

    if (plot_spinup)
      save_trajectory_png(out, file.path(plot_dir,
        sprintf("trajectory_%dpools_iter%d.png", length(state), i)))

    stab      <- check_stability(out, abs_floor = abs_floor)
    use       <- stab[stab$above_floor, ]
    max_drift <- if (nrow(use)) max(use$rel_drift, na.rm = TRUE) else 0
    worst     <- if (nrow(use)) use$pool[which.max(use$rel_drift)] else NA
    if (verbose)
      cat("Max annual-mean drift:", signif(max_drift, 3), "in", worst,
          "(", sum(!stab$above_floor), "negligible pools ignored )\n")

    if (is.finite(max_drift) && max_drift < tol) {
      if (verbose) cat("Seasonal limit cycle reached.\n")
      return(list(out = out, final_state = final_state(out), stability = stab,
                  max_drift = max_drift, converged = TRUE, iterations = i))
    }
    state   <- final_state(out)
    if(i == 3){
      n_years <- n_years * 5
    }
  }
  if (verbose) cat("Reached max iterations; worst drift", signif(max_drift, 3),
                   "in", worst, "- inspect the PNGs / raise max_iter or tol.\n")
  list(out = out, final_state = final_state(out), stability = stab,
       max_drift = max_drift, converged = FALSE, iterations = max_iter)
}

# Last row of an ODE output, as a named numeric state vector (drops time and
# the mass_balance_check column).
final_state <- function(out) {
  df   <- as.data.frame(out)
  cols <- setdiff(names(df), c("time", "mass_balance_check"))
  v    <- as.numeric(df[nrow(df), cols])
  setNames(v, cols)
}
