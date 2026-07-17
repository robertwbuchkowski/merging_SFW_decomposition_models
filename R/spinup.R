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
# sized for the number of panels (default 4 columns).
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
# (means are phase-independent, so an uneven `by` does not create spurious
# drift). Pools whose mean is below `abs_floor` are flagged negligible and
# excluded from the convergence decision (so a near-zero pool can't dominate).
# Returns a table sorted by relative drift.
# ------------------------------------------------------------
# ------------------------------------------------------------
# check_by(): the output step MUST divide the year, or the seasonal cycle is
# sampled at a drifting phase and the annual-mean stability test never settles.
# 365 = 5 x 73, so the only clean steps are 1, 5, 73 and 365. check_stability()
# interpolates to a daily grid so a bad `by` still gives a correct test, but it
# costs seasonal resolution, hence the warning.
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
    res <- data.frame(pool = cols, abs_drift = abs(as.numeric(abs_drift)), rel_drift = abs(as.numeric(rel_drift)),
                      above_floor = T)
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
  
  # t(sapply(1:365, parms$climate_forcing)) %>% as.data.frame() %>% tibble() %>% mutate(doy = 1:365)  %>% pivot_longer(!doy) %>% ggplot(aes(x = doy, y = value)) + geom_point() + facet_wrap(.~name, scales = "free_y")

  for (i in seq_len(max_iter)) {
    if (verbose) cat("\n--- Seasonal spin-up iteration", i, "(", round(n_years), "yr ) ---\n")
    times <- seq(0, 365 * n_years, by = by)
    out   <- deSolve::ode(y = state, times = times, func = model_fn, parms = parms)

    if (plot_spinup)
      save_trajectory_png(out, file.path(plot_dir,
        sprintf("trajectory_%dpools_iter%d_%s_%s.png", length(state), i, model, scenario)))

    stab      <- check_stability(out, abs_floor = abs_floor)
    use       <- stab[stab$above_floor, ]
    max_drift <- if (nrow(use)) max(use$rel_drift, na.rm = TRUE) else 0
    worst     <- if (nrow(use)) use$pool[which.max(use$rel_drift)] else NA
    if (verbose)
      cat("Max annual-mean drift:", signif(max_drift, 3), "in", worst,
          "(", sum(!stab$above_floor), "negligible pools ignored )\n")
    cat("Annual drift in all pools is: \n")
    print(stab)

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

# ============================================================
# LIMIT-CYCLE SPIN-UP BY SHOOTING / NEWTON
# ------------------------------------------------------------
# A periodically-forced ODE has a limit cycle y*(t) with period P (= 365 d).
# Plain forward integration ("run for N years until it stops changing") can only
# approach it at the rate of the SLOWEST eigenvalue of the year-map, so a system
# with an ~18-yr mode needs many centuries to settle. Instead we solve directly
# for the periodic orbit:
#
#   define the one-year (stroboscopic) map    Phi(y0) = state after 365 d,
#   the limit cycle start satisfies           G(y0) = Phi(y0) - y0 = 0,
#   solve G(y0)=0 by Newton:                   y0 <- y0 - (M - I)^{-1} G(y0),
#   where M = dPhi/dy0 is the monodromy matrix (finite-differenced by
#   integrating one year per perturbed coordinate).
#
# Convergence is Newton-quadratic and INDEPENDENT of the slow eigenvalue: a few
# one-year integrations, not hundreds of years. Stability of the found cycle is
# read off for free from the eigenvalues of M (all |lambda| < 1 => stable).
# ============================================================

# one-year map: integrate exactly one period and return the end state (ordered
# to match names(y0)). `by` only sets output density; period end is exact.
.year_map <- function(y0, parms, model_fn, period = 365, by = 5) {
  times <- sort(unique(c(seq(0, period, by = by), period)))
  out   <- deSolve::ode(y = y0, times = times, func = model_fn, parms = parms,
                        atol = 1e-9, rtol = 1e-9)
  fs <- final_state(out)              # drops time + mass_balance_check
  fs[names(y0)]
}

# ------------------------------------------------------------
# spinup_limit_cycle(): Newton shooting for the seasonal limit cycle.
#   y0          starting guess for the cycle's t=0 state (e.g. an equilibrium
#               from spinup_equilibrium(), or any reasonable state)
#   model_fn    the derivative function (obj$wrapped_model)
#   tol         stop when max|Phi(y0)-y0| / max(|y0|, floor) < tol
#   max_newton  max Newton iterations (each costs length(y0)+1 one-year solves)
#   damping     step is scaled by this if a full Newton step increases the
#               residual or drives a pool negative (simple line-search safeguard)
# Returns the cycle start state, the within-year trajectory, the residual, the
# monodromy spectral radius (|lambda|max), and convergence flag.
# ------------------------------------------------------------
spinup_limit_cycle <- function(y0, parms, model_fn,
                               period = 365, by = 5,
                               tol = 1e-8, max_newton = 20,
                               fd_rel = 1e-6, floor = 1e-6,
                               damping = 0.5, verbose = TRUE) {
  if (missing(model_fn) || is.null(model_fn))
    stop("spinup_limit_cycle(): supply model_fn (e.g. obj$wrapped_model).")
  parms$climate_forcing <- make_climate_forcing(parms)
  y  <- y0
  n  <- length(y)
  I  <- diag(n)

  resid_of <- function(v) .year_map(v, parms, model_fn, period, by) - v
  rel_norm <- function(g, v) max(abs(g) / pmax(abs(v), floor))

  G  <- resid_of(y)
  rho <- NA_real_
  for (it in seq_len(max_newton)) {
    rn <- rel_norm(G, y)
    if (verbose) cat(sprintf("  Newton %2d: max relative |Phi(y)-y| = %.3e\n", it, rn))
    if (is.finite(rn) && rn < tol) {
      if (verbose) cat("  Limit cycle found.\n")
      break
    }

    # monodromy M = dPhi/dy0 by forward differences (one one-year solve per column)
    Phi0 <- G + y                                  # Phi(y) = G + y
    M <- matrix(0, n, n)
    for (j in seq_len(n)) {
      h  <- fd_rel * max(abs(y[j]), floor)
      yj <- y; yj[j] <- yj[j] + h
      M[, j] <- (.year_map(yj, parms, model_fn, period, by) - Phi0) / h
    }
    rho <- tryCatch(max(Mod(eigen(M, only.values = TRUE)$values)), error = function(e) NA_real_)

    step <- tryCatch(solve(M - I, G), error = function(e) {           # (M-I) dy = G
      solve(M - I + 1e-8 * I, G)                                      # tiny ridge if singular
    })

    # damped line search: shrink step until residual drops and stays feasible
    lam <- 1
    repeat {
      y_new <- y - lam * step
      if (all(is.finite(y_new)) && all(y_new > -floor)) {
        G_new <- resid_of(y_new)
        if (rel_norm(G_new, y_new) < rn || lam < 1e-3) break
      }
      lam <- lam * damping
      if (lam < 1e-6) break
    }
    y <- pmax(y_new, 0)                             # clamp tiny negatives from FD
    G <- resid_of(y)
  }

  out <- deSolve::ode(y = y, times = sort(unique(c(seq(0, period, by = by), period))),
                      func = model_fn, parms = parms, atol = 1e-9, rtol = 1e-9)
  list(final_state = y, out = out,
       residual = rel_norm(G, y),
       monodromy_rho = rho,
       stable = is.finite(rho) && rho < 1 + 1e-6,
       converged = is.finite(rel_norm(G, y)) && rel_norm(G, y) < tol,
       iterations = it)
}

# ------------------------------------------------------------
# dynamic_spinup_newton(): finds the limit cycle by Newton shooting and returns
# the same shape as dynamic_spinup() (list with final_state / out / converged).
# Warm-starts from the constant-forcing equilibrium
# (cheap and close), which makes Newton converge in a handful of iterations.
# ------------------------------------------------------------
dynamic_spinup_newton <- function(obj, from = NULL, by = 5,
                                  tol = 1e-8, max_newton = 20, verbose = TRUE) {
  y0 <- if (!is.null(from)) {
    from[names(obj$working_state)]
  } else if (!is.null(obj$init_state_spin)) {
    obj$init_state_spin[names(obj$working_state)]         # constant-forcing equilibrium
  } else {
    obj$working_state
  }
  if (any(is.na(y0))) stop("dynamic_spinup_newton(): start state missing pools of the setup.")

  res <- spinup_limit_cycle(y0, obj$parms, obj$wrapped_model,
                            by = by, tol = tol, max_newton = max_newton, verbose = verbose)
  if (verbose) {
    cat(sprintf("Newton limit-cycle: residual %.2e, monodromy rho = %.4f (%s), %d iters\n",
                res$residual, res$monodromy_rho,
                if (isTRUE(res$stable)) "stable" else "UNSTABLE / check", res$iterations))
  }
  list(out = res$out, final_state = res$final_state,
       converged = res$converged, stable = res$stable,
       monodromy_rho = res$monodromy_rho, residual = res$residual,
       iterations = res$iterations, method = "newton_shooting")
}
