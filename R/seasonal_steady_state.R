#------------------------------------------------------------------------#
# Seasonal steady state (limit cycle) for repeated model evaluations ----
#------------------------------------------------------------------------#
# Used by the sensitivity, uncertainty and identifiability analyses
# (Scripts/4, 5, 7) to evaluate the model at its SEASONAL steady state:
# the annual limit cycle under the seasonal climate forcing
# (make_climate_forcing), rather than the equilibrium under constant mean
# forcing. Because the model is nonlinear the two differ (10-20% for some
# soil pools), and field data are collected under seasonal climate.
#
# The limit cycle is the fixed point of the one-year map Phi (integrate one
# year from y): Phi(y) = y. It is found by Newton's method on
# G(y) = Phi(y) - y, as in spinup_limit_cycle() (R/spinup.R), but made cheap
# enough for thousands of evaluations:
#   * the Jacobian dG/dy (n one-year solves by finite differences) is taken
#     from the NOMINAL run and only UPDATED between iterations (Broyden's
#     method: one one-year solve per iteration); it is recomputed by finite
#     differences only if an iteration fails to reduce the residual;
#   * the first guess is the constant-forcing equilibrium of the evaluated
#     parameters (cheap, and it fixes the solution branch, e.g. animals present
#     or not), scaled pool by pool by the nominal ratio cycle / equilibrium.
#
# Which forcing the analyses use is set by forcing_mode in
# R/sensitivity_settings.R ("seasonal" or "constant").

#------------------------------------------------------------------------#
# limit_cycle_solve() ----
#------------------------------------------------------------------------#
# Used in: called by steady_state()
# Broyden-Newton shooting for the annual limit cycle.
#   y0        first guess for the cycle state at t = 0 (1 January)
#   parms     model parameters (seasonal forcing is (re)built here)
#   model_fn  derivative function (obj$wrapped_model)
#   J0        starting Jacobian of G = Phi(y) - y (e.g. from the nominal
#             run); NULL = finite differences at y0
#   tol       stop when max_i |G_i| / max(|y_i|, floor) < tol
#   rtol, atol  ODE tolerances of the one-year solves (must be below tol; the
#             one-year solve dominates the cost: rtol 1e-7 is ~1.6x faster
#             than 1e-9)
#   max_forced  times a full Newton step is taken although it does not lower
#             the residual (Newton often overshoots far from the cycle)
#   shrink_max  the Newton step is shortened so that no pool shrinks by more
#             than this factor in one step (keeps pools positive without
#             clipping them to an absorbing 0)
# Returns list(state, converged, residual, J, J_fresh, n_solves, iterations).
limit_cycle_solve <- function(y0, parms, model_fn, J0 = NULL, tol = 1e-6,
                              max_iter = 40, max_forced = 3, fd_rel = 1e-5, floor = 1e-6,
                              rtol = 1e-7, atol = 1e-7, shrink_max = 10) {
  parms$climate_forcing <- make_climate_forcing(parms)
  pools   <- names(y0)
  n_solve <- 0L
  Phi <- function(v) {
    n_solve <<- n_solve + 1L
    out <- tryCatch(suppressWarnings(deSolve::ode(y = v, times = c(0, 365), func = model_fn,
                                                  parms = parms, rtol = rtol, atol = atol)),
                    error = function(e) NULL)
    if (is.null(out) || nrow(out) < 2) return(rep(NA_real_, length(v)))
    as.numeric(out[2, pools])
  }
  rel_norm <- function(g, v) max(abs(g) / pmax(abs(v), floor))
  fd_jac <- function(v, Pv) {
    J <- matrix(0, length(v), length(v))
    for (j in seq_along(v)) {
      h  <- fd_rel * max(abs(v[j]), floor)
      vj <- v; vj[j] <- vj[j] + h
      J[, j] <- (Phi(vj) - Pv) / h
    }
    J - diag(length(v))
  }

  y  <- y0
  Py <- Phi(y)
  G  <- Py - y
  if (any(!is.finite(G)))
    return(list(state = y, converged = FALSE, residual = NA_real_, J = J0,
                J_fresh = FALSE, n_solves = n_solve, iterations = 0L))
  rn <- rel_norm(G, y)
  fresh <- is.null(J0)
  J  <- if (fresh) fd_jac(y, Py) else J0
  it <- 0L; n_forced <- 0L
  while (rn >= tol && it < max_iter) {
    it   <- it + 1L
    step <- tryCatch(-solve(J, G), error = function(e)
      tryCatch(-solve(J + 1e-8 * diag(length(y)), G), error = function(e) NULL))
    ok <- FALSE; best <- NULL
    if (!is.null(step) && all(is.finite(step))) {
      # positivity (fraction to the boundary): shorten the WHOLE step so that
      # no pool shrinks more than shrink_max-fold. Clipping single pools at 0
      # (absorbing: B = 0 stays 0) or at a floor (changes the step direction)
      # both let Newton walk into a false, microbe-free cycle.
      dn  <- step < 0 & y > 0
      lam <- if (any(dn)) min(1, min(y[dn] * (1 - 1 / shrink_max) / -step[dn])) else 1
      for (ls in seq_len(if (fresh) 8 else 3)) {   # short search on a Broyden Jacobian
        y_new <- pmax(y + lam * step, 0)
        G_new <- Phi(y_new) - y_new
        if (all(is.finite(G_new))) {
          rn_new <- rel_norm(G_new, y_new)
          if (rn_new < rn) { ok <- TRUE; break }
          if (is.null(best)) best <- list(y = y_new, G = G_new, rn = rn_new)
        }
        lam <- lam / 2
      }
    }
    if (!ok) {
      if (!fresh) {                             # stale Broyden Jacobian: refresh it
        J <- fd_jac(y, G + y); fresh <- TRUE
        next
      }
      # fresh Jacobian and still no decrease of the max-norm residual: take the
      # full Newton step anyway (as spinup_limit_cycle() does), a few times at most
      n_forced <- n_forced + 1L
      if (is.null(best) || n_forced > max_forced) break
      y_new <- best$y; G_new <- best$G; rn_new <- best$rn
    }
    s <- y_new - y; dG <- G_new - G
    if (sum(s^2) > 0) J <- J + outer(as.numeric(dG - J %*% s), s) / sum(s^2)   # Broyden
    fresh <- FALSE
    y <- y_new; G <- G_new; rn <- rn_new
  }
  list(state = setNames(y, pools), converged = is.finite(rn) && rn < tol, residual = rn,
       J = J, J_fresh = fresh, n_solves = n_solve, iterations = it)
}

#------------------------------------------------------------------------#
# cycle_trajectory() ----
#------------------------------------------------------------------------#
# Used in: called by steady_state()
# One year along the limit cycle at daily resolution (t = 0..365) from the
# cycle state at t = 0. Returns a matrix: rows = day 0..365, columns = pools.
cycle_trajectory <- function(state, parms, model_fn, rtol = 1e-7, atol = 1e-7) {
  parms$climate_forcing <- make_climate_forcing(parms)
  out <- deSolve::ode(y = state, times = 0:365, func = model_fn, parms = parms,
                      rtol = rtol, atol = atol)
  as.matrix(out[, names(state), drop = FALSE])
}

#------------------------------------------------------------------------#
# steady_state() ----
#------------------------------------------------------------------------#
# Used in: called by animal_effect_totalC(), headline_effects(),
#   targets_for(), make_base_pair(), ident_observations()
# Needs (set in R/sensitivity_settings.R): forcing_mode, eq_max_time, lc_tol, lc_ode_tol,
#   lc_reuse_jacobian
# The steady state of one setup object under the chosen forcing.
#   mode = "constant"  equilibrium under constant mean forcing (runsteady),
#                      the previous behaviour
#   mode = "seasonal"  annual limit cycle under the seasonal forcing
#   ref         nominal reference from an earlier steady_state() call on the
#               same arm (obj$ss_ref): supplies the first-guess scaling and
#               the starting Jacobian (seasonal mode only)
#   warm_start  starting state for the constant-forcing equilibrium (e.g. the
#               baseline equilibrium when solving the treatment)
#   reuse_jacobian  TRUE: start from the nominal Jacobian ref$J (Broyden);
#               best for small perturbations (identifiability). FALSE: fresh
#               finite-difference Jacobian; cheaper on average for points far
#               from nominal (Morris, Latin hypercube)
#   fallback_years  if Newton fails from the scaled guess and from the
#               equilibrium, integrate this many years forward and try again
#   collapse_check_rel  a cycle whose annual mean of some pool is below this
#               fraction of that pool at the constant-forcing equilibrium is
#               re-solved from the equilibrium; the collapse-free cycle is kept
#               if one is found (a real collapse is kept and flagged downstream)
# Adds to obj:
#   ss_mean     the state used as "the steady state": annual mean of the
#               limit cycle (seasonal) or the equilibrium (constant)
#   ss_const    constant-forcing equilibrium (always computed; first guess)
#   ss_cycle    cycle state at t = 0 (1 January)            (seasonal only)
#   ss_traj     daily cycle trajectory, day 0..365          (seasonal only)
#   ss_ref      list(const, cycle, J) for later calls with ref = obj$ss_ref
#   spin_info   converged (TRUE/FALSE), method, n_solves, residual
steady_state <- function(obj, ref = NULL, warm_start = NULL, mode = forcing_mode,
                         tol = lc_tol, rtol = lc_ode_tol, atol = lc_ode_tol,
                         fallback_years = 20, reuse_jacobian = lc_reuse_jacobian,
                         collapse_check_rel = 1e-4) {
  mode <- match.arg(mode, c("seasonal", "constant"))
  if (mode == "constant") ref <- NULL      # constant forcing: same as before (cold start)
  ce <- spinup_equilibrium(obj, warm_start = if (!is.null(warm_start)) warm_start else ref$const,
                           max_time = eq_max_time, verbose = FALSE)
  obj$ss_const <- ce$init_state_spin
  obj$ss_mean  <- ce$init_state_spin
  obj$init_state_spin <- ce$init_state_spin
  obj$spin_info <- ce$spin_info
  obj$ss_ref   <- list(const = ce$init_state_spin, cycle = NULL, J = NULL)
  if (mode == "constant" || !isTRUE(ce$spin_info$converged)) return(obj)

  y0 <- ce$init_state_spin
  if (!is.null(ref$cycle) && !is.null(ref$const)) {
    sh <- intersect(names(y0), intersect(names(ref$cycle), names(ref$const)))
    r  <- ref$cycle[sh] / ref$const[sh]
    ok <- sh[is.finite(r) & r > 0 & ref$const[sh] > 0]
    # a reference in another regime (e.g. microbes collapsed) must not drag the
    # first guess there: the seasonal correction is limited to 0.1x .. 10x
    y0[ok] <- y0[ok] * pmin(pmax(r[ok], 0.1), 10)
  }
  J0 <- if (reuse_jacobian && !is.null(ref$J) && nrow(ref$J) == length(y0)) ref$J else NULL
  lc <- limit_cycle_solve(y0, obj$parms, obj$wrapped_model, J0 = J0, tol = tol,
                          rtol = rtol, atol = atol)
  n_sol <- lc$n_solves
  # a cycle in which a pool has (nearly) vanished although it is clearly
  # present at the constant-forcing equilibrium is suspect: near-zero
  # microbes regrow so slowly that the cycle condition is met in a false,
  # microbe-free state. Re-solve from the equilibrium itself (no scaling) and
  # keep that cycle if it is free of collapse.
  collapsed <- function(st) {
    if (is.null(st)) return(FALSE)
    m <- colMeans(cycle_trajectory(st, obj$parms, obj$wrapped_model, rtol = rtol,
                                   atol = atol)[1:365, , drop = FALSE])
    ce_v <- ce$init_state_spin[names(m)]
    any(ce_v > 1e-9 & m < collapse_check_rel * ce_v)
  }
  if (lc$converged && collapsed(lc$state)) {
    lc2 <- limit_cycle_solve(ce$init_state_spin, obj$parms, obj$wrapped_model, J0 = NULL,
                             tol = tol, rtol = rtol, atol = atol)
    n_sol <- n_sol + lc2$n_solves
    if (lc2$converged && !collapsed(lc2$state)) lc <- lc2
  }
  if (!lc$converged) {                       # retry 1: plain Newton from the equilibrium
    lc <- limit_cycle_solve(ce$init_state_spin, obj$parms, obj$wrapped_model, J0 = NULL,
                            tol = tol, rtol = rtol, atol = atol)
    n_sol <- n_sol + lc$n_solves
  }
  if (!lc$converged) {                       # retry 2: integrate forward, then Newton
    p <- obj$parms; p$climate_forcing <- make_climate_forcing(p)
    fw <- tryCatch(suppressWarnings(deSolve::ode(y = ce$init_state_spin,
                                times = 365 * (0:fallback_years), func = obj$wrapped_model,
                                parms = p, rtol = rtol, atol = atol)),
                   error = function(e) NULL)
    if (!is.null(fw) && nrow(fw) == fallback_years + 1 &&
        all(is.finite(fw[nrow(fw), names(y0)]))) {
      lc <- limit_cycle_solve(pmax(fw[nrow(fw), names(y0)], 0), obj$parms, obj$wrapped_model,
                              J0 = NULL, tol = tol, rtol = rtol, atol = atol)
      n_sol <- n_sol + lc$n_solves + fallback_years
    }
  }
  obj$spin_info <- list(method = "limit_cycle", converged = lc$converged,
                        residual = lc$residual, n_solves = n_sol)
  if (!lc$converged) return(obj)
  tr <- cycle_trajectory(lc$state, obj$parms, obj$wrapped_model, rtol = rtol, atol = atol)
  obj$ss_cycle <- lc$state
  obj$ss_traj  <- tr
  obj$ss_mean  <- colMeans(tr[1:365, , drop = FALSE])     # days 0..364
  obj$ss_ref   <- list(const = ce$init_state_spin, cycle = lc$state, J = lc$J)
  obj
}

#------------------------------------------------------------------------#
# cycle_state_on() ----
#------------------------------------------------------------------------#
# Used in: called by ident_observations()
# State on day `doy` (0 = 1 January, 0..364) of a solved seasonal steady
# state; for constant forcing the equilibrium itself.
cycle_state_on <- function(obj, doy = 0) {
  if (is.null(obj$ss_traj)) return(obj$ss_mean)
  obj$ss_traj[doy + 1, ]
}
