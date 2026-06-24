# ============================================================
# SPIN-UP & STABILITY HELPERS
# ------------------------------------------------------------
# Two complementary tools, both easy to inspect:
#
#   spinup_equilibrium()  FAST warm-start to steady state under CONSTANT
#                         (annual-mean) forcing, by forward time-integration
#                         (rootSolve::runsteady). This REPLACES stode() for
#                         the equilibrium spin-up: it does not factorize the
#                         steady-state Jacobian, so the near-singular /
#                         ill-conditioned Jacobians that make stode() fail
#                         (e.g. MIMICS "negative diagonal at pool 7") are not
#                         a problem -- the trajectory simply flows to the
#                         attractor.
#
#   spinup_until_stable() OPTIONAL refinement: continuous SEASONAL integration
#                         until the annual limit cycle stops drifting.
#
# Models are called through the wrapper and return mass_balance_check as an
# extra output column, so you can always watch conservation.
# ============================================================

# ------------------------------------------------------------
# runsteady_spinup(): core engine. Integrate y forward in time to steady
# state under whatever forcing is already in `parms`. Returns a rich,
# debuggable result rather than just the state.
#   y         named initial-state vector (the active/working state)
#   model_fn  the (wrapped) model function
#   parms     parameter list (must already contain $climate_forcing)
#   max_time  integration horizon in days (raise if not converged)
#   stol      steady-state tolerance on the relative rate of change
# Returns list(final_state, converged, time_to_steady, max_deriv, raw).
# ------------------------------------------------------------
runsteady_spinup <- function(y, model_fn, parms,
                             max_time = 1e7, stol = 1e-8, verbose = TRUE) {
  ss <- rootSolve::runsteady(
    y     = y,
    time  = c(0, max_time),
    func  = model_fn,
    parms = parms,
    stol  = stol
  )

  y_eq <- ss$y
  names(y_eq) <- names(y)

  converged <- isTRUE(attr(ss, "steady"))
  t_steady  <- attr(ss, "time")

  # max |derivative| at the returned state = how close to true steady state
  d_end    <- model_fn(t_steady, y_eq, parms)[[1]]
  max_drv  <- max(abs(d_end))

  if (verbose) {
    cat(sprintf("  runsteady: %s | reached t = %.3g d | max|dy/dt| = %.2e\n",
                if (converged) "CONVERGED" else "NOT converged (raise max_time/relax stol)",
                t_steady, max_drv))
  }
  list(final_state = y_eq, converged = converged,
       time_to_steady = t_steady, max_deriv = max_drv, raw = ss)
}

# ------------------------------------------------------------
# spinup_equilibrium(): run runsteady on a setup object (from setup_model /
# setup_scenario). Spins up under CONSTANT forcing, then restores seasonal
# forcing on the returned object. Stores results on the object so you can
# inspect them: obj$init_state_spin and obj$spin_info.
#   warm_start  optional named vector to start from (e.g. a baseline
#               equilibrium) -- shared pools override the object's init.
# ------------------------------------------------------------
spinup_equilibrium <- function(obj, warm_start = NULL,
                               max_time = 1e7, stol = 1e-8, verbose = TRUE) {

  y0 <- obj$working_state
  if (!is.null(warm_start)) {
    shared <- intersect(names(y0), names(warm_start))
    y0[shared] <- warm_start[shared]
    if (verbose && length(shared))
      cat(sprintf("  warm-started %d shared pool(s) from a previous equilibrium\n",
                  length(shared)))
  }

  # constant annual-mean forcing for the steady-state solve
  obj$parms$climate_forcing <- make_climate_forcing_equilibrium(obj$parms)
  res <- runsteady_spinup(y0, obj$wrapped_model, obj$parms,
                          max_time = max_time, stol = stol, verbose = verbose)
  # restore seasonal forcing for any later seasonal runs
  obj$parms$climate_forcing <- make_climate_forcing(obj$parms)

  obj$init_state_spin <- res$final_state
  obj$spin_info <- list(method = "runsteady", converged = res$converged,
                        time_to_steady = res$time_to_steady,
                        max_deriv = res$max_deriv)
  obj
}

# ------------------------------------------------------------
# spinup_until_stable(): OPTIONAL seasonal refinement. Continuous integration
# under seasonal forcing, lengthening until the year-over-year drift is small.
# Start it from an equilibrium state (e.g. spinup_equilibrium()$init_state_spin)
# for fast settling onto the limit cycle.
# ------------------------------------------------------------
spinup_until_stable <- function(init_state, parms, model_fn,
                                n_years = 100, by = 1,
                                max_iter = 10, tol = 1e-4, verbose = TRUE, plot_spinup = TRUE) {
  if (missing(model_fn) || is.null(model_fn))
    stop("spinup_until_stable(): supply model_fn (e.g. obj$wrapped_model).")
  state <- init_state
  parms$climate_forcing <- make_climate_forcing(parms)

  for (i in seq_len(max_iter)) {
    if (verbose) cat("\n--- Seasonal spin-up iteration", i,
                     "(", round(n_years), "yr ) ---\n")
    times <- seq(0, 365 * n_years, by = by)
    out   <- deSolve::ode(y = state, times = times, func = model_fn, parms = parms)

    if(plot_spinup){
      pdf(paste0("Data/spinup/tragectory_",length(state),"_", i,".pdf"))
      plot(out)
      dev.off()
    }
    
    stab <- check_stability(out)
    max_drift <- max(abs(stab$rel_drift), na.rm = TRUE)
    max_drift_pool <- stab$pool[which(stab$rel_drift == max_drift)]
    if (verbose) cat("Max year-over-year drift:", signif(max_drift, 3), "in",max_drift_pool, "\n")

    if (max_drift < tol) {
      if (verbose) cat("Seasonal limit cycle reached.\n")
      return(list(out = out, final_state = final_state(out),
                  converged = TRUE, iterations = i))
    }
    state   <- final_state(out)
    n_years <- n_years * 1.5
  }
  if (verbose) cat("Reached max iterations without full seasonal stability.\n")
  list(out = out, final_state = final_state(out), converged = FALSE,
       iterations = max_iter)
}

# Last row of an ODE output, as a named numeric state vector (drops time and
# the mass_balance_check column). Use to continue from a spun-up state.
final_state <- function(out) {
  df   <- as.data.frame(out)
  cols <- setdiff(names(df), c("time", "mass_balance_check"))
  v    <- as.numeric(df[nrow(df), cols])
  setNames(v, cols)
}

# ------------------------------------------------------------
# check_stability(): compare the final year to the previous year at the same
# phase. Small rel_drift everywhere => spun up onto a stable cycle.
# ------------------------------------------------------------
check_stability <- function(out, period = 365) {
  df   <- as.data.frame(out)
  tt   <- df[[1]]
  cols <- setdiff(names(df), c("time", "mass_balance_check"))
  i2   <- which.min(abs(tt - max(tt)))
  i1   <- which.min(abs(tt - (max(tt) - period)))
  v2   <- as.numeric(df[i2, cols]); v1 <- as.numeric(df[i1, cols])
  res  <- data.frame(pool = cols, prev_year = v1, last_year = v2,
                     rel_drift = abs(v2 - v1) / pmax(abs(v1), 1e-8))
  res[order(-res$rel_drift), ]
}
