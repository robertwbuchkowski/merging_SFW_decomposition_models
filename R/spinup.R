# ============================================================
# SPIN-UP & STABILITY HELPERS  (simple, continuous, debuggable)
# ------------------------------------------------------------
# No year-by-year Poincare iteration. Spin-up is just ONE continuous
# integration under seasonal forcing for a length you choose. If the run
# isn't at a stable repeating cycle yet, increase `n_years` and re-run.
#
# The models are called DIRECTLY by deSolve (no wrapper). They return
# mass_balance_check as an extra output column, so you can watch it.
# ============================================================

# ------------------------------------------------------------
# spinup_run(): one continuous ODE integration (seasonal forcing).
#   y0       named initial-state vector (full model state)
#   model_fn the model function, e.g. century_model
#   parms    parameter list (must already contain everything; this sets
#            parms$climate_forcing to the seasonal forcing)
#   n_years  how long to run -- INCREASE THIS if not yet stable
#   by       output step in days (1 = daily; 365 = yearly snapshots)
# Returns the full deSolve output matrix (so you can plot/inspect it).
# ------------------------------------------------------------
spinup_run <- function(y0, model_fn, parms, n_years = 200, by = 1,
                       method = "lsoda", ...) {
  parms$climate_forcing <- make_climate_forcing(parms)
  times <- seq(0, 365 * n_years, by = by)
  deSolve::ode(y = y0, times = times, func = model_fn, parms = parms,
               method = method, ...)
}

# Last row of an ODE output, as a named numeric state vector (drops time
# and the mass_balance_check column). Use to continue from a spun-up state.
final_state <- function(out) {
  df   <- as.data.frame(out)
  cols <- setdiff(names(df), c("time", "mass_balance_check"))
  v    <- as.numeric(df[nrow(df), cols])
  setNames(v, cols)
}

# ------------------------------------------------------------
# check_stability(): compare the final year to the previous year at the
# same phase. Small rel_drift everywhere => spun up. If not, increase
# n_years in spinup_run() and try again.
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

# ------------------------------------------------------------
# (OPTIONAL) spinup_equilibrium(): fast warm-start to the annual-MEAN
# steady state under CONSTANT forcing, using rootSolve::runsteady. Handy
# to get a good starting point before spinup_run(); not required.
# ------------------------------------------------------------
spinup_equilibrium <- function(y0, model_fn, parms, maxtime = 1e6, ...) {
  parms$climate_forcing <- make_climate_forcing_equilibrium(parms)
  ss <- rootSolve::runsteady(y = y0, times = c(0, maxtime),
                             func = model_fn, parms = parms, ...)
  setNames(pmax(ss$y, 0), names(y0))
}
