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
# spinup_until_stable(): one continuous ODE integration (seasonal forcing).
#   init_state       named initial-state vector
#   model_fn the model function, e.g. century_model
#   parms    parameter list (must already contain everything; this sets
#            parms$climate_forcing to the seasonal forcing)
#   n_years  how long to run -- INCREASE THIS if not yet stable
#   by       output step in days (1 = daily; 365 = yearly snapshots)
# Returns the full deSolve output matrix (so you can plot/inspect it).
# ------------------------------------------------------------

spinup_until_stable <- function(init_state, 
                                parms,
                                model_fn = run$wrapped_model,
                                n_years = 100,
                                by = 1,
                                max_iter = 10,
                                tol = 1e-4,
                                verbose = TRUE) {
  
  state <- init_state
  parms$climate_forcing <- make_climate_forcing(parms)
  
  for (i in seq_len(max_iter)) {
    
    if (verbose) cat("\n--- Spin-up iteration", i, "---\n")
    
    times <- seq(0, 365 * n_years, by = by)
    out <- deSolve::ode(y = state, times = times, func = model_fn, parms = parms)
    
    # check stability
    stab <- check_stability(out)
    if (verbose) print(stab)
    
    # here assuming check_stability returns numeric drifts per pool
    max_drift <- max(abs(stab$rel_drift), na.rm = TRUE)
    
    if (verbose) cat("Max drift:", max_drift, "\n")
    
    if (max_drift < tol) {
      if (verbose) cat("System stabilized.\n")
      return(list(out = out,
                  final_state = final_state(out),
                  converged = TRUE,
                  iterations = i))
    }
    
    # update state for next loop
    state <- final_state(out)
    # Increase the simulation length
    n_years <- n_years*1.5
  }
  
  if (verbose) cat("Reached max iterations without full stability.\n")
  
  return(list(out = out,
              final_state = final_state(out),
              converged = FALSE,
              iterations = max_iter))
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


