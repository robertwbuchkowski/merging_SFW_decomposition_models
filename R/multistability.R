#------------------------------------------------------------------------#
# Multistability test: initial conditions, clustering, per-scenario test ----
#------------------------------------------------------------------------#
# Helper functions for Scripts/6_multistability_test.R. They were moved here from the script so
# that every custom function is defined BEFORE the script's analysis code
# runs (sourced at the top of the script), which also lets
# debugonce_functions() flag them all for step-through debugging.
#
# These functions read settings and data that the script creates at top level
# (listed per function under 'Needs'). They are looked up when the function is
# CALLED, so source this file first and run the script's settings before use.

#------------------------------------------------------------------------#
# make_starts() ----
#------------------------------------------------------------------------#
# Used in: called by test_scenario()
# Needs (set in 6_multistability_test.R): abs_floor
# make_starts(): a matrix of initial states (rows = pools, cols = starts).
# Columns are: the default; an all-low and an all-high corner; a microbial-
# collapse corner (tiny MIC/B); a microbial-bloom corner; then random draws
# where every pool is multiplied by an independent log-uniform factor. This
# spreads seeds across basins, including the ones microbial models flip between.
make_starts <- function(y0, n_random, span = 2) {
  pools <- names(y0)
  cols  <- list(default = y0)

  cols$all_low  <- y0 * 10^(-span)
  cols$all_high <- y0 * 10^(+span)
  micdead <- y0; micdead[intersect(c("MIC","B"), pools)] <- abs_floor * 10
  cols$microbial_collapse <- micdead
  micbloom <- y0; micbloom[intersect(c("MIC","B"), pools)] <-
    y0[intersect(c("MIC","B"), pools)] * 10^span
  cols$microbial_bloom <- micbloom

  for (i in seq_len(n_random)) {
    fac <- 10^runif(length(y0), -span, span)
    cols[[paste0("rand", i)]] <- y0 * fac
  }
  M <- do.call(cbind, cols)
  rownames(M) <- pools
  M
}

#------------------------------------------------------------------------#
# same_state() ----
#------------------------------------------------------------------------#
# Used in: called by cluster_states()
# Needs (set in 6_multistability_test.R): abs_floor, rel_tol
# distinct_states(): cluster converged equilibria into unique states by the
# max RELATIVE pool difference (pools below abs_floor ignored). Greedy: each
# state joins the first cluster it matches, else starts a new one.
same_state <- function(a, b) {
  keep <- (abs(a) > abs_floor) | (abs(b) > abs_floor)
  if (!any(keep)) return(TRUE)
  denom <- pmax(abs(a[keep]), abs(b[keep]), abs_floor)
  max(abs(a[keep] - b[keep]) / denom) < rel_tol
}

#------------------------------------------------------------------------#
# cluster_states() ----
#------------------------------------------------------------------------#
# Used in: called by test_scenario()
cluster_states <- function(states) {          # list of named numeric vectors
  reps <- list(); members <- integer(0)
  for (i in seq_along(states)) {
    hit <- NA_integer_
    for (k in seq_along(reps))
      if (same_state(states[[i]], reps[[k]])) { hit <- k; break }
    if (is.na(hit)) { reps[[length(reps)+1]] <- states[[i]]; members <- c(members, 1L) }
    else members[hit] <- members[hit] + 1L
  }
  list(reps = reps, counts = members)
}

#------------------------------------------------------------------------#
# soil_pools_of() ----
#------------------------------------------------------------------------#
# Used in: called by test_scenario()
# run one scenario: seed many starts, keep converged states, cluster them.
soil_pools_of <- function(nm) setdiff(nm, c("Earthworm","Detritivore","RootHerb"))

#------------------------------------------------------------------------#
# test_scenario() ----
#------------------------------------------------------------------------#
# Used in: Scripts/6_multistability_test.R (section: run all scenarios)
# Needs (set in 6_multistability_test.R): model, n_starts, scen, span_dec, use_treatment
test_scenario <- function(scenario) {
  obj <- setup_scenario(model, scen, scenario, animals = use_treatment)
  y0  <- obj$working_state
  S   <- make_starts(y0, n_starts, span = span_dec)

  eqs <- list(); conv <- logical(ncol(S))
  for (j in seq_len(ncol(S))) {
    start <- S[, j]; names(start) <- rownames(S)
    o <- tryCatch(spinup_equilibrium(obj, warm_start = start, verbose = FALSE),
                  error = function(e) NULL)
    ok <- !is.null(o) && isTRUE(o$spin_info$converged) && all(is.finite(o$init_state_spin))
    conv[j] <- ok
    if (ok) eqs[[length(eqs)+1]] <- o$init_state_spin
  }

  cl <- cluster_states(eqs)
  n_states <- length(cl$reps)

  # per-state total soil + root C, to describe how different the states are
  state_tab <- map_dfr(seq_len(n_states), function(k) {
    s  <- cl$reps[[k]]; sp <- soil_pools_of(names(s))
    hdr <- tibble(scenario = scenario, state = k, n_starts = cl$counts[k],
                  total_soil_C = sum(s[sp]))
    pool_cols <- as_tibble(as.list(round(s, 4)))       # one column per pool
    bind_cols(hdr, pool_cols)
  })

  list(
    summary = tibble(
      scenario         = scenario,
      n_converged      = sum(conv),
      n_tested         = ncol(S),
      n_stable_states  = n_states,
      multistable      = n_states > 1,
      min_total_C      = if (n_states) min(state_tab$total_soil_C) else NA_real_,
      max_total_C      = if (n_states) max(state_tab$total_soil_C) else NA_real_,
      spread_pct       = if (n_states > 1)
        round(100 * (max(state_tab$total_soil_C) - min(state_tab$total_soil_C)) /
                max(abs(state_tab$total_soil_C)), 1) else 0),
    states = state_tab)
}
