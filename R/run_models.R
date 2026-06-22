# ============================================================
# run_models.R - run scenarios across MULTIPLE models in one call
# ------------------------------------------------------------
# Wraps run_scenario() so you can compare models without re-running blocks
# by hand. Requires setup.R, spinup.R, compare_functions.R sourced.
# ============================================================

# ------------------------------------------------------------
# run_models(): every scenario x every model.
#   models     character vector, e.g. c("century","millennial","MIMICS")
#   scen       scenario list from read_scenarios()
#   scenarios  which scenarios to run (default: all)
#   ...        passed to run_scenario() (method, spinup_treatment, stol, ...)
#
# Returns list(
#   results   nested list: results[[model]][[scenario]] = run_scenario() output
#   combined  one tidy data frame stacking every eq_compare with model+scenario
# )
# ------------------------------------------------------------
run_models <- function(models, scen, scenarios = names(scen), ...) {
  results <- setNames(vector("list", length(models)), models)

  for (m in models) {
    results[[m]] <- setNames(
      lapply(scenarios, function(s) {
        out <- tryCatch(run_scenario(s, scen = scen, model = m, ...),
                        error = function(e) {
                          message("  !! ", m, " / ", s, " failed: ", conditionMessage(e))
                          NULL
                        })
        out
      }),
      scenarios
    )
  }

  combined <- do.call(rbind, lapply(models, function(m) {
    do.call(rbind, lapply(scenarios, function(s) {
      r <- results[[m]][[s]]
      if (is.null(r) || is.null(r$eq_compare) || !nrow(r$eq_compare)) return(NULL)
      data.frame(model = m, scenario = s, r$eq_compare,
                 row.names = NULL, stringsAsFactors = FALSE)
    }))
  }))

  list(results = results, combined = combined)
}

# ------------------------------------------------------------
# eq_compare_list(): pull a NAMED list of eq_compare data frames (one per
# model) for a single scenario -- the input shape plot_eq_compare() expects.
# ------------------------------------------------------------
eq_compare_list <- function(run_out, scenario) {
  models <- names(run_out$results)
  out <- lapply(models, function(m) {
    r <- run_out$results[[m]][[scenario]]
    if (is.null(r)) NULL else r$eq_compare
  })
  names(out) <- models
  out[!vapply(out, is.null, logical(1))]
}

# ------------------------------------------------------------
# spin_summary(): quick convergence table across all model x scenario runs,
# so you can see at a glance what converged (and how long it took).
# ------------------------------------------------------------
spin_summary <- function(run_out) {
  rows <- list()
  for (m in names(run_out$results)) for (s in names(run_out$results[[m]])) {
    r <- run_out$results[[m]][[s]]; if (is.null(r)) next
    for (arm in c("baseline", "treatment")) {
      si <- r$spin[[arm]]
      rows[[length(rows) + 1]] <- data.frame(
        model = m, scenario = s, arm = arm,
        method = if (is.null(si$method)) NA else si$method,
        converged = if (is.null(si$converged)) NA else si$converged,
        max_deriv = if (is.null(si$max_deriv)) NA else si$max_deriv,
        time_to_steady = if (is.null(si$time_to_steady)) NA else si$time_to_steady,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
