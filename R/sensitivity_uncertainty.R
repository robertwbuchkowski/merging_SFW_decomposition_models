#------------------------------------------------------------------------#
# Sensitivity analysis: uncertainty propagation (Latin hypercube) ----
#------------------------------------------------------------------------#
# Helper functions for Scripts/5_uncertainty_animal_effect.R. They were moved here from the scripts so
# that every custom function is defined BEFORE the script's analysis code
# runs (sourced at the top of the script), which also lets
# debugonce_functions() flag them all for step-through debugging.
#
# These functions read settings and data that the script creates at top level
# (listed per function under 'Needs'). They are looked up when the function is
# CALLED, so source this file first and run the script's settings before use.

#------------------------------------------------------------------------#
# q_trunc_norm() ----
#------------------------------------------------------------------------#
# Used in: called by up_quantile()
# inverse-CDF samplers (u in (0,1))
q_trunc_norm <- function(u, mean, sd, lo, hi) {
  a <- pnorm(lo, mean, sd); b <- pnorm(hi, mean, sd)
  qnorm(a + u * (b - a), mean, sd)
}

#------------------------------------------------------------------------#
# up_distribution() ----
#------------------------------------------------------------------------#
# Used in: called by run_uncertainty_propagation_fn()
# Sampling distribution for one scenario x parameter: a normal centred on the
# model value with the standard deviation from param_sd_info() (reported SD,
# else CV = cv_default, reduced by a reported Min/Max), truncated to the
# analysed range (+/- 2 SD, reported Min/Max, physical bounds).
up_distribution <- function(scenario, p, d0) {
  s <- param_sd_info(scenario, p, d0)
  list(dist = "truncated normal", lo = s$lo, hi = s$hi, centre = s$centre,
       spread = s$sd, cv = s$cv, sd_source = s$source, bound_hit = s$bound_hit)
}

#------------------------------------------------------------------------#
# up_quantile() ----
#------------------------------------------------------------------------#
# Used in: called by run_uncertainty_propagation_fn()
# Inverse CDF of the sampling distribution from up_distribution() at u in (0,1).
up_quantile <- function(u, ds) q_trunc_norm(u, ds$centre, ds$spread, ds$lo, ds$hi)

#------------------------------------------------------------------------#
# lhs_unit() ----
#------------------------------------------------------------------------#
# Used in: called by run_uncertainty_propagation_fn()
# Latin hypercube in (0,1)^k: one stratified draw per row in each column
lhs_unit <- function(n, k) {
  vapply(seq_len(k), function(j) (sample.int(n) - runif(n)) / n, numeric(n))
}

#------------------------------------------------------------------------#
# baseline_values() ----
#------------------------------------------------------------------------#
# Used in: called by headline_effects(), targets_for()
# Needs (set in R/sensitivity_settings.R): animal_pools
# baseline values used by the realism filter ("TotalC" = all non-animal pools)
baseline_values <- function(eb, pools) {
  soil <- setdiff(names(eb), animal_pools)
  v <- c(TotalC = sum(eb[soil]), eb[soil])
  setNames(vapply(pools, function(p) if (p %in% names(v)) unname(v[[p]]) else NA_real_,
                  numeric(1)), paste0("base_", pools))
}

#------------------------------------------------------------------------#
# headline_effects() ----
#------------------------------------------------------------------------#
# Used in: called by make_headline_fun(), run_uncertainty_propagation_fn()
# Needs (set in R/sensitivity_settings.R): animal_pools, eq_max_time
# total and direct equilibrium effects + baseline values for one parameter vector
headline_effects <- function(base_pair, pv, pools) {
  pair <- base_pair
  for (nm in names(pv)) {
    pair$treatment <- set_param(pair$treatment, nm, pv[[nm]])
    pair$baseline  <- set_param(pair$baseline,  nm, pv[[nm]])
  }
  out <- c(total = NA_real_, direct = NA_real_, converged = 0,
           microbe_collapse = 0, kb_floor = 0, infeasible = 0, base_min_pool = NA_real_,
           setNames(rep(NA_real_, length(pools)), paste0("base_", pools)))
  if (nzchar(infeasible_soil(pair$baseline$parms))) { out[["infeasible"]] <- 1; return(out) }
  tryCatch({
    b  <- spinup_equilibrium(pair$baseline, max_time = eq_max_time, verbose = FALSE)
    eb <- b$init_state_spin
    tr <- spinup_equilibrium(pair$treatment, warm_start = eb, max_time = eq_max_time,
                             verbose = FALSE)
    td <- spinup_equilibrium(zero_indirect_effects(pair$treatment, verbose = FALSE),
                             warm_start = eb, max_time = eq_max_time, verbose = FALSE)
    soil <- setdiff(intersect(names(eb), names(tr$init_state_spin)), animal_pools)
    conv <- isTRUE(b$spin_info$converged) && isTRUE(tr$spin_info$converged) &&
            isTRUE(td$spin_info$converged)
    rg_t <- regime_flags(eb, tr$init_state_spin, pair$treatment$parms)
    rg_d <- regime_flags(eb, td$init_state_spin, pair$treatment$parms)
    c(total  = sum(tr$init_state_spin[soil]) - sum(eb[soil]),
      direct = sum(td$init_state_spin[soil]) - sum(eb[soil]),
      converged = as.numeric(conv),
      microbe_collapse = as.numeric(grepl("microbe_collapse", paste(rg_t, rg_d))),
      kb_floor         = as.numeric(grepl("kb_floor", rg_t)),
      infeasible       = 0,
      base_min_pool    = min(eb[setdiff(names(eb), animal_pools)]),  # realism: all pools > 0
      baseline_values(eb, pools))
  }, error = function(e) out)
}

#------------------------------------------------------------------------#
# make_headline_fun() ----
#------------------------------------------------------------------------#
# Used in: called by run_uncertainty_propagation_fn()
# factory so a parallel worker receives only what the evaluator needs
make_headline_fun <- function(base_pair, pools) {
  force(base_pair); force(pools)
  function(pv) headline_effects(base_pair, pv, pools)
}

#------------------------------------------------------------------------#
# targets_for() ----
#------------------------------------------------------------------------#
# Used in: called by run_uncertainty_propagation_fn()
# Needs (set in 5_uncertainty_animal_effect.R): up_filter_factor, up_filter_factor_by_pool, up_filter_pools
# Needs (set in R/sensitivity_settings.R): eq_max_time
# realism targets for one scenario: defaults from the all-default baseline,
# overridden / extended by rows in up_targets_file
targets_for <- function(scenario, base_pair, user_targets) {
  eb  <- spinup_equilibrium(base_pair$baseline, max_time = eq_max_time,
                            verbose = FALSE)$init_state_spin
  ref <- baseline_values(eb, up_filter_pools)
  fac <- ifelse(up_filter_pools %in% names(up_filter_factor_by_pool),
                up_filter_factor_by_pool[up_filter_pools], up_filter_factor)
  tg  <- data.frame(scenario = scenario, pool = up_filter_pools,
                    reference = unname(ref), lower = unname(ref) / fac,
                    upper = unname(ref) * fac, factor = unname(fac), source = "default",
                    stringsAsFactors = FALSE)
  tg <- tg[is.finite(tg$reference), , drop = FALSE]
  if (!is.null(user_targets)) {
    u <- user_targets[user_targets$scenario == scenario, , drop = FALSE]
    for (i in seq_len(nrow(u))) {
      j <- which(tg$pool == u$pool[i])
      if (length(j)) { tg$lower[j] <- u$lower[i]; tg$upper[j] <- u$upper[i]; tg$source[j] <- "user" }
      else tg <- rbind(tg, data.frame(scenario = scenario, pool = u$pool[i],
                                      reference = if ("reference" %in% names(u)) u$reference[i] else NA_real_,
                                      lower = u$lower[i], upper = u$upper[i], factor = NA_real_,
                                      source = "user"))
    }
  }
  tg
}

#------------------------------------------------------------------------#
# run_uncertainty_propagation_fn() ----
#------------------------------------------------------------------------#
# Used in: Scripts/5_uncertainty_animal_effect.R (section: UNCERTAINTY PROPAGATION TO THE HEADLINE ANIMAL EFFECT)
# Needs (set in 5_uncertainty_animal_effect.R): cl, scenarios, up_n, up_pool_floor, up_seed,
#   up_targets_file
run_uncertainty_propagation_fn <- function() {
  user_targets <- if (file.exists(up_targets_file)) {
    u <- read.csv(up_targets_file, stringsAsFactors = FALSE)
    if (!all(c("scenario", "pool", "lower", "upper") %in% names(u)))
      stop(up_targets_file, " needs columns scenario, pool, lower, upper")
    u$scenario <- rename_scenarios(u$scenario)
    cat("Baseline realism targets: using", up_targets_file, "\n"); u
  } else NULL

  samp_rows <- summ_rows <- diag_rows <- input_rows <- target_rows <- list()
  q <- function(v, p) if (length(v)) unname(stats::quantile(v, p, na.rm = TRUE)) else NA_real_

  for (si in seq_along(scenarios)) {
    scenario   <- scenarios[si]
    base_pair  <- make_base_pair(scenario)
    parms_here <- base_pair$treatment$parms
    ps         <- param_set_for(scenario, parms_here)

    # sampling distributions
    dists <- list()
    for (p in ps$pnames) {
      d0 <- parms_here[[p]]
      if (!is.finite(d0) || d0 == 0) next
      ds <- up_distribution(scenario, p, d0)
      if (is.finite(ds$lo) && is.finite(ds$hi) && ds$hi > ds$lo) dists[[p]] <- ds
    }
    pn <- names(dists); k <- length(pn)
    input_rows[[scenario]] <- data.frame(
      scenario = scenario, parameter = pn, parameter_label = pretty_param(pn),
      param_type = ps$ptype[pn], default = vapply(pn, function(p) parms_here[[p]], 0),
      distribution = vapply(dists, `[[`, "", "dist"),
      sd_source = vapply(dists, `[[`, "", "sd_source"),
      sd = vapply(dists, `[[`, 0, "spread"), cv = vapply(dists, `[[`, 0, "cv"),
      lo = vapply(dists, `[[`, 0, "lo"), hi = vapply(dists, `[[`, 0, "hi"),
      bound_hit = vapply(dists, `[[`, "", "bound_hit"),
      row.names = NULL)

    # realism targets
    tg <- targets_for(scenario, base_pair, user_targets)
    target_rows[[scenario]] <- tg
    pools <- tg$pool

    cat("Uncertainty propagation: ", scenario, " - ", k, " parameters, ",
        up_n, " LHS samples\n", sep = "")
    set.seed(up_seed + si)
    U <- lhs_unit(up_n, k)
    X <- vapply(seq_len(k), function(j) up_quantile(U[, j], dists[[pn[j]]]), numeric(up_n))
    if (up_n == 1) X <- matrix(X, nrow = 1)
    colnames(X) <- pn

    default <- headline_effects(base_pair, setNames(numeric(0), character(0)), pools)
    rows_list <- lapply(seq_len(up_n), function(i) setNames(X[i, ], pn))
    hf <- make_headline_fun(base_pair, pools)
    Y  <- if (is.null(cl)) lapply(rows_list, hf) else parallel::parLapply(cl, rows_list, hf)
    Y  <- do.call(rbind, Y)

    # classification
    inf  <- Y[, "infeasible"] == 1
    ok   <- Y[, "converged"] == 1 & is.finite(Y[, "total"]) & is.finite(Y[, "direct"])
    flag <- Y[, "microbe_collapse"] == 1 | Y[, "kb_floor"] == 1
    pass_pool <- vapply(seq_len(nrow(tg)), function(j) {
      v <- Y[, paste0("base_", tg$pool[j])]
      is.finite(v) & v >= tg$lower[j] & v <= tg$upper[j]
    }, logical(up_n))
    if (up_n == 1) pass_pool <- matrix(pass_pool, nrow = 1)
    # realism filter: every listed pool within its target AND every baseline
    # pool positive and non-zero (above up_pool_floor)
    pos_ok   <- is.finite(Y[, "base_min_pool"]) & Y[, "base_min_pool"] > up_pool_floor
    accepted <- ok & apply(pass_pool, 1, all) & pos_ok

    n_inf <- sum(inf); n_bad <- sum(!ok) - n_inf
    cat(sprintf("  %d/%d realistic samples converged; %d unrealistic soils removed\n",
                sum(ok), up_n - n_inf, n_inf))
    cat(sprintf("  baseline realism filter: %d of %d converged samples accepted (%.0f%%)\n",
                sum(accepted), sum(ok), 100 * sum(accepted) / max(sum(ok), 1)))
    if (sum(ok) > 0 && sum(accepted) / sum(ok) < 0.1)
      warning(sprintf("Uncertainty propagation %s: only %d of %d samples pass the baseline realism filter -- the parameter distributions may be wider than plausible, or the targets too strict.",
                      scenario, sum(accepted), sum(ok)))
    if (n_bad > 0)
      warning(sprintf("Uncertainty propagation %s: %d samples did NOT reach a stable state; excluded.",
                      scenario, n_bad))

    samp_rows[[scenario]] <- data.frame(scenario = scenario, sample = seq_len(up_n),
                                        accepted = accepted, Y, X, check.names = FALSE)

    rej <- c(setNames(colSums(!pass_pool[ok, , drop = FALSE]), paste0("n_reject_", tg$pool)),
             n_reject_nonpositive_pool = sum(!pos_ok[ok]))
    diag_rows[[scenario]] <- data.frame(
      scenario = scenario, n_samples = up_n, n_infeasible = n_inf,
      n_nonconverged = n_bad, n_converged = sum(ok),
      n_microbe_collapse = sum(ok & Y[, "microbe_collapse"] == 1),
      n_kb_floor = sum(ok & Y[, "kb_floor"] == 1),
      n_accepted = sum(accepted), acceptance_rate = sum(accepted) / max(sum(ok), 1),
      t(rej), check.names = FALSE)

    summ <- function(idx, label) {
      tot <- Y[idx, "total"]; dir <- Y[idx, "direct"]
      pct <- 100 * dir / tot; pct <- pct[is.finite(pct)]
      data.frame(
        scenario = scenario, subset = label, n = sum(idx),
        total_default = default[["total"]],
        total_median = q(tot, 0.5), total_mean = if (length(tot)) mean(tot) else NA_real_,
        total_q025 = q(tot, 0.025), total_q975 = q(tot, 0.975),
        p_total_positive = if (length(tot)) mean(tot > 0) else NA_real_,
        p_sign_as_default = if (length(tot)) mean(sign(tot) == sign(default[["total"]])) else NA_real_,
        direct_default = default[["direct"]],
        direct_median = q(dir, 0.5), direct_q025 = q(dir, 0.025), direct_q975 = q(dir, 0.975),
        pct_direct_default = 100 * default[["direct"]] / default[["total"]],
        pct_direct_median = q(pct, 0.5), pct_direct_q025 = q(pct, 0.025),
        pct_direct_q975 = q(pct, 0.975))
    }
    summ_rows[[scenario]] <- rbind(summ(ok, "all converged"),
                                   summ(ok & !flag, "unflagged"),
                                   summ(accepted, "accepted"))
  }
  list(samples = bind_rows(samp_rows), summary = bind_rows(summ_rows),
       diagnostics = bind_rows(diag_rows), inputs = bind_rows(input_rows),
       targets = bind_rows(target_rows))
}
