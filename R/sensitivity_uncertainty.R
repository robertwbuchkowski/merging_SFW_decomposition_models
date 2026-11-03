#------------------------------------------------------------------------#
# Sensitivity analysis: uncertainty propagation (Latin hypercube) ----
#------------------------------------------------------------------------#
# Helper functions for Scripts/4_uncertainty_animal_effect.R. They were moved here from the scripts so
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
# Steady-state animal effects on total soil + root C for one parameter vector
# (equilibria under constant mean forcing), plus the baseline values:
#   total        treatment - baseline, all pathways on
#   direct_dmii  treatment with the trait-mediated indirect effects switched
#                off (zero_indirect_effects()) - baseline: the DIRECT effect
#                PLUS the density-mediated indirect effects
#   tmii         TRAIT-MEDIATED INDIRECT effect = total - direct_dmii
headline_effects <- function(base_pair, pv, pools) {
  pair <- base_pair
  for (nm in names(pv)) {
    pair$treatment <- set_param(pair$treatment, nm, pv[[nm]])
    pair$baseline  <- set_param(pair$baseline,  nm, pv[[nm]])
  }
  out <- c(total = NA_real_, tmii = NA_real_, direct_dmii = NA_real_, converged = 0,
           microbe_collapse = 0, kb_floor = 0, infeasible = 0, base_min_pool = NA_real_,
           setNames(rep(NA_real_, length(pools)), paste0("base_", pools)))
  if (nzchar(infeasible_soil(pair$baseline$parms))) { out[["infeasible"]] <- 1; return(out) }
  tryCatch({
    # equilibria under constant mean forcing: baseline, treatment, and the
    # treatment with the trait-mediated indirect effects switched off
    eqm <- function(o, ws = NULL) spinup_equilibrium(o, warm_start = ws, max_time = eq_max_time,
                                                     verbose = FALSE)
    b  <- eqm(pair$baseline)
    eb <- b$init_state_spin
    tr <- eqm(pair$treatment, eb)
    td <- eqm(zero_indirect_effects(pair$treatment, verbose = FALSE), eb)
    tr$ss_mean <- tr$init_state_spin; td$ss_mean <- td$init_state_spin
    soil <- setdiff(intersect(names(eb), names(tr$ss_mean)), animal_pools)
    conv <- isTRUE(b$spin_info$converged) && isTRUE(tr$spin_info$converged) &&
            isTRUE(td$spin_info$converged)
    rg_t <- regime_flags(eb, tr$ss_mean, pair$treatment$parms)
    rg_d <- regime_flags(eb, td$ss_mean, pair$treatment$parms)
    tot <- sum(tr$ss_mean[soil]) - sum(eb[soil])
    dd  <- sum(td$ss_mean[soil]) - sum(eb[soil])
    c(total = tot, tmii = tot - dd, direct_dmii = dd,
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
# read_user_targets() ----
#------------------------------------------------------------------------#
# Used in: called by run_uncertainty_propagation_fn(), run_sobol_scenario()
# Needs (set in R/sensitivity_settings.R): up_targets_file
# User-edited realism targets (columns scenario, pool, lower, upper), or NULL.
read_user_targets <- function(file = up_targets_file) {
  if (!file.exists(file)) return(NULL)
  u <- read.csv(file, stringsAsFactors = FALSE)
  if (!all(c("scenario", "pool", "lower", "upper") %in% names(u)))
    stop(file, " needs columns scenario, pool, lower, upper")
  u$scenario <- rename_scenarios(u$scenario)
  cat("Baseline realism targets: using", file, "\n")
  u
}

#------------------------------------------------------------------------#
# targets_for() ----
#------------------------------------------------------------------------#
# Used in: called by run_uncertainty_propagation_fn()
# Needs (set in 4_uncertainty_animal_effect.R): up_filter_factor, up_filter_factor_by_pool, up_filter_pools
# Uses base_pair$nominal_b (make_base_pair(), nominal no-animal equilibrium)
# realism targets for one scenario: defaults from the all-default baseline,
# overridden / extended by rows in up_targets_file
targets_for <- function(scenario, base_pair, user_targets) {
  eb  <- base_pair$nominal_b         # nominal no-animal equilibrium (make_base_pair())
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
# Used in: Scripts/4_uncertainty_animal_effect.R (section: UNCERTAINTY PROPAGATION TO THE HEADLINE ANIMAL EFFECT)
# Needs (set in 4_uncertainty_animal_effect.R): cl, scenarios, up_n, up_pool_floor, up_seed,
#   up_rarefy_n, up_rarefy_step,
#   up_targets_file
run_uncertainty_propagation_fn <- function() {
  user_targets <- read_user_targets()

  samp_rows <- summ_rows <- diag_rows <- input_rows <- target_rows <- list()
  keep <- list()            # per scenario: outputs needed for the rarefied summary
  q <- function(v, p) if (length(v)) unname(stats::quantile(v, p, na.rm = TRUE)) else NA_real_
  # summary of the total and trait-mediated indirect (tmii) effects (and of
  # direct_dmii = direct + density-mediated indirect) over the samples in idx
  summ <- function(Y, default, scenario, idx, label) {
    tot <- Y[idx, "total"]; tm <- Y[idx, "tmii"]; dd <- Y[idx, "direct_dmii"]
    pct <- 100 * tm / tot; pct <- pct[is.finite(pct)]
    data.frame(
      scenario = scenario, subset = label, n = sum(idx),
      total_default = default[["total"]],
      total_median = q(tot, 0.5), total_mean = if (length(tot)) mean(tot) else NA_real_,
      total_q025 = q(tot, 0.025), total_q975 = q(tot, 0.975),
      p_total_positive = if (length(tot)) mean(tot > 0) else NA_real_,
      p_sign_as_default = if (length(tot)) mean(sign(tot) == sign(default[["total"]])) else NA_real_,
      tmii_default = default[["tmii"]],
      tmii_median = q(tm, 0.5), tmii_q025 = q(tm, 0.025), tmii_q975 = q(tm, 0.975),
      pct_tmii_default = 100 * default[["tmii"]] / default[["total"]],
      pct_tmii_median = q(pct, 0.5), pct_tmii_q025 = q(pct, 0.025),
      pct_tmii_q975 = q(pct, 0.975),
      direct_dmii_default = default[["direct_dmii"]],
      direct_dmii_median = q(dd, 0.5), direct_dmii_q025 = q(dd, 0.025),
      direct_dmii_q975 = q(dd, 0.975))
  }

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
    ok   <- Y[, "converged"] == 1 & is.finite(Y[, "total"]) & is.finite(Y[, "direct_dmii"])
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

    summ_rows[[scenario]] <- rbind(summ(Y, default, scenario, ok, "all converged"),
                                   summ(Y, default, scenario, ok & !flag, "unflagged"),
                                   summ(Y, default, scenario, accepted, "accepted"))
    keep[[scenario]] <- list(Y = Y, default = default, accepted = accepted, converged = ok)
  }

  # RAREFACTION: every scenario is summarised over the SAME number of
  # samples, drawn at random (without replacement), so the quantiles of the
  # scenarios rest on equal sample sizes. Done for two subsets:
  #   accepted       converged + baseline realism filter (main table)
  #   all converged  converged only: the same subsetting as the Sobol'
  #                  analysis (impossible soils and unconverged runs removed,
  #                  no realism filter) -- the supplementary table
  rarefy_subset <- function(set, label, col, seed_off) {
    n_set  <- vapply(keep, function(k) sum(k[[set]]), 0)
    n_rare <- rarefy_size(n_set, up_rarefy_n, up_rarefy_step)
    cat(sprintf("Rarefaction (%s): every scenario summarised over %d samples (%s)\n",
                label, n_rare, paste(names(n_set), n_set, sep = " = ", collapse = ", ")))
    for (si in seq_along(scenarios)) {
      scenario <- scenarios[si]; k <- keep[[scenario]]
      set.seed(up_seed + seed_off + si)
      idx  <- which(k[[set]])
      pick <- if (n_rare > 0) idx[sort(sample.int(length(idx), n_rare))] else integer(0)
      r <- seq_len(nrow(k$Y)) %in% pick
      samp_rows[[scenario]][[col]] <<- r
      summ_rows[[scenario]] <<- rbind(summ_rows[[scenario]],
                                      summ(k$Y, k$default, scenario, r, paste0(label, " (rarefied)")))
      diag_rows[[scenario]][[paste0("n_", col)]] <<- n_rare
    }
    n_rare
  }
  cat("\n")
  n_rare   <- rarefy_subset("accepted",  "accepted",      "rarefied",           1000)
  n_rare_c <- rarefy_subset("converged", "all converged", "rarefied_converged", 2000)

  diagnostics <- bind_rows(diag_rows)
  # % of ALL draws accepted (the acceptance_rate column is of converged draws)
  diagnostics$pct_accepted_of_draws  <- 100 * diagnostics$n_accepted / diagnostics$n_samples
  diagnostics$pct_converged_of_draws <- 100 * diagnostics$n_converged / diagnostics$n_samples
  summary <- bind_rows(summ_rows)
  idx <- match(summary$scenario, diagnostics$scenario)
  summary <- cbind(summary, diagnostics[idx, c("n_samples", "n_accepted", "pct_accepted_of_draws",
                                               "n_converged", "pct_converged_of_draws",
                                               "n_rarefied", "n_rarefied_converged")],
                   row.names = NULL)
  list(samples = bind_rows(samp_rows), summary = summary,
       diagnostics = diagnostics, inputs = bind_rows(input_rows),
       targets = bind_rows(target_rows), n_rarefied = n_rare,
       n_rarefied_converged = n_rare_c)
}

#------------------------------------------------------------------------#
# rarefy_size() ----
#------------------------------------------------------------------------#
# Used in: called by run_uncertainty_propagation_fn()
# The common number of samples every scenario is rarefied to (n_accepted =
# the number in the subset per scenario):
#   fixed  a number (up_rarefy_n): used as is (must not exceed the smallest
#          accepted count)
#   NULL   the largest multiple of `step` (e.g. 500 -> 1000, 1500, 2000)
#          not exceeding the smallest accepted count; if fewer than `step`
#          samples were accepted, the largest multiple of 100, else all of them.
rarefy_size <- function(n_accepted, fixed = NULL, step = 500) {
  n_min <- min(n_accepted)
  if (!is.null(fixed)) {
    if (fixed > n_min)
      stop("up_rarefy_n = ", fixed, " but the smallest number of accepted samples is ",
           n_min, ": lower up_rarefy_n or raise up_n.")
    return(as.integer(fixed))
  }
  n <- floor(n_min / step) * step
  if (n == 0) n <- floor(n_min / 100) * 100
  if (n == 0) {
    warning("Fewer than 100 accepted samples in some scenario: rarefied to ", n_min, ".")
    n <- n_min
  }
  as.integer(n)
}

#------------------------------------------------------------------------#
# animal_effect_paper_table() ----
#------------------------------------------------------------------------#
# Used in: Scripts/4_uncertainty_animal_effect.R
# The manuscript tables: per scenario, the total effect and the
# trait-mediated indirect effect (g C m-2) at the default parameter values
# and as "median (2.5% to 97.5% quantile)" over the rarefied draws of a
# subset, the trait-mediated indirect effect as % of the total (default
# parameters), and the numbers of draws kept and used. Values rounded to two
# significant figures, with a thousands separator (text, ready to paste).
#   subset  "accepted"      main table (converged + baseline realism filter)
#           "all converged" supplementary table (same subsetting as Sobol')
animal_effect_paper_table <- function(summary, subset = c("accepted", "all converged")) {
  subset <- match.arg(subset)
  f2 <- function(x) {
    s <- formatC(signif(x, 2), format = "fg", digits = 2, flag = "#", big.mark = ",")
    s <- sub("\\.$", "", trimws(s))
    s[!is.finite(x)] <- "NA"
    s
  }
  ci <- function(m, lo, hi) sprintf("%s (%s to %s)", f2(m), f2(lo), f2(hi))
  d <- summary[summary$subset == paste0(subset, " (rarefied)"), , drop = FALSE]
  lab <- if (exists("pretty_scenario")) unname(pretty_scenario(d$scenario)) else d$scenario
  if (subset == "accepted") {
    n_kept <- d$n_accepted;  pct_kept <- d$pct_accepted_of_draws;  kept <- "accepted"
  } else {
    n_kept <- d$n_converged; pct_kept <- d$pct_converged_of_draws; kept <- "converged"
  }
  out <- data.frame(
    Scenario = lab,
    n_kept   = n_kept,
    pct_kept = sprintf("%.0f", pct_kept),
    `N used` = d$n,
    `Total default` = f2(d$total_default),
    `Total median (95% interval)` = ci(d$total_median, d$total_q025, d$total_q975),
    `Trait-mediated indirect default` = f2(d$tmii_default),
    `Trait-mediated indirect median (95% interval)` = ci(d$tmii_median, d$tmii_q025, d$tmii_q975),
    `Trait-mediated indirect (% of Total)` = f2(d$pct_tmii_default),
    check.names = FALSE, stringsAsFactors = FALSE)
  names(out)[2:3] <- c(paste("N", kept), paste("%", kept))
  out
}
