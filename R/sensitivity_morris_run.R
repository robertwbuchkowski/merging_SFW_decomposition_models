#------------------------------------------------------------------------#
# Sensitivity analysis: Morris drivers, outputs and figures ----
#------------------------------------------------------------------------#
# Helper functions for the sensitivity scripts. They were moved here from the scripts so
# that every custom function is defined BEFORE the script's analysis code
# runs (sourced at the top of the script), which also lets
# debugonce_functions() flag them all for step-through debugging.
#
# These functions read settings and data that the script creates at top level
# (listed per function under 'Needs'). They are looked up when the function is
# CALLED, so source this file first and run the script's settings before use.

#------------------------------------------------------------------------#
# run_morris() ----
#------------------------------------------------------------------------#
# Used in: Scripts/4_sensitivity_animal_effects.R (section: RUN)
# Needs (set in 4_sensitivity_animal_effects.R): boot_B, boot_top_k, cl, morris_levels, morris_r, morris_seed, scenarios
# Needs (set in R/sensitivity_settings.R): forcing_mode
# run_morris(): Morris screening for every scenario under one range mode.
# Each scenario is seeded with morris_seed + its index, so the main and
# standardized runs share the same unit-hypercube trajectories.
run_morris <- function(mode = c("main", "standardized")) {
  mode <- match.arg(mode)
  morris_rows <- list()
  ee_all      <- list()
  for (si in seq_along(scenarios)) {
    scenario  <- scenarios[si]
    base_pair  <- make_base_pair(scenario)
    parms_here <- base_pair$treatment$parms

    ps     <- param_set_for(scenario, parms_here)
    pnames <- ps$pnames; ptype <- ps$ptype; a_here <- ps$animal

    lo <- hi <- setNames(numeric(length(pnames)), pnames)
    src <- hit <- setNames(character(length(pnames)), pnames)
    sdv <- cvv <- setNames(rep(NA_real_, length(pnames)), pnames)
    for (p in pnames) {
      d0 <- parms_here[[p]]
      if (!is.finite(d0) || d0 == 0) { lo[p] <- hi[p] <- NA; next }
      rg <- range_for(scenario, p, d0, mode)
      lo[p] <- rg$lo; hi[p] <- rg$hi; src[p] <- rg$source
      sdv[p] <- rg$sd; cvv[p] <- rg$cv; hit[p] <- rg$bound_hit
    }
    keep <- names(lo)[is.finite(lo) & is.finite(hi) & (hi > lo)]
    lo <- lo[keep]; hi <- hi[keep]; src <- src[keep]
    hits <- hit[keep][nzchar(hit[keep])]
    if (length(hits))
      cat(sprintf("  %d parameters reach a bound: %s\n", length(hits),
                  paste(sprintf("%s (%s)", names(hits), hits), collapse = ", ")))

    cat("Morris [", mode, "]: ", scenario, " - ", length(keep), " parameters, ",
        morris_r, " trajectories, ", forcing_mode, " forcing\n", sep = "")
    set.seed(morris_seed + si)
    # impossible soils are caught BEFORE any simulation: trajectories with a
    # point whose mean moisture >= porosity (or other infeasible soil) are redrawn
    accept_fun <- function(X) trajectory_feasible(X, lo, hi, parms_here)
    trajs <- morris_trajectories(length(keep), morris_r, levels = morris_levels,
                                 accept = accept_fun)
    n_redrawn <- attr(trajs, "n_redrawn")
    if (n_redrawn > 0)
      cat(sprintf("  %d trajectories redrawn before running (contained an impossible soil)\n",
                  n_redrawn))
    ee    <- morris_run(trajs, lo, hi,
                        eval_fun = make_effect_fun(base_pair),
                        verbose = FALSE, cl = cl)
    n_bad <- attr(ee, "n_nonconverged"); n_tot <- attr(ee, "n_eval")
    # impossible soils are skipped before solving: count them separately
    sr    <- attr(ee, "skip_reasons")
    n_inf <- if (length(sr)) sum(sr[grepl("^infeasible soil", names(sr))]) else 0L
    n_bad <- n_bad - n_inf
    if (n_inf > 0)
      message(sprintf("Morris [%s] %s: %d of %d points are unrealistic soils (porosity <= 0, or soil moisture < 0 or >= porosity) and were removed.",
                      mode, scenario, n_inf, n_tot))
    rc    <- attr(ee, "regime_counts")
    n_mic <- if ("microbe_collapse" %in% names(rc)) rc[["microbe_collapse"]] else 0L
    n_kb  <- if ("kb_floor" %in% names(rc)) rc[["kb_floor"]] else 0L
    if (n_mic + n_kb > 0)
      message(sprintf("Morris [%s] %s: flagged evaluations - microbe_collapse %d, kb_floor %d (of %d).",
                      mode, scenario, n_mic, n_kb, n_tot))
    if (n_bad > 0)
      warning(sprintf("Morris [%s] %s: %d of %d equilibrium evaluations did NOT reach a stable state (%.1f%%).",
                      mode, scenario, n_bad, n_tot, 100 * n_bad / max(n_tot, 1)))
    out_rng <- keep[!(vapply(keep, function(p) parms_here[[p]], numeric(1)) >= lo &
                      vapply(keep, function(p) parms_here[[p]], numeric(1)) <= hi)]
    if (length(out_rng))
      warning(sprintf("Morris [%s] %s: model value lies OUTSIDE its range for: %s",
                      mode, scenario, paste(out_rng, collapse = ", ")))
    cat(sprintf("  stable-state check: %d/%d realistic evaluations converged%s; %d unrealistic soils removed\n",
                n_tot - n_inf - n_bad, n_tot - n_inf, if (n_bad > 0) "  <-- SEE WARNING" else "", n_inf))

    sm <- morris_summary(ee)
    bt <- morris_bootstrap(ee, B = boot_B, top_k = boot_top_k,
                           seed = morris_seed + 1000 + si)
    if (!is.null(bt)) sm <- merge(sm, bt, by = "parameter", all.x = TRUE, sort = FALSE)
    # same metrics using only steps where neither end is in a flagged regime
    ee_clean <- ee[!nzchar(ee$regime), , drop = FALSE]
    sm_clean <- if (nrow(ee_clean)) morris_summary(ee_clean) else
      data.frame(parameter = character(), mu_star = numeric(), sigma = numeric(), n_ee = integer())
    names(sm_clean)[names(sm_clean) != "parameter"] <-
      paste0(names(sm_clean)[names(sm_clean) != "parameter"], "_clean")
    sm <- merge(sm, sm_clean, by = "parameter", all.x = TRUE, sort = FALSE)
    # bootstrap CI of mu* on the clean steps too (used by the convergence check)
    bt_c <- if (nrow(ee_clean)) morris_bootstrap(ee_clean, B = boot_B, top_k = boot_top_k,
                                                 seed = morris_seed + 2000 + si) else NULL
    if (!is.null(bt_c)) {
      bt_c <- data.frame(parameter = bt_c$parameter, mu_star_clean_lo = bt_c$mu_star_lo,
                         mu_star_clean_hi = bt_c$mu_star_hi)
      sm <- merge(sm, bt_c, by = "parameter", all.x = TRUE, sort = FALSE)
    } else { sm$mu_star_clean_lo <- NA_real_; sm$mu_star_clean_hi <- NA_real_ }
    sm$n_ee_flagged <- vapply(sm$parameter, function(p)
      sum(ee$parameter == p & nzchar(ee$regime)), integer(1))
    ee$scenario <- scenario
    ee_all[[scenario]] <- ee
    sm$scenario       <- scenario
    sm$range_source   <- src[sm$parameter]
    sm$sd_used        <- sdv[sm$parameter]
    sm$cv_used        <- cvv[sm$parameter]
    sm$bound_hit      <- hit[sm$parameter]
    sm$range_lo       <- lo[sm$parameter]
    sm$range_hi       <- hi[sm$parameter]
    sm$default        <- vapply(sm$parameter, function(p) parms_here[[p]], numeric(1))
    # flag ranges that do not contain the value actually used in the model (e.g. a
    # fitted animal parameter whose reported SD/Min-Max was built around the
    # sheet value, not the fitted value)
    sm$default_in_range <- sm$default >= sm$range_lo & sm$default <= sm$range_hi
    sm$param_type     <- ptype[sm$parameter]
    sm$is_animal      <- sm$parameter %in% a_here
    # elementary effects lost because an end point was removed or unsolved
    sm$n_ee_removed <- vapply(sm$parameter, function(p)
      sum(ee$parameter == p & !is.finite(ee$ee)), integer(1))
    sm$n_eval_infeasible <- n_inf
    sm$n_traj_redrawn    <- n_redrawn
    sm$n_nonconverged <- n_bad
    sm$n_evaluations  <- n_tot
    sm$n_eval_microbe_collapse <- n_mic
    sm$n_eval_kb_floor         <- n_kb
    morris_rows[[scenario]] <- sm
  }

  tbl <- bind_rows(morris_rows) %>%
    mutate(parameter_label = pretty_param(parameter)) %>%
    group_by(scenario) %>% arrange(scenario, desc(mu_star)) %>%
    mutate(rank = row_number()) %>%
    mutate(rank_clean = rank(-mu_star_clean, ties.method = "first", na.last = "keep")) %>%
    ungroup() %>%
    select(scenario, parameter, parameter_label, param_type, is_animal,
           range_source, sd_used, cv_used, default, range_lo, range_hi, bound_hit,
           default_in_range,
           mu_star, mu_star_lo, mu_star_hi, sigma, sigma_lo, sigma_hi,
           n_ee, rank, rank_median, rank_lo, rank_hi, p_top_k,
           mu_star_clean, mu_star_clean_lo, mu_star_clean_hi, sigma_clean, n_ee_clean,
           rank_clean, n_ee_flagged,
           n_ee_removed, n_eval_infeasible, n_traj_redrawn, n_nonconverged, n_evaluations,
           n_eval_microbe_collapse, n_eval_kb_floor)
  attr(tbl, "ee") <- bind_rows(ee_all)          # raw elementary effects
  tbl
}

#------------------------------------------------------------------------#
# plot_morris() ----
#------------------------------------------------------------------------#
# Used in: called by write_morris_outputs()
# FIGURE helper: traditional Morris mu* vs sigma map, one facet per scenario.
# Colour = how the range was set; shape = parameter type.
plot_morris <- function(tbl, colours, labels, legend_name = "Range source") {
  ggplot(tbl, aes(mu_star, sigma)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
    geom_point(aes(colour = range_source, shape = param_type), alpha = 0.85, size = 2) +
    ggrepel::geom_text_repel(aes(label = parameter_label), size = 2.2, max.overlaps = 14) +
    facet_wrap(~scenario, scales = "free", labeller = scenario_labeller) +
    scale_colour_manual(values = colours, labels = labels, name = legend_name, drop = FALSE) +
    scale_shape_manual(values = c(general = 16, `scenario-specific` = 15, animal = 17),
                       name = "Parameter type") +
    labs(x = if (identical(effect_output, "pct")) "mu* (mean |elementary effect|, % of no-animal soil C)" else
           expression(mu*"* (mean |elementary effect|, g C m"^-2*")"),
         y = expression(sigma*" (sd of elementary effect)")) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
}

#------------------------------------------------------------------------#
# tidytext_reorder() ----
#------------------------------------------------------------------------#
# Used in: called by write_morris_outputs()
# order a label within each facet (suffix trick, stripped again on the axis)
tidytext_reorder <- function(x, by, within) {
  key <- paste(x, within, sep = "___")
  factor(key, levels = unique(key[order(within, by)]))
}

#------------------------------------------------------------------------#
# morris_convergence() ----
#------------------------------------------------------------------------#
# Used in: called by write_morris_outputs()
# Convergence of the Morris screening for each scenario, following Sarrazin,
# Pianosi & Wagener (2016, Environ. Model. Softw. 79: 135-152), from the
# bootstrap 95% CIs of mu*:
#   INDEX convergence      the top_k parameters' mu* are precise: 95% CI width
#                          / mu* <= ci_tol for every one of them.
#   SCREENING convergence  the set of clearly influential parameters is stable.
#                          A top_k parameter is SECURE when its lower CI is above
#                          the upper CI of every parameter outside the top_k.
#                          Parameters whose CIs overlap across the top_k boundary
#                          are TIED: their mu* are near-equal, so which of them
#                          is "k-th" is arbitrary. A tie is not a failure; it is
#                          reported as a group (more trajectories only narrow
#                          the CIs, they cannot separate parameters whose mu*
#                          really are about equal).
# Status per scenario:
#   converged       index converged and no ties at the top_k boundary
#   tied            index converged, but some parameters are tied at the
#                   boundary (report the tied parameters as a group)
#   screening only  the top_k set is secure, but its mu* CIs are wider than
#                   ci_tol x mu* (the influential set is settled, the values
#                   of mu* are not yet precise)
#   unconverged     CIs too wide AND the top_k set is not secure: increase morris_r
# effects = "all" uses mu_star / mu_star_lo / mu_star_hi; "clean" uses the
# same metrics computed from unflagged steps only (*_clean columns).
morris_convergence <- function(tbl, top_k = 5, ci_tol = 0.2,
                               effects = c("all", "clean")) {
  effects <- match.arg(effects)
  cols <- if (effects == "all") c("mu_star", "mu_star_lo", "mu_star_hi") else
    c("mu_star_clean", "mu_star_clean_lo", "mu_star_clean_hi")
  one <- function(d) {
    d <- data.frame(p = d$parameter, mu = d[[cols[1]]], lo = d[[cols[2]]],
                    hi = d[[cols[3]]], stringsAsFactors = FALSE)
    d <- d[is.finite(d$mu) & is.finite(d$lo) & is.finite(d$hi), , drop = FALSE]
    d <- d[order(-d$mu), , drop = FALSE]
    k <- min(top_k, nrow(d) - 1)
    if (k < 1) return(data.frame(n_parameters = nrow(d), status = "not evaluated",
                                 index_converged = NA, screening_secure = NA,
                                 max_rel_ci_width = NA_real_, max_ci_width_rel_top = NA_real_,
                                 n_secure = NA_integer_, n_tied = NA_integer_,
                                 secure_top = "", tied = "", stringsAsFactors = FALSE))
    top <- d[seq_len(k), ]; rest <- d[-seq_len(k), ]
    rel_w  <- (top$hi - top$lo) / top$mu
    secure <- top$lo > max(rest$hi)
    tied   <- c(top$p[!secure], rest$p[rest$hi >= min(top$lo)])
    idx_ok <- all(is.finite(rel_w)) && max(rel_w) <= ci_tol + 1e-12
    data.frame(
      n_parameters         = nrow(d),
      status               = if (idx_ok) { if (length(tied)) "tied" else "converged" } else {
                               if (length(tied)) "unconverged" else "screening only" },
      index_converged      = idx_ok,
      screening_secure     = !length(tied),
      max_rel_ci_width     = max(rel_w),                         # index convergence metric
      max_ci_width_rel_top = max(top$hi - top$lo) / top$mu[1],   # Sarrazin et al. scaling
      n_secure             = sum(secure),
      n_tied               = length(tied),
      secure_top           = paste(top$p[secure], collapse = ", "),
      tied                 = paste(tied, collapse = ", "),
      stringsAsFactors = FALSE)
  }
  out <- lapply(split(tbl, tbl$scenario), one)
  bind_rows(out, .id = "scenario") %>%
    mutate(effects = effects, top_k = top_k, ci_tol = ci_tol, .after = scenario)
}

#------------------------------------------------------------------------#
# write_morris_outputs() ----
#------------------------------------------------------------------------#
# Used in: Scripts/4_sensitivity_animal_effects.R (section: RUN)
# Needs (set in 4_sensitivity_animal_effects.R): boot_top_k, conv_ci_tol, fig_dir, res_dir
write_morris_outputs <- function(tbl, tag, colours, labels) {
  csv <- file.path(res_dir, paste0("animal_effect_morris", tag, ".csv"))
  tbl$effect_output <- effect_output          # output screened (read by Scripts/9)
  write_csv(tbl, csv)
  p_all <- plot_morris(tbl, colours, labels)
  ggsave(file.path(fig_dir, paste0("animal_effect_morris", tag, ".png")), p_all,
         width = 13, height = 9, dpi = 150)
  p_an <- plot_morris(tbl %>% filter(is_animal), colours, labels)
  ggsave(file.path(fig_dir, paste0("animal_effect_morris", tag, "_animal.png")), p_an,
         width = 13, height = 9, dpi = 150)
  # raw elementary effects (lets the bootstrap / summaries be recomputed later)
  ee <- attr(tbl, "ee")
  if (!is.null(ee))
    write_csv(ee, file.path(res_dir, paste0("animal_effect_morris", tag, "_elementary_effects.csv")))

  # bootstrap figure: mu* with 95% CI for the top 15 parameters per scenario
  top <- tbl %>% group_by(scenario) %>% slice_min(rank, n = 15) %>% ungroup() %>%
    mutate(lab = tidytext_reorder(parameter_label, mu_star, scenario))
  p_ci <- ggplot(top, aes(mu_star, lab, colour = range_source, shape = param_type)) +
    geom_errorbar(aes(xmin = mu_star_lo, xmax = mu_star_hi), width = 0.25, orientation = "y") +
    geom_point(size = 2) +
    facet_wrap(~scenario, scales = "free", labeller = scenario_labeller) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_colour_manual(values = colours, labels = labels, name = "Range source", drop = FALSE) +
    scale_shape_manual(values = c(general = 16, `scenario-specific` = 15, animal = 17),
                       name = "Parameter type") +
    labs(x = if (identical(effect_output, "pct")) "mu* with 95% bootstrap CI (% of no-animal soil C)" else
           expression(mu*"* with "*95*"% bootstrap CI (g C m"^-2*")"), y = NULL) +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, paste0("animal_effect_morris", tag, "_bootstrap.png")), p_ci,
         width = 13, height = 9, dpi = 150)

  # convergence per scenario (index + screening convergence; ties reported),
  # on all elementary effects and on the unflagged ("clean") ones
  conv <- bind_rows(morris_convergence(tbl, boot_top_k, conv_ci_tol, "all"),
                    morris_convergence(tbl, boot_top_k, conv_ci_tol, "clean"))
  conv_csv <- file.path(res_dir, paste0("animal_effect_morris", tag, "_convergence.csv"))
  write_csv(conv, conv_csv)
  cat("\nMorris convergence (", if (nzchar(tag)) sub("^_", "", tag) else "main", "): top-",
      boot_top_k, ", index tolerance = ", conv_ci_tol, " (95% CI width / mu*)\n", sep = "")
  short <- function(x, n = 50) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
  print(as.data.frame(conv %>% transmute(scenario, effects, status,
                                         max_rel_ci_width = round(max_rel_ci_width, 3),
                                         n_secure, n_tied, tied = short(tied))),
        row.names = FALSE)
  bad <- conv %>% filter(effects == "all", status == "unconverged")
  if (nrow(bad))
    warning("Morris screening not converged (top-", boot_top_k, " CIs wider than ", conv_ci_tol,
            " x mu* and overlapping the rest) for: ", paste(bad$scenario, collapse = ", "),
            " -- increase morris_r.")
  scr <- conv %>% filter(effects == "all", status == "screening only")
  if (nrow(scr))
    message("Top-", boot_top_k, " set settled but mu* CIs wider than ", conv_ci_tol,
            " x mu* (screening converged, index not): ", paste(scr$scenario, collapse = ", "))
  tie <- conv %>% filter(effects == "all", status == "tied")
  if (nrow(tie))
    message("Converged but tied at the top-", boot_top_k, " boundary (near-equal mu*, ",
            "report as a group): ", paste(tie$scenario, collapse = ", "))
  bad_c <- conv %>% filter(effects == "clean", status == "unconverged")
  if (nrow(bad_c))
    message("Clean-step mu* (flagged regimes removed) not converged for: ",
            paste(bad_c$scenario, collapse = ", "))

  if (interactive()) print(tbl %>% group_by(scenario) %>% slice_max(mu_star, n = 5))
  cat("Wrote: ", csv, ", ", conv_csv, " and figures animal_effect_morris", tag, "[ _animal].png\n", sep = "")
}

#------------------------------------------------------------------------#
# plot_morris_paired() ----
#------------------------------------------------------------------------#
# Used in: Scripts/4_sensitivity_animal_effects.R (section: Paired Morris ranking: current knowledge vs standardized)
# Needs (set in R/sensitivity_settings.R): cv_default
# Per scenario, the top parameters with bootstrap 95% intervals for BOTH range
# definitions. mu* is scaled to the most influential parameter within each
# scenario x analysis (raw mu* depends on range width, so it is not comparable
# between analyses). "*" after a label = the bootstrap rank intervals of the two
# analyses do not overlap (rank shift beyond sampling noise). Uses this
# session's runs, or the saved CSVs when the Morris runs were switched off.
plot_morris_paired <- function(main_tbl, std_tbl, n_top = 8) {
  prep <- function(d, tag) d %>%
    transmute(scenario, parameter, parameter_label, analysis = tag,
              mu_star, mu_star_lo, mu_star_hi, rank, rank_lo, rank_hi)
  mm <- bind_rows(prep(main_tbl, "main"), prep(std_tbl, "std")) %>%
    group_by(scenario, analysis) %>%
    mutate(s_max = max(mu_star, na.rm = TRUE), rel = mu_star / s_max,
           rel_lo = mu_star_lo / s_max, rel_hi = mu_star_hi / s_max) %>%
    ungroup()
  keep <- mm %>% filter(rank <= n_top) %>% distinct(scenario, parameter)
  pd   <- mm %>% semi_join(keep, by = c("scenario", "parameter"))
  shift <- pd %>% select(scenario, parameter, analysis, rank_lo, rank_hi) %>%
    pivot_wider(names_from = analysis, values_from = c(rank_lo, rank_hi)) %>%
    mutate(star = rank_hi_main < rank_lo_std | rank_hi_std < rank_lo_main) %>%
    select(scenario, parameter, star)
  ord <- pd %>% filter(analysis == "main") %>% select(scenario, parameter, rel_main = rel)
  pd <- pd %>% left_join(ord, by = c("scenario", "parameter")) %>%
    left_join(shift, by = c("scenario", "parameter")) %>%
    mutate(label = paste0(parameter_label, ifelse(!is.na(star) & star, " *", "")),
           key   = paste(label, scenario, sep = "___"))
  pd$key <- factor(pd$key, levels = unique(pd$key[order(pd$scenario, pd$rel_main)]))

  dodge <- c(main = 0.18, std = -0.18)       # current knowledge above, standardized below
  ggplot(pd, aes(y = key, colour = analysis)) +
    lapply(names(dodge), function(a) list(
      geom_errorbar(data = filter(pd, analysis == a),
                    aes(xmin = rel_lo, xmax = rel_hi), width = 0.25,
                    orientation = "y", position = position_nudge(y = dodge[[a]])),
      geom_point(data = filter(pd, analysis == a), aes(x = rel), size = 2,
                 position = position_nudge(y = dodge[[a]])))) +
    facet_wrap(~scenario, scales = "free_y", labeller = scenario_labeller) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_colour_manual(values = c(main = "#2166ac", std = "#b2182b"),
                        labels = c(main = "Current knowledge",
                                   std = paste0("Standardized (CV = ", cv_default, ")")),
                        name = NULL) +
    labs(x = expression(mu*"* relative to the most influential parameter (95% bootstrap CI)"),
         y = NULL,
         caption = paste0("Top ", n_top, " parameters from either analysis per scenario, ordered by ",
                          "current-knowledge mu*. * = rank shift beyond bootstrap noise.")) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}
