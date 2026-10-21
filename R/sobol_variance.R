#------------------------------------------------------------------------#
# Variance-based sensitivity (Sobol') after Morris screening ----
#------------------------------------------------------------------------#
# Helper functions for Scripts/9_sobol_variance_partition.R. Two-stage global
# sensitivity analysis (Saltelli et al. 2008; Campolongo et al. 2007):
#   1. SCREENING (Morris, Scripts/4): per scenario, keep the parameters whose
#      mu* is at least screen_rel_cut x the largest mu* in the MAIN or the
#      STANDARDIZED Morris analysis (or whose bootstrap CI reaches that cut,
#      i.e. tied with it), plus every animal parameter (morris_shortlist()).
#   2. VARIANCE PARTITION (Sobol'): first-order (S) and total-effect (ST)
#      indices of the shortlisted parameters for the animal effect on total
#      soil + root C -- as a PERCENTAGE of the no-animal soil + root C
#      (primary, "pct") and in g C m-2 ("abs"), from the same runs. The
#      absolute effect scales with the size of the baseline soil, which varies
#      by orders of magnitude over the parameter distributions (runaway
#      carbon at implausible parameter combinations), so its variance is
#      dominated by a few extreme runs; the relative effect is far better
#      behaved and is the primary output. A third output, "rank", is the
#      rank-transformed % effect: indices robust to the remaining heavy tails.
#      Parameters drawn from the MAIN (knowledge-based)
#      distributions of the uncertainty propagation (Scripts/5). All other
#      parameters vary TOGETHER as one group ("all other parameters"): its total
#      index is the share of output variance that would be lost by fixing them,
#      i.e. a direct check of the screening (screening_ok if ST_rest is below
#      rest_st_max).
# Estimators (Saltelli et al. 2010, Comput. Phys. Commun. 181: 259-270):
#   two independent sample matrices A and B (N x k) and, for every factor j, A
#   with the columns of factor j taken from B (A_B^j); N (k + 2) model runs.
#     S_j  = mean(f(B) (f(A_B^j) - f(A))) / V                (Saltelli 2010)
#     ST_j = mean((f(A) - f(A_B^j))^2) / (2 V)                 (Jansen 1999)
#   with V the variance of f over A and B; 95% CIs by bootstrapping the N rows.
# The REALISM FILTER (Scripts/5) is NOT applied to the Sobol' design (Sobol'
# indices assume independent inputs; conditioning on a realistic baseline
# would make them dependent). It is reported separately: which samples have
# a realistic baseline and how the effect differs among them (option A).
# PAWN (Pianosi & Wagener 2015, Environ. Model. Softw. 67: 1-11), a
# distribution-based index that does not rely on the variance (robust to the
# heavy tails caused by regime shifts), is computed from the same A and B
# runs for EVERY parameter, as a cross-check of both stages.
#
# Settings are read from the script (listed per function under 'Needs').

#------------------------------------------------------------------------#
# morris_shortlist() ----
#------------------------------------------------------------------------#
# Used in: Scripts/9_sobol_variance_partition.R (section: Screening)
# Needs (set in Scripts/9): screen_rel_cut, include_all_animal, max_screened
# Shortlist for one scenario from the Morris tables (main and, if available,
# standardized; NULL is allowed). A parameter is kept when, in either
# analysis, mu* / max(mu*) >= screen_rel_cut ("main" / "standardized"), or
# its upper 95% bootstrap bound reaches the cut ("tie"), or it is an animal
# parameter (include_all_animal). At most max_screened NON-animal parameters
# are kept (largest relative mu* first); the rest join the group of fixed
# parameters, whose influence the Sobol' group index then measures.
# Returns one row per parameter with the relative mu* and the reason.
morris_shortlist <- function(main_tbl, std_tbl, scenario) {
  rel_of <- function(d, tag) {
    if (is.null(d)) return(NULL)
    d <- d[d$scenario == scenario & is.finite(d$mu_star), , drop = FALSE]
    if (!nrow(d)) return(NULL)
    m <- max(d$mu_star)
    data.frame(parameter = d$parameter, is_animal = d$is_animal,
               rel = d$mu_star / m, rel_hi = d$mu_star_hi / m, analysis = tag,
               stringsAsFactors = FALSE)
  }
  r <- rbind(rel_of(main_tbl, "main"), rel_of(std_tbl, "standardized"))
  if (is.null(r)) stop("No Morris results for ", scenario)
  pars <- unique(r$parameter)
  get <- function(p, tag, col) { v <- r[[col]][r$parameter == p & r$analysis == tag]
                                 if (length(v)) v[1] else NA_real_ }
  out <- data.frame(
    parameter = pars,
    is_animal = vapply(pars, function(p) any(r$is_animal[r$parameter == p]), logical(1)),
    rel_main  = vapply(pars, get, numeric(1), tag = "main", col = "rel"),
    rel_std   = vapply(pars, get, numeric(1), tag = "standardized", col = "rel"),
    relhi_main = vapply(pars, get, numeric(1), tag = "main", col = "rel_hi"),
    relhi_std  = vapply(pars, get, numeric(1), tag = "standardized", col = "rel_hi"),
    stringsAsFactors = FALSE, row.names = NULL)
  ge <- function(x) !is.na(x) & x >= screen_rel_cut
  out$reason <- ifelse(ge(out$rel_main) & ge(out$rel_std), "main + standardized",
                ifelse(ge(out$rel_main), "main",
                ifelse(ge(out$rel_std), "standardized",
                ifelse(ge(out$relhi_main) | ge(out$relhi_std), "tie (CI reaches cut)", ""))))
  out$rel_max  <- pmax(out$rel_main, out$rel_std, na.rm = TRUE)
  out$selected <- nzchar(out$reason)
  # cap the number of non-animal parameters
  na_sel <- which(out$selected & !out$is_animal)
  if (length(na_sel) > max_screened) {
    drop <- na_sel[order(-out$rel_max[na_sel])][-seq_len(max_screened)]
    out$selected[drop] <- FALSE
    out$reason[drop]   <- paste0(out$reason[drop], " (over max_screened: in fixed group)")
  }
  if (include_all_animal) {
    a <- out$is_animal & !out$selected
    out$selected[a] <- TRUE
    out$reason[a]   <- "animal parameter"
  }
  out$scenario <- scenario
  out[order(!out$selected, -out$rel_max), c("scenario", setdiff(names(out), "scenario"))]
}

#------------------------------------------------------------------------#
# soil_values_feasible() ----
#------------------------------------------------------------------------#
# Used in: called by sobol_design()
# Needs (set in R/sensitivity_settings.R): soil_feasibility_params
# TRUE when a full parameter vector `pv` describes a physically possible soil
# (only the soil parameters matter; checked without solving the model).
soil_values_feasible <- function(pv, base_parms) {
  sp <- intersect(names(pv), soil_feasibility_params)
  if (!length(sp)) return(TRUE)
  p <- base_parms
  for (nm in sp) p[[nm]] <- pv[[nm]]
  !nzchar(infeasible_soil(apply_relative_moisture(derive_fn(p))))
}

#------------------------------------------------------------------------#
# sobol_design() ----
#------------------------------------------------------------------------#
# Used in: called by run_sobol_scenario()
# The Saltelli design in PARAMETER units.
#   dists    named list of sampling distributions (up_distribution())
#   factors  named list: factor -> parameter names (single parameters, plus
#            the group of all other parameters)
#   N        number of base rows
# A and B are Latin-hypercube samples of the unit cube (independent), mapped
# through the inverse CDFs of the distributions. Rows for which ANY of A, B or
# A_B^j would be an impossible soil are redrawn BEFORE any simulation (only
# BD, MAtheta and the moisture amplitude decide this), so every run is a
# physically possible soil; the count is returned as n_redrawn.
# Returns list(A, B, AB = list per factor, n_redrawn).
sobol_design <- function(dists, factors, N, base_parms, seed) {
  set.seed(seed)
  pn <- names(dists); k <- length(pn)
  to_par <- function(U) {
    X <- vapply(seq_len(k), function(j) up_quantile(U[, j], dists[[pn[j]]]), numeric(nrow(U)))
    if (nrow(U) == 1) X <- matrix(X, nrow = 1)
    colnames(X) <- pn; X
  }
  A <- B <- matrix(numeric(0), 0, k, dimnames = list(NULL, pn))
  n_drawn <- 0L
  while (nrow(A) < N) {
    need <- N - nrow(A)
    U  <- lhs_unit(max(need, 10), 2 * k)
    Xa <- to_par(U[, seq_len(k), drop = FALSE]); Xb <- to_par(U[, k + seq_len(k), drop = FALSE])
    n_drawn <- n_drawn + nrow(Xa)
    ok <- vapply(seq_len(nrow(Xa)), function(i) {
      a <- Xa[i, ]; b <- Xb[i, ]
      if (!soil_values_feasible(a, base_parms) || !soil_values_feasible(b, base_parms)) return(FALSE)
      for (f in factors) { ab <- a; ab[f] <- b[f]; if (!soil_values_feasible(ab, base_parms)) return(FALSE) }
      TRUE
    }, logical(1))
    A <- rbind(A, Xa[ok, , drop = FALSE]); B <- rbind(B, Xb[ok, , drop = FALSE])
    if (n_drawn > 50 * N) stop("sobol_design: fewer than 2% of draws are possible soils")
  }
  A <- A[seq_len(N), , drop = FALSE]; B <- B[seq_len(N), , drop = FALSE]
  AB <- lapply(factors, function(f) { m <- A; m[, f] <- B[, f]; m })
  list(A = A, B = B, AB = AB, n_redrawn = n_drawn - N)
}

#------------------------------------------------------------------------#
# sobol_model_output() / make_sobol_fun() ----
#------------------------------------------------------------------------#
# Used in: called by run_sobol_scenario() (via make_sobol_fun on the workers)
# Needs (set in R/sensitivity_settings.R): animal_pools, forcing_mode
# One model run: the steady-state TOTAL animal effect on total soil + root C
# (treatment - baseline, as in Scripts/4 and 5) and the baseline values the
# realism filter needs. effect is NA when either steady state is not reached.
sobol_model_output <- function(base_pair, pv, pools) {
  pair <- base_pair
  for (nm in names(pv)) {
    pair$treatment <- set_param(pair$treatment, nm, pv[[nm]])
    pair$baseline  <- set_param(pair$baseline,  nm, pv[[nm]])
  }
  out <- c(effect = NA_real_, converged = 0, microbe_collapse = 0, kb_floor = 0,
           infeasible = 0, base_min_pool = NA_real_,
           setNames(rep(NA_real_, length(pools)), paste0("base_", pools)))
  if (nzchar(infeasible_soil(pair$baseline$parms))) { out[["infeasible"]] <- 1; return(out) }
  tryCatch({
    b  <- steady_state(pair$baseline, ref = pair$ref_b)
    tr <- steady_state(pair$treatment, ref = pair$ref_t, warm_start = b$ss_const)
    eb <- b$ss_mean
    soil <- setdiff(intersect(names(eb), names(tr$ss_mean)), animal_pools)
    conv <- isTRUE(b$spin_info$converged) && isTRUE(tr$spin_info$converged)
    rg   <- regime_flags(eb, tr$ss_mean, pair$treatment$parms)
    c(effect = if (conv) sum(tr$ss_mean[soil]) - sum(eb[soil]) else NA_real_,
      converged = as.numeric(conv),
      microbe_collapse = as.numeric(grepl("microbe_collapse", rg)),
      kb_floor = as.numeric(grepl("kb_floor", rg)), infeasible = 0,
      base_min_pool = min(eb[soil]), baseline_values(eb, pools))
  }, error = function(e) out)
}

make_sobol_fun <- function(base_pair, pools) {
  force(base_pair); force(pools)
  function(pv) sobol_model_output(base_pair, pv, pools)
}

#------------------------------------------------------------------------#
# sobol_indices() ----
#------------------------------------------------------------------------#
# Used in: called by run_sobol_scenario()
# First-order and total-effect indices with bootstrap percentile CIs.
#   yA, yB  model output for A and B (length N); yAB  N x (number of factors)
# Rows where any of f(A), f(B), f(A_B^j) is missing (steady state not
# reached) are dropped from ALL estimators (complete cases), so every index
# uses the same rows; their number is returned (n_used).
sobol_indices <- function(yA, yB, yAB, B_boot = 1000, conf = 0.95, seed = NULL) {
  ok <- is.finite(yA) & is.finite(yB) & apply(is.finite(yAB), 1, all)
  yA <- yA[ok]; yB <- yB[ok]; yAB <- yAB[ok, , drop = FALSE]
  n  <- length(yA)
  est <- function(i) {
    a <- yA[i]; b <- yB[i]; ab <- yAB[i, , drop = FALSE]
    V <- stats::var(c(a, b))
    rbind(S  = colMeans(b * (ab - a)) / V,
          ST = colMeans((a - ab)^2) / (2 * V))
  }
  if (n < 10) return(NULL)
  pt <- est(seq_len(n))
  if (!is.null(seed)) set.seed(seed)
  bs <- replicate(B_boot, est(sample.int(n, n, replace = TRUE)))   # 2 x k x B
  al <- (1 - conf) / 2
  q  <- function(r, p) apply(bs[r, , , drop = FALSE], 2, stats::quantile, probs = p, na.rm = TRUE)
  data.frame(factor = colnames(yAB),
             S = pt["S", ], S_lo = q("S", al), S_hi = q("S", 1 - al),
             ST = pt["ST", ], ST_lo = q("ST", al), ST_hi = q("ST", 1 - al),
             n_used = n, variance = stats::var(c(yA, yB)),
             stringsAsFactors = FALSE, row.names = NULL)
}

#------------------------------------------------------------------------#
# pawn_indices() ----
#------------------------------------------------------------------------#
# Used in: called by run_sobol_scenario()
# PAWN, given-data version (Pianosi & Wagener 2018, Environ. Model. Softw.
# 108: 197-207): for each parameter, the output sample is split into n_cond
# equal-count classes of that parameter, and the Kolmogorov-Smirnov distance
# between each conditional and the unconditional output distribution is
# summarised by its median (and max). Values above ks_crit (two-sample KS
# critical value, alpha = 0.05) indicate an influence beyond sampling noise.
#   X  runs x parameters (all parameters); y  output (NA rows dropped)
pawn_indices <- function(X, y, n_cond = 10) {
  ok <- is.finite(y); X <- X[ok, , drop = FALSE]; y <- y[ok]
  n  <- length(y)
  Fy <- stats::ecdf(y); grid <- sort(unique(y))
  rows <- lapply(colnames(X), function(p) {
    br  <- unique(stats::quantile(X[, p], probs = seq(0, 1, length.out = n_cond + 1)))
    cls <- cut(X[, p], br, include.lowest = TRUE)
    ks  <- vapply(split(y, cls), function(yc)
      if (length(yc) < 5) NA_real_ else max(abs(stats::ecdf(yc)(grid) - Fy(grid))), numeric(1))
    m <- n / n_cond
    data.frame(parameter = p, pawn_median = stats::median(ks, na.rm = TRUE),
               pawn_max = max(ks, na.rm = TRUE),
               ks_crit = 1.36 * sqrt((n + m) / (n * m)), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#------------------------------------------------------------------------#
# run_sobol_scenario() ----
#------------------------------------------------------------------------#
# Used in: Scripts/9_sobol_variance_partition.R (section: Sobol' analysis)
# Needs (set in Scripts/9): sobol_screen, sobol_N, sobol_seed, sobol_boot_B, rest_st_max,
#   pawn_n_cond, cl, morris_main, morris_std; (R/sensitivity_settings.R)
#   up_filter_pools, up_pool_floor, ... (realism filter)
# Screening, Sobol' design, model runs, indices, PAWN and the realism report
# for one scenario. Returns a list of tables.
run_sobol_scenario <- function(scenario, si = 1) {
  base_pair  <- make_base_pair(scenario)
  parms_here <- base_pair$treatment$parms
  ps         <- param_set_for(scenario, parms_here)

  # sampling distributions (main analysis, same as Scripts/5)
  dists <- list()
  for (p in ps$pnames) {
    d0 <- parms_here[[p]]
    if (!is.finite(d0) || d0 == 0) next
    ds <- up_distribution(scenario, p, d0)
    if (is.finite(ds$lo) && is.finite(ds$hi) && ds$hi > ds$lo) dists[[p]] <- ds
  }
  pn <- names(dists)

  # screening
  sl <- morris_shortlist(morris_main, morris_std, scenario)
  sel  <- if (sobol_screen) intersect(sl$parameter[sl$selected], pn) else pn   # FALSE: all
  if (!sobol_screen) sl$selected <- sl$parameter %in% pn
  rest <- setdiff(pn, sel)
  factors <- c(setNames(as.list(sel), sel),
               if (length(rest)) list(all_other_parameters = rest))
  cat(sprintf("Sobol' [%s]: %d shortlisted (%d animal) + group of %d other parameters; N = %d -> %d runs\n",
              scenario, length(sel), sum(sel %in% ps$animal), length(rest), sobol_N,
              sobol_N * (length(factors) + 2)))

  # design and runs
  des <- sobol_design(dists, factors, sobol_N, parms_here, seed = sobol_seed + si)
  if (des$n_redrawn > 0)
    cat(sprintf("  %d design rows redrawn before running (impossible soil in A, B or A_B^j)\n",
                des$n_redrawn))
  user_targets <- read_user_targets()
  tg    <- targets_for(scenario, base_pair, user_targets)
  pools <- union("TotalC", tg$pool)          # TotalC needed for the relative effect
  all_X <- rbind(des$A, des$B, do.call(rbind, des$AB))
  rows  <- lapply(seq_len(nrow(all_X)), function(i) all_X[i, ])
  fun   <- make_sobol_fun(base_pair, pools)
  t0 <- Sys.time()
  Y  <- if (is.null(cl)) lapply(rows, fun) else parallel::parLapplyLB(cl, rows, fun)
  Y  <- do.call(rbind, Y)
  cat(sprintf("  %d runs in %.1f min\n", nrow(Y), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  N  <- sobol_N; nf <- length(factors)
  Y  <- cbind(Y, effect_pct = 100 * Y[, "effect"] / Y[, "base_TotalC"])
  # rank-transformed % effect: Sobol' indices of the ranks are robust to the
  # heavy tails (they measure influence on the ORDER of the outcomes)
  Y  <- cbind(Y, effect_rank = rank(Y[, "effect_pct"], na.last = "keep"))
  out_col <- c(pct = "effect_pct", rank = "effect_rank", abs = "effect")

  # Sobol' indices and PAWN for each output (same runs)
  si_tab <- pawn <- NULL
  for (o in names(out_col)) {
    y  <- Y[, out_col[[o]]]
    yA <- y[seq_len(N)]; yB <- y[N + seq_len(N)]
    yAB <- matrix(y[2 * N + seq_len(N * nf)], N, nf, dimnames = list(NULL, names(factors)))
    st <- sobol_indices(yA, yB, yAB, B_boot = sobol_boot_B, seed = sobol_seed + 100 + si)
    if (is.null(st)) stop("Sobol' [", scenario, "]: fewer than 10 complete rows")
    st$output <- o
    si_tab <- rbind(si_tab, st)
    pw <- pawn_indices(rbind(des$A, des$B), c(yA, yB), n_cond = pawn_n_cond)
    pw$output <- o
    pawn <- rbind(pawn, pw)
  }
  si_tab$scenario  <- scenario
  si_tab$type      <- ifelse(si_tab$factor == "all_other_parameters", "group", "parameter")
  si_tab$is_animal <- si_tab$factor %in% ps$animal
  si_tab$label     <- ifelse(si_tab$type == "group",
                             sprintf("All other parameters (%d, fixed by screening)", length(rest)),
                             pretty_param(si_tab$factor))
  si_tab$n_params  <- vapply(factors[si_tab$factor], length, integer(1))
  pawn$scenario    <- scenario
  pawn$shortlisted <- pawn$parameter %in% sel
  pawn$label       <- pretty_param(pawn$parameter)

  # realism filter, reported separately (option A): A and B rows
  YAB <- Y[seq_len(2 * N), , drop = FALSE]
  pass_pool <- vapply(seq_len(nrow(tg)), function(j) {
    v <- YAB[, paste0("base_", tg$pool[j])]
    is.finite(v) & v >= tg$lower[j] & v <= tg$upper[j]
  }, logical(2 * N))
  if (2 * N == 1) pass_pool <- matrix(pass_pool, nrow = 1)
  realistic <- YAB[, "converged"] == 1 & apply(pass_pool, 1, all) &
               is.finite(YAB[, "base_min_pool"]) & YAB[, "base_min_pool"] > up_pool_floor
  q <- function(v, p) if (sum(is.finite(v))) unname(stats::quantile(v, p, na.rm = TRUE)) else NA_real_
  kurt <- function(v) { v <- v[is.finite(v)]; mean((v - mean(v))^4) / stats::var(v)^2 }

  summ <- do.call(rbind, lapply(names(out_col), function(o) {
  eff <- YAB[, out_col[[o]]]
  sio <- si_tab[si_tab$output == o, , drop = FALSE]
  rest_row <- sio[sio$type == "group", , drop = FALSE]
  data.frame(
    scenario = scenario, output = o, forcing = forcing_mode, N = N, n_factors = nf,
    n_shortlisted = length(sel), n_animal = sum(sel %in% ps$animal), n_fixed_group = length(rest),
    n_runs = nrow(Y), n_design_redrawn = des$n_redrawn,
    n_not_converged = sum(Y[, "converged"] == 0 & Y[, "infeasible"] == 0),
    n_microbe_collapse = sum(Y[, "microbe_collapse"] == 1),
    n_kb_floor = sum(Y[, "kb_floor"] == 1),
    n_rows_used = sio$n_used[1], output_variance = sio$variance[1],
    output_kurtosis = kurt(eff),
    sum_S = sum(sio$S), interaction_share = 1 - sum(sio$S),
    ST_rest = if (nrow(rest_row)) rest_row$ST else NA_real_,
    ST_rest_hi = if (nrow(rest_row)) rest_row$ST_hi else NA_real_,
    screening_ok = if (nrow(rest_row)) rest_row$ST_hi <= rest_st_max else NA,   # NA: no screening
    effect_median_all = q(eff, 0.5), effect_q025_all = q(eff, 0.025), effect_q975_all = q(eff, 0.975),
    n_realistic = sum(realistic), realistic_rate = mean(realistic),
    effect_median_realistic = q(eff[realistic], 0.5),
    effect_q025_realistic = q(eff[realistic], 0.025),
    effect_q975_realistic = q(eff[realistic], 0.975),
    stringsAsFactors = FALSE)
  }))

  samples <- data.frame(scenario = scenario,
                        block = rep(c("A", "B", paste0("AB_", names(factors))), each = N),
                        row = rep(seq_len(N), nf + 2), Y, check.names = FALSE)
  samples$realistic <- NA
  samples$realistic[seq_len(2 * N)] <- realistic

  list(shortlist = sl, indices = si_tab, pawn = pawn, summary = summ,
       samples = samples, design = des, targets = tg)
}

#------------------------------------------------------------------------#
# plot_sobol_indices() ----
#------------------------------------------------------------------------#
# Used in: Scripts/9_sobol_variance_partition.R (section: Figures)
# First-order (S, dark) and total-effect (ST, light) indices with 95%
# bootstrap CIs per factor (the n_top largest ST), one panel per scenario, ordered by ST. The group
# of parameters fixed by the screening is drawn in grey: a small ST there
# means the screening lost little.
plot_sobol_indices <- function(ind, summ, output_use = "pct", n_top = 20) {
  ind  <- ind[ind$output == output_use, , drop = FALSE]
  summ <- summ[summ$output == output_use, , drop = FALSE]
  # the n_top factors with the largest ST per scenario (the fixed group always shown)
  ind <- ind %>% group_by(scenario) %>%
    filter(type == "group" | rank(-ST, ties.method = "first") <= n_top) %>% ungroup()
  d <- ind %>%
    mutate(key = paste(label, scenario, sep = "___")) %>%
    group_by(scenario) %>% arrange(ST, .by_group = TRUE) %>% ungroup()
  d$key <- factor(d$key, levels = unique(d$key))
  long <- bind_rows(
    d %>% transmute(scenario, key, type, is_animal, index = "Total effect (ST)", value = ST, lo = ST_lo, hi = ST_hi),
    d %>% transmute(scenario, key, type, is_animal, index = "First order (S)", value = S, lo = S_lo, hi = S_hi)) %>%
    mutate(index = factor(index, levels = c("Total effect (ST)", "First order (S)")),
           what = ifelse(type == "group", "group", ifelse(is_animal, "animal", "other")))
  lab <- summ %>% transmute(scenario, txt = if_else(is.na(ST_rest_hi),
    sprintf("sum S = %.2f; all parameters analysed (no screening)", sum_S),
    sprintf("sum S = %.2f; screening %s (ST rest <= %.2f)", sum_S,
            ifelse(screening_ok, "OK", "NOT OK"), ST_rest_hi)))
  dodge <- position_dodge(width = 0.7)
  ggplot(long, aes(value, key, fill = interaction(what, index), group = index)) +
    geom_col(position = dodge, width = 0.65) +
    geom_errorbar(aes(xmin = lo, xmax = hi), position = dodge, width = 0.25,
                  orientation = "y", colour = "grey30", linewidth = 0.3) +
    geom_text(data = lab, aes(x = Inf, y = -Inf, label = txt), inherit.aes = FALSE,
              hjust = 1.02, vjust = -0.6, size = 2.6, colour = "grey30") +
    facet_wrap(~scenario, scales = "free_y", labeller = scenario_labeller) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_fill_manual(values = c(`other.Total effect (ST)` = "#9ecae1", `other.First order (S)` = "#2171b5",
                                 `animal.Total effect (ST)` = "#fdae6b", `animal.First order (S)` = "#d94801",
                                 `group.Total effect (ST)` = "grey80", `group.First order (S)` = "grey45"),
                      breaks = c("other.First order (S)", "other.Total effect (ST)",
                                 "animal.First order (S)", "animal.Total effect (ST)",
                                 "group.First order (S)", "group.Total effect (ST)"),
                      labels = c("Soil/plant parameter: S", "Soil/plant parameter: ST",
                                 "Animal parameter: S", "Animal parameter: ST",
                                 "Fixed group: S", "Fixed group: ST"), name = NULL) +
    coord_cartesian(xlim = c(min(0, min(long$lo, na.rm = TRUE)), NA)) +
    labs(x = "Share of output variance (Sobol' index, 95% bootstrap CI)", y = NULL,
         title = paste0("Variance of the animal effect on total soil + root C",
                        switch(output_use, pct = " (% of the no-animal value)",
                               rank = " (ranks of the % effect)", abs = " (g C m-2)")),
         subtitle = "Dark = first order (alone); light = total (incl. interactions)") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold")) +
    guides(fill = guide_legend(nrow = 2))
}

#------------------------------------------------------------------------#
# plot_pawn() ----
#------------------------------------------------------------------------#
# Used in: Scripts/9_sobol_variance_partition.R (section: Figures)
# PAWN (median KS) for every parameter, shortlisted ones highlighted; the
# dashed line is the KS critical value (influence beyond sampling noise).
# A non-shortlisted parameter far above the line would question the screening.
plot_pawn <- function(pw, n_top = 20, output_use = "pct") {
  d <- pw %>% filter(output == output_use) %>% group_by(scenario) %>% slice_max(pawn_median, n = n_top, with_ties = FALSE) %>%
    ungroup() %>% mutate(key = paste(label, scenario, sep = "___"))
  d$key <- factor(d$key, levels = unique(d$key[order(d$scenario, d$pawn_median)]))
  crit <- d %>% distinct(scenario, ks_crit)
  ggplot(d, aes(pawn_median, key, colour = shortlisted)) +
    geom_vline(data = crit, aes(xintercept = ks_crit), linetype = "dashed", colour = "grey50") +
    geom_segment(aes(x = 0, xend = pawn_median, yend = key), linewidth = 0.6) +
    geom_point(size = 2) +
    facet_wrap(~scenario, scales = "free_y", labeller = scenario_labeller) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_colour_manual(values = c(`TRUE` = "#2171b5", `FALSE` = "grey55"),
                        labels = c(`TRUE` = "Shortlisted (Sobol')", `FALSE` = "Fixed by screening"),
                        name = NULL) +
    labs(x = "PAWN index (median KS distance)", y = NULL,
         title = "Distribution-based sensitivity (PAWN), all parameters",
         subtitle = paste0("Top ", n_top, " per scenario; dashed = KS critical value (alpha = 0.05)")) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))
}

#------------------------------------------------------------------------#
# plot_sobol_realism() ----
#------------------------------------------------------------------------#
# Used in: Scripts/9_sobol_variance_partition.R (section: Figures)
# The animal effect over all A and B samples vs the samples with a realistic
# baseline (realism filter reported separately, option A).
plot_sobol_realism <- function(samples, summ, output_use = "pct") {
  summ <- summ[summ$output == output_use, , drop = FALSE]
  samples$value <- samples[[if (output_use == "abs") "effect" else "effect_pct"]]
  d <- samples %>% filter(block %in% c("A", "B"), is.finite(value))
  d <- bind_rows(d %>% mutate(set = "All samples (Sobol' design)"),
                 d %>% filter(realistic %in% TRUE) %>% mutate(set = "Realistic baseline only"))
  lab <- summ %>% transmute(scenario, txt = sprintf("%.0f%% realistic", 100 * realistic_rate))
  ggplot(d, aes(set, value, fill = set)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_violin(colour = NA, alpha = 0.6, scale = "width") +
    stat_summary(fun = median, geom = "point", shape = 21, size = 2.2, fill = "white") +
    geom_text(data = lab, aes(x = 1.5, y = Inf, label = txt), inherit.aes = FALSE,
              vjust = 1.5, size = 2.8, colour = "grey30") +
    facet_wrap(~scenario, scales = "free_y", labeller = scenario_labeller) +
    scale_fill_manual(values = c("#9ecae1", "#2171b5"), name = NULL) +
    scale_x_discrete(labels = c("All", "Realistic")) +
    labs(x = NULL, y = if (output_use != "abs") "Animal effect on total soil + root C (% of no-animal value)" else
           expression("Animal effect on total soil + root C (g C m"^-2*")"),
         title = "Effect of the realism filter (reported separately from the Sobol' indices)") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
}
