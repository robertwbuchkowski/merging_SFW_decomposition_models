#------------------------------------------------------------------------#
# Parameter identifiability: artificial data, sensitivities, collinearity ----
#------------------------------------------------------------------------#
# Helper functions for Scripts/7_parameter_identifiability.R. Local
# identifiability analysis after Brun, Reichert & Kuensch (2001, Water Resour.
# Res. 37: 1015-1030) as implemented in the FME package (Soetaert & Petzoldt
# 2010, J. Stat. Softw. 33(3)):
#   1. The model generates ARTIFICIAL observations for three data types
#      (ident_observations()):
#        ss_noanimal  steady-state pools without the animal
#        ss_animal    steady-state pools with the animal (incl. animal biomass)
#        ts_y1..yN    pools (incl. animal biomass) once per year for N years
#                     after the animal is introduced into the no-animal soil
#   2. FME::sensFun() gives the normalized local sensitivity of every
#      observation to every parameter: S_ij = (dy_i / dp_j) * p_j / y_i, i.e.
#      the relative change in the observation per relative change in the
#      parameter (equivalent to assuming the same relative measurement error
#      for every observation).
#   3. For each data design (a combination of data types) and parameter set
#      (all parameters, or the animal parameters only with the rest assumed
#      known), FME::collin() gives the collinearity index
#      gamma = 1 / sqrt(min eigenvalue of S_n' S_n), S_n = S with unit-length
#      columns. gamma ~ how much a change in one parameter can be compensated
#      by the others; a parameter set is called identifiable when gamma is
#      below ident_threshold (Brun et al. 2001 use 10-15; 20 is also common)
#      and every parameter carries information (RMS sensitivity >= min_sens).
#
# Settings are passed explicitly in a list (see ident_settings in Scripts/7),
# so the functions also run on parallel workers. Functions that build the
# scenario objects read the same top-level objects as the Morris scripts
# (listed per function under 'Needs').

#------------------------------------------------------------------------#
# set_params() ----
#------------------------------------------------------------------------#
# Used in: called by ident_observations()
# Needs (set in Scripts/7): derive_fn
# Set several parameters at once, then rebuild derived parameters and the
# seasonal forcing ONCE (set_param() does that after every single change).
set_params <- function(obj, pv) {
  for (nm in names(pv)) obj$parms[[nm]] <- pv[[nm]]
  obj$parms                 <- derive_fn(obj$parms)
  obj$parms$climate_forcing <- make_climate_forcing(obj$parms)
  obj
}

#------------------------------------------------------------------------#
# ident_observations() ----
#------------------------------------------------------------------------#
# Used in: called by make_ident_fun(), ident_reference()
# Needs (set in R/sensitivity_settings.R): animal_pools
# The artificial data set for one parameter vector `pv` (named numeric):
#   ss_noanimal:<pool>  no-animal equilibrium (constant mean forcing)
#   ss_animal:<pool>    with-animal equilibrium (constant mean forcing)
#   ts_y<k>:<pool>      k = 1..ts_years after introducing the animal into the
#                       no-animal soil (seeded at ts_seed_mod x its input
#                       biomass). ts_forcing = "seasonal": the run starts from
#                       the no-animal SEASONAL limit cycle (Newton shooting) and
#                       is sampled on the same calendar day each year;
#                       "constant": starts from the no-animal equilibrium and
#                       uses constant mean forcing (fast; no seasonal signal).
# `ref` (from ident_reference()) supplies nominal states used as warm starts,
# so a small parameter change stays on the same solution branch.
# Stops with an error if an equilibrium or limit cycle is not reached (a
# finite-difference sensitivity from an unconverged solve is meaningless).
ident_observations <- function(pv, base_pair, settings, ref = NULL) {
  s <- settings
  pair <- base_pair
  pair$baseline  <- set_params(pair$baseline,  pv)
  pair$treatment <- set_params(pair$treatment, pv)
  why <- infeasible_soil(pair$baseline$parms)
  if (nzchar(why)) stop(why)

  b <- spinup_equilibrium(pair$baseline, warm_start = ref$ss_noanimal,
                          max_time = s$eq_max_time, stol = s$eq_stol, verbose = FALSE)
  t <- spinup_equilibrium(pair$treatment,
                          warm_start = if (is.null(ref)) b$init_state_spin else ref$ss_animal,
                          max_time = s$eq_max_time, stol = s$eq_stol, verbose = FALSE)
  if (!isTRUE(b$spin_info$converged) || !isTRUE(t$spin_info$converged))
    stop("equilibrium not reached")
  ss0 <- b$init_state_spin; ss1 <- t$init_state_spin

  # starting state of the introduction experiment (no-animal soil)
  if (s$ts_forcing == "seasonal") {
    lc <- dynamic_spinup_newton(b, from = if (is.null(ref)) ss0 else ref$lc_noanimal,
                                by = 73, tol = s$lc_tol, verbose = FALSE)
    if (!isTRUE(lc$converged))
      stop(sprintf("no-animal limit cycle not reached (residual %.1e > lc_tol %.0e)",
                   lc$residual, s$lc_tol))
    start <- lc$final_state
  } else start <- ss0

  ws      <- pair$treatment$working_state
  shared  <- intersect(names(ws), names(start))
  ws[shared] <- start[shared]
  animals <- intersect(animal_pools, names(ws))
  ws[animals] <- pair$treatment$working_state[animals] * s$ts_seed_mod
  parms <- pair$treatment$parms
  parms$climate_forcing <- if (s$ts_forcing == "seasonal") make_climate_forcing(parms) else
    make_climate_forcing_equilibrium(parms)
  out <- deSolve::ode(y = ws, times = 365 * (0:s$ts_years), func = pair$treatment$wrapped_model,
                      parms = parms, atol = 1e-9, rtol = 1e-9)
  out <- as.data.frame(out)
  if (nrow(out) < s$ts_years + 1) stop("introduction run failed")
  pools <- names(ws)
  ts <- unlist(lapply(seq_len(s$ts_years), function(k)
    setNames(as.numeric(out[k + 1, pools]), paste0("ts_y", k, ":", pools))))

  y <- c(setNames(ss0, paste0("ss_noanimal:", names(ss0))),
         setNames(ss1, paste0("ss_animal:",   names(ss1))), ts)
  if (s$ts_forcing == "seasonal") attr(y, "lc") <- start   # nominal cycle for warm starts
  y
}

#------------------------------------------------------------------------#
# ident_reference() ----
#------------------------------------------------------------------------#
# Used in: called by ident_sensitivity()
# Nominal run (all parameters at their model values): the reference
# observations and the nominal states used as warm starts by every perturbed run.
ident_reference <- function(pv, base_pair, settings) {
  y <- ident_observations(pv, base_pair, settings, ref = NULL)
  grab <- function(prefix) {
    v <- y[startsWith(names(y), prefix)]
    setNames(as.numeric(v), sub(prefix, "", names(v), fixed = TRUE))
  }
  list(y = c(y), ss_noanimal = grab("ss_noanimal:"), ss_animal = grab("ss_animal:"),
       lc_noanimal = attr(y, "lc"))
}

#------------------------------------------------------------------------#
# make_ident_fun() ----
#------------------------------------------------------------------------#
# Used in: called by ident_sensitivity()
# The model function handed to FME::sensFun(): parameter vector -> artificial
# observations. Built by a factory so its environment holds only what it needs.
# sensFun() works on MULTIPLIERS m (nominal 1; parameter = pv0 x m) rather
# than on the parameters themselves: sensFun() floors its perturbation at
# `tiny` in ABSOLUTE terms, which for small parameters (e.g. 1e-4) would be a
# huge relative step. With multipliers the step is exactly fd_rel x value, and
# since dy/dm = p dy/dp the normalized sensitivities are unchanged.
make_ident_fun <- function(base_pair, settings, ref, pv0) {
  force(base_pair); force(settings); force(ref); force(pv0)
  function(m) {
    m <- m[names(pv0)]
    y <- tryCatch(ident_observations(pv0 * m, base_pair, settings, ref),
                  error = function(e) stop("perturbed run failed (",
                    paste(names(m)[m != 1], collapse = ", "), "): ", conditionMessage(e),
                    call. = FALSE))
    attr(y, "lc") <- NULL
    y
  }
}

#------------------------------------------------------------------------#
# ident_sensitivity() ----
#------------------------------------------------------------------------#
# Used in: Scripts/7_parameter_identifiability.R (section: Sensitivities)
# Needs (set in Scripts/7): model, scen, fitted_params, unc_all, animal_unc,
#   scenspec_unc, derive_fn (via make_base_pair(), param_set_for())
# Local sensitivities of all artificial observations to all screened
# parameters of one scenario (the same parameter set as the Morris screening;
# parameters whose model value is 0 are skipped because a relative
# perturbation is undefined). Observations whose nominal value is below
# settings$obs_floor (numerically zero) or not in settings$obs_pools are dropped.
# Returns list(scenario, sens = sensFun object, ptype, animal, obs, time_s).
ident_sensitivity <- function(scenario, settings) {
  t0 <- Sys.time()
  bp    <- make_base_pair(scenario)
  parms <- bp$treatment$parms
  ps    <- param_set_for(scenario, parms)
  val   <- vapply(ps$pnames, function(p) as.numeric(parms[[p]]), numeric(1))
  pv    <- val[is.finite(val) & val != 0]
  ref   <- ident_reference(pv, bp, settings)
  y     <- ref$y
  pool  <- sub("^[^:]*:", "", names(y))
  keep  <- is.finite(y) & abs(y) > settings$obs_floor
  if (!is.null(settings$obs_pools)) keep <- keep & pool %in% settings$obs_pools
  obs   <- names(y)[keep]
  m  <- setNames(rep(1, length(pv)), names(pv))          # multipliers, see make_ident_fun()
  sF <- FME::sensFun(make_ident_fun(bp, settings, ref, pv), parms = m, sensvar = obs,
                     map = NULL, tiny = settings$fd_rel)
  attr(sF, "pars") <- pv                                  # report the actual values
  list(scenario = scenario, sens = sF, ptype = ps$ptype[names(pv)],
       animal = intersect(ps$animal, names(pv)),
       skipped_zero = setdiff(ps$pnames, names(pv)),
       obs = data.frame(observation = obs, value = as.numeric(y[obs]),
                        stringsAsFactors = FALSE),
       time_s = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

#------------------------------------------------------------------------#
# ident_rows() ----
#------------------------------------------------------------------------#
# Used in: called by ident_collinearity()
# Rows of a sensFun object that belong to the data types of one design.
# "ts" selects every ts_y<k> observation.
ident_rows <- function(var, types) {
  prefix <- sub(":.*$", "", var)
  prefix[grepl("^ts_y[0-9]+$", prefix)] <- "ts"
  prefix %in% types
}

#------------------------------------------------------------------------#
# collin_index() ----
#------------------------------------------------------------------------#
# Used in: called by ident_collinearity()
# Collinearity index of one parameter subset via FME::collin().
# sub_sens: sensFun rows of one design (columns x, var and parameters).
collin_index <- function(sub_sens, pars) {
  if (length(pars) < 2) return(1)
  # only the subset's columns go to collin(): its zero-sensitivity check fails
  # when several parameters of the full set have zero sensitivity
  suppressWarnings(FME::collin(sub_sens[, c("x", "var", pars)],
                               parset = seq_along(pars))$collinearity[1])
}

#------------------------------------------------------------------------#
# ident_collinearity() ----
#------------------------------------------------------------------------#
# Used in: called by run_identifiability()
# Identifiability of one parameter set under one data design.
#   sens        sensFun object (all observations, all parameters)
#   params      parameters to identify (the others are assumed known)
#   types       data types in the design (ss_noanimal, ss_animal, ts)
#   threshold   collinearity index below which a set is identifiable
#   min_sens    minimum RMS normalized sensitivity (FME's L2) for a parameter
#               to carry information in this design
#   exhaustive_max  test every subset when at most this many parameters carry
#               information (FME::collin with all N); otherwise greedy only
# Identifiable subset (greedy, Brun et al. 2001): parameters carrying
# information are taken in order of decreasing sensitivity and added when the
# subset's collinearity index stays <= threshold. The exhaustive search (small
# sets) gives the largest identifiable subset size and its best subset.
# Returns list(summary = one row, params = one row per parameter,
#              pairs = most collinear parameter pairs).
ident_collinearity <- function(sens, params, types, threshold = 20, min_sens = 1e-3,
                               exhaustive_max = 12, n_pairs = 10) {
  rows <- ident_rows(sens$var, types)
  sub  <- as.data.frame(sens)[rows, c("x", "var", params), drop = FALSE]
  S    <- as.matrix(sub[, params, drop = FALSE])
  L1   <- colMeans(abs(S))
  L2   <- sqrt(colMeans(S^2))
  info <- params[is.finite(L2) & L2 >= min_sens]
  info <- info[order(-L2[info])]

  # greedy identifiable subset, most sensitive parameters first
  ident <- character(0)
  for (p in info) {
    if (!length(ident) || collin_index(sub, c(ident, p)) <= threshold) ident <- c(ident, p)
  }
  gamma_all <- if (length(info) >= 2) collin_index(sub, info) else NA_real_

  # exhaustive search over all subsets of the informative parameters
  ex_n <- NA_integer_; ex_best <- ""; ex_gamma <- NA_real_
  if (length(info) >= 2 && length(info) <= exhaustive_max) {
    cl_all <- suppressWarnings(FME::collin(sub[, c("x", "var", info)], maxcomb = 1e6))
    ok <- cl_all[cl_all$collinearity <= threshold, , drop = FALSE]
    if (nrow(ok)) {
      ex_n <- max(ok$N)
      best <- ok[ok$N == ex_n, , drop = FALSE]
      best <- best[which.min(best$collinearity), , drop = FALSE]
      ex_best  <- paste(info[unlist(best[1, info]) == 1], collapse = ", ")
      ex_gamma <- best$collinearity
    } else { ex_n <- 1L; ex_best <- info[1]; ex_gamma <- 1 }
  } else if (length(info) == 1) { ex_n <- 1L; ex_best <- info; ex_gamma <- 1 }

  # most collinear pairs among the informative parameters
  pairs <- NULL
  if (length(info) >= 2) {
    cp <- suppressWarnings(FME::collin(sub[, c("x", "var", info)], N = 2, maxcomb = 1e6))
    cp <- cp[order(-cp$collinearity), , drop = FALSE][seq_len(min(n_pairs, nrow(cp))), , drop = FALSE]
    pairs <- data.frame(
      param_1 = apply(cp[, info, drop = FALSE], 1, function(r) info[r == 1][1]),
      param_2 = apply(cp[, info, drop = FALSE], 1, function(r) info[r == 1][2]),
      collinearity = cp$collinearity, stringsAsFactors = FALSE)
  }

  status <- ifelse(!params %in% info, "no information",
                   ifelse(params %in% ident, "identifiable", "collinear"))
  list(
    summary = data.frame(
      n_obs = nrow(S), n_params = length(params), n_informative = length(info),
      collinearity_all = gamma_all,
      all_identifiable = length(info) == length(params) &&
                         (length(info) < 2 || gamma_all <= threshold),
      n_identifiable_greedy = length(ident),
      collinearity_greedy = if (length(ident) >= 2) collin_index(sub, ident) else 1,
      n_identifiable_exhaustive = ex_n, collinearity_exhaustive = ex_gamma,
      identifiable_greedy = paste(ident, collapse = ", "),
      identifiable_exhaustive = ex_best,
      no_information = paste(setdiff(params, info), collapse = ", "),
      stringsAsFactors = FALSE),
    params = data.frame(parameter = params, L1 = L1[params], L2 = L2[params],
                        status = status, greedy_order = match(params, ident),
                        stringsAsFactors = FALSE, row.names = NULL),
    pairs = pairs)
}

#------------------------------------------------------------------------#
# run_identifiability() ----
#------------------------------------------------------------------------#
# Used in: Scripts/7_parameter_identifiability.R (section: Collinearity)
# Every scenario x parameter set x data design. `sens_list` is the list of
# ident_sensitivity() results; `designs` a named list of data-type vectors;
# `param_sets` a subset of c("all", "animal").
# Returns list(summary, params, pairs) of tidy data frames.
run_identifiability <- function(sens_list, designs, param_sets = c("all", "animal"),
                                threshold = 20, min_sens = 1e-3, exhaustive_max = 12) {
  summ <- prm <- prs <- list()
  for (sr in sens_list) {
    for (set in param_sets) {
      params <- if (set == "all") names(sr$ptype) else sr$animal
      if (!length(params)) next
      for (d in names(designs)) {
        r <- ident_collinearity(sr$sens, params, designs[[d]], threshold = threshold,
                                min_sens = min_sens, exhaustive_max = exhaustive_max)
        key <- data.frame(scenario = sr$scenario, param_set = set, design = d,
                          data_types = paste(designs[[d]], collapse = " + "),
                          stringsAsFactors = FALSE)
        summ[[length(summ) + 1]] <- cbind(key, r$summary)
        prm[[length(prm) + 1]]   <- cbind(key[rep(1, nrow(r$params)), ], r$params,
                                          param_type = sr$ptype[r$params$parameter],
                                          row.names = NULL)
        if (!is.null(r$pairs))
          prs[[length(prs) + 1]] <- cbind(key[rep(1, nrow(r$pairs)), ], r$pairs, row.names = NULL)
      }
    }
  }
  list(summary = bind_rows(summ), params = bind_rows(prm), pairs = bind_rows(prs))
}

#------------------------------------------------------------------------#
# plot_ident_summary() ----
#------------------------------------------------------------------------#
# Used in: Scripts/7_parameter_identifiability.R (section: Figures)
# Number of parameters that carry information and that are identifiable
# (greedy subset) out of the total, per data design; rows = parameter set,
# columns = scenario.
plot_ident_summary <- function(summ, design_levels) {
  d <- summ %>%
    transmute(scenario, param_set, design = factor(design, levels = rev(design_levels)),
              `Total` = n_params, `Carry information` = n_informative,
              `Identifiable` = n_identifiable_greedy) %>%
    pivot_longer(c(Total, `Carry information`, Identifiable), names_to = "what",
                 values_to = "n") %>%
    mutate(what = factor(what, levels = c("Total", "Carry information", "Identifiable")),
           param_set = factor(ifelse(param_set == "all", "All parameters",
                                     "Animal parameters only"),
                              levels = c("All parameters", "Animal parameters only")))
  ggplot(d, aes(n, design, fill = what)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    facet_grid(param_set ~ scenario, scales = "free_x",
               labeller = labeller(scenario = scenario_labeller)) +
    scale_fill_manual(values = c(Total = "grey80", `Carry information` = "#92c5de",
                                 Identifiable = "#2166ac"), name = NULL) +
    labs(x = "Number of parameters", y = "Artificial data used") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank())
}

#------------------------------------------------------------------------#
# plot_ident_params() ----
#------------------------------------------------------------------------#
# Used in: Scripts/7_parameter_identifiability.R (section: Figures)
# Status of every parameter (identifiable / collinear / no information) under
# every data design, one panel per scenario, most sensitive parameters on top.
# The number in an identifiable cell is the order in which the greedy search
# added the parameter (1 = most sensitive in that design).
plot_ident_params <- function(prm, design_levels, set = "animal") {
  d <- prm %>% filter(param_set == set) %>%
    group_by(scenario, parameter) %>% mutate(L2_max = max(L2, na.rm = TRUE)) %>% ungroup() %>%
    mutate(design = factor(design, levels = design_levels),
           label  = pretty_param(parameter),
           key    = paste(label, scenario, sep = "___"),
           status = factor(status, levels = c("identifiable", "collinear", "no information")))
  # most sensitive parameters at the top of each panel
  d$key <- factor(d$key, levels = unique(d$key[order(d$scenario, d$L2_max)]))
  ggplot(d, aes(design, key, fill = status)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(is.na(greedy_order), "", greedy_order)), size = 2.2,
              colour = "white") +
    facet_wrap(~scenario, scales = "free_y", labeller = scenario_labeller) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_fill_manual(values = c(identifiable = "#2166ac", collinear = "#f4a582",
                                 `no information` = "grey85"), name = NULL, drop = FALSE) +
    labs(x = "Artificial data used", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom",
          panel.grid = element_blank())
}
