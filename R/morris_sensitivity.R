#------------------------------------------------------------------------#
# MORRIS ELEMENTARY-EFFECTS SENSITIVITY ----
#------------------------------------------------------------------------#
# MORRIS ELEMENTARY-EFFECTS SENSITIVITY (method of Morris, 1991; Campolongo 2007)
# A screening method that varies parameters in COMBINATION rather than one at a
# time. It builds r random "trajectories" through parameter space; along each
# trajectory exactly one parameter changes between consecutive model runs, by a
# fixed step. The per-step change in the model output, divided by the step, is
# an "elementary effect" (EE) for that parameter. Because the moves start from
# scattered points across the whole space, the EEs capture how a parameter acts
# in the presence of the others.
#
# Sensitivity measures per parameter (over its r elementary effects):
#   mu_star mean |EE|             (TOTAL sensitivity, incl. interactions)  <- rank on this
#   sigma   sd of EE              (interactions / non-linearity; high => the
#                                  parameter's effect depends on the others)
#
# This file is model-agnostic: give it parameter bounds and a function that maps
# a named parameter vector to a scalar output. Scripts supply the bounds (from a
# +/- buffer in script 4, or from the Excel SD/Min-Max in script 5) and the
# output function (the equilibrium animal effect).

#------------------------------------------------------------------------#
# morris_trajectories() ----
#------------------------------------------------------------------------#
# morris_trajectories(): build r trajectories over the unit hypercube using the
# standard Morris construction with a p-level grid and step delta = p/(2(p-1)).
# Returns a list of (k+1) x k matrices in [0,1]; each row is a design point and
# consecutive rows differ in exactly one coordinate by +/- delta.
morris_trajectories <- function(k, r, levels = 4L) {
  p     <- as.integer(levels)
  delta <- p / (2 * (p - 1))                 # canonical Morris step (in [0,1])
  grid  <- seq(0, 1 - delta, length.out = p / 2)   # feasible starting levels
  if (length(grid) < 1) grid <- 0

  make_one <- function() {
    xstar <- sample(grid, k, replace = TRUE)         # random base point (1 x k)
    Dstar <- diag(sample(c(-1, 1), k, replace = TRUE))   # random +/- directions
    Pstar <- diag(k)[sample(k), , drop = FALSE]          # random order of moves
    J     <- matrix(1, k + 1, k)
    B     <- lower.tri(matrix(1, k + 1, k + 1))[, seq_len(k)] * 1  # strictly-lower incidence
    # Morris (1991) orientation matrix:
    Bstar <- (J * matrix(xstar, k + 1, k, byrow = TRUE) +
              (delta / 2) * ((2 * B - J) %*% Dstar + J)) %*% Pstar
    # clamp into [0,1] against floating error
    pmin(pmax(Bstar, 0), 1)
  }
  replicate(r, make_one(), simplify = FALSE)
}

#------------------------------------------------------------------------#
# morris_run() ----
#------------------------------------------------------------------------#
# morris_run(): evaluate one Morris design and return the elementary effects.
#   traj_list  list of trajectories in [0,1] (from morris_trajectories)
#   lo, hi     named numeric bounds mapping [0,1] -> real parameter values
#   eval_fun   function(named_param_vector) -> scalar output. It may attach
#              attributes: "converged" (FALSE = not a stable equilibrium; the
#              point is dropped) and "regime" (a non-empty string flags a
#              special regime, e.g. "microbe_collapse"; the point is kept).
#   cl         optional cluster from morris_cluster(): trajectories are then
#              evaluated in parallel (they are independent).
# Returns a data frame: parameter, trajectory, ee, regime ("" unless either
# end of the step was flagged). Attributes: n_eval, n_nonconverged and
# regime_counts (number of evaluations in each flagged regime).
# Factory for the per-trajectory evaluator. Built at top level so that its
# environment holds only what a parallel worker needs (not the cluster itself).
.make_eval_traj <- function(pnames, span, lo, eval_fun) {
  function(X) {
    real <- sweep(sweep(X, 2, span, `*`), 2, lo, `+`)   # map to real values
    n <- nrow(real)
    y <- rep(NA_real_, n); conv <- logical(n); reg <- character(n)
    for (i in seq_len(n)) {
      v <- tryCatch(eval_fun(setNames(as.numeric(real[i, ]), pnames)),
                    error = function(e) structure(NA_real_, converged = FALSE))
      cv <- attr(v, "converged"); rg <- attr(v, "regime")
      conv[i] <- if (is.null(cv)) is.finite(as.numeric(v)) else isTRUE(cv)
      y[i]    <- if (conv[i]) as.numeric(v) else NA_real_
      reg[i]  <- if (is.null(rg) || is.na(rg)) "" else rg
    }
    list(y = y, conv = conv, regime = reg)
  }
}

morris_run <- function(traj_list, lo, hi, eval_fun, verbose = TRUE, cl = NULL) {
  pnames    <- names(lo)
  eval_traj <- .make_eval_traj(pnames, hi - lo, lo, eval_fun)
  res <- if (is.null(cl)) lapply(traj_list, eval_traj) else
    parallel::parLapply(cl, traj_list, eval_traj)

  ee_rows <- list()
  for (ti in seq_along(traj_list)) {
    X <- traj_list[[ti]]; r <- res[[ti]]
    for (i in seq_len(nrow(X) - 1)) {           # consecutive points differ in one param
      dx <- X[i + 1, ] - X[i, ]
      j  <- which(abs(dx) > 1e-9)
      if (length(j) != 1) next
      rg <- unique(setdiff(unlist(strsplit(c(r$regime[i], r$regime[i + 1]), "\\+")), ""))
      ee_rows[[length(ee_rows) + 1]] <- data.frame(
        parameter = pnames[j], trajectory = ti,
        ee = (r$y[i + 1] - r$y[i]) / dx[j],      # EE on the [0,1] scale
        regime = paste(sort(rg), collapse = "+"),
        stringsAsFactors = FALSE)
    }
    if (verbose) cat(sprintf("  trajectory %d/%d done\n", ti, length(traj_list)))
  }
  out <- do.call(rbind, ee_rows)

  all_conv <- unlist(lapply(res, `[[`, "conv"))
  all_reg  <- unlist(lapply(res, `[[`, "regime"))
  tags     <- unlist(strsplit(all_reg[nzchar(all_reg)], "\\+"))
  attr(out, "n_eval")         <- length(all_conv)
  attr(out, "n_nonconverged") <- sum(!all_conv)
  attr(out, "regime_counts")  <- table(factor(tags))
  out
}

#------------------------------------------------------------------------#
# morris_cluster() ----
#------------------------------------------------------------------------#
# morris_cluster(): start a parallel cluster whose workers can run the model.
# Works on Windows, macOS and Linux (PSOCK sockets). Each worker sets the
# project root as its working directory, loads `packages`, sources
# `source_files`, then receives the objects named in `export` from `envir`.
# Returns NULL when n_cores <= 1 (serial). Stop it with parallel::stopCluster().
morris_cluster <- function(n_cores, source_files, packages = c("deSolve", "rootSolve", "yaml"),
                           export = character(0), envir = globalenv()) {
  if (is.null(n_cores) || n_cores <= 1) return(NULL)
  cl <- parallel::makeCluster(n_cores)
  wd <- getwd()
  parallel::clusterCall(cl, function(wd, pkgs, files) {
    setwd(wd)
    for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
    for (f in files) source(f)
    NULL
  }, wd, packages, source_files)
  if (length(export)) parallel::clusterExport(cl, export, envir = envir)
  cl
}

#------------------------------------------------------------------------#
# morris_summary() ----
#------------------------------------------------------------------------#
# morris_summary(): the traditional Morris metrics per parameter from the
# elementary effects:
#   mu_star mean |EE|  -> TOTAL sensitivity (rank on this)
#   sigma   sd of EE    -> interactions / non-linearity
morris_summary <- function(ee_df) {
  ee_df <- ee_df[is.finite(ee_df$ee), , drop = FALSE]
  agg <- lapply(split(ee_df$ee, ee_df$parameter), function(v) {
    data.frame(mu_star = mean(abs(v)),
               sigma = if (length(v) > 1) sd(v) else 0,
               n_ee = length(v))
  })
  out <- do.call(rbind, Map(function(nm, d) cbind(parameter = nm, d),
                            names(agg), agg))
  rownames(out) <- NULL
  out[order(-out$mu_star), , drop = FALSE]
}

#------------------------------------------------------------------------#
# morris_bootstrap() ----
#------------------------------------------------------------------------#
# morris_bootstrap(): convergence check for the Morris ranking. Resamples the r
# trajectories with replacement B times and recomputes mu_star and sigma, giving
# percentile confidence intervals and the stability of each parameter's rank.
#   ee_df  elementary effects from morris_run() (parameter, trajectory, ee)
#   B      bootstrap replicates
#   top_k  report P(parameter ranks in the top_k by mu_star)
#   conf   interval level
# Returns one row per parameter: mu_star_lo/hi, sigma_lo/hi, rank_median,
# rank_lo/hi, p_top_k. Narrow mu_star intervals and tight rank intervals mean
# the number of trajectories was sufficient.
# Implementation: per-trajectory sums of EE, |EE|, EE^2 and counts are combined
# with multinomial resampling weights, so each replicate is a matrix product.
morris_bootstrap <- function(ee_df, B = 1000, top_k = 5, conf = 0.95, seed = NULL) {
  ee_df <- ee_df[is.finite(ee_df$ee), , drop = FALSE]
  if (!nrow(ee_df)) return(NULL)
  if (!is.null(seed)) set.seed(seed)

  traj   <- sort(unique(ee_df$trajectory))
  params <- sort(unique(ee_df$parameter))
  ti <- match(ee_df$trajectory, traj); pj <- match(ee_df$parameter, params)
  nt <- length(traj); np <- length(params)

  acc <- function(v) { m <- matrix(0, nt, np); for (i in seq_along(v)) m[ti[i], pj[i]] <- m[ti[i], pj[i]] + v[i]; m }
  S_abs <- acc(abs(ee_df$ee)); S <- acc(ee_df$ee); Q <- acc(ee_df$ee^2)
  N     <- acc(rep(1, nrow(ee_df)))

  # multinomial resampling weights: W[b, t] = times trajectory t was drawn
  W <- t(vapply(seq_len(B), function(b) tabulate(sample.int(nt, nt, replace = TRUE), nt),
                numeric(nt)))
  if (nt == 1) W <- matrix(W, ncol = 1)

  n_b    <- W %*% N
  mustar <- (W %*% S_abs) / n_b
  mean_b <- (W %*% S) / n_b
  var_b  <- ((W %*% Q) / n_b - mean_b^2) * n_b / pmax(n_b - 1, 1)
  sigma  <- sqrt(pmax(var_b, 0))
  colnames(mustar) <- colnames(sigma) <- params

  # rank within each replicate (1 = most sensitive); NA-safe
  mustar_r <- mustar; mustar_r[!is.finite(mustar_r)] <- -Inf
  ranks <- t(apply(mustar_r, 1, function(x) rank(-x, ties.method = "average")))
  if (np == 1) ranks <- matrix(ranks, ncol = 1)
  colnames(ranks) <- params

  a  <- (1 - conf) / 2
  qf <- function(m, p) apply(m, 2, stats::quantile, probs = p, na.rm = TRUE)
  data.frame(
    parameter   = params,
    mu_star_lo  = qf(mustar, a), mu_star_hi = qf(mustar, 1 - a),
    sigma_lo    = qf(sigma, a),  sigma_hi   = qf(sigma, 1 - a),
    rank_median = qf(ranks, 0.5),
    rank_lo     = qf(ranks, a),  rank_hi    = qf(ranks, 1 - a),
    p_top_k     = colMeans(ranks <= top_k),
    stringsAsFactors = FALSE, row.names = NULL)
}
