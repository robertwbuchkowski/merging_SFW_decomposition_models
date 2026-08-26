# ============================================================
# MORRIS ELEMENTARY-EFFECTS SENSITIVITY (method of Morris, 1991; Campolongo 2007)
# ------------------------------------------------------------
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
# ============================================================

# ------------------------------------------------------------
# morris_trajectories(): build r trajectories over the unit hypercube using the
# standard Morris construction with a p-level grid and step delta = p/(2(p-1)).
# Returns a list of (k+1) x k matrices in [0,1]; each row is a design point and
# consecutive rows differ in exactly one coordinate by +/- delta.
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# morris_run(): evaluate one Morris design and return the elementary effects.
#   traj_list  list of trajectories in [0,1] (from morris_trajectories)
#   lo, hi     named numeric bounds mapping [0,1] -> real parameter values
#   eval_fun   function(named_param_vector) -> scalar output (NA allowed)
# Returns a data frame of elementary effects: parameter, trajectory, ee.
# ------------------------------------------------------------
morris_run <- function(traj_list, lo, hi, eval_fun, verbose = TRUE) {
  pnames <- names(lo)
  k      <- length(pnames)
  span   <- hi - lo
  ee_rows <- list()

  for (ti in seq_along(traj_list)) {
    X   <- traj_list[[ti]]                    # (k+1) x k in [0,1]
    real <- sweep(X, 2, span, `*`)
    real <- sweep(real, 2, lo, `+`)           # map to real parameter values
    colnames(real) <- pnames

    y <- rep(NA_real_, nrow(real))
    for (i in seq_len(nrow(real))) {
      pv <- setNames(as.numeric(real[i, ]), pnames)
      y[i] <- tryCatch(eval_fun(pv), error = function(e) NA_real_)
    }

    # each consecutive pair differs in exactly one parameter
    for (i in seq_len(nrow(X) - 1)) {
      dx  <- X[i + 1, ] - X[i, ]
      j   <- which(abs(dx) > 1e-9)
      if (length(j) != 1) next                 # skip degenerate steps
      dy  <- y[i + 1] - y[i]
      ee  <- dy / dx[j]                         # EE on the [0,1] scale
      ee_rows[[length(ee_rows) + 1]] <- data.frame(
        parameter = pnames[j], trajectory = ti, ee = ee,
        stringsAsFactors = FALSE)
    }
    if (verbose) cat(sprintf("  trajectory %d/%d done\n", ti, length(traj_list)))
  }
  do.call(rbind, ee_rows)
}

# ------------------------------------------------------------
# morris_summary(): the traditional Morris metrics per parameter from the
# elementary effects:
#   mu_star mean |EE|  -> TOTAL sensitivity (rank on this)
#   sigma   sd of EE    -> interactions / non-linearity
# ------------------------------------------------------------
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
