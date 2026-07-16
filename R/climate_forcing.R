# Climate forcing function:
#
# All quantities are SMOOTH, continuous functions of `time` (no day-of-year
# table lookups), so there are no artificial daily "kinks" for the ODE solver's
# adaptive step controller to trip over. Temp/theta are sinusoids; the three
# input-timing weights are closed-form functions of time, each divided by a
# precomputed annual normaliser so they integrate to 1 over the year.

# ------------------------------------------------------------
# .forcing_norms(): precompute (once per parameter set) the annual normalising
# constants for the growing-season and leaf-litter weights, by integrating
# their UNNORMALISED shapes over a fine within-year grid. Attaching these to
# `parms` avoids recomputing them on every derivative evaluation.
# ------------------------------------------------------------
.grow_shape <- function(doy, p) {
  Temp  <- p$MAT + p$T_amp * sin(2 * pi * (doy - 110) / 365)
  theta <- if (p$N_theta_peaks == 2) p$MAtheta + p$theta_amp * cos(4 * pi * (doy - 110) / 365)
           else                      p$MAtheta + p$theta_amp * sin(2 * pi * (doy - 110) / 365)
  fT  <- 1 / (1 + exp(-p$k_root_dormancy * (Temp - p$root_dormancy_temp)))
  fth <- pmin(1, pmax(0, theta) / p$theta_opt)
  fT * fth                                     # unnormalised growing-season activity, smooth in doy
}
# wrapped Gaussian litterfall pulse (smooth, periodic), unnormalised.
# Sum the Gaussian and its wrapped copies; works for scalar OR vector `doy`
# (Reduce accumulates a numeric of the same length as `doy`, so there is no
# sapply/rowSums shape ambiguity when doy has length 1).
.litter_shape <- function(doy, p, K = 2) {
  ks <- (-K):K
  Reduce(`+`, lapply(ks, function(k)
    dnorm(doy - p$litter_peak_doy - k * 365, sd = p$litter_width_d)),
    accumulate = FALSE)
}
forcing_norms <- function(parms) {
  p  <- as.list(parms)
  gg <- 1:365
  list(grow   = mean(.grow_shape(gg, p)) * 365,      # integral over the year of the grow shape
       litter = mean(.litter_shape(gg, p)) * 365)    # integral over the year of the litter pulse
}

climate_forcing_function <- function(time, parms) {
  p <- as.list(parms)

  # Day-of-year as a CONTINUOUS value in [1, 366) (no floor -> no daily steps)
  doy <- (time %% 365) + 1

  Temp  <- p$MAT + p$T_amp * sin(2 * pi * (doy - 110) / 365)
  if (p$N_theta_peaks == 2) {
    theta <- p$MAtheta + p$theta_amp * cos(4 * pi * (doy - 110) / 365)
  } else if (p$N_theta_peaks == 1) {
    theta <- p$MAtheta + p$theta_amp * sin(2 * pi * (doy - 110) / 365)
  } else stop("N_theta_peaks must be 1 or 2.")

  # annual normalisers (precomputed on parms if available, else compute now)
  nrm <- if (!is.null(p$.forcing_norms)) p$.forcing_norms else forcing_norms(p)

  # --------------------------------------------------
  # ALLOCATION-SPECIFIC INPUT TIMING, evaluated continuously at `doy`. Each
  # weight is a rate per day that integrates to 1 over the year:
  #   root_input_weight   growing-season activity / (its annual integral)
  #   wood_input_weight   uniform, 1/365
  #   leaf_litter_weight  fall Gaussian pulse + small growing-season trickle,
  #                       mixed by leaf_litter_summer_frac
  # --------------------------------------------------
  grow_w <- if (nrm$grow > 0) .grow_shape(doy, p) / nrm$grow else 1/365
  wood_w <- 1/365
  sfrac  <- if (!is.null(p$leaf_litter_summer_frac)) p$leaf_litter_summer_frac else 0.2
  litter_pulse <- if (nrm$litter > 0) .litter_shape(doy, p) / nrm$litter else 1/365
  leaf_w <- (1 - sfrac) * litter_pulse + sfrac * grow_w      # already integrates to 1

  c(
    Temp = Temp,
    theta = theta,
    root_input_weight  = unname(grow_w),
    wood_input_weight  = unname(wood_w),
    leaf_litter_weight = unname(leaf_w)
  )
}


# ------------------------------------------------------------
# Seasonal forcing (a function of time)
# ------------------------------------------------------------
make_climate_forcing <- function(parms) {
  p <- as.list(parms)
  p$.forcing_norms <- forcing_norms(p)          # precompute annual normalisers once
  function(time) {
    forcing <- climate_forcing_function(time = time, parms = p)
    stopifnot(is.numeric(forcing), !any(is.na(forcing)))
    forcing
  }
}

# ------------------------------------------------------------
# Equilibrium (constant) forcing: annual-mean climate, with uniform daily
# weights so NPP_daily = NPP_annual / 365 and litterfall is spread evenly.
# Used for stode()/runsteady() warm-starts.
# ------------------------------------------------------------
make_climate_forcing_equilibrium <- function(parms) {
  function(time) {
    forcing <- c(
      Temp = parms$MAT,
      theta = parms$MAtheta,
      root_input_weight  = 1/365,
      wood_input_weight  = 1/365,
      leaf_litter_weight = 1/365
    )
    stopifnot(is.numeric(forcing), !any(is.na(forcing)))
    forcing
  }
}
