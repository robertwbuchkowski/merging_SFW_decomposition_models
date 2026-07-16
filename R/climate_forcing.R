# Climate forcing function:

climate_forcing_function <- function(time, parms) {
  parms <- as.list(parms)

  # Day-of-year (1..365): no leap years:
  doy <- (time %% 365) + 1

  Temp  <- parms$MAT     + parms$T_amp     * sin(2 * pi * (doy - 110) / 365)
  
  if(parms$N_theta_peaks == 2){
    theta <- parms$MAtheta + parms$theta_amp * cos(4 * pi * (doy - 110) / 365)
  }else{
    if(parms$N_theta_peaks == 1){
      theta <- parms$MAtheta + parms$theta_amp * sin(2 * pi * (doy - 110) / 365)
    }else{
      stop("N_theta_peaks must be 1 or 2.")
    }
  }
  

  # --------------------------------------------------
  # Seasonal litterfall (shoots & leaves): normalized to sum to 1 over a year
  # --------------------------------------------------
  wrapped_pdf <- function(t, mu, sigma, L = 365, K = 1) {
    ks <- (-K):K
    rowSums(sapply(ks, function(k)
      dnorm(t - mu - k * L, sd = sigma)))
  }

  pdf_vals  <- wrapped_pdf(1:365, parms$litter_peak_doy, parms$litter_width_d)
  prob_vals <- pdf_vals / sum(pdf_vals)

  # --------------------------------------------------
  # ALLOCATION-SPECIFIC INPUT TIMING. Instead of one seasonal NPP weight, each
  # destination gets its OWN within-year weight (each sums to 1 over the year,
  # so annual input = allocation x annual NPP regardless of shape):
  #   root_input_weight   growing-season only (0 in dormancy) -> the activity
  #                       index, normalized.
  #   wood_input_weight   even over the whole year (uniform) -> CWD gets a
  #                       steady drip (dead wood does not track phenology).
  #   leaf_litter_weight  a big autumn peak PLUS a small growing-season trickle
  #                       (green-leaf turnover): a blend of the fall litterfall
  #                       pulse and the growing-season weight, mixed by
  #                       leaf_litter_summer_frac (0 = all autumn, 1 = all
  #                       growing-season).
  # --------------------------------------------------
  doy_all   <- 1:365
  Temp_all  <- parms$MAT     + parms$T_amp     * sin(2 * pi * (doy_all - 110) / 365)
  if(parms$N_theta_peaks == 2){
    theta_all <- parms$MAtheta + parms$theta_amp * cos(4 * pi * (doy_all - 110) / 365)
  }else{
    if(parms$N_theta_peaks == 1){
      theta_all <- parms$MAtheta + parms$theta_amp * sin(2 * pi * (doy_all - 110) / 365)
    }else{
      stop("N_theta_peaks must be 1 or 2.")
    }
  }
  fT_all    <- 1 / (1 + exp(-parms$k_root_dormancy * (Temp_all - parms$root_dormancy_temp)))
  fth_all   <- pmin(1, pmax(0, theta_all) / parms$theta_opt)
  act_all   <- fT_all #* fth_all
  grow_w    <- if (sum(act_all) > 0) act_all / sum(act_all) else rep(1/365, 365)   # growing season

  wood_w    <- rep(1/365, 365)                                                     # uniform

  sfrac     <- if (!is.null(parms$leaf_litter_summer_frac)) parms$leaf_litter_summer_frac else 0.2
  leaf_w    <- (1 - sfrac) * prob_vals + sfrac * grow_w                            # fall peak + summer trickle
  leaf_w    <- leaf_w / sum(leaf_w)                                                # (already ~1; normalize to be safe)

  c(
    Temp = Temp,
    theta = theta,
    root_input_weight  = grow_w[doy],
    wood_input_weight  = wood_w[doy],
    leaf_litter_weight = leaf_w[doy]
  )
}


# ------------------------------------------------------------
# Seasonal forcing (a function of time)
# ------------------------------------------------------------
make_climate_forcing <- function(parms) {
  function(time) {
    forcing <- climate_forcing_function(time = time, parms = parms)
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
