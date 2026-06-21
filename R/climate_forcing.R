# Climate forcing function:

climate_forcing_function <- function(time, parms) {
  parms <- as.list(parms)

  # Day-of-year (1..365): no leap years:
  doy <- (time %% 365) + 1

  Temp  <- parms$MAT     + parms$T_amp     * sin(2 * pi * (doy - 110) / 365)
  theta <- parms$MAtheta + parms$theta_amp * cos(4 * pi * (doy - 110) / 365)

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
  # Seasonal NPP weight: the SAME growing-season activity that drives root/wood
  # dormancy (logistic temperature activity x moisture limitation), normalized
  # to sum to 1 over the year. Multiplying an ANNUAL NPP parameter by this
  # weight spreads it over the year and integrates back to the annual total.
  # --------------------------------------------------
  doy_all   <- 1:365
  Temp_all  <- parms$MAT     + parms$T_amp     * sin(2 * pi * (doy_all - 110) / 365)
  theta_all <- parms$MAtheta + parms$theta_amp * cos(4 * pi * (doy_all - 110) / 365)
  fT_all    <- 1 / (1 + exp(-parms$k_root_dormancy * (Temp_all - parms$root_dormancy_temp)))
  fth_all   <- pmin(1, pmax(0, theta_all) / parms$theta_opt)
  act_all   <- fT_all * fth_all
  npp_w     <- if (sum(act_all) > 0) act_all / sum(act_all) else rep(1/365, 365)

  c(
    Temp = Temp,
    theta = theta,
    litterfall_prob_val = prob_vals[doy],
    npp_weight_val = npp_w[doy]
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
      litterfall_prob_val = 1/365,
      npp_weight_val = 1/365
    )
    stopifnot(is.numeric(forcing), !any(is.na(forcing)))
    forcing
  }
}
