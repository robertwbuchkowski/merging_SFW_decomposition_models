# Climate forcing function:

climate_forcing_function <- function(time, parms) {
  parms <- as.list(parms)
  
  # Day-of-year (1..365): no leap years:
  doy <- (time %% 365) + 1
  
  Temp <- parms$MAT + parms$T_amp * sin(2 * pi * (doy - 110) / 365)
  theta <- parms$MAtheta + parms$theta_amp * cos(4 * pi * (doy - 110) / 365)
  
  # --------------------------------------------------
  # Seasonal litterfall (shoots & leaves)
  # --------------------------------------------------
  wrapped_pdf <- function(t, mu, sigma, L = 365, K = 1) {
    ks <- (-K):K
    rowSums(sapply(ks, function(k)
      dnorm(t - mu - k * L, sd = sigma)))
  }
  
  pdf_vals  <- wrapped_pdf(1:365, parms$litter_peak_doy, parms$litter_width_d)
  prob_vals <- pdf_vals / sum(pdf_vals)
  
  c(
    Temp = Temp,
    theta = theta,
    litterfall_prob_val = prob_vals[doy]
  )
}


# ------------------------------------------------------------
# make_climate_forcing.R
# Construct a vegetation forcing function from tree model
# ------------------------------------------------------------

make_climate_forcing <- function(parms) {
  
  # Return a function of time
  function(time) {
    
    forcing <- climate_forcing_function(
      time  = time,
      parms = parms
    )
    
    # Defensive checks (optional but recommended)
    stopifnot(
      is.numeric(forcing),
      !any(is.na(forcing))
    )
    
    forcing
  }
}

# ------------------------------------------------------------
# make_climate_forcing.R
# Construct a vegetation forcing function from tree model
# ------------------------------------------------------------

make_climate_forcing_equilibrium <- function(parms) {
  
  # Return a function of time
  function(time) {
    
    forcing <- c(
      Temp = parms$MAT,
      theta = parms$MAtheta,
      litterfall_prob_val = 1/365
    )
    
    # Defensive checks (optional but recommended)
    stopifnot(
      is.numeric(forcing),
      !any(is.na(forcing))
    )
    
    forcing
  }
}

