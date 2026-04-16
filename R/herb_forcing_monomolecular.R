# ------------------------------------------------------------
# Forcing function: Monomolecular herbaceous biomass + litter/exudates
# ------------------------------------------------------------
herb_forcing_monomolecular <- function(time, parms) {
  parms <- as.list(parms)
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # ---- Required parameters ----
  req <- c(
    "TBHmax",          # max total herbaceous biomass
    "k_monomol",        # monomolecular growth rate (1/yr)
    "B0",               # initial biomass
    "a_root_herb",      # root fraction (shoot = 1 - a_root_herb)
    "k_mort_root",
    "k_mort_leaf",
    "k_exudate",
    "k_litterfall_ann"  # controls annual shoot turnover
  )
  miss <- setdiff(req, names(parms))
  if (length(miss)) stop("Missing required parms: ", paste(miss, collapse = ", "))
  
  TBHmax     <- parms$TBHmax
  k_monomol  <- parms$k_monomol
  B0         <- parms$B0
  t0         <- parms$t0 %||% 0
  
  a_root     <- parms$a_root_herb
  a_shoot    <- 1 - a_root
  
  k_mort_root <- parms$k_mort_root
  k_mort_leaf <- parms$k_mort_leaf
  k_exudate   <- parms$k_exudate
  k_litterfall_herb_ann <- parms$k_litterfall_herb_ann
  
  # ---- Seasonality parameters ----
  year_length_d    <- parms$year_length_d %||% 365
  litter_peak_doy  <- parms$litter_peak_doy %||% 288
  litter_width_d   <- parms$litter_width_d  %||% 30
  
  # ---- Time handling ----
  tau <- max(0, time / 365 - t0)
  doy <- (time %% year_length_d) + 1
  
  # ---- Climate forcing ----
  Temp_function <- function(t) {
    parms$MAT + parms$T_amp * sin(2 * pi * (t - 110) / 365)
  }
  
  theta_function <- function(t) {
    parms$MAtheta + parms$theta_amp * cos(4 * pi * (t - 110) / 365)
  }
  
  Temp  <- Temp_function(doy)
  theta <- theta_function(doy)
  
  # ---- Herbaceous biomass growth (monomolecular) ----
  B_herb <- TBHmax - (TBHmax - B0) * exp(-k_monomol * tau)
  
  B_shoot <- a_shoot * B_herb
  B_root  <- a_root  * B_herb
  
  # ---- Seasonal litterfall weighting (normalized) ----
  wrapped_pdf <- function(t, mu, sigma, L = 365, K = 1) {
    ks <- (-K):K
    rowSums(sapply(ks, function(k) dnorm(t - mu - k * L, sd = sigma)))
  }
  
  pdf_vals  <- wrapped_pdf(1:365, mu = litter_peak_doy, sigma = litter_width_d)
  prob_vals <- pdf_vals / sum(pdf_vals)
  
  # ---- Fluxes ----
  # Aboveground herbs fully senesce annually (controlled by k_litterfall_ann)
  litterfall <- k_litterfall_herb_ann * prob_vals[doy] * B_shoot
  
  leaf_mortality <- k_mort_leaf * B_shoot
  root_mortality <- k_mort_root * B_root
  exudates       <- k_exudate   * B_root
  
  total_losses <- litterfall + root_mortality + exudates
  
  c(
    B_tree        = 0,
    B_herb       = B_herb,
    B_root        = B_root,
    litterfall_f  = litterfall,
    root_mort     = root_mortality,
    wood_mort     = 0,
    leaf_mort     = leaf_mortality,
    exudates_f    = exudates,
    losses_herb   = total_losses,
    Temp          = Temp,
    theta         = theta
  )
}
