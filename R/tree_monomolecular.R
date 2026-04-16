# Forcing function: Monomolecular tree biomass + litter/exudate fluxes
# ------------------------------------------------------------
tree_forcing_monomolecular <- function(time, parms) {
  parms <- as.list(parms)
  
  # Helper: default if NULL
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # ---- Required parameters ----
  req <- c(
    "TBTmax","k_monomol","B0",
    "a_leaf","a_wood","a_root",
    "k_litterfall_ann","k_mort_leaf","k_mort_wood","k_mort_root","k_exudate"
  )
  miss <- setdiff(req, names(parms))
  if (length(miss)) stop("Missing required parms: ", paste(miss, collapse = ", "))
  
  TBTmax         <- parms$TBTmax
  k_monomol      <- parms$k_monomol
  B0             <- parms$B0
  t0             <- parms$t0 %||% 0
  
  a_leaf         <- parms$a_leaf
  a_wood         <- parms$a_wood
  a_root         <- parms$a_root
  a_root_herb    <- parms$a_root_herb
  
  # Annual-mean litterfall rate (1/yr): controls annual total
  k_litterfall_ann <- parms$k_litterfall_ann
  k_litterfall_herb_ann <- parms$k_litterfall_herb_ann
  
  k_mort_leaf    <- parms$k_mort_leaf
  k_mort_wood    <- parms$k_mort_wood
  k_mort_root    <- parms$k_mort_root
  k_exudate      <- parms$k_exudate
  
  # ---- Optional litterfall seasonality parameters (sane defaults) ----
  Lfs <- parms$year_length_d %||% 365          # keep simple; default 365
  # Peak day-of-year for fall litterfall (default ~ Oct 15 = 288)
  litter_peak_doy <- parms$litter_peak_doy %||% 288
  # Width (sd) of litterfall pulse in days; smaller = more concentrated in fall
  litter_width_d  <- parms$litter_width_d  %||% 30
  
  # ---- Optional removals (default 0) ----
  leaf_harvest    <- parms$leaf_harvest    %||% 0
  wood_harvest    <- parms$wood_harvest    %||% 0
  root_harvest    <- parms$root_harvest    %||% 0
  herbivory_leaf  <- parms$herbivory_leaf  %||% 0
  herbivory_wood  <- parms$herbivory_wood  %||% 0
  herbivory_root  <- parms$herbivory_root  %||% 0
  
  # ---- Time handling ----
  # time is in days; convert to years for growth equation (k_monomol in 1/yr)
  tau <- max(0, time / 365 - t0)
  
  # Day-of-year (1..365). Keeps it simple (ignores leap years).
  doy <- (time %% 365) + 1
  
  # Model temperature:
  Temp_function <- function(t) {
    parms$MAT + parms$T_amp * sin(2 * pi * (t - 110) / 365)
  }
  
  theta_function <- function(t) {
    parms$MAtheta + parms$theta_amp * cos(4 * pi * (t - 110) / 365)
  }
  
  theta = theta_function(doy)
  Temp = Temp_function(doy)
  
  # ---- Monomolecular growth of total tree biomass ----
  B_tree <- TBTmax - (TBTmax - B0) * exp(-k_monomol * tau)
  
  # ---- Allocation ----
  if (abs(a_leaf + a_wood + a_root - 1) > 1e-8) {
    stop("a_leaf + a_wood + a_root must sum to 1")
  }
  B_leaf <- a_leaf * B_tree
  B_wood <- a_wood * B_tree
  B_root <- a_root * B_tree
  
  # Herbaceous plant NPP during stand regeneration:
  # Based on Fig. 1 from Crow et al. 1991 in Ecol. App.
  # assuming a linear decline in proportion herb based on their data.
  # This is aboveground biomass
  
  B_herb = ifelse(tau == 0,
                  0,
                  (0.67*exp(-0.2*tau))*B_tree/(1 - (0.67*exp(-0.2*tau)))
  )
  
  # ---- Seasonal litterfall weighting (mean ~ 1 over year) ----
  wrapped_pdf <- function(t, mu, sigma, Lfss = 365, K = 1) {
    ks <- (-K):K
    rowSums( sapply(ks, function(k) dnorm(t - mu - k*Lfss, sd = sigma)) )
  }
  
  pdf_vals <- wrapped_pdf(1:365, mu = 288, sigma = 30)
  prob_vals <- pdf_vals / sum(pdf_vals)
  
  # ---- Fluxes ----
  litterfall     <- k_litterfall_ann * prob_vals[doy] * B_leaf + 
    # Herbaceous plants are 100% lost:
    k_litterfall_herb_ann * prob_vals[doy] * B_herb
  
  leaf_mortality <- k_mort_leaf * B_leaf
  wood_mortality <- k_mort_wood * B_wood
  root_mortality <- k_mort_root * B_root + k_mort_root* B_herb * (1- a_root_herb)/a_root_herb
  exudates       <- k_exudate   * B_root + k_exudate* B_herb * (1- a_root_herb)/a_root_herb
  
  # Losses from tree biomass (diagnostic)
  losses_tree <- litterfall + leaf_mortality + wood_mortality + root_mortality + exudates +
    leaf_harvest + wood_harvest + root_harvest +
    herbivory_leaf + herbivory_wood + herbivory_root
  
  c(
    B_tree      = B_tree,
    B_herb      = B_herb,
    litterfall_f  = litterfall,
    leaf_mort   = leaf_mortality,
    wood_mort   = wood_mortality,
    root_mort   = root_mortality,
    exudates_f    = exudates,
    losses_tree = losses_tree,
    Temp = Temp,
    theta = theta
  )
}
