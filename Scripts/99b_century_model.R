# Century model simulations using Ambramoff code:

# Load in the parameters:

parms_century <- list(
  w1 = 33.61,
  w2 = 8.42,
  t1 = 18.08,
  t2 = 12.99,
  t3 = 25.92,
  t4 = 0.038,
  c1 = 0.83,
  c2 = 0.67,
  k_active = 0.02,
  k_slow = 0.0005,
  k_passive = 0.000091,
  slow_to_active = 0.42,
  slow_to_passive = 0.02,
  passive_to_active = 0.45,
  active_to_passive = 0.004,
  k_strlitter = 0.1,
  k_metlitter = 0.045,
  metlitter_to_active = 0.45,
  strlitter_to_active = 0.5,
  strlitter_to_slow = 0.7,
  LigFrac = 0.2,
  
  # OLD PARAMS:
  force_eqm = 0,
  B0 = 0,
  # Drivers and allocation
  k_monomol = 0.014,      # Saturation control for the monomolecular model, from the reference cited in the model for now.
  TBTmax = 10000,      # Maximum biomass estimated from the portions. Could get higher, but need to explore.
  a_leaf = 0.75*0.1,
  a_wood = 0.75*0.9,
  a_root = 0.25,
  a_root_herb = 0.5,
  # Plant turnover (CBM-inspired)
  k_litterfall_ann = 0.35,   # yr^-1 (foliar litter; HW may approach ~1 yr^-1)
  litter_peak_doy = 288, # Fall
  litter_width_d = 30, # width of peak
  k_mort_leaf  = 0.02/365,   # yr^-1 (non-litterfall leaf mortality)
  k_mort_wood  = 0.025/365,   # yr^-1 (wood mortality to CWD)
  k_mort_root  = 0.08/365,   # yr^-1 (aggregate root mortality)
  k_exudate    = 0.01/365,   # yr^-1 (placeholder; tune to data)
  root_to_organic = 0.7,
  
  # Forcings:
  MAT = 6, # From the occupancy model sites
  T_amp = 26, # Set to get summer and winter temperatures correct
  MAtheta = 0.30, # m^3 m^-3; moderately moist soil
  theta_amp = 0.15,
  
  pct_claysilt = 70,      # %; moderate clay+silt content; from my data.
  
  field_cap = 0.39
)



# ------------------------------------------------------------
# Updated Century (From: Millennial v2 ODE)
# ------------------------------------------------------------
#Author: Robert Buchkowski after Rose Abramoff, after Parton et al. (1987)
#Date: Jan 25, 2026; Sep 11, 2021

#This function contains the system of equations for the Century model, which was developed in Parton et al. (1987). 
#The equation numbers correspond to those in Abramoff et al. (2021) Appendix B, with parameters defined in Table A2.

century_ode <- function(time, state, parms){
  
  with(c(state, parms), {
    # Replace tree forcing with equilibrium if needed:
    if(force_eqm){
      # ----------------------------
      # Tree forcing (no plant states)
      # ----------------------------
      tf <- tree_forcing_monomolecular(time, parms)
      
      tf <- tree_forcing_monomolecular(36500000, parms)
      
      TotalBiomassTree <- tf["B_tree"]
      litterfall       <- tf["litterfall_f"]
      leaf_mortality   <- tf["leaf_mort"]
      wood_mortality   <- tf["wood_mort"]
      root_mortality   <- tf["root_mort"]
      exudates         <- tf["exudates_f"]
      # ----------------------------
      # Forcings at current time (Fi, T, theta)
      # ----------------------------
      T_t              <- tf["Temp"] # °C
      theta_t          <- tf["theta"] # m^3 m^-3
    }else{
      # ----------------------------
      # Tree forcing (no plant states)
      # ----------------------------
      TotalBiomassTree <- B_tree(time)
      litterfall       <- litterfall_f(time)
      leaf_mortality   <- leaf_mort(time)
      wood_mortality   <- wood_mort(time)
      root_mortality   <- root_mort(time)
      exudates         <- exudates_f(time)
      # ----------------------------
      # Forcings at current time (Fi, T, theta)
      # ----------------------------
      T_t              <- Temp(time) # °C
      theta_t          <- theta(time) # m^3 m^-3
    }
    
    # ----------------------------------#
    # Abiotic scalars for decomposition:
    # ----------------------------------#
    
    #Equation B1
    t_scalar <- (t2 + (t3 / pi) * atan(pi * t4 * (T_t - t1))) /
      (t2 + (t3 / pi) * atan(pi * t4 *(30.0 - t1)))
    
    #Equation B2
    w_scalar <- 1.0 / (1.0 + w1 * exp(-w2 * theta_t/field_cap))
    
    #Equation B3
    f_TEX = c1 - c2*pct_claysilt*0.01
    
    #Equation B4
    f_StrLitter = StrLitter * k_strlitter * t_scalar * w_scalar * exp(-3*LigFrac)
    
    #Equation B5
    f_MetLitter = MetLitter * k_metlitter * t_scalar * w_scalar  
    
    #Equation B6
    f_ACTIVE <- ACTIVE * k_active * t_scalar * w_scalar * f_TEX
    
    #Equation B7 
    f_SLOW <- SLOW * k_slow * t_scalar * w_scalar
    
    #Equation B8
    f_PASSIVE <- PASSIVE * k_passive * t_scalar * w_scalar
    
    #Equation B9
    dStrLitter = wood_mortality + root_to_organic*root_mortality - f_StrLitter
    
    #Equation B10
    dMetLitter = litterfall + leaf_mortality + exudates - f_MetLitter
    
    #Equation B11
    dACTIVE <- (1 - root_to_organic)*root_mortality + (1-LigFrac) * strlitter_to_active * f_StrLitter + metlitter_to_active * f_MetLitter  + f_SLOW * slow_to_active + f_PASSIVE * passive_to_active - f_ACTIVE
    
    #Equation B12
    dSLOW <-  LigFrac * strlitter_to_slow * f_StrLitter + f_ACTIVE * (1-f_TEX-active_to_passive) - f_SLOW
    
    #Equation B13
    dPASSIVE <- f_ACTIVE * active_to_passive + f_SLOW * slow_to_passive - f_PASSIVE

    # ---------------------------
    # Return list for deSolve
    # ---------------------------
    list(
      c(dStrLitter, dMetLitter, dACTIVE, dSLOW, dPASSIVE)
    )
  })
}


