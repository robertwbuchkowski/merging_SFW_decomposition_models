
# ------------------------------------------------------------
# Updated Century
# ------------------------------------------------------------
#Author: Robert Buchkowski after Rose Abramoff, after Parton et al. (1987)
#Date: Jan 25, 2026; Sep 11, 2021

#This function contains the system of equations for the Century model, which was developed in Parton et al. (1987). 
#The equation numbers correspond to those in Abramoff et al. (2021) Appendix B, with parameters defined in Table A2.

century_model <- function(time, state, parms){
  
  with(c(state, parms), {
    
    # ----------------------------
    # ---- Get climate forcing ----
    # ----------------------------
    forcing <- climate_forcing(time)
    T_t              <- forcing["Temp"] # °C 
    theta_t          <- forcing["theta"] # m^3 m^-3
    litterfall_prob_val        <- forcing["litterfall_prob_val"] # probability of litterfall
    
    # --------------------------------------------------
    # Climate scalars for plant productivity
    # --------------------------------------------------
    f_T     <- Q10 ^ ((T_t - Tref) / 10)
    f_theta <- pmin(1, theta_t / theta_opt)
    
    # --------------------------------------------------
    # Herbaceous plant carbon fluxes
    # --------------------------------------------------
    if(C_leaf_herb > 0){
      GPP_herb <- GPPmax_herb * f_T * f_theta
      Ra_herb <- maint_resp * (C_leaf_herb + C_root_herb) +
        growth_resp * GPP_herb
    }else{
      GPP_herb <- 0
      Ra_herb <- 0
    }
    
    NPP_herb <- pmax(0, GPP_herb - Ra_herb)
    
    leaf_growth_herb <- (1 - a_root_herb) * NPP_herb
    root_growth_herb  <- a_root_herb * NPP_herb
    
    # --------------------------------------------------
    # Tree plant carbon fluxes
    # --------------------------------------------------
    if(C_leaf_tree > 0){
      GPP_tree <- GPPmax_tree * f_T * f_theta
      
      Ra_tree <- maint_resp * (C_leaf_tree + C_wood_tree + C_root_tree) +
        growth_resp * GPP_tree
    }else{
      GPP_tree <- 0
      Ra_tree <- 0
    }
    
    NPP_tree <- pmax(0, GPP_tree - Ra_tree)
    
    if(abs(a_leaf_tree + a_wood_tree + a_root_tree - 1) > 1e-6) {stop("Tree allocation fractions must sum to 1")}
    
    leaf_growth_tree <- a_leaf_tree * NPP_tree
    wood_growth_tree <- a_wood_tree * NPP_tree
    root_growth_tree  <- a_root_tree * NPP_tree
    
    # --------------------------------------------------
    # Continuous plant losses
    # --------------------------------------------------
    
    litterfall_tree     <- k_litterfall_ann * litterfall_prob_val * C_leaf_tree
    litterfall_herb     <- k_litterfall_herb_ann * litterfall_prob_val * C_leaf_herb
    
    leaf_mortality_tree <- k_mort_leaf_tree * C_leaf_tree
    leaf_mortality_herb <- k_mort_leaf_herb * C_leaf_herb
    
    # Root and wood winter dormancy:
    act <- winter_root_act_prop + 
      (1 - winter_root_act_prop) /
      (1 + exp(-k_root_dormancy * (T_t - root_dormancy_temp)))
    
    wood_mortality_tree <- k_mort_wood_tree * C_wood_tree*act
    
    root_mortality_herb <- k_mort_root_herb * C_root_herb*act
    root_mortality_tree <- k_mort_root_tree * C_root_tree*act
    
    
    exudates_herb       <- (k_exudate_intercept + RootHerb*k_exudate_slope)* C_root_herb*act
    exudates_tree       <- k_exudate_tree* C_root_tree*act
    
    # ----------------------------
    # Earthworm rates
    # ----------------------------
 
    Fed_earthworm_MetLitter = c_earthworm_litter*MetLitter*Earthworm
    
    Fed_earthworm_PASSIVE = c_earthworm_soil*Earthworm*PASSIVE
    
    Fed_earthworm_ACTIVE = c_earthworm_soil*Earthworm*ACTIVE
    
    Fed_earthworm_SLOW = c_earthworm_soil*Earthworm*SLOW

    Waste_earthworm_ACTIVE = prop_feaces_earthworm_LMWC*((1-a_earthworm)*(Fed_earthworm_litter) + (1-a_earthworm_soil)*(Fed_earthworm_PASSIVE + Fed_earthworm_ACTIVE + Fed_earthworm_SLOW))
    
    Carcass_earthworm_P = d_earthworm*Earthworm^2
    
    Waste_earthworm_A = (1-prop_feaces_earthworm_LMWC)*((1-a_earthworm)*(Fed_earthworm_litter) + (1-a_earthworm_soil)*(Fed_earthworm_PASSIVE + Fed_earthworm_ACTIVE + Fed_earthworm_SLOW))
    
    Respiration_earthworm = (1-p_earthworm)*(a_earthworm*(Fed_earthworm_MetLitter) + a_earthworm_soil*(Fed_earthworm_PASSIVE + Fed_earthworm_ACTIVE + Fed_earthworm_SLOW)) + E_earthworm*Earthworm
    
    # ----------------------------
    # Detritiviory rates
    # ----------------------------
    
    Fed_det_MetLitter = c_detritivores*MetLitter*Detritivore
    
    Carcass_det_SLOW = d_detritivores*Detritivore^2
    
    Waste_det_SLOW = (1-a_detritivores)*(Fed_det_MetLitter)
    
    Respiration_detritivore = (1-p_detritivores)*a_detritivores*(Fed_det_MetLitter)
    
    # ----------------------------
    # Predator rates
    # ----------------------------
    
    Fed_detpred_det = c_detpredator*Detritivore*DetPredator
    
    Carcass_detpred_om = d_detpredator*DetPredator^2
    
    Waste_detpred_om = (1-a_detpredator)*c_detpredator*Detritivore*DetPredator
    
    Respiration_detpred = (1-p_detpredator)*a_detpredator*c_detpredator*Detritivore*DetPredator
    
    # --------------------------------------------------
    # Root herbivores
    # --------------------------------------------------
    Fed_rootherb_herb = c_rootherb*C_root_herb*RootHerb
    
    Carcass_rootherb_P = d_rootherb*RootHerb^2
    
    Waste_rootherb_P = (1-a_rootherb)*Fed_rootherb_herb
    
    Respiration_rootherb = (1-p_rootherb)*a_rootherb*Fed_rootherb_herb
    
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
    dMetLitter = litterfall + leaf_mortality + exudates - f_MetLitter  - 
      Fed_det_MetLitter -
      Fed_earthworm_MetLitter
    
    #Equation B11
    dACTIVE <- (1 - root_to_organic)*root_mortality + (1-LigFrac) * strlitter_to_active * f_StrLitter + metlitter_to_active * f_MetLitter  + f_SLOW * slow_to_active + f_PASSIVE * passive_to_active - f_ACTIVE  -
      Fed_earthworm_ACTIVE
    
    #Equation B12
    dSLOW <-  LigFrac * strlitter_to_slow * f_StrLitter + f_ACTIVE * (1-f_TEX-active_to_passive) - f_SLOW -
      Fed_earthworm_SLOW -
      Carcass_det_SLOW -
      Waste_det_SLOW
    
    #Equation B13
    dPASSIVE <- f_ACTIVE * active_to_passive + f_SLOW * slow_to_passive - f_PASSIVE -
      Fed_earthworm_PASSIVE 

    # ---------------------------
    # Return list for deSolve
    # ---------------------------
    list(
      c(dStrLitter, dMetLitter, dACTIVE, dSLOW, dPASSIVE)
    )
  })
}


