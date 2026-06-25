# ------------------------------------------------------------
# Updated Century  (with plants + animals + mass balance)
# ------------------------------------------------------------
# Author: Robert Buchkowski after Rose Abramoff, after Parton et al. (1987)
# Date:   Jan 25, 2026; Sep 11, 2021
# Updated: animal food web + closed mass balance added, abiotic
#          (Century soil) equations and plant equations unchanged.
#
# This function contains the system of equations for the Century model,
# developed in Parton et al. (1987). The abiotic equation numbers correspond
# to those in Abramoff et al. (2021) Appendix B, with parameters in Table A2.
#
# ------------------------------------------------------------------------
# MAPPING NOTE (Millennial -> Century pools for animal-derived carbon)
# ------------------------------------------------------------------------
# The Millennial animal sub-model deposits faeces / necromass into pools
# (L = LMWC, A = aggregate, P = POM, Organic) that do not exist in Century.
# Those destinations are remapped onto Century pools as follows. These are
# modelling choices that affect dynamics but NOT mass balance, so adjust
# them freely:
#   - labile faeces (LMWC analogue)        -> ACTIVE
#   - stabilised faeces (aggregate analogue) -> SLOW
#   - all necromass / carcasses            -> SLOW
# Animals feed from:
#   - earthworm: MetLitter (litter) + ACTIVE, SLOW, PASSIVE (soil)
#   - detritivore: MetLitter
#   - predator: detritivores
#   - root herbivore: herbaceous fine roots (C_root_herb)
# ------------------------------------------------------------------------

century_model <- function(time, state, parms){
  
  with(c(state, parms), {
    
    # ----------------------------
    # ---- Get climate forcing ----
    # ----------------------------
    forcing <- climate_forcing(time)
    T_t              <- forcing["Temp"] # °C
    theta_t          <- forcing["theta"] # m^3 m^-3
    litterfall_prob_val        <- forcing["litterfall_prob_val"] # probability of litterfall
    npp_weight_val             <- forcing["npp_weight_val"] # within-year NPP weight (sums to 1 over the year)
    
    # --------------------------------------------------
    # Shared growing-season activity (temperature + moisture), in [0, 1].
    # This single activity index drives BOTH the seasonal allocation of NPP
    # and the winter dormancy of roots/wood, so the two follow the same logic.
    #   f_T_act : logistic temperature activity (warm -> 1, cold -> 0)
    #   f_theta : moisture limitation (theta / theta_opt, capped at 1)
    # --------------------------------------------------
    f_T_act  <- 1 / (1 + exp(-k_root_dormancy * (T_t - root_dormancy_temp)))
    f_theta  <- pmin(1, pmax(0, theta_t) / theta_opt)
    activity <- f_T_act * f_theta

    # --------------------------------------------------
    # NPP input (the system input). NPP_herb / NPP_tree are ANNUAL parameters
    # (g C m-2 yr-1, from the scenarios file). npp_weight_val (from the climate
    # forcing) sums to 1 over the year and is proportional to the same activity
    # above, so the annual total delivered equals the parameter. No leaves ->
    # no NPP, which keeps a group off when its pools are zeroed.
    # --------------------------------------------------
    NPP_herb_ann <- NPP_herb
    NPP_tree_ann <- NPP_tree

    if (C_leaf_herb > 0) {
      NPP_herb <- NPP_herb_ann * npp_weight_val
    } else {
      NPP_herb <- 0
    }

    if (C_leaf_tree > 0) {
      NPP_tree <- NPP_tree_ann * npp_weight_val
    } else {
      NPP_tree <- 0
    }

    leaf_growth_herb <- (1 - a_root_herb) * NPP_herb
    root_growth_herb <- a_root_herb * NPP_herb

    if(abs(a_leaf_tree + a_wood_tree + a_root_tree - 1) > 1e-6) {stop("Tree allocation fractions must sum to 1")}

    leaf_growth_tree <- a_leaf_tree * NPP_tree
    wood_growth_tree <- a_wood_tree * NPP_tree
    root_growth_tree <- a_root_tree * NPP_tree

    # --------------------------------------------------
    # Continuous plant losses
    # --------------------------------------------------
    litterfall_tree     <- k_litterfall_ann * litterfall_prob_val * C_leaf_tree
    litterfall_herb     <- k_litterfall_herb_ann * litterfall_prob_val * C_leaf_herb

    leaf_mortality_tree <- k_mort_leaf_tree * C_leaf_tree
    leaf_mortality_herb <- k_mort_leaf_herb * C_leaf_herb

    # Root and wood winter dormancy: same activity index that shapes NPP,
    # with a winter floor (winter_root_act_prop) so some activity persists.
    act <- winter_root_act_prop + (1 - winter_root_act_prop) * activity

    wood_mortality_tree <- k_mort_wood_tree * C_wood_tree*act

    root_mortality_herb <- k_mort_root_herb * C_root_herb*act
    root_mortality_tree <- k_mort_root_tree * C_root_tree*act


    exudates_herb       <- (k_exudate_intercept + RootHerb*k_exudate_slope)* C_root_herb*act
    exudates_tree       <- k_exudate_tree* C_root_tree*act
    
    # --------------------------------------------------
    # Combined plant -> soil litter inputs
    # (aggregate herb + tree so the abiotic equations have single terms)
    # --------------------------------------------------
    litterfall      <- litterfall_tree     + litterfall_herb
    leaf_mortality  <- leaf_mortality_tree + leaf_mortality_herb
    wood_mortality  <- wood_mortality_tree
    root_mortality  <- root_mortality_tree + root_mortality_herb
    exudates        <- exudates_tree       + exudates_herb
    
    # ----------------------------
    # Earthworm rates
    #   litter source: MetLitter ; soil sources: ACTIVE, SLOW, PASSIVE
    # ----------------------------
    
    Fed_earthworm_MetLitter = c_earthworm_litter*MetLitter*Earthworm
    
    Fed_earthworm_ACTIVE  = c_earthworm_soil*Earthworm*ACTIVE
    Fed_earthworm_SLOW    = c_earthworm_soil*Earthworm*SLOW
    Fed_earthworm_PASSIVE = 0 #c_earthworm_soil*Earthworm*PASSIVE
    
    # Assimilated carbon (litter assimilation vs. soil assimilation efficiencies)
    Assim_earthworm = a_earthworm*Fed_earthworm_MetLitter +
      a_earthworm_soil*(Fed_earthworm_ACTIVE + Fed_earthworm_SLOW + Fed_earthworm_PASSIVE)
    
    # Egested (unassimilated) carbon -> faeces
    Egest_earthworm = (1 - a_earthworm)*Fed_earthworm_MetLitter +
      (1 - a_earthworm_soil)*(Fed_earthworm_ACTIVE + Fed_earthworm_SLOW + Fed_earthworm_PASSIVE)
    
    # Faeces split: labile -> ACTIVE, stabilised -> SLOW
    Waste_earthworm_ACTIVE = prop_feaces_earthworm_LMWC      * Egest_earthworm
    Waste_earthworm_SLOW   = (1 - prop_feaces_earthworm_LMWC) * Egest_earthworm
    
    # Necromass -> SLOW
    Carcass_earthworm_SLOW = d_earthworm*Earthworm^2
    
    # Respiration (growth/assimilation inefficiency + maintenance E)
    Respiration_earthworm = (1 - p_earthworm)*Assim_earthworm + E_earthworm*Earthworm
    
    # ----------------------------
    # Detritivory rates (feed on MetLitter)
    # ----------------------------
    
    Fed_det_MetLitter = c_detritivores*MetLitter*Detritivore
    
    Carcass_det_SLOW = d_detritivores*Detritivore^2                 # necromass -> SLOW
    
    Waste_det_SLOW = (1 - a_detritivores)*Fed_det_MetLitter         # faeces   -> SLOW
    
    Respiration_detritivore = (1 - p_detritivores)*a_detritivores*Fed_det_MetLitter
    
    # ----------------------------
    # Predator rates (feed on detritivores)
    # ----------------------------
    
    Fed_detpred_det = c_detpredator*Detritivore*DetPredator
    
    Carcass_detpred_SLOW = d_detpredator*DetPredator^2             # necromass -> SLOW
    
    Waste_detpred_SLOW = (1 - a_detpredator)*Fed_detpred_det       # faeces   -> SLOW
    
    Respiration_detpred = (1 - p_detpredator)*a_detpredator*Fed_detpred_det
    
    # --------------------------------------------------
    # Root herbivores (feed on herbaceous fine roots)
    # --------------------------------------------------
    Fed_rootherb_herb = c_rootherb*C_root_herb*RootHerb
    
    Carcass_rootherb_SLOW = d_rootherb*RootHerb^2                  # necromass -> SLOW
    
    Waste_rootherb_SLOW = (1 - a_rootherb)*Fed_rootherb_herb       # faeces   -> SLOW
    
    Respiration_rootherb = (1 - p_rootherb)*a_rootherb*Fed_rootherb_herb
    
    # ----------------------------------#
    # Abiotic scalars for decomposition: (UNCHANGED)
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
    
    # ---- Animal effects on fluxes ------
    
    # Detritivores increase fragmentation:
    f_MetLitter = pmax(0, f_MetLitter + f_MetLitter * slope_pint_det_k_frag_litter * Detritivore)
    
    # Earthworms slow down PASSIVE to ACTIVE transfer
    f_PASSIVE = pmax(0, f_PASSIVE + Earthworm * k_b_slope_pint * f_PASSIVE)
    
    
    # --------------------------------------------------
    # Plant differential equations:
    # --------------------------------------------------
    
    dC_leaf_herb <- leaf_growth_herb - litterfall_herb - leaf_mortality_herb
    dC_root_herb <- root_growth_herb - root_mortality_herb - exudates_herb - Fed_rootherb_herb
    
    dC_leaf_tree <- leaf_growth_tree - litterfall_tree - leaf_mortality_tree
    dC_wood_tree <- wood_growth_tree - wood_mortality_tree
    dC_root_tree <- root_growth_tree - root_mortality_tree - exudates_tree
    
    # --------------------------------------------------
    # Animal differential equations:
    # --------------------------------------------------
    
    # Earthworms:
    dEarthworm <- p_earthworm*Assim_earthworm -
      Carcass_earthworm_SLOW -
      E_earthworm*Earthworm
    
    # Detritivore:
    dDetritivore <- p_detritivores*a_detritivores*Fed_det_MetLitter -
      Carcass_det_SLOW -
      Fed_detpred_det
    
    # DetPredator:
    dDetPredator <- p_detpredator*a_detpredator*Fed_detpred_det -
      Carcass_detpred_SLOW
    
    # Root herbivores:
    dRootHerb <- p_rootherb*a_rootherb*Fed_rootherb_herb -
      Carcass_rootherb_SLOW
    
    # --------------------------------------------------
    # Soil (Century) differential equations:
    #   abiotic terms UNCHANGED; animal feeding subtracted,
    #   animal faeces / necromass added.
    # --------------------------------------------------
    
    #Equation B9
    dStrLitter = wood_mortality + root_to_organic*root_mortality - f_StrLitter
    
    #Equation B10
    dMetLitter = litterfall + leaf_mortality + exudates - f_MetLitter -
      Fed_det_MetLitter -
      Fed_earthworm_MetLitter
    
    #Equation B11
    dACTIVE <- (1 - root_to_organic)*root_mortality +
      (1 - LigFrac) * strlitter_to_active * f_StrLitter +
      metlitter_to_active * f_MetLitter +
      f_SLOW * slow_to_active + f_PASSIVE * passive_to_active - f_ACTIVE -
      Fed_earthworm_ACTIVE +
      Waste_earthworm_ACTIVE
    
    #Equation B12
    dSLOW <- LigFrac * strlitter_to_slow * f_StrLitter +
      f_ACTIVE * (1 - f_TEX - active_to_passive) - f_SLOW -
      Fed_earthworm_SLOW +
      Waste_earthworm_SLOW + Carcass_earthworm_SLOW +
      Carcass_det_SLOW + Waste_det_SLOW +
      Carcass_detpred_SLOW + Waste_detpred_SLOW +
      Carcass_rootherb_SLOW + Waste_rootherb_SLOW
    
    #Equation B13
    dPASSIVE <- f_ACTIVE * active_to_passive + f_SLOW * slow_to_passive - f_PASSIVE -
      Fed_earthworm_PASSIVE
    
    # --------------------#
    # MASS BALANCE CHECK:
    # --------------------#
    # Heterotrophic (decomposition) respiration = the fraction of each
    # decomposition flux not transferred to a receiving pool.
    Resp_soil <-
      f_StrLitter * (1 - (1 - LigFrac) * strlitter_to_active - LigFrac * strlitter_to_slow) +
      f_MetLitter * (1 - metlitter_to_active) +
      f_ACTIVE    * f_TEX +
      f_SLOW      * (1 - slow_to_active - slow_to_passive) +
      f_PASSIVE   * (1 - passive_to_active)
    
    mass_balance_check <- (
      # change in all stocks
      dC_leaf_herb + dC_root_herb + dC_leaf_tree + dC_wood_tree + dC_root_tree +
        dEarthworm + dDetritivore + dDetPredator + dRootHerb +
        dStrLitter + dMetLitter + dACTIVE + dSLOW + dPASSIVE
    ) + (
      # add respiration losses back
      Respiration_earthworm +
        Respiration_detritivore +
        Respiration_detpred +
        Respiration_rootherb +
        Resp_soil
    ) - (
      # subtract external inputs (no leaching term in this Century version)
      NPP_herb + NPP_tree
    )
    
    # browser()
    
    # ---------------------------
    # Return list for deSolve
    # ---------------------------
    list(
      c(
        # Plant pools:
        dC_leaf_herb,
        dC_root_herb,
        dC_leaf_tree,
        dC_wood_tree,
        dC_root_tree,
        
        # Animal pools:
        dEarthworm,
        dDetritivore,
        dDetPredator,
        dRootHerb,
        
        # Soil (Century) pools:
        dStrLitter,
        dMetLitter,
        dACTIVE,
        dSLOW,
        dPASSIVE
      ),
      mass_balance_check = unname(mass_balance_check)
    )
  })
}