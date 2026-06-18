# ------------------------------------------------------------
# Updated MIMICS model ODE  (internal plants + soil animals + mass balance)
# ------------------------------------------------------------

# Original license for MIMICS code:
# The MIT License (MIT)
#
# Copyright (c) 2015 will wieder
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, ITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Modified by Robert Buchkowski
# Updated: internal plant pools (identical dynamics to the Millennial and
#          Century models) + soil-animal food web (earthworm, detritivore,
#          predator, root herbivore) with closed carbon mass balance.
#          The MIMICS abiotic equations are UNCHANGED. Plant litter is fed
#          into MIMICS's existing I[1]/I[2] inputs.
#
# ------------------------------------------------------------------------
# STATE VECTOR ORDER (deSolve matches by position):
#   Plants : C_leaf_herb, C_root_herb, C_leaf_tree, C_wood_tree, C_root_tree
#   Animals: Earthworm, Detritivore, DetPredator, RootHerb
#   MIMICS : LIT_1, LIT_2, MIC_1, MIC_2, SOM_1, SOM_2, SOM_3
#
# MAPPING NOTE (animal-derived carbon -> MIMICS pools; adjust freely):
#   Earthworm feeds on : LIT_1, LIT_2 (litter) + SOM_1, SOM_2, SOM_3 (soil)
#   Detritivore feeds on: LIT_1 + MIC_1, MIC_2 + SOM_3
#   Predator feeds on   : Detritivore
#   Root herbivore eats : living herbaceous fine roots (C_root_herb)
#   labile faeces            -> SOM_3 (available SOM)
#   earthworm casts (stable) -> SOM_1 (physically protected SOM)
#   all necromass / carcasses-> SOM_3 (available SOM)
#
# NOTE: the plant block is identical to century_model.R /
#       millennial_model_wplant. The external input to the whole system is
#       NPP_herb + NPP_tree; plant litter (I[1]+I[2]) is now an internal
#       plant -> soil transfer rather than an external forcing.
# ------------------------------------------------------------------------

MIMICS_model <- function(time, state, parms) {
  
  with(as.list(c(state, parms)), {
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
    
    # --------------------------------------------------
    # Combined plant -> soil litter inputs (aggregate herb + tree)
    # --------------------------------------------------
    litterfall      <- litterfall_tree     + litterfall_herb
    leaf_mortality  <- leaf_mortality_tree + leaf_mortality_herb
    wood_mortality  <- wood_mortality_tree
    root_mortality  <- root_mortality_tree + root_mortality_herb
    exudates        <- exudates_tree       + exudates_herb
    
    # ------------------------------#
    # MIMICS Model code:  -- abiotic UNCHANGED
    # ------------------------------#
    
    I <- c(
      litterfall + leaf_mortality + exudates,
      wood_mortality + root_mortality
    )
    
    # Calculate the proportion of litter that is metabolic
    if((I[1] + I[2]) > 0){
      fMET = I[1]/(I[1] + I[2])
    }else{
      fMET = 0.5
    }
    
    # ---- Temperature-dependent parameters ----
    Vmax <- exp(T_t * parms$Vslope + parms$Vint) * parms$aV
    Km   <- exp(parms$Kslope * T_t + parms$Kint) * parms$aK
    
    tao <- c(
      parms$tao_lit_coeff[1] * exp(parms$tao_lit_exp[1] * fMET),
      parms$tao_lit_coeff[2] * exp(parms$tao_lit_exp[2] * fMET)
    ) * Tao_MOD1
    
    
    fCHEM <- c(
      parms$fCHEM_coeff[1] * exp(-parms$fCHEM_exp[1] * fMET),
      parms$fCHEM_coeff[2] * exp(-parms$fCHEM_exp[2] * fMET)
    )
    
    fAVAI <- 1 - (fPHYS + fCHEM)
    
    # ---- Final kinetic matrices ----
    VMAX <- Vmax * parms$MOD1
    KM   <- Km / MOD2
    
    # ------------------------------#
    # MIMICS fluxes:  -- UNCHANGED
    # ------------------------------#
    
    LITmin = c(NA, NA)
    MICtrn = c(NA, NA, NA)
    SOMmin = c(NA, NA)
    
    #Flows to and from MIC_1
    LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + LIT_1)   #MIC_1 decomp of MET lit
    LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + LIT_2)   #MIC_1 decomp of STRUC lit
    MICtrn[1] = MIC_1 * tao[1]  * fPHYS[1]                  #MIC_1 turnover to PHYSICAL SOM
    MICtrn[2] = MIC_1 * tao[1]  * fCHEM[1]                  #MIC_1 turnover to CHEMICAL SOM
    MICtrn[3] = MIC_1 * tao[1]  * fAVAI[1]                  #MIC_1 turnover to AVAILABLE SOM
    SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3)   #decomp of SOMa by MIC_1
    
    #Flows to and from MIC_2
    LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + LIT_1)   #decomp of MET litter
    LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + LIT_2)   #decomp of SRUCTURAL litter
    MICtrn[4] = MIC_2 * tao[2]  * fPHYS[2]                  #MIC_2 turnover to PHYSICAL  SOM
    MICtrn[5] = MIC_2 * tao[2]  * fCHEM[2]                  #MIC_2 turnover to CHEMICAL  SOM
    MICtrn[6] = MIC_2 * tao[2]  * fAVAI[2]                  #MIC_2 turnover to AVAILABLE SOM
    SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)   #decomp of SOMa by MIC_2
    
    DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)		#desorbtion of PHYS to AVAIL (function of fCLAY)
    OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) +
                   (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2)))  #oxidation of C to A
    
    # ------------------------------#
    # Soil-animal food web
    # ------------------------------#
    
    # ---- Earthworm ----
    #   litter sources: LIT_1, LIT_2 (assimilation a_earthworm)
    #   soil sources  : SOM_1, SOM_2, SOM_3 (assimilation a_earthworm_soil)
    Fed_earthworm_LIT1 = c_earthworm_litter * LIT_1 * Earthworm
    Fed_earthworm_LIT2 = c_earthworm_litter * LIT_2 * Earthworm
    Fed_earthworm_SOM1 = c_earthworm_soil   * SOM_1 * Earthworm
    Fed_earthworm_SOM2 = c_earthworm_soil   * SOM_2 * Earthworm
    Fed_earthworm_SOM3 = c_earthworm_soil   * SOM_3 * Earthworm
    
    Fed_earthworm_litter_tot = Fed_earthworm_LIT1 + Fed_earthworm_LIT2
    Fed_earthworm_soil_tot   = Fed_earthworm_SOM1 + Fed_earthworm_SOM2 + Fed_earthworm_SOM3
    
    Assim_earthworm = a_earthworm      * Fed_earthworm_litter_tot +
      a_earthworm_soil * Fed_earthworm_soil_tot
    Egest_earthworm = (1 - a_earthworm)      * Fed_earthworm_litter_tot +
      (1 - a_earthworm_soil) * Fed_earthworm_soil_tot
    
    Waste_earthworm_SOM3 = prop_feaces_earthworm_LMWC       * Egest_earthworm  # labile faeces -> available SOM
    Waste_earthworm_SOM1 = (1 - prop_feaces_earthworm_LMWC) * Egest_earthworm  # casts        -> physical SOM
    Carcass_earthworm_SOM3 = d_earthworm * Earthworm^2                          # necromass    -> available SOM
    Respiration_earthworm  = (1 - p_earthworm) * Assim_earthworm + E_earthworm * Earthworm
    
    # ---- Detritivore ----
    #   sources: LIT_1 (litter), MIC_1, MIC_2 (microbial grazing), SOM_3 (available SOM)
    Fed_det_LIT1 = c_detritivores * LIT_1 * Detritivore
    Fed_det_MIC1 = c_detritivores * MIC_1 * Detritivore
    Fed_det_MIC2 = c_detritivores * MIC_2 * Detritivore
    Fed_det_SOM3 = c_detritivores * SOM_3 * Detritivore
    Fed_det_tot  = Fed_det_LIT1 + Fed_det_MIC1 + Fed_det_MIC2 + Fed_det_SOM3
    
    Carcass_det_SOM3 = d_detritivores * Detritivore^2
    Waste_det_SOM3   = (1 - a_detritivores) * Fed_det_tot
    Respiration_detritivore = (1 - p_detritivores) * a_detritivores * Fed_det_tot
    
    # ---- Predator ----
    Fed_detpred_det = c_detpredator * Detritivore * DetPredator
    Carcass_detpred_SOM3 = d_detpredator * DetPredator^2
    Waste_detpred_SOM3   = (1 - a_detpredator) * Fed_detpred_det
    Respiration_detpred  = (1 - p_detpredator) * a_detpredator * Fed_detpred_det
    
    # ---- Root herbivore (feeds on living herbaceous fine roots) ----
    Fed_rootherb_herb = c_rootherb * C_root_herb * RootHerb
    Carcass_rootherb_SOM3 = d_rootherb * RootHerb^2
    Waste_rootherb_SOM3   = (1 - a_rootherb) * Fed_rootherb_herb
    Respiration_rootherb  = (1 - p_rootherb) * a_rootherb * Fed_rootherb_herb
    
    # --------------------------------------------------
    # Plant differential equations: (identical to century / millennial)
    # --------------------------------------------------
    
    dC_leaf_herb <- leaf_growth_herb - litterfall_herb - leaf_mortality_herb
    dC_root_herb <- root_growth_herb - root_mortality_herb - exudates_herb - Fed_rootherb_herb
    
    dC_leaf_tree <- leaf_growth_tree - litterfall_tree - leaf_mortality_tree
    dC_wood_tree <- wood_growth_tree - wood_mortality_tree
    dC_root_tree <- root_growth_tree - root_mortality_tree - exudates_tree
    
    # --------------------------------------------------
    # Animal differential equations:
    # --------------------------------------------------
    
    dEarthworm <- p_earthworm * Assim_earthworm -
      Carcass_earthworm_SOM3 -
      E_earthworm * Earthworm
    
    dDetritivore <- p_detritivores * a_detritivores * Fed_det_tot -
      Carcass_det_SOM3 -
      Fed_detpred_det
    
    dDetPredator <- p_detpredator * a_detpredator * Fed_detpred_det -
      Carcass_detpred_SOM3
    
    dRootHerb <- p_rootherb * a_rootherb * Fed_rootherb_herb -
      Carcass_rootherb_SOM3
    
    # ------------------------------#
    # MIMICS ODEs (abiotic terms unchanged; animal terms added):
    # ------------------------------#
    
    dLIT_1 = I[1]*(1-FI[1]) - LITmin[1] - LITmin[3] -
      Fed_earthworm_LIT1 - Fed_det_LIT1
    
    dLIT_2 = I[2]*(1-FI[2]) - LITmin[2] - LITmin[4] -
      Fed_earthworm_LIT2
    
    dMIC_1 = CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3]) -
      Fed_det_MIC1
    
    dMIC_2 = CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6]) -
      Fed_det_MIC2
    
    dSOM_1 = I[1]*FI[1] + MICtrn[1] + MICtrn[4] - DEsorb -
      Fed_earthworm_SOM1 +
      Waste_earthworm_SOM1
    
    dSOM_2 = I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT -
      Fed_earthworm_SOM2
    
    dSOM_3  = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2] -
      Fed_earthworm_SOM3 - Fed_det_SOM3 +
      Waste_earthworm_SOM3 + Carcass_earthworm_SOM3 +
      Carcass_det_SOM3 + Waste_det_SOM3 +
      Carcass_detpred_SOM3 + Waste_detpred_SOM3 +
      Carcass_rootherb_SOM3 + Waste_rootherb_SOM3
    
    # --------------------#
    # MASS BALANCE CHECK:
    # --------------------#
    # MIMICS heterotrophic respiration = the (1-CUE) fraction of every
    # microbial uptake flux that is not retained as biomass.
    Resp_MIMICS <- (1 - CUE[1]) * (LITmin[1] + SOMmin[1]) +
      (1 - CUE[2]) *  LITmin[2] +
      (1 - CUE[3]) * (LITmin[3] + SOMmin[2]) +
      (1 - CUE[4]) *  LITmin[4]
    
    mass_balance_check <- (
      # change in all stocks
      dC_leaf_herb + dC_root_herb + dC_leaf_tree + dC_wood_tree + dC_root_tree +
        dEarthworm + dDetritivore + dDetPredator + dRootHerb +
        dLIT_1 + dLIT_2 + dMIC_1 + dMIC_2 + dSOM_1 + dSOM_2 + dSOM_3
    ) + (
      # add respiration losses back
      Resp_MIMICS +
        Respiration_earthworm +
        Respiration_detritivore +
        Respiration_detpred +
        Respiration_rootherb
    ) - (
      # subtract external inputs (NPP; no leaching in MIMICS)
      NPP_herb + NPP_tree
    )
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
        
        # MIMICS pools:
        dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3
      ),
      mass_balance_check = unname(mass_balance_check)
    )
  })
}