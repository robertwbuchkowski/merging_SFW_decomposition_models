# -------------------------------------------------------
# Updated Millennial v2 ODE
# --------------------------------------------------------

# Original Millennial model license:
# MIT License
# 
# Copyright (c) 2021 rabramoff
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#   
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


millennial_model_wplant <- function(time, state, parms){
  
  with(c(state, parms), {

    # ----------------------------------------------------------------
    # Feeding-rate adjustment factors: one knob per animal scales ALL of
    # that animal's feeding coefficients together (default 1 = no change).
    # The animal fitting tunes these (adj_*) instead of the individual c_*
    # rates, so an animal that feeds on several pools keeps its relative
    # food preferences while its overall feeding rate goes up or down.
    # ----------------------------------------------------------------
    c_earthworm_litter <- adj_earthworm    * c_earthworm_litter
    c_earthworm_soil   <- adj_earthworm    * c_earthworm_soil
    c_earthworm_om     <- adj_earthworm    * c_earthworm_om
    c_detritivores     <- adj_detritivores * c_detritivores
    c_detpredator      <- adj_detpredator  * c_detpredator
    c_rootherb         <- adj_rootherb     * c_rootherb
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
    
    # ----------------------------
    # Earthworm rates
    # ----------------------------
    
    Fed_earthworm_litter = c_earthworm_litter*Litter*Earthworm
    
    Fed_earthworm_om = c_earthworm_om*Organic*Earthworm
    
    Fed_earthworm_M = c_earthworm_soil*Earthworm*M
    
    Fed_earthworm_P = c_earthworm_soil*Earthworm*P
    
    Fed_earthworm_L = c_earthworm_soil*Earthworm*L
    
    Fed_earthworm_A = c_earthworm_soil*Earthworm*A
    
    
    Waste_earthworm_L = prop_feaces_earthworm_LMWC*((1-a_earthworm)*(Fed_earthworm_litter + Fed_earthworm_om) + (1-a_earthworm_soil)*(Fed_earthworm_M + Fed_earthworm_P + Fed_earthworm_L + Fed_earthworm_A))
    
    Carcass_earthworm_P = d_earthworm*Earthworm^2
    
    Waste_earthworm_A = (1-prop_feaces_earthworm_LMWC)*((1-a_earthworm)*(Fed_earthworm_litter + Fed_earthworm_om) + (1-a_earthworm_soil)*(Fed_earthworm_M + Fed_earthworm_P + Fed_earthworm_L + Fed_earthworm_A))
    
    Respiration_earthworm = (1-p_earthworm)*(a_earthworm*(Fed_earthworm_litter + Fed_earthworm_om) + a_earthworm_soil*(Fed_earthworm_M + Fed_earthworm_P + Fed_earthworm_L + Fed_earthworm_A)) + E_earthworm*Earthworm
    
    # ----------------------------
    # Detritiviory rates
    # ----------------------------
    
    Fed_det_mic = c_detritivores*MIC*Detritivore
    
    Fed_det_om = c_detritivores*Organic*Detritivore
    
    Fed_det_lit = c_detritivores*Litter*Detritivore
    
    Carcass_det_om = d_detritivores*Detritivore^2
    
    Waste_det_om = (1-a_detritivores)*(Fed_det_mic + Fed_det_om + Fed_det_lit)
    
    Respiration_detritivore = (1-p_detritivores)*a_detritivores*(Fed_det_mic + Fed_det_om + Fed_det_lit)
    
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
    
    # ----------------------------
    # Fragmentation and physical transfer to organic and mineral soil
    # ----------------------------
    fragmentation_litter   <- (k_frag_litter + k_frag_litter*slope_pint_det_k_frag_litter*Detritivore)  * Litter   # -> Organic
    fragmentation_CWD      <- k_frag_CWD     * CWD      # -> Organic
    fragmentation_organic  <- (k_frag_organic + k_frag_organic*slope_pint_det_k_frag_organic*Detritivore) * Organic  # -> POM
    
    # ----------------------------
    # Sorption capacity Qmax (Eq. 11)
    # ----------------------------
    Qmax <- depth * BD * pct_claysilt * p_c
    
    # ----------------------------
    # Binding affinity for L sorption (Eq. 10)
    # ----------------------------
    K_lm <- exp(-p1*pH - p2) * K_ld
    
    # ----------------------------
    # Temperature functions (Arrhenius; Eqs. 3 & 14)
    # ----------------------------
    V_pl <- alpha_pl * exp(-Ea_pl / (Rgas*(T_t + 273.15)))
    V_lb <- alpha_lb * exp(-Ea_lb / (Rgas*(T_t + 273.15)))
    
    V_ol <- alpha_ol * exp(-Ea_pl / (Rgas*(T_t + 273.15)))
    V_ob <- alpha_ob * exp(-Ea_pl / (Rgas*(T_t + 273.15)))
    
    # ----------------------------
    # Moisture sensitivity
    # ----------------------------
    # Diffusion limitation (Eq. 4)
    S_wD <- (theta_t / phi_por)^0.5
    
    # Biological moisture scalar (Eq. 15)
    oxygen_term <- k_a_min + (1 - k_a_min) * ((max(phi_por - theta_t, 0)) / phi_por)^0.5
    S_wB <- exp(lambda_mat * psi_matric) * oxygen_term * S_wD
    
    # ----------------------------
    # Carbon use efficiency (Eqs. 21–22)
    # ----------------------------
    CUE <- CUE_ref - CUE_T * (T_t - T_ref)
    CUE <- max(0, min(1, CUE))
    
    # ----------------------------
    # KINETICS: depolymerization (F_pl) & uptake (F_lb)
    # ----------------------------
    if(kinetics == 1){
      # Reverse MM for depolymerization; forward MM for uptake
      F_pl <- V_pl * S_wD * P * B / (K_pl + B)
      F_lb <- V_lb * S_wB * B * L / (K_lb + L)
      
      # Degradation in the organic horizon:
      F_Litter_DOM  <- V_ol * S_wD * Litter  * MIC / (K_ol + MIC)
      F_CWD_DOM     <- V_ol * S_wD * CWD     * MIC / (K_ol + MIC)
      F_Organic_DOM <- V_ol * S_wD * Organic * MIC / (K_ol + MIC)
      F_DOM_MIC     <- V_ob * S_wB * MIC * DOM / (K_ob + DOM)
      
    } else if(kinetics == 2){
      # Equilibrium Chemistry Approximation (Eqs. 2a, 13a)
      F_pl <- V_pl * S_wD * P * B / (K_pl + B + P)
      F_lb <- V_lb * S_wB * B * L / (K_lb + L + B)
      
      # Degradation in the organic horizon:
      F_Litter_DOM  <- V_ol * S_wD * Litter  * MIC / (K_ol + MIC + Litter)
      F_CWD_DOM     <- V_ol * S_wD * CWD     * MIC / (K_ol + MIC + CWD)
      F_Organic_DOM <- V_ol * S_wD * Organic * MIC / (K_ol + MIC + Organic)
      F_DOM_MIC     <- V_ob * S_wB * MIC * DOM / (K_ob + DOM + MIC)
      
    } else if(kinetics == 3){
      # Linear kinetics (Eqs. 2b, 13b)
      F_pl <- V_pl * S_wD * P * B / K_pl
      F_lb <- V_lb * S_wB * B * L / K_lb
      
      # Degradation in the organic horizon:
      F_Litter_DOM  <- V_ol * S_wD * Litter  * MIC / K_ol
      F_CWD_DOM     <- V_ol * S_wD * CWD     * MIC / K_ol
      F_Organic_DOM <- V_ol * S_wD * Organic * MIC / K_ol
      F_DOM_MIC     <- V_ob * S_wB * MIC * DOM / K_ob
      
    } else {
      stop("Unknown kinetics option. Use 1, 2, or 3")
    }
    
    # ----------------------------
    # Aggregation fluxes
    # ----------------------------
    
    k_b_cur = pmax(0.001, k_b + Earthworm * k_b_slope_pint * k_b)

    F_pa <- k_pa * S_wD * P   # Eq. 5
    F_a  <- k_b_cur  * S_wD * A   # Eq. 6
    F_ma <- k_ma * S_wD * M   # Eq. 18
    
    # ----------------------------
    # Sorption/desorption
    # ----------------------------
    sat_term <- max(0, 1 - (M / max(Qmax, .Machine$double.eps)))
    F_lm <- S_wD * K_lm * L * sat_term        # Eq. 9
    F_ld <- K_ld * (M / max(Qmax, .Machine$double.eps))  # Eq. 12
    
    # ----------------------------
    # Leaching
    # ----------------------------
    F_l <- k_l * S_wD * L             # Eq. 8
    F_l_organic <- k_l_o * S_wD * DOM # organic horizon leaching
    
    # ----------------------------
    # Microbial mortality + respiration partitioning
    # ----------------------------
    F_bm <- k_bd * B^2                # Eq. 16
    F_bg <- F_lb * CUE
    F_mr <- F_lb * (1 - CUE)
    
    # Organic horizon microbial mortality/resp
    F_MIC_mortality  <- k_MICd * MIC^2
    F_MIC_respiration <- F_DOM_MIC * (1 - CUE)
    
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
    dEarthworm <- 
      p_earthworm*(
        a_earthworm*(Fed_earthworm_litter + Fed_earthworm_om) + 
          a_earthworm_soil*(Fed_earthworm_M + Fed_earthworm_P + Fed_earthworm_L + Fed_earthworm_A)) - 
      d_earthworm*Earthworm^2 - 
      E_earthworm*Earthworm
    
    # Detritivore:
    dDetritivore <- p_detritivores*a_detritivores*(Fed_det_mic + Fed_det_om + Fed_det_lit) - Carcass_det_om - Fed_detpred_det
    
    # DetPredator:
    dDetPredator <- p_detpredator*a_detpredator*(Fed_detpred_det) - Carcass_detpred_om
    
    # Root herbivores:
    dRootHerb <- p_rootherb*a_rootherb*Fed_rootherb_herb - Carcass_rootherb_P
    
    # -------------------------------
    # Organic horizon pools
    # -------------------------------
    
    # Detritus pools
    dLitter  <- litterfall_herb + leaf_mortality_herb + litterfall_tree + leaf_mortality_tree - F_Litter_DOM - fragmentation_litter - Fed_earthworm_litter - Fed_det_lit
    
    dCWD     <- wood_mortality_tree - F_CWD_DOM - fragmentation_CWD
    
    dOrganic <- fragmentation_litter + fragmentation_CWD +
      root_to_organic * (root_mortality_tree + root_mortality_herb) -
      F_Organic_DOM - fragmentation_organic - 
      Fed_earthworm_om - 
      Fed_det_om + Carcass_det_om + Waste_det_om + 
      Carcass_detpred_om + Waste_detpred_om
    
    dDOM <- F_Litter_DOM + F_CWD_DOM + F_Organic_DOM + F_MIC_mortality - F_DOM_MIC - F_l_organic
    
    dMIC <- F_DOM_MIC - F_MIC_mortality - F_MIC_respiration - Fed_det_mic
    
    # Transfer to mineral:
    Fi_t_part <- (1 - root_to_organic) * (root_mortality_tree + root_mortality_herb) +
      fragmentation_organic
    
    Fi_t_dissolved <- F_l_organic + exudates_tree + exudates_herb
    
    Fi_t <- Fi_t_part + Fi_t_dissolved
    
    # Guard against Fi_t = 0
    Fi_safe <- max(Fi_t, .Machine$double.eps)
    p_i <- Fi_t_part / Fi_safe
    
    # -------------------------
    # Millennial model pools:
    # -------------------------
    
    # Eq. 1
    dP <- p_i * Fi_t + p_a * F_a - F_pa - F_pl - 
      Fed_earthworm_P + Carcass_earthworm_P + 
      Carcass_rootherb_P+ Waste_rootherb_P
    
    # Eq. 7
    dL <- Fi_t * (1 - p_i) - F_l + F_pl - F_lm - F_lb + (1 - p_b) * F_bm + F_ld - Fed_earthworm_L + Waste_earthworm_L
    
    # Eq. 17
    dA <- F_ma + F_pa - F_a - Fed_earthworm_A + Waste_earthworm_A
    
    # Eq. 19
    dM <- F_lm - F_ld + p_b * F_bm - F_ma + F_a * (1 - p_a) - Fed_earthworm_M
    
    # Eq. 20
    dB <- F_lb - F_bm - F_mr
    
    # --------------------#
    # MASS BALANCE CHECK:
    # --------------------#
    mass_balance_check <- (
      dC_leaf_herb + dC_root_herb + dC_leaf_tree + dC_wood_tree + dC_root_tree +
        dEarthworm + dDetritivore + dDetPredator + dRootHerb +
        dLitter + dCWD + dOrganic + dDOM + dMIC +
        dP + dL + dA + dM + dB
    ) + (
      # add respiration losses back
        Respiration_earthworm +
        Respiration_detritivore +
        Respiration_detpred +
        Respiration_rootherb +
        F_mr + F_MIC_respiration
    ) + (
      # add leaching losses back
      F_l
    ) - (
      # subtract external inputs
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
        
        # Organic horizons:
        dLitter, 
        dCWD, 
        dOrganic, 
        dDOM, 
        dMIC, 
        
        # Mineral horizons:
        dP, 
        dL, 
        dA, 
        dM, 
        dB),mass_balance_check = unname(mass_balance_check
    ))
  })
}