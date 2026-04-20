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
    GPP_herb <- GPPmax_herb * f_T * f_theta
    
    Ra_herb <- maint_resp * (C_leaf_herb + C_root_herb) +
      growth_resp * GPP_herb
    
    NPP_herb <- pmax(0, GPP_herb - Ra_herb)
    
    leaf_growth_herb <- (1 - a_root_herb) * NPP_herb
    root_growth_herb  <- a_root_herb * NPP_herb
    
    # --------------------------------------------------
    # Tree plant carbon fluxes
    # --------------------------------------------------
    GPP_tree <- GPPmax_tree * f_T * f_theta
    
    Ra_tree <- maint_resp * (C_leaf_tree + C_wood_tree + C_root_tree) +
      growth_resp * GPP_tree
    
    NPP_tree <- pmax(0, GPP_tree - Ra_tree)
    
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
    if(T_t < root_dormancy_temp){
      wood_mortality_tree <- k_mort_wood_tree * C_wood_tree
      
      root_mortality_herb <- k_mort_root_herb * C_root_herb
      root_mortality_tree <- k_mort_root_tree * C_root_tree
      
      
      exudates_herb       <- k_exudate_herb* C_root_herb
      exudates_tree       <- k_exudate_tree* C_root_tree
    }else{
      wood_mortality_tree <- k_mort_wood_tree * C_wood_tree*winter_root_act_prop
      
      root_mortality_herb <- k_mort_root_herb * C_root_herb*winter_root_act_prop
      root_mortality_tree <- k_mort_root_tree * C_root_tree*winter_root_act_prop
      
      
      exudates_herb       <- k_exudate_herb* C_root_herb*winter_root_act_prop
      exudates_tree       <- k_exudate_tree* C_root_tree*winter_root_act_prop
    }
    
    # ----------------------------
    # Fragmentation and physical transfer to organic and mineral soil
    # ----------------------------
    fragmentation_litter   <- k_frag_litter  * Litter   # -> Organic
    fragmentation_CWD      <- k_frag_CWD     * CWD      # -> Organic
    fragmentation_organic  <- k_frag_organic * Organic  # -> POM
    
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
    F_pa <- k_pa * S_wD * P   # Eq. 5
    F_a  <- k_b  * S_wD * A   # Eq. 6
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
    dC_root_herb <- root_growth_herb - root_mortality_herb - exudates_herb
    
    dC_leaf_tree <- leaf_growth_tree - litterfall_tree - leaf_mortality_tree
    dC_wood_tree <- wood_growth_tree - wood_mortality_tree
    dC_root_tree <- root_growth_tree - root_mortality_tree - exudates_tree
    
    
    # -------------------------------
    # Organic horizon pools
    # -------------------------------
    
    # Detritus pools
    dLitter  <- litterfall_herb + leaf_mortality_herb + litterfall_tree + leaf_mortality_tree - F_Litter_DOM - fragmentation_litter
    
    dCWD     <- wood_mortality_tree - F_CWD_DOM - fragmentation_CWD
    
    dOrganic <- fragmentation_litter + fragmentation_CWD +
      root_to_organic * (root_mortality_tree + root_mortality_herb) +
      F_Organic_DOM - fragmentation_organic
    
    dDOM <- F_Litter_DOM + F_CWD_DOM + F_Organic_DOM + F_MIC_mortality - F_DOM_MIC - F_l_organic
    
    dMIC <- F_DOM_MIC - F_MIC_mortality - F_MIC_respiration
    
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
    dP <- p_i * Fi_t + p_a * F_a - F_pa - F_pl
    
    # Eq. 7
    dL <- Fi_t * (1 - p_i) - F_l + F_pl - F_lm - F_lb + (1 - p_b) * F_bm + F_ld
    
    # Eq. 17
    dA <- F_ma + F_pa - F_a
    
    # Eq. 19
    dM <- F_lm - F_ld + p_b * F_bm - F_ma + F_a * (1 - p_a)
    
    # Eq. 20
    dB <- F_lb - F_bm - F_mr
    
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
        dB)
    )
  })
}