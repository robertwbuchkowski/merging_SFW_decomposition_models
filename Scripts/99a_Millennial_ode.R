# ------------------------------------------------------------
# Updated Millennial v2 ODE (no Leaf/Wood/Root plant states)
# ------------------------------------------------------------

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


millennial_v2_ode <- function(time, state, parms){
  
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
  
    # Other optional external biological fluxes used later
    detritivory_litter  <- if (exists("detritivory_litter"))  detritivory_litter  else 0
    detritivory_CWD     <- if (exists("detritivory_CWD"))     detritivory_CWD     else 0
    detritivory_organic <- if (exists("detritivory_organic")) detritivory_organic else 0
    
    faeces  <- if (exists("faeces"))  faeces  else 0
    carcass <- if (exists("carcass")) carcass else 0
    
    extra_mineral_input_test <- if (exists("extra_mineral_input_test")) extra_mineral_input_test else 0
    
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
    
    # -------------------------------
    # Above-mineral pools (no plant pools)
    # -------------------------------
    
    # Save net detritus inputs from vegetation (to modeled detritus/soil system)
    net_det_inputs <- litterfall + leaf_mortality + wood_mortality + root_mortality + exudates
    
    # Detritus pools
    dLitter  <- litterfall + leaf_mortality - F_Litter_DOM - fragmentation_litter - detritivory_litter
    dCWD     <- wood_mortality - F_CWD_DOM - fragmentation_CWD - detritivory_CWD
    dOrganic <- fragmentation_litter + fragmentation_CWD +
      root_to_organic * root_mortality +
      faeces_to_organic * faeces +
      carcass_to_organic * carcass -
      F_Organic_DOM - fragmentation_organic - detritivory_organic
    
    dDOM <- F_Litter_DOM + F_CWD_DOM + F_Organic_DOM + F_MIC_mortality - F_DOM_MIC - F_l_organic
    dMIC <- F_DOM_MIC - F_MIC_mortality - F_MIC_respiration
    
    # Transfer to mineral:
    Fi_t_part <- (1 - root_to_organic) * root_mortality +
      (1 - faeces_to_organic) * faeces +
      (1 - carcass_to_organic) * carcass +
      fragmentation_organic
    
    Fi_t_dissolved <- F_l_organic + exudates
    
    Fi_t <- Fi_t_part + Fi_t_dissolved + extra_mineral_input_test
    
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
    # Check system mass-balance (for modeled pools only)
    # ---------------------------
    dState <- dLitter + dCWD + dOrganic + dDOM + dMIC + dP + dL + dA + dM + dB
    
    # Inputs to modeled system: vegetation detritus fluxes + any explicit extra mineral input
    Inputs <- net_det_inputs + extra_mineral_input_test
    
    # Outputs from modeled system: respiratory CO2 + leaching losses + any explicit herbivory/harvest terms (if you want them counted)
    Outputs <- F_mr + F_MIC_respiration + F_l +
      leaf_harvest + wood_harvest + root_harvest +
      herbivory_leaf + herbivory_wood + herbivory_root
    
    # ---------------------------
    # Return list for deSolve
    # ---------------------------
    list(
      c(dLitter, dCWD, dOrganic, dDOM, dMIC, dP, dL, dA, dM, dB),
      c(
        # core diagnostics you already output
        F_pl=F_pl, F_lb=F_lb, F_pa=F_pa, F_a=F_a, F_ma=F_ma,
        F_lm=F_lm, F_ld=F_ld, F_l=F_l, F_bm=F_bm, F_bg=F_bg, F_mr=F_mr,
        Qmax=Qmax, K_lm=K_lm, S_wD=S_wD, S_wB=S_wB, CUE=CUE, T=T_t, theta=theta_t,
        Fi=Fi_t, p_i=p_i,
        
        # organic horizon diagnostics
        F_MIC_respiration=F_MIC_respiration,
        F_MIC_mortality=F_MIC_mortality,
        F_l_organic=F_l_organic,
        
        # new tree forcing diagnostics
        B_tree = TotalBiomassTree,
        net_det_inputs = net_det_inputs,
        
        # balance diagnostics
        dCO2 = F_mr + F_MIC_respiration,
        dState = dState, Inputs = Inputs, Outputs = Outputs
      )
    )
  })
}
