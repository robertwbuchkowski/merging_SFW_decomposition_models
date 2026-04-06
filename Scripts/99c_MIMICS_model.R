# ------------------------------------------------------------
# Updated MIMICS model ODE
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

MIMICS_v1_ode <- function(time, state, parms){
  
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
    
    # ----------------------------#
    # Fragmentation and physical transfer to organic and mineral soil
    # ----------------------------#
    fragmentation_litter   <- k_frag_litter  * Litter   # -> Organic
    fragmentation_CWD      <- k_frag_CWD     * CWD      # -> Organic
    fragmentation_organic  <- k_frag_organic * Organic  # -> POM
    
    # ------------------------------#
    # MIMICS Model code:
    # ------------------------------#
    
    Vslope   <- array(0.063,dim=6)
    Vint     <- 5.47
    aV       <- 8e-6
    Vmax     <- exp(TSOI * Vslope + Vint) * aV
    
    Kslope   <- array(NA,dim=6)
    Kslope[1]<- 0.017 #META LIT to MIC_1
    Kslope[2]<- 0.027 #STRU LIT to MIC_1 
    Kslope[3]<- 0.017 #AVAI SOM to MIC_1 
    Kslope[4]<- 0.017 #META LIT to MIC_2
    Kslope[5]<- 0.027 #STRU LIT to MIC_2
    Kslope[6]<- 0.017 #AVAI SOM to MIC_2
    Kint     <- 3.19
    aK       <- 10
    Km       <- exp(Kslope * TSOI + Kint) * aK
    
    CUE        <- c(0.55, 0.25, 0.75, 0.35)  #for LITm and LITs entering MICr and MICK, respectively
    #ANPP strongly correlated with MAP
    Tao_MOD1 <- sqrt(ANPP[s]/100)  #basicaily standardize against NWT
    tao      <- c(5.2e-4*exp(0.3*fMET), 2.4e-4*exp(0.1*fMET))	
    tao      <- tao * Tao_MOD1
    
    #------NEW Parameters--------------
    fPHYS    <- c(0.3 * exp(1.3*fCLAY), 0.2 * exp(0.8*fCLAY)) 	#fraction to SOMp
    fCHEM    <- c(0.1 * exp(-3*fMET)  , 0.3 * exp(-3*fMET)  ) 	#fraction to SOMc
    fAVAI    <- 1- (fPHYS + fCHEM)
    desorb   <- 9e-4 * exp(-3*(sqrt(fCLAY))) #if modified by MIC!
    desorb   <- 3e-4 * exp(-4*(sqrt(fCLAY))) #if stand alone rate constant
    desorb   <- 1.5e-5 * exp(-1.5*(fCLAY))      #CHANGED FOR GLOBAL RUN!!!   
    
    k        <- 2.0    #2.0			#REDUCED FROM 3 TO 1, REDUCES TEXTURE EFFECTS ON SOMa decay
    a        <- 2.0    #2.2			#increased from 4.0 to 4.5
    
    cMAX     <- 1.4                    #ORIG 1.4 Maximum CHEM SOM scalar w/   0% Clay 
    cMIN     <- 1.2                    #ORIG 1.4 Minimum CHEM SOM scalar w/ 100% Clay 
    cSLOPE   <- cMIN - cMAX            #Slope of linear function of cSCALAR for CHEM SOM  
    
    pSCALAR  <- a * exp(-k*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
    
    #------------!!MODIFIERS AS IN MIMICS2_b!!---------------
    MOD1     <- c(10, 2, 10, 3, 3, 2) 
    MOD2     <- c( 8, 2 ,4 * pSCALAR, 2, 4, 6 * pSCALAR) 	
    
    VMAX     <- Vmax * MOD1 
    KM       <- Km / MOD2
    KO       <- c(4,4)      #scalar modifies Km of Oxidat	
    I        <- array(NA, dim=2)              #Litter inputs to MET/STR
    I[1]     <- (EST_LIT / depth) * fMET      #partitioned to layers
    I[2]     <- (EST_LIT / depth) * (1-fMET)
    
    
    
    
    
    
    
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
    #can make fluxes from CHEM a function of microbial biomass size?
    
    dLIT_1 = I[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
    dMIC_1 = CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])
    dSOM_1 = I[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb 
    
    dLIT_2 = I[2] * (1-FI[2]) - LITmin[2] - LITmin[4]
    dMIC_2 = CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  
    dSOM_2 = I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT
    
    dSOM_3  = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]
    
    list(c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3))
    
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
      c(dLitter, dCWD, dOrganic, dDOM, dMIC, dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3)
    )
  })
}
