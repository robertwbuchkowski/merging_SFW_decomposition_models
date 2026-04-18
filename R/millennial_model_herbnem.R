millennial_model_herbnem <- function(time, state, parms){
  
  with(c(state, parms), {
    # --------------------------------------------------
    # Time + climate
    # --------------------------------------------------
    doy <- (time %% 365) + 1
    
    Temp <- MAT + T_amp * sin(2 * pi * (doy - 110) / 365)
    theta <- MAtheta + theta_amp * cos(4 * pi * (doy - 110) / 365)
    
    # --------------------------------------------------
    # Climate scalars for plant productivity
    # --------------------------------------------------
    f_T     <- Q10 ^ ((Temp - Tref) / 10)
    f_theta <- pmin(1, theta / theta_opt)
    
    # --------------------------------------------------
    # Herbaceous plant carbon fluxes
    # --------------------------------------------------
    GPP_herb <- GPPmax_herb * f_T * f_theta
    
    Ra_herb <- maint_resp * (C_shoot + C_root) +
      growth_resp * GPP_herb
    
    NPP_herb <- pmax(0, GPP_herb - Ra_herb)
    
    shoot_growth <- (1 - a_root_herb) * NPP_herb
    root_growth  <- a_root_herb * NPP_herb
    
    # --------------------------------------------------
    # Seasonal litterfall (shoots)
    # --------------------------------------------------
    wrapped_pdf <- function(t, mu, sigma, L = 365, K = 1) {
      ks <- (-K):K
      rowSums(sapply(ks, function(k)
        dnorm(t - mu - k * L, sd = sigma)))
    }
    
    pdf_vals  <- wrapped_pdf(1:365, litter_peak_doy, litter_width_d)
    prob_vals <- pdf_vals / sum(pdf_vals)
    
    litterfall <- k_litterfall_herb_ann * prob_vals[doy] * C_shoot
    
    # --------------------------------------------------
    # Continuous plant losses
    # --------------------------------------------------
    leaf_mortality <- k_mort_leaf * C_shoot
    
    if(Temp < root_dormancy_temp){
      root_mortality <- k_mort_root * C_root*winter_root_act_prop
      exudates       <- (k_exudate_intercept + RootHerb*k_exudate_slope) * C_root*winter_root_act_prop
    }else{
      root_mortality <- k_mort_root * C_root
      exudates       <- (k_exudate_intercept + RootHerb*k_exudate_slope) * C_root
    }
    
    
    
    # --------------------------------------------------
    # Root herbivores
    # --------------------------------------------------
    ConsumpRootHerb = c_rootherb*C_root*RootHerb
    DeathRootHerb = d_rootherb*RootHerb^2
    FaecesRootHerb = (1-a_rootherb)*c_rootherb*C_root*RootHerb
    
    # --------------------------------------------------
    # Herbaceous pool ODEs
    # --------------------------------------------------
    dC_shoot <- shoot_growth - litterfall - leaf_mortality
    dC_root  <- root_growth  - root_mortality - exudates - ConsumpRootHerb
    
    # ==================================================
    # ========== ORIGINAL MILLENNIAL CORE ===============
    # ==================================================
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
    V_pl <- alpha_pl * exp(-Ea_pl / (Rgas*(Temp + 273.15)))
    V_lb <- alpha_lb * exp(-Ea_lb / (Rgas*(Temp + 273.15)))
    
    V_ol <- alpha_ol * exp(-Ea_pl / (Rgas*(Temp + 273.15)))
    V_ob <- alpha_ob * exp(-Ea_pl / (Rgas*(Temp + 273.15)))
    
    # ----------------------------
    # Moisture sensitivity
    # ----------------------------
    # Diffusion limitation (Eq. 4)
    S_wD <- (theta / phi_por)^0.5
    
    # Biological moisture scalar (Eq. 15)
    oxygen_term <- k_a_min + (1 - k_a_min) * ((max(phi_por - theta, 0)) / phi_por)^0.5
    S_wB <- exp(lambda_mat * psi_matric) * oxygen_term * S_wD
    
    # ----------------------------
    # Carbon use efficiency (Eqs. 21–22)
    # ----------------------------
    CUE <- CUE_ref - CUE_T * (Temp - T_ref)
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
    
    # Root herbivores:
    dRootHerb <- p_rootherb*a_rootherb*ConsumpRootHerb - DeathRootHerb
      
    # --------------------------------------------------
    # Return ODEs + diagnostics
    # --------------------------------------------------
    list(
      c(
        dC_shoot, dC_root,
        dLitter, dCWD, dOrganic, dDOM, dMIC,
        dP, dL, dA, dM, dB,
        dRootHerb
      )
    )
  })
}