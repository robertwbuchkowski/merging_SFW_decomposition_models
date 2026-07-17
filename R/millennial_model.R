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


.POOLS_MILLENNIAL <- c(
  "C_root_herb", "C_root_tree",
  "Earthworm", "Detritivore", "DetPredator", "RootHerb",
  "Litter", "CWD", "Organic", "DOM", "MIC",
  "P", "L", "A", "M", "B"
)

millennial_model_wplant <- function(time, state, parms){

  # Read the parameters this model uses directly from the parms list.
  BD                           <- parms[["BD"]]
  CUE_T                        <- parms[["CUE_T"]]
  CUE_ref                      <- parms[["CUE_ref"]]
  E_earthworm                  <- parms[["E_earthworm"]]
  Ea_lb                        <- parms[["Ea_lb"]]
  Ea_pl                        <- parms[["Ea_pl"]]
  K_lb                         <- parms[["K_lb"]]
  K_ld                         <- parms[["K_ld"]]
  K_ob                         <- parms[["K_ob"]]
  K_ol                         <- parms[["K_ol"]]
  K_pl                         <- parms[["K_pl"]]
  NPP_herb                     <- parms[["NPP_herb"]]
  NPP_tree                     <- parms[["NPP_tree"]]
  Rgas                         <- parms[["Rgas"]]
  T_ref                        <- parms[["T_ref"]]
  a_detpredator                <- parms[["a_detpredator"]]
  a_detritivores               <- parms[["a_detritivores"]]
  a_earthworm                  <- parms[["a_earthworm"]]
  a_earthworm_soil             <- parms[["a_earthworm_soil"]]
  a_leaf_tree                  <- parms[["a_leaf_tree"]]
  a_root_herb                  <- parms[["a_root_herb"]]
  a_root_tree                  <- parms[["a_root_tree"]]
  a_rootherb                   <- parms[["a_rootherb"]]
  a_wood_tree                  <- parms[["a_wood_tree"]]
  adj_detpredator              <- parms[["adj_detpredator"]]
  adj_detritivores             <- parms[["adj_detritivores"]]
  adj_earthworm                <- parms[["adj_earthworm"]]
  adj_rootherb                 <- parms[["adj_rootherb"]]
  alpha_lb                     <- parms[["alpha_lb"]]
  alpha_ob                     <- parms[["alpha_ob"]]
  alpha_ol                     <- parms[["alpha_ol"]]
  alpha_pl                     <- parms[["alpha_pl"]]
  c_detpredator                <- parms[["c_detpredator"]]
  c_detritivores               <- parms[["c_detritivores"]]
  c_earthworm_litter           <- parms[["c_earthworm_litter"]]
  c_earthworm_om               <- parms[["c_earthworm_om"]]
  c_earthworm_soil             <- parms[["c_earthworm_soil"]]
  c_rootherb                   <- parms[["c_rootherb"]]
  climate_forcing              <- parms[["climate_forcing"]]
  d_detpredator                <- parms[["d_detpredator"]]
  d_detritivores               <- parms[["d_detritivores"]]
  d_earthworm                  <- parms[["d_earthworm"]]
  d_rootherb                   <- parms[["d_rootherb"]]
  depth                        <- parms[["depth"]]
  k_MICd                       <- parms[["k_MICd"]]
  k_a_min                      <- parms[["k_a_min"]]
  k_b                          <- parms[["k_b"]]
  k_b_slope_pint               <- parms[["k_b_slope_pint"]]
  k_bd                         <- parms[["k_bd"]]
  k_exudate_intercept          <- parms[["k_exudate_intercept"]]
  k_exudate_slope              <- parms[["k_exudate_slope"]]
  k_exudate_tree               <- parms[["k_exudate_tree"]]
  k_frag_CWD                   <- parms[["k_frag_CWD"]]
  k_frag_litter                <- parms[["k_frag_litter"]]
  k_frag_organic               <- parms[["k_frag_organic"]]
  k_l                          <- parms[["k_l"]]
  k_l_o                        <- parms[["k_l_o"]]
  k_ma                         <- parms[["k_ma"]]
  k_mort_root_herb             <- parms[["k_mort_root_herb"]]
  k_mort_root_tree             <- parms[["k_mort_root_tree"]]
  k_pa                         <- parms[["k_pa"]]
  k_root_dormancy              <- parms[["k_root_dormancy"]]
  kinetics                     <- parms[["kinetics"]]
  lambda_mat                   <- parms[["lambda_mat"]]
  p1                           <- parms[["p1"]]
  p2                           <- parms[["p2"]]
  pH                           <- parms[["pH"]]
  p_a                          <- parms[["p_a"]]
  p_b                          <- parms[["p_b"]]
  p_c                          <- parms[["p_c"]]
  p_detpredator                <- parms[["p_detpredator"]]
  p_detritivores               <- parms[["p_detritivores"]]
  p_earthworm                  <- parms[["p_earthworm"]]
  p_rootherb                   <- parms[["p_rootherb"]]
  pct_claysilt                 <- parms[["pct_claysilt"]]
  phi_por                      <- parms[["phi_por"]]
  prop_feaces_earthworm_LMWC   <- parms[["prop_feaces_earthworm_LMWC"]]
  psi_matric                   <- parms[["psi_matric"]]
  root_dormancy_temp           <- parms[["root_dormancy_temp"]]
  root_to_organic              <- parms[["root_to_organic"]]
  slope_pint_det_k_frag_litter <- parms[["slope_pint_det_k_frag_litter"]]
  slope_pint_det_k_frag_organic <- parms[["slope_pint_det_k_frag_organic"]]
  theta_opt                    <- parms[["theta_opt"]]
  winter_root_act_prop         <- parms[["winter_root_act_prop"]]

  # Direct state reads. .ns is the set of pools actually passed in; any
  # pool NOT in .ns defaults to 0, so this SAME function works either on
  # the full state (via the wrapper) or on a scenario-specific reduced
  # state that omits inactive pools entirely.
  .ns <- names(state)
  C_root_herb <- if ("C_root_herb" %in% .ns) state[["C_root_herb"]] else 0
  C_root_tree <- if ("C_root_tree" %in% .ns) state[["C_root_tree"]] else 0
  Earthworm  <- if ("Earthworm" %in% .ns) state[["Earthworm"]] else 0
  Detritivore <- if ("Detritivore" %in% .ns) state[["Detritivore"]] else 0
  DetPredator <- if ("DetPredator" %in% .ns) state[["DetPredator"]] else 0
  RootHerb   <- if ("RootHerb" %in% .ns) state[["RootHerb"]] else 0
  Litter     <- if ("Litter" %in% .ns) state[["Litter"]] else 0
  CWD        <- if ("CWD" %in% .ns) state[["CWD"]] else 0
  Organic    <- if ("Organic" %in% .ns) state[["Organic"]] else 0
  DOM        <- if ("DOM" %in% .ns) state[["DOM"]] else 0
  MIC        <- if ("MIC" %in% .ns) state[["MIC"]] else 0
  P          <- if ("P" %in% .ns) state[["P"]] else 0
  L          <- if ("L" %in% .ns) state[["L"]] else 0
  A          <- if ("A" %in% .ns) state[["A"]] else 0
  M          <- if ("M" %in% .ns) state[["M"]] else 0
  B          <- if ("B" %in% .ns) state[["B"]] else 0


    # ----------------------------------------------------------------
    # Feeding-rate adjustment factors: one knob per animal scales ALL of
    # that animal's feeding coefficients together (default 1 = no change).
    # The animal fitting tunes these (adj_*) rather than the individual c_*
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
    root_input_weight  <- forcing["root_input_weight"]   # growing-season only (sums to 1/yr)
    wood_input_weight  <- forcing["wood_input_weight"]   # uniform over the year (sums to 1/yr)
    leaf_litter_weight <- forcing["leaf_litter_weight"]  # autumn peak + summer trickle (sums to 1/yr)

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
    # INPUTS FROM ANNUAL NPP + ALLOCATION + TIMING FORCINGS.
    # NPP_herb / NPP_tree are ANNUAL parameters (g C m-2 yr-1). We work directly
    # from the annual value and route each allocation to its destination with
    # its OWN within-year timing weight (each sums to 1 over the year, so the
    # annual input equals allocation x annual NPP). A group is switched off by
    # zeroing its ROOT pool (the remaining plant pool), which zeroes its NPP.
    #   leaves (1 - a_root_herb)[herb] + a_leaf_tree[tree] -> Litter, autumn peak
    #                                                         + summer trickle
    #   wood   a_wood_tree[tree]                            -> CWD,    uniform
    #   roots  a_root_herb[herb], a_root_tree[tree]         -> root pools,
    #                                                         growing season only
    # --------------------------------------------------
    NPP_herb_ann <- if (C_root_herb > 0) NPP_herb else 0
    NPP_tree_ann <- if (C_root_tree > 0) NPP_tree else 0

    if(abs(a_leaf_tree + a_wood_tree + a_root_tree - 1) > 1e-6) {stop("Tree allocation fractions must sum to 1")}

    # leaves -> Litter (autumn peak + small summer trickle)
    leaf_litter_input <- ((1 - a_root_herb) * NPP_herb_ann + a_leaf_tree * NPP_tree_ann) * leaf_litter_weight

    # wood -> CWD (evenly over the year)
    cwd_input         <- (a_wood_tree * NPP_tree_ann) * wood_input_weight

    # roots -> explicit root pools (growing season only)
    root_growth_herb  <- a_root_herb * NPP_herb_ann * root_input_weight
    root_growth_tree  <- a_root_tree * NPP_tree_ann * root_input_weight

    # total plant C entering the tracked pools this instant (for mass balance)
    total_plant_input <- leaf_litter_input + cwd_input + root_growth_herb + root_growth_tree

    # --------------------------------------------------
    # Root winter dormancy: activity index with a winter floor.
    # --------------------------------------------------
    act <- winter_root_act_prop + (1 - winter_root_act_prop) * activity

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

    dC_root_herb <- root_growth_herb - root_mortality_herb - exudates_herb - Fed_rootherb_herb

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
    dLitter  <- leaf_litter_input - F_Litter_DOM - fragmentation_litter - Fed_earthworm_litter - Fed_det_lit

    dCWD     <- cwd_input - F_CWD_DOM - fragmentation_CWD

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
      dC_root_herb + dC_root_tree +
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
      # subtract external inputs (annual NPP delivered via the timing forcings)
      total_plant_input
    )

    # browser()

    # ---------------------------
    # Return list for deSolve
    # ---------------------------
    .dvec <- c(
        # Plant pools (roots only):
        dC_root_herb,
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
        dB)
    names(.dvec) <- .POOLS_MILLENNIAL
    list(.dvec[names(state)], mass_balance_check = unname(mass_balance_check))
}