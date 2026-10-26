#------------------------------------------------------------------------#
# Updated Millennial v2 ODE ----
#------------------------------------------------------------------------#

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


#------------------------------------------------------------------------#
# Simplified pool structure ----
#------------------------------------------------------------------------#
# The separate organic-horizon pools of the earlier version were merged into
# the mineral-soil pools of Millennial v2:
#   Organic (fragmented litter)  -> P (particulate organic C, POM)
#   DOM (dissolved organic C)    -> L (low molecular weight C, LMWC)
#   MIC (organic-horizon microbes) -> B (microbial biomass)
# Litter and CWD therefore fragment straight into P, and soil microbes (B)
# depolymerize litter and CWD into L with the same enzyme kinetics they use
# on POM (alpha_pl, Ea_pl, K_pl). All dissolved C (litter/CWD
# depolymerization, root exudates, microbial necromass) enters L.

.POOLS_MILLENNIAL <- c(
  "C_root_herb", "C_root_tree",
  "Earthworm", "Detritivore", "RootHerb",
  "Litter", "CWD",
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
  K_pl                         <- parms[["K_pl"]]
  NPP_herb                     <- parms[["NPP_herb"]]
  NPP_tree                     <- parms[["NPP_tree"]]
  Rgas                         <- parms[["Rgas"]]
  T_ref                        <- parms[["T_ref"]]
  a_detritivores               <- parms[["a_detritivores"]]
  a_earthworm                  <- parms[["a_earthworm"]]
  a_earthworm_soil             <- parms[["a_earthworm_soil"]]
  a_leaf_tree                  <- parms[["a_leaf_tree"]]
  a_root_herb                  <- parms[["a_root_herb"]]
  a_root_tree                  <- parms[["a_root_tree"]]
  a_rootherb                   <- parms[["a_rootherb"]]
  a_wood_tree                  <- parms[["a_wood_tree"]]
  adj_detritivores             <- parms[["adj_detritivores"]]
  adj_earthworm                <- parms[["adj_earthworm"]]
  adj_rootherb                 <- parms[["adj_rootherb"]]
  alpha_lb                     <- parms[["alpha_lb"]]
  alpha_pl                     <- parms[["alpha_pl"]]
  c_detritivores               <- parms[["c_detritivores"]]
  c_earthworm_litter           <- parms[["c_earthworm_litter"]]
  c_earthworm_soil             <- parms[["c_earthworm_soil"]]
  c_rootherb                   <- parms[["c_rootherb"]]
  climate_forcing              <- parms[["climate_forcing"]]
  d_detritivores               <- parms[["d_detritivores"]]
  d_earthworm                  <- parms[["d_earthworm"]]
  d_rootherb                   <- parms[["d_rootherb"]]
  depth                        <- parms[["depth"]]
  k_a_min                      <- parms[["k_a_min"]]
  k_b                          <- parms[["k_b"]]
  k_b_slope_pint               <- parms[["k_b_slope_pint"]]
  k_bd                         <- parms[["k_bd"]]
  k_exudate_intercept          <- parms[["k_exudate_intercept"]]
  k_exudate_slope              <- parms[["k_exudate_slope"]]
  k_exudate_tree               <- parms[["k_exudate_tree"]]
  k_frag_CWD                   <- parms[["k_frag_CWD"]]
  k_frag_litter                <- parms[["k_frag_litter"]]
  k_l                          <- parms[["k_l"]]
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
  p_detritivores               <- parms[["p_detritivores"]]
  p_earthworm                  <- parms[["p_earthworm"]]
  p_rootherb                   <- parms[["p_rootherb"]]
  pct_claysilt                 <- parms[["pct_claysilt"]]
  phi_por                      <- parms[["phi_por"]]
  prop_feaces_earthworm_LMWC   <- parms[["prop_feaces_earthworm_LMWC"]]
  psi_matric                   <- parms[["psi_matric"]]
  root_dormancy_temp           <- parms[["root_dormancy_temp"]]
  slope_pint_det_k_frag_litter <- parms[["slope_pint_det_k_frag_litter"]]
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
  RootHerb   <- if ("RootHerb" %in% .ns) state[["RootHerb"]] else 0
  Litter     <- if ("Litter" %in% .ns) state[["Litter"]] else 0
  CWD        <- if ("CWD" %in% .ns) state[["CWD"]] else 0
  P          <- if ("P" %in% .ns) state[["P"]] else 0
  L          <- if ("L" %in% .ns) state[["L"]] else 0
  A          <- if ("A" %in% .ns) state[["A"]] else 0
  M          <- if ("M" %in% .ns) state[["M"]] else 0
  B          <- if ("B" %in% .ns) state[["B"]] else 0


    #------------------------------------------------------------------------#
    # Feeding-rate adjustment factors ----
    #------------------------------------------------------------------------#
    # Feeding-rate adjustment factors: one knob per animal scales ALL of
    # that animal's feeding coefficients together (default 1 = no change).
    # The animal fitting tunes these (adj_*) rather than the individual c_*
    # rates, so an animal that feeds on several pools keeps its relative
    # food preferences while its overall feeding rate goes up or down.
    c_earthworm_litter <- adj_earthworm    * c_earthworm_litter
    c_earthworm_soil   <- adj_earthworm    * c_earthworm_soil
    c_detritivores     <- adj_detritivores * c_detritivores
    c_rootherb         <- adj_rootherb     * c_rootherb
    #------------------------------------------------------------------------#
    # Get climate forcing ----
    #------------------------------------------------------------------------#
    forcing <- climate_forcing(time)
    T_t              <- forcing["Temp"] # °C
    theta_t          <- forcing["theta"] # m^3 m^-3
    root_input_weight  <- forcing["root_input_weight"]   # growing-season only (sums to 1/yr)
    wood_input_weight  <- forcing["wood_input_weight"]   # uniform over the year (sums to 1/yr)
    leaf_litter_weight <- forcing["leaf_litter_weight"]  # autumn peak + summer trickle (sums to 1/yr)

    #------------------------------------------------------------------------#
    # Shared growing-season activity (temperature + moisture) ----
    #------------------------------------------------------------------------#
    # Shared growing-season activity (temperature + moisture), in [0, 1].
    # This single activity index drives BOTH the seasonal allocation of NPP
    # and the winter dormancy of roots/wood, so the two follow the same logic.
    #   f_T_act : logistic temperature activity (warm -> 1, cold -> 0)
    #   f_theta : moisture limitation (theta / theta_opt, capped at 1)
    f_T_act  <- 1 / (1 + exp(-k_root_dormancy * (T_t - root_dormancy_temp)))
    f_theta  <- pmin(1, pmax(0, theta_t) / theta_opt)
    activity <- f_T_act * f_theta

    #------------------------------------------------------------------------#
    # Inputs from annual NPP + allocation + timing forcings ----
    #------------------------------------------------------------------------#
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

    #------------------------------------------------------------------------#
    # Root winter dormancy ----
    #------------------------------------------------------------------------#
    # Root winter dormancy: activity index with a winter floor.
    act <- winter_root_act_prop + (1 - winter_root_act_prop) * activity

    root_mortality_herb <- k_mort_root_herb * C_root_herb*act   # -> P
    root_mortality_tree <- k_mort_root_tree * C_root_tree*act   # -> P

    exudates_herb       <- (k_exudate_intercept + RootHerb*k_exudate_slope)* C_root_herb*act  # -> L
    exudates_tree       <- k_exudate_tree* C_root_tree*act                                  # -> L

    #------------------------------------------------------------------------#
    # Earthworm rates ----
    #------------------------------------------------------------------------#
    # Earthworms eat surface litter (c_earthworm_litter, assimilation
    # a_earthworm) and mineral soil: P, L, A, M (c_earthworm_soil, assimilation
    # a_earthworm_soil). Faeces go to L (prop_feaces_earthworm_LMWC) and A.

    Fed_earthworm_litter = c_earthworm_litter*Litter*Earthworm

    Fed_earthworm_M = c_earthworm_soil*Earthworm*M

    Fed_earthworm_P = c_earthworm_soil*Earthworm*P

    Fed_earthworm_L = c_earthworm_soil*Earthworm*L

    Fed_earthworm_A = c_earthworm_soil*Earthworm*A

    Fed_earthworm_soil = Fed_earthworm_M + Fed_earthworm_P + Fed_earthworm_L + Fed_earthworm_A

    Waste_earthworm = (1-a_earthworm)*Fed_earthworm_litter + (1-a_earthworm_soil)*Fed_earthworm_soil

    Waste_earthworm_L = prop_feaces_earthworm_LMWC*Waste_earthworm

    Waste_earthworm_A = (1-prop_feaces_earthworm_LMWC)*Waste_earthworm

    Carcass_earthworm_P = d_earthworm*Earthworm^2

    Respiration_earthworm = (1-p_earthworm)*(a_earthworm*Fed_earthworm_litter + a_earthworm_soil*Fed_earthworm_soil) + E_earthworm*Earthworm

    #------------------------------------------------------------------------#
    # Detritivory rates ----
    #------------------------------------------------------------------------#
    # Detritivores eat litter, particulate organic C (P) and microbes (B);
    # faeces and carcasses go to P.

    Fed_det_mic = c_detritivores*B*Detritivore

    Fed_det_P = c_detritivores*P*Detritivore

    Fed_det_lit = c_detritivores*Litter*Detritivore

    Fed_det = Fed_det_mic + Fed_det_P + Fed_det_lit

    Carcass_det_P = d_detritivores*Detritivore^2

    Waste_det_P = (1-a_detritivores)*Fed_det

    Respiration_detritivore = (1-p_detritivores)*a_detritivores*Fed_det

    #------------------------------------------------------------------------#
    # Root herbivores ----
    #------------------------------------------------------------------------#
    Fed_rootherb_herb = c_rootherb*C_root_herb*RootHerb

    Carcass_rootherb_P = d_rootherb*RootHerb^2

    Waste_rootherb_P = (1-a_rootherb)*Fed_rootherb_herb

    Respiration_rootherb = (1-p_rootherb)*a_rootherb*Fed_rootherb_herb

    #------------------------------------------------------------------------#
    # Fragmentation of surface detritus into POM ----
    #------------------------------------------------------------------------#
    fragmentation_litter   <- (k_frag_litter + k_frag_litter*slope_pint_det_k_frag_litter*Detritivore)  * Litter   # -> P
    fragmentation_CWD      <- k_frag_CWD     * CWD      # -> P

    #------------------------------------------------------------------------#
    # Sorption capacity Qmax (Eq. 11) ----
    #------------------------------------------------------------------------#
    Qmax <- depth * BD * pct_claysilt * p_c

    #------------------------------------------------------------------------#
    # Binding affinity for L sorption (Eq. 10) ----
    #------------------------------------------------------------------------#
    K_lm <- exp(-p1*pH - p2) * K_ld

    #------------------------------------------------------------------------#
    # Temperature functions (Arrhenius; Eqs. 3 & 14) ----
    #------------------------------------------------------------------------#
    V_pl <- alpha_pl * exp(-Ea_pl / (Rgas*(T_t + 273.15)))
    V_lb <- alpha_lb * exp(-Ea_lb / (Rgas*(T_t + 273.15)))

    #------------------------------------------------------------------------#
    # Moisture sensitivity ----
    #------------------------------------------------------------------------#
    # Diffusion limitation (Eq. 4)
    S_wD <- (theta_t / phi_por)^0.5

    # Biological moisture scalar (Eq. 15)
    oxygen_term <- k_a_min + (1 - k_a_min) * ((max(phi_por - theta_t, 0)) / phi_por)^0.5
    S_wB <- exp(lambda_mat * psi_matric) * oxygen_term * S_wD

    #------------------------------------------------------------------------#
    # Carbon use efficiency (Eqs. 21–22) ----
    #------------------------------------------------------------------------#
    CUE <- CUE_ref - CUE_T * (T_t - T_ref)
    CUE <- max(0, min(1, CUE))

    #------------------------------------------------------------------------#
    # KINETICS: depolymerization (F_pl, litter, CWD) & uptake (F_lb) ----
    #------------------------------------------------------------------------#
    # Microbial enzymes depolymerize POM (P), litter and CWD into L with the
    # same kinetics (V_pl, K_pl), each substrate on its own.
    if(kinetics == 1){
      # Reverse MM for depolymerization; forward MM for uptake
      F_pl      <- V_pl * S_wD * P      * B / (K_pl + B)
      F_Litter_L <- V_pl * S_wD * Litter * B / (K_pl + B)
      F_CWD_L    <- V_pl * S_wD * CWD    * B / (K_pl + B)
      F_lb <- V_lb * S_wB * B * L / (K_lb + L)

    } else if(kinetics == 2){
      # Equilibrium Chemistry Approximation (Eqs. 2a, 13a)
      F_pl      <- V_pl * S_wD * P      * B / (K_pl + B + P)
      F_Litter_L <- V_pl * S_wD * Litter * B / (K_pl + B + Litter)
      F_CWD_L    <- V_pl * S_wD * CWD    * B / (K_pl + B + CWD)
      F_lb <- V_lb * S_wB * B * L / (K_lb + L + B)

    } else if(kinetics == 3){
      # Linear kinetics (Eqs. 2b, 13b)
      F_pl      <- V_pl * S_wD * P      * B / K_pl
      F_Litter_L <- V_pl * S_wD * Litter * B / K_pl
      F_CWD_L    <- V_pl * S_wD * CWD    * B / K_pl
      F_lb <- V_lb * S_wB * B * L / K_lb

    } else {
      stop("Unknown kinetics option. Use 1, 2, or 3")
    }

    #------------------------------------------------------------------------#
    # Aggregation fluxes ----
    #------------------------------------------------------------------------#

    k_b_cur = pmax(0.001, k_b + Earthworm * k_b_slope_pint * k_b)

    F_pa <- k_pa * S_wD * P   # Eq. 5
    F_a  <- k_b_cur  * S_wD * A   # Eq. 6
    F_ma <- k_ma * S_wD * M   # Eq. 18

    #------------------------------------------------------------------------#
    # Sorption/desorption ----
    #------------------------------------------------------------------------#
    sat_term <- max(0, 1 - (M / max(Qmax, .Machine$double.eps)))
    F_lm <- S_wD * K_lm * L * sat_term        # Eq. 9
    F_ld <- K_ld * (M / max(Qmax, .Machine$double.eps))  # Eq. 12

    #------------------------------------------------------------------------#
    # Leaching ----
    #------------------------------------------------------------------------#
    F_l <- k_l * S_wD * L             # Eq. 8 (loss from the system)

    #------------------------------------------------------------------------#
    # Microbial mortality + respiration partitioning ----
    #------------------------------------------------------------------------#
    F_bm <- k_bd * B^2                # Eq. 16
    F_bg <- F_lb * CUE
    F_mr <- F_lb * (1 - CUE)

    #------------------------------------------------------------------------#
    # Plant differential equations ----
    #------------------------------------------------------------------------#

    dC_root_herb <- root_growth_herb - root_mortality_herb - exudates_herb - Fed_rootherb_herb

    dC_root_tree <- root_growth_tree - root_mortality_tree - exudates_tree

    #------------------------------------------------------------------------#
    # Animal differential equations ----
    #------------------------------------------------------------------------#

    # Earthworms:
    dEarthworm <-
      p_earthworm*(a_earthworm*Fed_earthworm_litter + a_earthworm_soil*Fed_earthworm_soil) -
      Carcass_earthworm_P -
      E_earthworm*Earthworm

    # Detritivore:
    dDetritivore <- p_detritivores*a_detritivores*Fed_det - Carcass_det_P

    # Root herbivores:
    dRootHerb <- p_rootherb*a_rootherb*Fed_rootherb_herb - Carcass_rootherb_P

    #------------------------------------------------------------------------#
    # Surface detritus pools ----
    #------------------------------------------------------------------------#

    dLitter  <- leaf_litter_input - F_Litter_L - fragmentation_litter - Fed_earthworm_litter - Fed_det_lit

    dCWD     <- cwd_input - F_CWD_L - fragmentation_CWD

    #------------------------------------------------------------------------#
    # Millennial model pools ----
    #------------------------------------------------------------------------#

    # Eq. 1 (POM): fragmented litter and CWD, dead roots, animal carcasses,
    # detritivore and root-herbivore faeces
    dP <- fragmentation_litter + fragmentation_CWD +
      root_mortality_tree + root_mortality_herb +
      p_a * F_a - F_pa - F_pl -
      Fed_earthworm_P + Carcass_earthworm_P -
      Fed_det_P + Carcass_det_P + Waste_det_P +
      Carcass_rootherb_P + Waste_rootherb_P

    # Eq. 7 (LMWC): depolymerized POM, litter and CWD, root exudates
    dL <- F_Litter_L + F_CWD_L + exudates_tree + exudates_herb -
      F_l + F_pl - F_lm - F_lb + (1 - p_b) * F_bm + F_ld - Fed_earthworm_L + Waste_earthworm_L

    # Eq. 17
    dA <- F_ma + F_pa - F_a - Fed_earthworm_A + Waste_earthworm_A

    # Eq. 19
    dM <- F_lm - F_ld + p_b * F_bm - F_ma + F_a * (1 - p_a) - Fed_earthworm_M

    # Eq. 20
    dB <- F_lb - F_bm - F_mr - Fed_det_mic

    # --------------------#
    # MASS BALANCE CHECK:
    # --------------------#
    mass_balance_check <- (
      dC_root_herb + dC_root_tree +
        dEarthworm + dDetritivore + dRootHerb +
        dLitter + dCWD +
        dP + dL + dA + dM + dB
    ) + (
      # add respiration losses back
        Respiration_earthworm +
        Respiration_detritivore +
        Respiration_rootherb +
        F_mr
    ) + (
      # add leaching losses back
      F_l
    ) - (
      # subtract external inputs (annual NPP delivered via the timing forcings)
      total_plant_input
    )

    #------------------------------------------------------------------------#
    # Return list for deSolve ----
    #------------------------------------------------------------------------#
    .dvec <- c(
        # Plant pools (roots only):
        dC_root_herb,
        dC_root_tree,

        # Animal pools:
        dEarthworm,
        dDetritivore,
        dRootHerb,

        # Surface detritus:
        dLitter,
        dCWD,

        # Mineral soil:
        dP,
        dL,
        dA,
        dM,
        dB)
    names(.dvec) <- .POOLS_MILLENNIAL
    list(.dvec[names(state)], mass_balance_check = unname(mass_balance_check))
}
