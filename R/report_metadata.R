# ============================================================
# REPORT METADATA  -  units and descriptions used by the scenario reports.
#
# The model integrates carbon pools in g C m-2; every derivative (dPool/dt) is
# in g C m-2 day-1. Times are days; one model year = 365 days. The tables below
# give a human-readable label + unit + one-line description for each pool,
# parameter, and forcing so a co-author can vet the numbers without the code.
#
# USAGE
#   source("R/report_metadata.R")
#   pool_units[["Litter"]]              # "g C m-2"
#   param_info("k_frag_litter")         # list(label, unit, desc)
#   forcing_info                        # data.frame of the forcing functions
# ============================================================

# ---- global unit strings -------------------------------------------------
UNIT_POOL       <- "g C m-2"
UNIT_FLUX       <- "g C m-2 day-1"
UNIT_RATE_DAY   <- "day-1"
UNIT_NPP_ANNUAL <- "g C m-2 yr-1"

# ---- pools: pretty label + unit ------------------------------------------
# (labels mirror R/compare_functions.R name_lookup; all pools are g C m-2)
pool_label <- c(
  C_root_herb = "Herbaceous Root C", C_root_tree = "Tree Root C",
  Earthworm = "Earthworms", Detritivore = "Detritivores",
  DetPredator = "Detritivore Predators", RootHerb = "Root Herbivores",
  Litter = "Litter", CWD = "Coarse Woody Debris", Organic = "Organic Matter",
  DOM = "Dissolved Organic Matter", MIC = "Microbial Biomass (Organic horizon)",
  P = "Particulate OC (POC)", L = "Low-MW C (LMWC)", A = "Aggregate C",
  M = "Mineral-Associated OC (MAOC)", B = "Microbial Biomass (Mineral soil)")

pool_units <- setNames(rep(UNIT_POOL, length(pool_label)), names(pool_label))

# which pools are animals (biomass) vs the soil/plant carbon they act on
animal_pools <- c("Earthworm", "Detritivore", "DetPredator", "RootHerb")
soil_pools   <- setdiff(names(pool_label), animal_pools)

# ---- parameters: label + unit + description ------------------------------
# One row per parameter that a reader is likely to inspect. Unknown params fall
# back to a generic entry via param_info().
param_meta <- list(
  # environmental drivers
  MAT          = list("Mean annual temperature", "deg C", "Constant temperature used for the equilibrium; annual mean of the seasonal forcing."),
  MAtheta      = list("Mean annual soil moisture", "m3 m-3", "Constant volumetric water content at equilibrium; annual mean of the seasonal forcing."),
  T_amp        = list("Temperature seasonal amplitude", "deg C", "Half-range of the seasonal temperature sinusoid (dynamic runs only)."),
  theta_amp    = list("Moisture seasonal amplitude", "m3 m-3", "Half-range of the seasonal moisture cycle (dynamic runs only)."),
  theta_opt    = list("Moisture optimum", "m3 m-3", "Soil moisture at which growing-season plant activity is maximal."),
  pH           = list("Soil pH", "-", "Enters the LMWC->MAOC sorption affinity."),
  # productivity + allocation
  NPP_herb     = list("Herbaceous NPP", UNIT_NPP_ANNUAL, "Annual herbaceous net primary productivity (system input)."),
  NPP_tree     = list("Tree NPP", UNIT_NPP_ANNUAL, "Annual tree net primary productivity (system input)."),
  a_leaf_tree  = list("Tree leaf allocation", "fraction", "Fraction of tree NPP to leaves (-> Litter). Leaf+wood+root = 1."),
  a_wood_tree  = list("Tree wood allocation", "fraction", "Fraction of tree NPP to wood (-> CWD)."),
  a_root_tree  = list("Tree root allocation", "fraction", "Fraction of tree NPP to roots (-> Tree Root C pool)."),
  a_root_herb  = list("Herb root allocation", "fraction", "Fraction of herb NPP to roots; the rest becomes leaf litter."),
  root_to_organic = list("Root death -> Organic fraction", "fraction", "Share of dead root C routed to Organic Matter (rest to POC)."),
  # plant turnover / exudation
  k_mort_root_herb = list("Herb root mortality", UNIT_RATE_DAY, "First-order herbaceous root turnover (x activity)."),
  k_mort_root_tree = list("Tree root mortality", UNIT_RATE_DAY, "First-order tree root turnover (x activity)."),
  k_exudate_intercept = list("Herb root exudation rate", UNIT_RATE_DAY, "Baseline root-exudate C flux per unit herb root (x activity)."),
  k_exudate_tree   = list("Tree root exudation rate", UNIT_RATE_DAY, "Baseline root-exudate C flux per unit tree root (x activity)."),
  k_root_dormancy  = list("Root dormancy steepness", "deg C-1", "Logistic steepness of temperature control on plant activity."),
  root_dormancy_temp = list("Root dormancy temperature", "deg C", "Temperature midpoint of the activity logistic."),
  winter_root_act_prop = list("Winter activity floor", "fraction", "Minimum plant activity retained in dormancy."),
  # fragmentation / detrital turnover
  k_frag_litter  = list("Litter fragmentation", UNIT_RATE_DAY, "First-order litter -> Organic fragmentation; modified by detritivores."),
  k_frag_CWD     = list("CWD fragmentation", UNIT_RATE_DAY, "First-order coarse woody debris -> Organic fragmentation."),
  k_frag_organic = list("Organic fragmentation", UNIT_RATE_DAY, "First-order Organic -> POC fragmentation; modified by detritivores."),
  # microbial + enzyme kinetics
  CUE_ref      = list("Carbon-use efficiency (ref)", "fraction", "Microbial CUE at the reference temperature T_ref."),
  CUE_T        = list("CUE temperature slope", "deg C-1", "Linear decline of CUE with temperature above T_ref."),
  T_ref        = list("Reference temperature", "deg C", "Temperature at which CUE = CUE_ref."),
  k_MICd       = list("Microbial turnover (organic)", "g C-1 m2 day-1", "Density-dependent microbial (organic-horizon) death."),
  k_bd         = list("Microbial turnover (mineral)", "g C-1 m2 day-1", "Density-dependent microbial (mineral-soil) death."),
  alpha_pl     = list("POC->LMWC pre-exponential", "g C m-2 day-1", "Arrhenius pre-exponential for depolymerization."),
  alpha_lb     = list("LMWC->microbe pre-exponential", "g C m-2 day-1", "Arrhenius pre-exponential for LMWC uptake."),
  alpha_ol     = list("Organic->DOM pre-exponential", "g C m-2 day-1", "Arrhenius pre-exponential for organic breakdown."),
  alpha_ob     = list("DOM->microbe pre-exponential", "g C m-2 day-1", "Arrhenius pre-exponential for DOM uptake."),
  Ea_pl        = list("POC->LMWC activation energy", "J mol-1", "Temperature sensitivity of depolymerization."),
  Ea_lb        = list("LMWC->microbe activation energy", "J mol-1", "Temperature sensitivity of LMWC uptake."),
  K_pl         = list("POC->LMWC half-saturation", UNIT_POOL, "Michaelis constant for depolymerization."),
  K_lb         = list("LMWC->microbe half-saturation", UNIT_POOL, "Michaelis constant for LMWC uptake."),
  K_ol         = list("Organic->DOM half-saturation", UNIT_POOL, "Michaelis constant for organic breakdown."),
  K_ob         = list("DOM->microbe half-saturation", UNIT_POOL, "Michaelis constant for DOM uptake."),
  K_ld         = list("LMWC->MAOC max sorption", UNIT_RATE_DAY, "Scales Langmuir sorption of LMWC to minerals."),
  # aggregate / mineral cycle
  k_pa         = list("POC->aggregate", UNIT_RATE_DAY, "Occlusion of POC into aggregates."),
  k_ma         = list("MAOC->aggregate", UNIT_RATE_DAY, "Occlusion of MAOC into aggregates."),
  k_b          = list("Aggregate breakdown", UNIT_RATE_DAY, "Aggregate disruption back to POC/MAOC; modified by earthworms."),
  k_l          = list("LMWC leaching/loss", UNIT_RATE_DAY, "First-order LMWC loss (leaching)."),
  k_l_o        = list("DOM->POC turnover", UNIT_RATE_DAY, "Dissolved organic matter routed to particulate."),
  k_a_min      = list("Minimum aeration factor", "fraction", "Floor on the oxygen/aeration limitation of microbial activity."),
  lambda_mat   = list("Matric-potential sensitivity", "MPa-1", "Sensitivity of microbial activity to matric potential."),
  psi_matric   = list("Matric potential", "MPa", "Soil matric potential used in the moisture stress term."),
  p1           = list("pH sorption coefficient 1", "-", "pH term in the LMWC->MAOC affinity."),
  p2           = list("pH sorption coefficient 2", "-", "Intercept term in the LMWC->MAOC affinity."),
  p_a          = list("Aggregate->POC partition", "fraction", "Fraction of disrupted aggregate C returned to POC (rest to MAOC)."),
  p_b          = list("Microbial necromass->MAOC", "fraction", "Fraction of mineral-microbial necromass routed to MAOC."),
  p_c          = list("Clay+silt protection coefficient", "-", "Scales maximum MAOC sorption capacity (Qmax)."),
  # soil physical
  depth        = list("Soil depth", "m", "Depth of the modelled soil layer; sets stocks per area and Qmax."),
  BD           = list("Bulk density", "kg m-3", "Soil bulk density; sets Qmax and porosity."),
  rho_p        = list("Particle density", "kg m-3", "Mineral particle density; sets porosity."),
  pct_claysilt = list("Clay + silt", "percent", "Fine-fraction texture; sets MAOC capacity (Qmax)."),
  # animal effect levers (indirect pathways)
  slope_pint_det_k_frag_litter  = list("Detritivore -> litter fragmentation", "g C-1 m2", "Increase in k_frag_litter per unit detritivore biomass."),
  slope_pint_det_k_frag_organic = list("Detritivore -> organic fragmentation", "g C-1 m2", "Increase in k_frag_organic per unit detritivore biomass."),
  k_b_slope_pint  = list("Earthworm -> aggregate breakdown", "g C-1 m2", "Change in k_b per unit earthworm biomass."),
  k_exudate_slope = list("Root herbivore -> exudation", "g C-1 m2 day-1", "Increase in root exudation per unit root-herbivore biomass.")
)

# resolve one parameter's metadata (label / unit / desc), with a safe fallback
param_info <- function(name) {
  m <- param_meta[[name]]
  if (is.null(m)) return(list(label = name, unit = "-", desc = "(no description on file)"))
  list(label = m[[1]], unit = m[[2]], desc = m[[3]])
}

# ---- forcing functions ---------------------------------------------------
# The climate forcing returns Temp, theta and three input-timing weights. Each
# weight integrates to 1 over the year, so (annual NPP x allocation x weight)
# has units g C m-2 day-1.
forcing_info <- data.frame(
  forcing = c("Temp", "theta", "root_input_weight", "wood_input_weight", "leaf_litter_weight"),
  unit    = c("deg C", "m3 m-3", "proportion", "proportion", "proportion"),
  description = c(
    "Soil temperature driving Arrhenius kinetics and plant activity.",
    "Volumetric soil moisture driving moisture-stress and plant activity.",
    "Within-year timing of NPP allocated to roots (growing-season weighted).",
    "Within-year timing of NPP allocated to wood -> CWD (uniform over the year).",
    "Within-year timing of NPP allocated to leaves -> Litter (autumn peak + summer trickle)."),
  stringsAsFactors = FALSE)
