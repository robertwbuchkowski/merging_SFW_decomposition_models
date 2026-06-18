# Initial state for the Century model (g C m^-2)
# ORDER MUST MATCH century_model() return: plants, animals, soil.
# Plant + animal initial values mirror init_millennial_state();
# soil values are the original Century pool sizes.
init_century_state <- function() {
  c(
    # ---- Plant pools ----
    C_leaf_herb = 200,
    C_root_herb = 200,
    C_leaf_tree = 750,
    C_wood_tree = 6750,
    C_root_tree = 2500,

    # ---- Animal pools ----
    Earthworm   = 0.48,
    Detritivore = 0.1,
    DetPredator = 0.01,
    RootHerb    = 0.1,

    # ---- Soil (Century) pools ----
    StrLitter   = 1000,   # structural litter
    MetLitter   = 210.1,  # metabolic litter
    ACTIVE      = 2450,   # active SOM
    SLOW        = 3000,   # slow SOM
    PASSIVE     = 1000    # passive SOM
  )
}
