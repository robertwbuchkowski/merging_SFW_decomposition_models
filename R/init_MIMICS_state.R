# Initial state for the MIMICS model (g C m^-2)
# ORDER MUST MATCH MIMICS_model() return: plants, animals, MIMICS pools.
init_MIMICS_state <- function() {
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

    # ---- MIMICS pools ----
    LIT_1 = 200,
    LIT_2 = 3000,
    MIC_1 = 20,
    MIC_2 = 20,
    SOM_1 = 400,
    SOM_2 = 1000,
    SOM_3 = 3000
  )
}
