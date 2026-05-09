init_millennial_state <- function(Detritivore = F, Earthworm = F, HerbNem = F) {
  # Units: g C m^-2
  c(
    C_leaf_herb = 200,
    C_root_herb = 200,
    C_leaf_tree = 750,
    C_wood_tree = 6750,
    C_root_tree = 2500,
    
    Earthworm = 0.48,
    Detritivore = 0.1,
    RootHerb = 0.001,
    
    Litter  = 200,
    CWD     = 1000,
    Organic = 2000,
    DOM     = 10,
    MIC     = 0.1,
    P       = 400,
    L       = 10,
    A       = 1000,
    M       = 3000,
    B       = 40
  )
}
