# ------------------------------------------------------------
# derive_century_parms.R
# Century needs no derived quantities (all parameters are direct),
# but this gives a uniform setup pipeline across the three models
# and validates allocation, mirroring derive_millennial_parms().
# ------------------------------------------------------------
derive_century_parms <- function(parms) {

  if (abs(parms$a_leaf_tree + parms$a_wood_tree + parms$a_root_tree - 1) > 1e-6)
    stop("Tree biomass allocation must sum to 1.")

  return(parms)
}
