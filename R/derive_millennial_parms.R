derive_millennial_parms <- function(parms) {

  # Use a tolerance instead of exact float equality: 0.075 + 0.675 + 0.25
  # is not exactly 1 in floating point, which made `!= 1` spuriously stop().
  if (abs(parms$a_leaf_tree + parms$a_wood_tree + parms$a_root_tree - 1) > 1e-6)
    stop("Tree biomass allocation must sum to 1.")

  parms$phi_por <- 1 - parms$BD / parms$rho_p

  return(parms)
}
