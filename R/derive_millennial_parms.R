derive_millennial_parms <- function(parms) {
  
  if(parms$a_leaf_tree + parms$a_wood_tree + parms$a_root_tree != 1) stop("Tree biomass allocation must sum to 1.")
  
  parms$phi_por <- 1 - parms$BD / parms$rho_p
  
  return(parms)
}