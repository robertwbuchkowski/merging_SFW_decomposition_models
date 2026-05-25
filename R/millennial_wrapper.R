# ============================================
# STATE CONFIGURATION SYSTEM
# ============================================

# Full list of state variables (order must match your model output)
full_names <- c(
  "C_leaf_herb", "C_root_herb",
  "C_leaf_tree", "C_wood_tree", "C_root_tree",
  "Earthworm", "Detritivore", "DetPredator","RootHerb",
  "Litter", "CWD", "Organic", "DOM", "MIC",
  "P", "L", "A", "M", "B"
)

# Group state variables by module
state_groups <- list(
  herb        = c("C_leaf_herb", "C_root_herb"),
  tree        = c("C_leaf_tree", "C_wood_tree", "C_root_tree"),
  earthworm   = c("Earthworm"),
  detritivore = c("Detritivore"),
  detpredator = c("DetPredator"),
  rootherb    = c("RootHerb")
)

# --------------------------------------------
# Get active state names based on config
# --------------------------------------------
get_active_names <- function(config){
  inactive <- unlist(state_groups[!unlist(config)])
  setdiff(full_names, inactive)
}

# --------------------------------------------
# Build reduced initial conditions
# --------------------------------------------
build_initial_state <- function(init_full, config){
  active <- get_active_names(config)
  init_full[active]
}

# --------------------------------------------
# Wrapper for reduced ODE system
# --------------------------------------------
millennial_wrapper <- function(time, state_reduced, parms){
  
  config <- parms$config
  
  # Recreate FULL state vector (missing vars = 0)
  full_state <- setNames(
    rep(0, length(full_names)),
    full_names
  )
  
  full_state[names(state_reduced)] <- state_reduced
  
  # Run full model
  d_full <- millennial_model_wplant(time, full_state, parms)[[1]]
  names(d_full) <- full_names
  
  # Return only active states
  active_names <- names(state_reduced)
  
  return(list(d_full[active_names]))
}