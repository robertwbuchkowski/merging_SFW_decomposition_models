# --------------------------------------------
# Generic wrapper for any full ODE model
# --------------------------------------------
make_model_wrapper <- function(model_fun, full_names, state_groups = NULL) {
  
  # Return an actual wrapper function
  function(time, state_reduced, parms) {
    
    config <- parms$config
    
    # ---- Determine active states ----
    if (!is.null(state_groups) && !is.null(config)) {
      inactive <- unlist(state_groups[!unlist(config)])
      active_names <- setdiff(full_names, inactive)
    } else {
      active_names <- names(state_reduced)
    }
    
    # ---- Reconstruct full state vector ----
    full_state <- setNames(
      rep(0, length(full_names)),
      full_names
    )
    
    full_state[names(state_reduced)] <- state_reduced
    
    # ---- Run the full model ----
    d_full <- model_fun(time, full_state, parms)
    d_full1 <- d_full[[1]]
    names(d_full1) <- full_names
    
    # ---- Return only active states ----
    return(list(d_full1[active_names], mass_balance_check = d_full[[2]]))
  }
}