# ------------------------------------------------------------
# make_tree_forcing.R
# Construct a vegetation forcing function from tree model
# ------------------------------------------------------------

make_tree_forcing <- function(parms) {
  
  # Return a function of time
  function(time) {
    
    forcing <- tree_forcing_monomolecular(
      time  = time,
      parms = parms
    )
    
    # Defensive checks (optional but recommended)
    stopifnot(
      is.numeric(forcing),
      !any(is.na(forcing))
    )
    
    forcing
  }
}

# ------------------------------------------------------------
# make_tree_forcing_equilibrium.R
# Construct constant tree forcing from climatological means
# ------------------------------------------------------------
make_tree_forcing_equilibrium <- function(parms,
                                          t_eq = 365 * 300,
                                          avg_days = 365) {
  
  # ---- 1. Times over which to average (tree already equilibrated) ----
  times_avg <- seq(t_eq, t_eq + avg_days - 1, by = 1)
  
  # ---- 2. Evaluate tree forcing over averaging window ----
  tf_mat <- vapply(
    times_avg,
    function(t) tree_forcing_monomolecular(t, parms),
    FUN.VALUE = tree_forcing_monomolecular(t_eq, parms)
  )
  # tf_mat: matrix with rows = forcing variables, cols = time
  
  # ---- 3. Take mean of ALL outputs ----
  tf_mean <- rowMeans(tf_mat)
  
  # ---- 4. Return constant forcing function ----
  function(time) {
    tf_mean
  }
}


# ------------------------------------------------------------
# make_herb_forcing.R
# Construct a vegetation forcing function from tree model
# ------------------------------------------------------------

make_herb_forcing <- function(parms) {
  
  # Return a function of time
  function(time) {
    
    forcing <- herb_forcing_monomolecular(
      time  = time,
      parms = parms
    )
    
    # Defensive checks (optional but recommended)
    stopifnot(
      is.numeric(forcing),
      !any(is.na(forcing))
    )
    
    forcing
  }
}

# ------------------------------------------------------------
# make_tree_forcing_equilibrium.R
# Construct constant tree forcing from climatological means
# ------------------------------------------------------------
make_herb_forcing_equilibrium <- function(parms,
                                          t_eq = 365 * 300,
                                          avg_days = 365) {
  
  # ---- 1. Times over which to average (tree already equilibrated) ----
  times_avg <- seq(t_eq, t_eq + avg_days - 1, by = 1)
  
  # ---- 2. Evaluate tree forcing over averaging window ----
  tf_mat <- vapply(
    times_avg,
    function(t) herb_forcing_monomolecular(t, parms),
    FUN.VALUE = herb_forcing_monomolecular(t_eq, parms)
  )
  # tf_mat: matrix with rows = forcing variables, cols = time
  
  # ---- 3. Take mean of ALL outputs ----
  tf_mean <- rowMeans(tf_mat)
  
  # ---- 4. Return constant forcing function ----
  function(time) {
    tf_mean
  }
}
