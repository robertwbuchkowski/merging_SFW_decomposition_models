single_parameter_sensitivity_millennial <- function(
    param_name,
    multipliers = c(0.5, 0.75, 1, 1.25, 1.5),
    y0 = init_millennial_state(),
    plot_results = TRUE,
    habitat_type = "forest") {
  
  # ---------------------------------------------------
  # Load required packages
  # ---------------------------------------------------
  library(deSolve)
  library(yaml)
  
  # ---------------------------------------------------
  # Load model & helper code
  # (assumes project structure you established)
  # ---------------------------------------------------
  source("R/tree_monomolecular.R")
  source("R/make_tree_forcing.R")
  source("R/millennial_model.R")
  source("R/load_config.R")
  source("R/derive_millennial_parms.R")
  
  # --------------------------------------------------
  # Load baseline parameters
  # --------------------------------------------------
  parms_base <- load_config("tree_monomolecular")
  parms_base <- modifyList(parms_base, yaml::read_yaml("config/millennial.yml"))
  parms_base <- derive_millennial_parms(parms_base)
  
  # Safety check
  if (!param_name %in% names(parms_base)) {
    stop("Parameter not found in MILLENNIAL parameters: ", param_name)
  }
  
  baseline_value <- parms_base[[param_name]]
  
  # -------------------------------------------------
  # Prepare results storage
  # -------------------------------------------------
  
  state_names <- names(init_millennial_state())
  
  results <- data.frame(
    multiplier = multipliers,
    parameter_value = NA_real_,
    Total = NA_real_,
    matrix(NA_real_,
           nrow = length(multipliers),
           ncol = length(state_names),
           dimnames = list(NULL, state_names)),
    stringsAsFactors = FALSE
  )
  
  cat("Running sensitivity for parameter:", param_name, "\n")
  
  # --------------------------------------------------
  # Sensitivity loop
  # --------------------------------------------------
  for (i in seq_along(multipliers)) {
    
    parms_test <- parms_base
    
    # Modify one parameter
    parms_test[[param_name]] <- baseline_value * multipliers[i]
    
    if(param_name == "a_root"){
      flo = 1- parms_test$a_root
      
      lprop = parms_test$a_leaf/(parms_test$a_leaf + parms_test$a_wood)
      
      parms_test$a_leaf = flo*lprop
      parms_test$a_wood = flo*(1-lprop)
      
    }
    
    # Recompute derived parameters
    parms_test <- derive_millennial_parms(parms_test)
    
    # Build equilibrium forcing
    if(habitat_type == "forest"){
      parms_test$tree_forcing <- make_tree_forcing_equilibrium(parms_test)
    }else{
      parms_test$tree_forcing <- make_herb_forcing_equilibrium(parms_test)
    }
    
    # Solve equilibrium
    eq <- rootSolve::stode(
      y     = y0,
      func  = millennial_model,
      parms = parms_test
    )
    
    # Store parameter value
    results$parameter_value[i] <- parms_test[[param_name]]
    
    # Store individual state variables
    results[i, state_names] <- eq$y[state_names]
    
    # Store the total C
    results[i, "Total"] <- sum(eq$y)
    
    cat(
      "  ×", multipliers[i], "\n"
    )
  }
  
  # -----------------------------------------------
  # Optional plot
  # -----------------------------------------------
  if (plot_results) {
    
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    
    # Identify state variable columns
    state_names <- names(y0)
    
    # Convert results from wide → long format
    results_long <- results %>%
      pivot_longer(
        cols = all_of(state_names),
        names_to = "state",
        values_to = "value"
      )
    
    # Plot
    print(ggplot(results_long,
                 aes(x = parameter_value,
                     y = value)) +
            geom_line() +
            geom_point() +
            facet_wrap(.~state, scales = "free_y") + 
            labs(
              x = paste(param_name, "value"),
              y = "Equilibrium pool size (g C m⁻²)",
              title = "Single-parameter sensitivity (equilibrium)"
            ) +
            theme_minimal())
  }
  
  return(results)
}
