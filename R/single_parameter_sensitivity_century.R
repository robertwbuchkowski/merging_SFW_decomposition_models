single_parameter_sensitivity_century <- function(
    param_name,
    multipliers = c(0.5, 0.75, 1, 1.25, 1.5),
    y0 = init_century_state(),
    plot_results = TRUE
) {
  
  # ---------------------------------------------------
  # Load required packages
  # ---------------------------------------------------
  library(deSolve)
  library(yaml)
  
  # ---------------------------------------------------
  # Load model & helper code
  # ---------------------------------------------------
  source("R/tree_monomolecular.R")
  source("R/make_tree_forcing.R")
  source("R/century_model.R")
  source("R/load_config.R")
  source("R/init_century_state.R")
  
  # --------------------------------------------------
  # Load baseline parameters
  # --------------------------------------------------
  parms_base <- load_config("tree_monomolecular")
  parms_base <- modifyList(
    parms_base,
    yaml::read_yaml("config/century.yml")
  )
  
  # --------------------------------------------------
  # Safety check
  # --------------------------------------------------
  if (!param_name %in% names(parms_base)) {
    stop("Parameter not found in CENTURY parameters: ", param_name)
  }
  
  baseline_value <- parms_base[[param_name]]
  state_names <- names(y0)
  
  # -------------------------------------------------
  # Prepare results storage
  # -------------------------------------------------
  results <- data.frame(
    multiplier = multipliers,
    parameter_value = NA_real_,
    Total = NA_real_,
    matrix(
      NA_real_,
      nrow = length(multipliers),
      ncol = length(state_names),
      dimnames = list(NULL, state_names)
    ),
    stringsAsFactors = FALSE
  )
  
  cat("Running CENTURY sensitivity for parameter:", param_name, "\n")
  
  # --------------------------------------------------
  # Sensitivity loop
  # --------------------------------------------------
  for (i in seq_along(multipliers)) {
    
    parms_test <- parms_base
    parms_test[[param_name]] <- baseline_value * multipliers[i]
    
    parms_test$tree_forcing <- make_tree_forcing_equilibrium(parms_test)
    
    eq <- rootSolve::stode(
      y     = y0,
      func  = century_model,
      parms = parms_test
    )
    
    results$parameter_value[i] <- parms_test[[param_name]]
    results$Total[i] <- sum(eq$y)
    results[i, state_names] <- eq$y[state_names]
    
    cat("  ×", multipliers[i], "completed\n")
  }
  
  # -----------------------------------------------
  # Optional plot
  # -----------------------------------------------
  if (plot_results) {
    
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    
    results_long <- results %>%
      pivot_longer(
        cols = all_of(state_names),
        names_to = "state",
        values_to = "value"
      )
    
    print(
      ggplot(results_long,
             aes(x = parameter_value, y = value)) +
        geom_line() +
        geom_point() +
        facet_wrap(~state, scales = "free_y") +
        labs(
          x = paste(param_name, "value"),
          y = "Equilibrium pool size (g C m⁻²)",
          title = "CENTURY single-parameter sensitivity (equilibrium)"
        ) +
        theme_minimal()
    )
  }
  
  return(results)
}