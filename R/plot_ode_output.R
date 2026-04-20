plot_ode_output <- function(
    ode_out,
    variable_cols = NULL,
    time_col = 1,
    start_time = NULL
) {
  
  # --------------------------------------------------
  # Load required packages
  # --------------------------------------------------
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  # --------------------------------------------------
  # Convert to data frame
  # --------------------------------------------------
  df <- as.data.frame(ode_out)
  
  # Extract time variable
  time_name <- names(df)[time_col]
  names(df)[time_col] <- "time"
  
  # --------------------------------------------------
  # Identify state variables
  # --------------------------------------------------
  if(is.null(variable_cols)){
    state_names <- names(df)[-time_col]
  }else{
    state_names <- variable_cols
  }
  
  
  
  # --------------------------------------------------
  # Get starting values
  # --------------------------------------------------
  start_vals <- df[1, state_names, drop = FALSE] %>%
    pivot_longer(
      cols = everything(),
      names_to = "state",
      values_to = "start_value"
    )
  
  # --------------------------------------------------
  # Convert to long format for ggplot
  # --------------------------------------------------
  df_long <- df %>%
    pivot_longer(
      cols = all_of(state_names),
      names_to = "state",
      values_to = "value"
    )
  
  # --------------------------------------------------
  # Plot
  # --------------------------------------------------
  p <- ggplot(df_long, aes(x = time, y = value)) +
    geom_line() +
    geom_hline(
      data = start_vals,
      aes(yintercept = start_value),
      linetype = "dashed",
      color = "grey40"
    ) +
    facet_wrap(~state, scales = "free_y") +
    labs(
      x = "Time",
      y = "State variable value",
      title = "ODE model output over time"
    ) +
    theme_minimal()
  
  return(p)
}
