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
  
  if(all(abs(df$mass_balance_check) < 1e-6)){
    state_names <- subset(state_names, state_names != "mass_balance_check")
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

# ------------------------------------------------------------
# plot_followup_comparison(): overlay TWO deSolve runs on the same pools/time
# axis -- built for comparing the continued no-animal baseline (control)
# against the add-animals follow-up, but works for any two comparable runs.
#   control_out, treatment_out   deSolve `ode()` output matrices/data frames
#   control_label, treatment_label   legend labels for the two runs
# ------------------------------------------------------------
plot_followup_comparison <- function(control_out, treatment_out,
                                     control_label = "Continued baseline (no animals)",
                                     treatment_label = "Animals added",
                                     variable_cols = NULL, title = NULL, by = NULL) {
  library(ggplot2); library(tidyr); library(dplyr)

  prep <- function(ode_out, run_label) {
    df <- as.data.frame(ode_out)
    names(df)[1] <- "time"
    cols <- if (is.null(variable_cols)) setdiff(names(df), "time") else variable_cols
    cols <- setdiff(cols, "mass_balance_check")
    df %>%
      select(time, all_of(cols)) %>%
      pivot_longer(-time, names_to = "state", values_to = "value") %>%
      mutate(run = run_label)
  }

  df_long <- bind_rows(prep(control_out, control_label),
                       prep(treatment_out, treatment_label))
  
  if(!is.null(by)) df_long <- df_long %>% filter(time %in% seq(0,max(df_long$time), by =by))

  ggplot(df_long, aes(x = time, y = value, color = run)) +
    geom_line() +
    facet_wrap(~state, scales = "free_y") +
    labs(x = "Time", y = "State variable value", color = NULL,
        title = if (is.null(title)) "Follow-up: baseline continued vs. animals added" else title) +
    theme_minimal() +
    theme(legend.position = "top")
}

# ------------------------------------------------------------
# plot_followup_add(): convenience wrapper -- loads the saved
# "continue_baseline" and "add" follow-ups for a model/scenario (as written by
# save_followup() in Scripts/followup_analysis.R) and plots them together.
# ------------------------------------------------------------
plot_followup_add <- function(model, scenario, dir = "Data/followup",by = NULL, ...) {
  control   <- load_followup(model, scenario, "continue_baseline", dir = dir)
  treatment <- load_followup(model, scenario, "add", dir = dir)
  plot_followup_comparison(control$out, treatment$out,
                           title = sprintf("%s / %s: baseline continued vs. animals added",
                                           model, scenario), by = by,...)
}
