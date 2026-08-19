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

# ============================================================
# STACKED CHANGE-OVER-TIME PLOT (added animals)
# ------------------------------------------------------------
# For each scenario, the animal EFFECT over the follow-up = the added-animals
# run minus the time-matched continued-baseline (no-animal) run, per pool.
# Pools are relabelled and ordered from the shared pool_names vector (in
# R/compare_functions.R) and stacked so you see the net change and how each
# pool contributes through time. One facet per scenario. Annual steps (by = 365).
# ============================================================

# ------------------------------------------------------------
# followup_change_long(): change (added - continued baseline) per pool over
# time, for ONE scenario, at annual (or `by`-day) steps. Returns a tidy frame
# with raw `state` and the relabelled `pretty` (from pool_names). Animal pools
# are dropped so the stack shows only the soil/plant carbon they affect.
# ------------------------------------------------------------
followup_change_long <- function(model, scenario, dir = "Data/followup",
                                 by = 365) {
  add  <- load_followup(model, scenario, "add", dir = dir)$out
  ctrl <- load_followup(model, scenario, "continue_baseline", dir = dir)$out
  a <- as.data.frame(add);  names(a)[1] <- "time"
  b <- as.data.frame(ctrl); names(b)[1] <- "time"

  drop <- c("time", "mass_balance_check", animal_pool_names)
  keep <- setdiff(intersect(names(a), names(b)), drop)
  tt   <- intersect(a$time, b$time)
  if (!is.null(by)) tt <- tt[tt %in% seq(0, max(tt), by = by)]
  a <- a[a$time %in% tt, , drop = FALSE]
  b <- b[b$time %in% tt, , drop = FALSE]

  ch <- a[, keep, drop = FALSE] - b[match(a$time, b$time), keep, drop = FALSE]
  data.frame(
    time     = rep(a$time, times = length(keep)),
    state    = rep(keep, each = nrow(a)),
    change   = as.vector(as.matrix(ch)),
    pretty   = rep(relabel_pools(keep), each = nrow(a)),
    scenario = scenario,
    stringsAsFactors = FALSE)
}

# ------------------------------------------------------------
# plot_followup_stacked(): the faceted stacked graphic. Stacks the per-pool
# change (positive above 0, negative below), one facet per scenario.
# `scenarios` defaults to every scenario with saved runs.
# ------------------------------------------------------------
plot_followup_stacked <- function(model,
                                  scenarios = NULL,
                                  dir = "Data/followup",
                                  by = 365) {
  library(ggplot2); library(dplyr)

  if (is.null(scenarios)) {
    fs <- list.files(dir, pattern = paste0("^", model, "_.*_add\\.rds$"))
    scenarios <- sub(paste0("^", model, "_(.*)_add\\.rds$"), "\\1", fs)
    scenarios <- scenarios[file.exists(file.path(dir,
                     sprintf("%s_%s_continue_baseline.rds", model, scenarios)))]
  }
  if (!length(scenarios)) stop("No scenarios with saved add + continue_baseline runs in ", dir)

  dat <- bind_rows(lapply(scenarios, function(s)
    followup_change_long(model, s, dir = dir, by = by)))

  # sum change within each pool x time x scenario for the stack, ordered by pool_names
  stacked <- dat %>%
    group_by(scenario, pretty, time) %>%
    summarise(change = sum(change, na.rm = TRUE), .groups = "drop") %>%
    mutate(group = factor(pretty, levels = pool_order),
           years = time / 365)

  # net line (total change across all grouped pools)
  net <- stacked %>% group_by(scenario, years) %>%
    summarise(change = sum(change), .groups = "drop")

  ggplot(stacked, aes(x = years, y = change, fill = group)) +
    geom_area(alpha = 0.9, colour = "white", linewidth = 0.1) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    geom_line(data = net, aes(x = years, y = change),
              inherit.aes = FALSE, linewidth = 0.6, linetype = "dashed") +
    facet_wrap(~scenario, scales = "free_y") +
    scale_fill_viridis_d() +
    labs(x = "Years after animals added", y = expression("Animal Effect (g C m"^-2*")"),
         fill = "Pool") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")
}

# ============================================================
# FOLLOW-UP TRAJECTORY + EQUILIBRIUM END-STATE (single combined figure)
# ------------------------------------------------------------
# Shows the stacked animal-effect trajectory over the follow-up AND the
# equilibrium animal effect (from Scripts/2_spinup_dynamic.R ->
# Results/animal_eq_effect.csv) as a stacked bar at the right, past a broken
# x-axis. Pool fills match between the trajectory and the equilibrium bar; a
# dashed horizontal line marks the NET equilibrium change per scenario. This
# replaces the need for a separate equilibrium-effect figure.
#
#   eq_csv   the equilibrium effect table written by script 2. Columns used:
#            name, difference (absolute change, g C m-2), scenario, type.
#   eq_type  which effect to show at equilibrium ("total" or "direct").
# ============================================================
read_eq_change <- function(eq_csv = "Results/animal_eq_effect.csv",
                           eq_type = "total") {
  if (!file.exists(eq_csv)) {
    warning("equilibrium effect file not found: ", eq_csv,
            " -- run Scripts/2_spinup_dynamic.R. Equilibrium bar omitted.")
    return(NULL)
  }
  e <- utils::read.csv(eq_csv, stringsAsFactors = FALSE)
  e <- e[e$type == eq_type & !is.na(e$difference) &
         !(e$name %in% animal_pool_names), , drop = FALSE]
  data.frame(scenario = e$scenario,
             pretty   = relabel_pools(e$name),
             change   = e$difference,
             stringsAsFactors = FALSE)
}

plot_followup_with_eq <- function(model,
                                  scenarios = NULL,
                                  dir = "Data/followup",
                                  by = 365,
                                  eq_csv = "Results/animal_eq_effect.csv",
                                  eq_type = "total",
                                  gap_frac = 0.06,   # x-gap width as frac of the time span
                                  bar_frac = 0.14) { # eq-bar width as frac of the time span
  library(ggplot2); library(dplyr)

  if (is.null(scenarios)) {
    fs <- list.files(dir, pattern = paste0("^", model, "_.*_add\\.rds$"))
    scenarios <- sub(paste0("^", model, "_(.*)_add\\.rds$"), "\\1", fs)
    scenarios <- scenarios[file.exists(file.path(dir,
                     sprintf("%s_%s_continue_baseline.rds", model, scenarios)))]
  }
  if (!length(scenarios)) stop("No scenarios with saved add + continue_baseline runs in ", dir)

  # --- trajectory (stacked area), reuse the shared grouping/relabelling ---
  traj <- bind_rows(lapply(scenarios, function(s)
    followup_change_long(model, s, dir = dir, by = by))) %>%
    group_by(scenario, pretty, time) %>%
    summarise(change = sum(change, na.rm = TRUE), .groups = "drop") %>%
    mutate(years = time / 365)

  tmax <- max(traj$years)
  gap  <- gap_frac * tmax
  barw <- bar_frac * tmax
  x_bar_centre <- tmax + gap + barw / 2          # where the equilibrium bar sits

  # --- equilibrium end-state (stacked bar) ---
  eq <- read_eq_change(eq_csv, eq_type)
  if (!is.null(eq)) eq <- eq %>% filter(scenario %in% scenarios)

  # consistent factor levels / fill order across both layers
  lv <- pool_order
  lv <- lv[lv %in% unique(traj$pretty)]
  traj$group <- factor(traj$pretty, levels = lv)

  # net trajectory + net equilibrium (dashed reference line)
  net_traj <- traj %>% group_by(scenario, years) %>%
    summarise(change = sum(change), .groups = "drop")
  net_eq <- if (!is.null(eq)) eq %>% group_by(scenario) %>%
    summarise(net = sum(change), .groups = "drop") else NULL

  p <- ggplot(traj, aes(years, change, fill = group)) +
    geom_area(alpha = 0.9, colour = "white", linewidth = 0.1) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    geom_line(data = net_traj, aes(years, change),
              inherit.aes = FALSE, linewidth = 0.6, linetype = "dashed")

  if (!is.null(eq) && nrow(eq)) {
    eq$group <- factor(relabel_pools(eq$pretty), levels = lv)
    eq$xmin  <- x_bar_centre - barw / 2
    eq$xmax  <- x_bar_centre + barw / 2
    # stack the equilibrium bar: cumulative positive and negative separately
    eq <- eq %>% group_by(scenario) %>%
      arrange(group, .by_group = TRUE) %>%
      mutate(pos = pmax(change, 0), neg = pmin(change, 0),
             ytop_pos = cumsum(pos), ybot_pos = ytop_pos - pos,
             ybot_neg = cumsum(neg), ytop_neg = ybot_neg - neg) %>%
      ungroup()

    p <- p +
      # positive segments
      geom_rect(data = eq %>% filter(change > 0),
                aes(xmin = xmin, xmax = xmax, ymin = ybot_pos, ymax = ytop_pos, fill = group),
                inherit.aes = FALSE, colour = "white", linewidth = 0.1) +
      # negative segments
      geom_rect(data = eq %>% filter(change < 0),
                aes(xmin = xmin, xmax = xmax, ymin = ytop_neg, ymax = ybot_neg, fill = group),
                inherit.aes = FALSE, colour = "white", linewidth = 0.1)

    if (!is.null(net_eq))
      p <- p +
        geom_segment(data = net_eq,
                     aes(x = tmax + gap, xend = x_bar_centre + barw / 2,
                         y = net, yend = net),
                     inherit.aes = FALSE, linetype = "dashed", linewidth = 0.5,
                     colour = "grey20")

    # broken-axis cue: label the equilibrium bar, mark the gap
    brk_lab <- data.frame(x = x_bar_centre, y = -Inf, lab = "Equilibrium")
    p <- p +
      geom_vline(xintercept = tmax + gap / 2, linetype = "dotted",
                 colour = "grey60", linewidth = 0.3) +
      scale_x_continuous(
        breaks = c(pretty(c(0, tmax)), x_bar_centre),
        labels = c(as.character(pretty(c(0, tmax))), "eq"),
        expand = expansion(mult = c(0.01, 0.02)))
  }

  p + facet_wrap(~scenario, scales = "free_y") +
    scale_fill_viridis_d(drop = FALSE) +
    labs(x = "Years after animals added  (eq = equilibrium)",
         y = expression("Animal Effect (g C m"^-2*")"), fill = "Pool") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")
}
