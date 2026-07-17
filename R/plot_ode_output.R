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
# STACKED, GROUPED CHANGE-OVER-TIME PLOT (added animals)
# ------------------------------------------------------------
# For each scenario, the animal EFFECT over the follow-up = the added-animals
# run minus the time-matched continued-baseline (no-animal) run, per pool.
# Pools are relabelled (via name_lookup) and COMBINED into a few groups, then
# stacked so you see the net change and how each group contributes through time.
# One facet per scenario. Annual time steps by default (by = 365).
# ============================================================

# default relabelling key (same names as Scripts/spinup_dynamic.R)
followup_name_lookup <- c(
  C_leaf_herb = "Herbaceous Leaf C", C_root_herb = "Herbaceous Root C",
  C_leaf_tree = "Tree Leaf C", C_wood_tree = "Tree Wood C", C_root_tree = "Tree Root C",
  Earthworm = "Earthworms", Litter = "Litter", CWD = "Coarse Woody Debris",
  Organic = "Organic Matter", DOM = "Dissolved Organic Matter",
  MIC = "Microbial Biomass (Organic horizon)", P = "POC", L = "LWMC",
  A = "Aggregate C", M = "MAOC", B = "Microbial Biomass (Mineral Soil)",
  Detritivore = "Detritivores", RootHerb = "Root Herbivores")

# default grouping of pools into stacked categories (by RAW pool name)
followup_pool_groups <- list(
  "Litter & CWD"        = c("Litter", "CWD"),
  "Particulate & DOM"   = c("P", "DOM"),
  "Organic & leached"   = c("Organic", "L"),
  "Mineral-associated"  = c("A", "M"),
  "Microbial biomass"   = c("MIC", "B"),
  "Plant roots"         = c("C_root_herb", "C_root_tree"))

# ------------------------------------------------------------
# followup_change_long(): change (added - continued baseline) per pool over
# time, for ONE scenario, at annual (or `by`-day) steps. Returns a tidy frame
# with raw `state`, relabelled `pretty`, and the assigned `group`.
# ------------------------------------------------------------
followup_change_long <- function(model, scenario, dir = "Data/followup",
                                 by = 365,
                                 name_lookup = followup_name_lookup,
                                 groups = followup_pool_groups) {
  add  <- load_followup(model, scenario, "add", dir = dir)$out
  ctrl <- load_followup(model, scenario, "continue_baseline", dir = dir)$out
  a <- as.data.frame(add);  names(a)[1] <- "time"
  b <- as.data.frame(ctrl); names(b)[1] <- "time"

  keep <- setdiff(intersect(names(a), names(b)), c("time", "mass_balance_check"))
  tt   <- intersect(a$time, b$time)
  if (!is.null(by)) tt <- tt[tt %in% seq(0, max(tt), by = by)]
  a <- a[a$time %in% tt, , drop = FALSE]
  b <- b[b$time %in% tt, , drop = FALSE]

  # map each pool to a group; pools not in any group are dropped from the stack
  pool2group <- unlist(lapply(names(groups), function(g)
    setNames(rep(g, length(groups[[g]])), groups[[g]])))

  ch <- a[, keep, drop = FALSE] - b[match(a$time, b$time), keep, drop = FALSE]
  long <- data.frame(
    time     = rep(a$time, times = length(keep)),
    state    = rep(keep, each = nrow(a)),
    change   = as.vector(as.matrix(ch)),
    stringsAsFactors = FALSE)
  long$pretty   <- ifelse(long$state %in% names(name_lookup),
                          name_lookup[long$state], long$state)
  long$group    <- pool2group[long$state]
  long$scenario <- scenario
  long[!is.na(long$group), , drop = FALSE]
}

# ------------------------------------------------------------
# plot_followup_stacked(): the faceted stacked graphic. Aggregates the per-pool
# change into the groups, stacks them (positive above 0, negative below), one
# facet per scenario. `scenarios` defaults to every scenario with saved runs.
# ------------------------------------------------------------
plot_followup_stacked <- function(model,
                                  scenarios = NULL,
                                  dir = "Data/followup",
                                  by = 365,
                                  name_lookup = followup_name_lookup,
                                  groups = followup_pool_groups,
                                  title = NULL) {
  library(ggplot2); library(dplyr)

  if (is.null(scenarios)) {
    fs <- list.files(dir, pattern = paste0("^", model, "_.*_add\\.rds$"))
    scenarios <- sub(paste0("^", model, "_(.*)_add\\.rds$"), "\\1", fs)
    scenarios <- scenarios[file.exists(file.path(dir,
                     sprintf("%s_%s_continue_baseline.rds", model, scenarios)))]
  }
  if (!length(scenarios)) stop("No scenarios with saved add + continue_baseline runs in ", dir)

  dat <- bind_rows(lapply(scenarios, function(s)
    followup_change_long(model, s, dir = dir, by = by,
                         name_lookup = name_lookup, groups = groups)))

  # sum change within each group x time x scenario for the stack
  stacked <- dat %>%
    group_by(scenario, group, time) %>%
    summarise(change = sum(change, na.rm = TRUE), .groups = "drop") %>%
    mutate(group = factor(group, levels = names(groups)),
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
    scale_fill_brewer(palette = "Set2") +
    labs(title = if (is.null(title))
           "Animal effect over the follow-up (added - continued baseline)" else title,
         subtitle = "Stacked C change by pool group; dashed line = net change",
         x = "Years after animals added", y = "Change in C (added - baseline)",
         fill = "Pool group") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")
}
