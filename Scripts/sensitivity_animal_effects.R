# ============================================================
# SENSITIVITY OF ANIMAL EFFECTS TO MODEL PARAMETERS (at equilibrium)
# ------------------------------------------------------------
# Like Scripts/spinup_dynamic.R but WITHOUT the long dynamic spin-up: for each
# scenario it only computes the equilibrium animal effect, and it repeats that
# calculation while sweeping chosen model parameters, one at a time, over a
# grid. For every parameter value it records:
#   TOTAL  effect = treatment (animals, all pathways) - baseline (no animals)
#   DIRECT effect = treatment with indirect (_pint / exudate) slopes zeroed
#                   - baseline
# (both as % change vs the no-animal baseline equilibrium, per pool), then plots
# how those effects respond to each parameter.
#
# Nothing is spun up dynamically and nothing is saved to Data/spinup -- this is
# a fast, equilibrium-only sensitivity sweep. Run from the project root.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R"); source("R/derive_millennial_parms.R");
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")

model   <- "millennial"
scen    <- read_scenarios("Data/scenarios.xlsx")
scen$MitePredator <- NULL

use_fitted_params <- TRUE
fitted_params <- if (use_fitted_params && file.exists("Results/fitted_animal_params.csv"))
  load_fitted_params("Results/fitted_animal_params.csv") else NULL

fig_dir <- "Results/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
res_dir <- "Results";         dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# animal biomass pools (excluded from the "effect on pools" output)
animal_pools <- c("Earthworm", "Detritivore", "DetPredator", "RootHerb")
derive_fn    <- match.fun(model_table[[model]]$derive)

# ------------------------------------------------------------
# set_param(): change ONE parameter on a setup object and rebuild everything
# that depends on it (derived params + climate forcing), so a sweep over a
# parameter that feeds phi_por / the forcing uses consistent values.
# ------------------------------------------------------------
set_param <- function(obj, param, value) {
  obj$parms[[param]]        <- value
  obj$parms                 <- derive_fn(obj$parms)
  obj$parms$climate_forcing <- make_climate_forcing(obj$parms)
  obj
}

# ------------------------------------------------------------
# animal_effect_eq(): the core calculation (no dynamic spin-up). Given a
# treatment/baseline pair, returns a long data frame of the TOTAL and DIRECT
# animal effect on every shared pool, at equilibrium.
# ------------------------------------------------------------
animal_effect_eq <- function(pair) {
  pair$baseline  <- spinup_equilibrium(pair$baseline, verbose = FALSE)
  pair$treatment <- spinup_equilibrium(pair$treatment,
                                       warm_start = pair$baseline$init_state_spin,
                                       verbose = FALSE)
  eq_b <- pair$baseline$init_state_spin
  eq_t <- pair$treatment$init_state_spin

  # direct-only: zero the indirect (_pint + exudate) slopes, re-equilibrate
  pair_d <- zero_indirect_effects(pair, verbose = FALSE)
  eq_d   <- spinup_equilibrium(pair_d$treatment,
                               warm_start = eq_b, verbose = FALSE)$init_state_spin

  tot <- compare_vectors(eq_t, eq_b); tot$type <- "total"
  dir <- compare_vectors(eq_d, eq_b); dir$type <- "direct"
  out <- rbind(tot, dir)
  out[!is.na(out$baseline) & !(out$name %in% animal_pools), ]     # shared, non-animal pools
}

# ------------------------------------------------------------
# sweep_param(): rebuild the scenario pair, apply fitted params, then vary ONE
# parameter over `values` and recompute the equilibrium animal effect each time.
# ------------------------------------------------------------
sweep_param <- function(scenario, param, values) {
  base_pair <- setup_scenario_pair(model, scen, scenario)
  if (!is.null(fitted_params))
    base_pair$treatment <- apply_fitted_params(base_pair$treatment, fitted_params,
                                               model, scenario, verbose = FALSE)
  default <- base_pair$treatment$parms[[param]]

  rows <- list()
  for (v in values) {
    pair <- base_pair
    pair$treatment <- set_param(pair$treatment, param, v)
    pair$baseline  <- set_param(pair$baseline,  param, v)
    eff <- tryCatch(animal_effect_eq(pair), error = function(e) {
      message("  ", scenario, " ", param, "=", signif(v, 4), " failed: ", conditionMessage(e)); NULL })
    if (is.null(eff)) next
    eff$param    <- param
    eff$value    <- v
    eff$default  <- default
    eff$scenario <- scenario
    rows[[length(rows) + 1]] <- eff
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# ============================================================
# CHOOSE what to sweep.
#   sweep_params  the model parameters to test (must exist in parms)
#   n_points      grid points per parameter
#   buffer        +/- fractional range around each parameter's default (linear)
#   scenarios     which animal scenarios to run
# ============================================================
sweep_params <- c("k_frag_litter", "k_frag_organic", "k_l_o", "k_l",
                  "k_b", "k_pa", "k_ma", "pct_claysilt",
                  "k_MICd", "k_bd", "root_to_organic", "a_root_herb",
                  "MAT", "MAtheta")
n_points  <- 7
buffer    <- 0.5
scenarios <- names(scen)

# run the full sweep: scenario x parameter
all_rows <- list()
for (scenario in scenarios) {
  pair0 <- setup_scenario_pair(model, scen, scenario)
  avail <- intersect(sweep_params, names(pair0$treatment$parms))
  for (param in avail) {
    d0   <- pair0$treatment$parms[[param]]
    vals <- fit_param_grid(d0, buffer = buffer, n = n_points, scale = "linear")
    cat("sweeping", scenario, "/", param, "around", signif(d0, 4), "\n")
    s <- sweep_param(scenario, param, vals)
    if (!is.null(s)) all_rows[[paste(scenario, param)]] <- s
  }
}
sens <- bind_rows(all_rows)
write_csv(sens, file.path(res_dir, "animal_effect_sensitivity.csv"))

# ------------------------------------------------------------
# SENSITIVITY SUMMARY: for each scenario x parameter x effect-type, how much
# does the animal effect (summed |% change| across pools) move across the swept
# range? This is one clean overview number per (parameter, effect type).
# ------------------------------------------------------------
summary_tbl <- sens %>%
  group_by(scenario, param, type, value, default) %>%
  summarise(total_abs_effect = sum(abs(percent_change), na.rm = TRUE), .groups = "drop") %>%
  group_by(scenario, param, type) %>%
  summarise(effect_range = max(total_abs_effect) - min(total_abs_effect),
            effect_at_default = total_abs_effect[which.min(abs(value - default))],
            .groups = "drop")
write_csv(summary_tbl, file.path(res_dir, "animal_effect_sensitivity_summary.csv"))

# ------------------------------------------------------------
# PLOT 1 - overview: sensitivity of the whole-soil animal effect to each
# parameter. x = parameter relative to its default; y = summed |effect| across
# pools; colour = parameter; solid = total, dashed = direct; facet by scenario.
# ------------------------------------------------------------
overview <- sens %>%
  group_by(scenario, param, type, value, default) %>%
  summarise(total_abs_effect = sum(abs(difference), na.rm = TRUE), .groups = "drop") %>%
  mutate(rel_param = value / default)

p_overview <- ggplot(overview,
                     aes(rel_param, total_abs_effect, colour = param, linetype = type)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey70") +
  geom_line(linewidth = 0.7) + geom_point(size = 0.8) +
  facet_wrap(~scenario, scales = "free_y") +
  scale_colour_viridis_d() +
  scale_linetype_manual(values = c(total = "solid", direct = "dashed")) +
  labs(x = "Parameter value (relative to default)", y = "Animal effect on total C (g C m^-2)",
       colour = "Parameter", linetype = "Effect") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "animal_effect_sensitivity_overview.png"),
       p_overview, width = 12, height = 8, dpi = 150)
print(p_overview)

# ------------------------------------------------------------
# PLOT 2 - detail for ONE scenario: the animal effect on EACH pool as a
# parameter varies (one panel per parameter, one line per pool, total effect).
# Change `focus_scenario` to inspect others.
# ------------------------------------------------------------
for(iii in 1:4){
  focus_scenario <- scenarios[iii]
  detail <- sens %>% filter(scenario == focus_scenario, type == "total") %>%
    mutate(rel_param = value / default)
  
  p_detail <- ggplot(detail, aes(rel_param, percent_change, colour = name)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_line(linewidth = 0.6) +
    facet_wrap(.~param, scales = "free") +
    labs(title = paste0("Animal effect per pool vs parameter -- ", focus_scenario, " (total effect)"),
         x = "Parameter value (relative to default)", y = "Effect on pool (%)",
         colour = "Pool") +
    theme_minimal(base_size = 10) + theme(legend.position = "right")
  ggsave(file.path(fig_dir, paste0("animal_effect_sensitivity_detail_", focus_scenario, ".png")),
         p_detail, width = 12, height = 8, dpi = 150)
  print(p_detail)
}

cat("\nWrote:\n  ", file.path(res_dir, "animal_effect_sensitivity.csv"),
    "\n  ", file.path(res_dir, "animal_effect_sensitivity_summary.csv"),
    "\n  ", file.path(fig_dir, "animal_effect_sensitivity_overview.png"),
    "\n  ", file.path(fig_dir, paste0("animal_effect_sensitivity_detail_", focus_scenario, ".png")), "\n")
