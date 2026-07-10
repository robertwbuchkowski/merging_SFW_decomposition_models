# ============================================================
# GRAPHICS - millennial model: the effect of each animal scenario
#   (A) AT EQUILIBRIUM  : animal (treatment) minus no-animal (baseline)
#                         equilibrium, per soil/plant pool, per scenario.
#   (B) DURING FOLLOW-UP: the ~100-yr seasonal trajectory after animals are
#                         added, vs. the continued no-animal baseline.
#
# Also demonstrates the NO-INDIRECT-EFFECTS variant (all _pint parameters = 0,
# so animals act only through direct feeding). Set indirect_effects <- FALSE
# to reproduce every figure with indirect effects switched off.
#
# Run from project root. (A) is self-contained (fast equilibria). (B) needs the
# saved spin-ups + follow-up runs from Scripts/spinup_dynamic.R and
# Scripts/followup_analysis.R; it falls back to a message if they are missing.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")

model  <- "millennial"
scen   <- read_scenarios("Data/scenarios.xlsx")
scen$MitePredator <- NULL

indirect_effects <- TRUE          # FALSE => zero all _pint params (direct feeding only)
use_fitted_params <- TRUE         # apply calibrated adj_*/effect params if available
fig_dir <- "Results/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
suffix  <- if (indirect_effects) "" else "_noindirect"

fitted_params <- if (use_fitted_params && file.exists("Results/fitted_animal_params.csv"))
  load_fitted_params("Results/fitted_animal_params.csv") else NULL

# ============================================================
# (A) EQUILIBRIUM EFFECT OF EACH ANIMAL SCENARIO
# ============================================================
eq_rows <- list()
for (scenario in names(scen)) {
  pair <- setup_scenario_pair(model, scen, scenario)
  if (!indirect_effects) pair <- zero_indirect_effects(pair)      # _pint = 0
  if (!is.null(fitted_params))
    pair$treatment <- apply_fitted_params(pair$treatment, fitted_params, model, scenario)

  pair$baseline  <- spinup_equilibrium(pair$baseline,  verbose = FALSE)
  pair$treatment <- spinup_equilibrium(pair$treatment,
                                       warm_start = pair$baseline$init_state_spin, verbose = FALSE)

  cmp <- compare_vectors(pair$treatment$init_state_spin, pair$baseline$init_state_spin)
  cmp$scenario <- scenario
  eq_rows[[scenario]] <- cmp
  cat("equilibrium done:", scenario, "\n")
}
eq <- bind_rows(eq_rows)

animal_pools <- c("Earthworm","Detritivore","DetPredator","RootHerb")
eq_soil <- eq %>% filter(!is.na(baseline), !name %in% animal_pools)   # shared pools only

# --- (A) figure: % change in each pool caused by the animals, faceted by pool
pA <- ggplot(eq_soil, aes(x = scenario, y = percent_change, fill = scenario)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  facet_wrap(~name, scales = "free_y") +
  labs(title = paste0("Millennial: equilibrium effect of each animal scenario",
                      if (!indirect_effects) " (no indirect effects)" else ""),
       subtitle = "Treatment (with animals) - baseline (no animals), % change per pool",
       x = NULL, y = "% change vs no-animal baseline") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none")

ggsave(file.path(fig_dir, paste0("millennial_equilibrium_effect", suffix, ".png")),
       pA, width = 11, height = 7, dpi = 150)
write_csv(eq_soil, file.path(fig_dir, paste0("millennial_equilibrium_effect", suffix, ".csv")))
print(pA)

# ============================================================
# (B) FOLLOW-UP: animals added vs continued baseline
# Uses saved runs from followup_analysis.R if present; otherwise computes them
# here from the saved spin-ups (Scripts/spinup_dynamic.R).
# ============================================================
followup_one <- function(scenario, n_years = 100, by = 30) {
  add_file  <- sprintf("Data/followup/%s_%s_add.rds", model, scenario)
  ctrl_file <- sprintf("Data/followup/%s_%s_continue_baseline.rds", model, scenario)
  if (file.exists(add_file) && file.exists(ctrl_file))
    return(list(add = load_followup(model, scenario, "add"),
                control = load_followup(model, scenario, "continue_baseline")))

  base_file <- sprintf("Data/spinup/%s_%s_baseline.rds", model, scenario)
  if (!file.exists(base_file)) return(NULL)                     # nothing to work from
  base_saved <- load_spinup(base_file)

  treatment_setup <- setup_scenario(model, scen, scenario, animals = TRUE)
  baseline_setup  <- setup_scenario(model, scen, scenario, animals = FALSE)
  if (!indirect_effects) {
    treatment_setup <- zero_indirect_effects(treatment_setup)
    baseline_setup  <- zero_indirect_effects(baseline_setup, verbose = FALSE)
  }
  if (!is.null(fitted_params))
    treatment_setup <- apply_fitted_params(treatment_setup, fitted_params, model, scenario)

  list(add     = followup_add_animals(base_saved, treatment_setup, n_years = n_years, by = by),
       control = followup_continue_baseline(base_saved, baseline_setup, n_years = n_years, by = by))
}

for (scenario in names(scen)) {
  fu <- tryCatch(followup_one(scenario), error = function(e) {
    message("follow-up skipped (", scenario, "): ", conditionMessage(e)); NULL })
  if (is.null(fu)) { message("no saved spin-up for ", scenario,
                             " -- run spinup_dynamic.R first."); next }

  pB <- plot_followup_comparison(
    fu$control$out, fu$add$out,
    control_label   = "Continued baseline (no animals)",
    treatment_label = "Animals added",
    title = paste0("Millennial / ", scenario, ": follow-up",
                   if (!indirect_effects) " (no indirect effects)" else ""))
  ggsave(file.path(fig_dir, paste0("millennial_followup_", scenario, suffix, ".png")),
         pB, width = 11, height = 7, dpi = 150)
  print(pB)
  cat("follow-up figure done:", scenario, "\n")
}

cat("\nFigures written to", fig_dir, "\n")
