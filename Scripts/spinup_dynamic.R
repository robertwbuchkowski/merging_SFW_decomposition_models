# ============================================================
# SPIN-UP  -  run once per scenario, then reuse the saved states.
# Pipeline: equilibrium spin-up -> apply saved fitted animal params ->
# record the equilibrium animal effect (total + direct) -> seasonal spin-up to
# the limit cycle -> save the stable state to Data/spinup/.
#
# HOW TO USE
#   * Set the toggles below (models, which arms to spin up, Newton vs forward).
#   * Run from the project root; this is the slow step, so it is kept separate.
#   * Outputs: Results/animal_eq_effect.csv, Plots/total_effect.png, and the
#     saved limit-cycle states under Data/spinup/ (consumed by
#     Scripts/followup_analysis.R).
# Requires Results/fitted_animal_params.csv from Scripts/fit_all_animals.R when
# use_fitted_params = TRUE.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")

scen   <- read_scenarios("Data/scenarios.xlsx")
models <- c("millennial")

# ---- toggles -------------------------------------------------------------
use_fitted_params <- TRUE    # apply saved fitted params (from fit_all_animals.R)
do_spinup         <- TRUE    # run the seasonal dynamic spin-up (the slow part)
do_treatment      <- FALSE   # also spin up the treatment arm now (else baseline only)
use_newton        <- TRUE    # TRUE = Newton shooting (fast, exact limit cycle);
                             # FALSE = forward-integration spin-up

scen$MitePredator <- NULL    # drop scenarios you are not running

# fitted animal parameters, keyed by model x scenario x param
fitted_params <- if (use_fitted_params)
  load_fitted_params("Results/fitted_animal_params.csv") else NULL

animal_eq_effect <- list()
for (model in models) {
  for (scenario in names(scen)) {

    pair <- setup_scenario_pair(model, scen, scenario)

    # baseline equilibrium: the no-animal reference every effect is measured against
    pair$baseline <- spinup_equilibrium(pair$baseline)

    # apply the calibrated animal feeding/effect parameters to the treatment arm
    if (use_fitted_params) {
      pair$treatment <- apply_fitted_params(pair$treatment, fitted_params, model, scenario)
    }

    # treatment equilibrium + a quick print of the animal biomass reached
    if (is.null(pair$treatment$init_state_spin))
      pair$treatment <- spinup_equilibrium(pair$treatment,
                                           warm_start = pair$baseline$init_state_spin)
    eq_t    <- pair$treatment$init_state_spin
    animals <- intersect(c("Earthworm", "Detritivore", "DetPredator", "RootHerb"),
                         names(eq_t))
    cat("\nEquilibrium animal biomass (", model, "/", scenario, "):\n", sep = "")
    print(round(eq_t[animals], 4))

    # DIRECT effect: re-equilibrate with the indirect (_pint / exudate) slopes
    # zeroed, so the animal acts through direct feeding only
    pair2     <- zero_indirect_effects(pair)
    eq_direct <- spinup_equilibrium(pair2$treatment,
                                    warm_start = pair2$baseline$init_state_spin)$init_state_spin

    # store total (all pathways) and direct effect vs the baseline equilibrium
    animal_eq_effect[[length(animal_eq_effect) + 1]] <- rbind(
      cbind(compare_vectors(eq_t,      pair$baseline$init_state_spin),
            model = model, scenario = scenario, type = "total"),
      cbind(compare_vectors(eq_direct, pair$baseline$init_state_spin),
            model = model, scenario = scenario, type = "direct"))
    rm(pair2, eq_direct)

    if (do_spinup) {
      # baseline: seasonal spin-up to the limit cycle, then save
      cat("\n--- Baseline dynamic spin-up ---\n")
      dyn_b <- if (use_newton) dynamic_spinup_newton(pair$baseline)
               else dynamic_spinup(pair$baseline, n_years = 600, by = 1, tol = 1e-4)
      cat("baseline converged:", dyn_b$converged, "\n")
      save_spinup(pair$baseline, dyn_b$final_state, scenario, "baseline")

      # treatment: optional here (can be run later in its own session)
      if (do_treatment) {
        cat("\n--- Treatment dynamic spin-up ---\n")
        if (is.null(pair$treatment$init_state_spin))
          pair$treatment <- spinup_equilibrium(pair$treatment,
                                               warm_start = pair$baseline$init_state_spin)
        dyn_t <- if (use_newton) dynamic_spinup_newton(pair$treatment)
                 else dynamic_spinup(pair$treatment, n_years = 600, by = 1, tol = 1e-4)
        cat("treatment converged:", dyn_t$converged, "\n")
        save_spinup(pair$treatment, dyn_t$final_state, scenario, "treatment")
      }

      cat("\nStable states saved under Data/spinup/. Use Scripts/followup_analysis.R next.\n")
    }

  }
}

animal_eq_effect <- do.call("rbind", animal_eq_effect)
write_csv(animal_eq_effect, "Results/animal_eq_effect.csv")

# ---- equilibrium animal-effect figure ------------------------------------
# Pools to show (soil + herbaceous root C); name_lookup / plot_order come from
# R/compare_functions.R.

dir.create("Plots", showWarnings = FALSE)
png("Plots/total_effect.png", width = 4, height = 8, units = "in", res = 600)
animal_eq_effect %>% filter(!is.na(baseline)) %>%
  mutate(pretty_name = name_lookup[name]) %>%
  mutate(pretty_name = factor(pretty_name, levels = plot_order)) %>%
  filter(pretty_name %in% keep_plot) %>%
  ggplot(aes(x = pretty_name, y = percent_change, fill = type)) +
  geom_col(position = "dodge") +
  facet_wrap(. ~ scenario, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Animal Effect (%)") + xlab("") +
  scale_fill_manual(name = "Type",values = c("black", "blue"))
dev.off()
