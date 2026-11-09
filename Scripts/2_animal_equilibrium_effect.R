#------------------------------------------------------------------------#
# Animal effects at equilibrium ----
#------------------------------------------------------------------------#
# EQUILIBRIUM ANIMAL EFFECTS at the default parameter values, per scenario,
# under constant mean climate forcing (the same steady state the animals are
# fitted at in Scripts/1):
#   baseline equilibrium (no animals) -> apply the fitted animal parameters ->
#   treatment equilibrium (with animals) -> treatment equilibrium with the
#   trait-mediated indirect effects switched off -> effects per pool:
#     total        all pathways on
#     tmii         TRAIT-MEDIATED INDIRECT effect = total - direct_dmii
#     direct_dmii  trait-mediated indirect effects off (direct + density-mediated)
#   -> carbon flux networks with and without animals.
# The seasonal spin-up used by the invasion runs is in Scripts/3.
#
# Outputs:
#   Results/animal_eq_effect.csv       effects per pool and scenario
#   Results/animal_effect_summary.csv  effects on total soil + root C, per scenario
#   Results/figures/flux_network_<scenario>.png (+ .csv)  carbon fluxes with and
#     without animals, and the change
#   Results/figures/flux_network_change_all_scenarios.png  the change panel of
#     every scenario in one figure
#   Results/figures/flux_network_wiring.png  how the model is wired (no values)
# Run from the project root. Requires Results/fitted_animal_params.csv from
# Scripts/1_fit_all_animals.R when use_fitted_params = TRUE.
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")
source("R/millennial_model.R"); source("R/derive_millennial_parms.R")
source("R/init_millennial_state.R"); source("R/flux_network.R")
source("R/run_utils.R")
debugonce_functions()   # only active when options(sfw.debugonce = TRUE); see Scripts/0_run_all.R

scen   <- read_scenarios("Data/scenarios.xlsx")
models <- c("millennial")

#------------------------------------------------------------------------#
# toggles ----
#------------------------------------------------------------------------#
use_fitted_params <- TRUE    # apply saved fitted params (from fit_all_animals.R)

# fitted animal parameters, keyed by model x scenario x param
fitted_params <- if (use_fitted_params)
  load_fitted_params("Results/fitted_animal_params.csv") else NULL

animal_eq_effect <- list()
flux_nets <- list()       # flux networks of every scenario, for the combined figure
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
    animals <- intersect(c("Earthworm", "Detritivore", "RootHerb"),
                         names(eq_t))
    cat("\nEquilibrium animal biomass (", model, "/", scenario, "):\n", sep = "")
    print(round(eq_t[animals], 4))

    # re-equilibrate with the trait-mediated indirect effects switched off
    # (_pint / exudate slopes zeroed): the animal then acts through direct
    # feeding AND density-mediated indirect effects only ("direct_dmii")
    pair2     <- zero_indirect_effects(pair)
    eq_direct <- spinup_equilibrium(pair2$treatment,
                                    warm_start = pair2$baseline$init_state_spin)$init_state_spin

    # effects vs the baseline equilibrium, per pool:
    #   total        all pathways on
    #   direct_dmii  trait-mediated indirect effects off (direct + density-mediated)
    #   tmii         TRAIT-MEDIATED INDIRECT effect = total - direct_dmii
    eff_total  <- compare_vectors(eq_t,      pair$baseline$init_state_spin)
    eff_dd     <- compare_vectors(eq_direct, pair$baseline$init_state_spin)
    eff_tmii   <- eff_total
    stopifnot(identical(eff_total$name, eff_dd$name))
    eff_tmii$treatment  <- NA_real_          # a difference of effects, not a state
    eff_tmii$difference <- eff_total$difference - eff_dd$difference
    eff_tmii$percent_change <- 100 * eff_tmii$difference / eff_total$baseline
    animal_eq_effect[[length(animal_eq_effect) + 1]] <- rbind(
      cbind(eff_total, model = model, scenario = scenario, type = "total"),
      cbind(eff_tmii,  model = model, scenario = scenario, type = "tmii"),
      cbind(eff_dd,    model = model, scenario = scenario, type = "direct_dmii"))
    rm(pair2, eq_direct, eff_total, eff_dd, eff_tmii)

    # carbon flux network at both equilibria + the change
    # (Results/figures/flux_network_*.png, flux table alongside as .csv)
    net <- flux_network_scenario(model, scenario, scen, pair = pair)
    save_flux_network(net, sprintf("Results/figures/flux_network_%s.png", scenario))
    flux_nets[[scenario]] <- net
    rm(net)
  }
}

animal_eq_effect <- do.call("rbind", animal_eq_effect)

# change in carbon fluxes (with - without animals) of every scenario in one
# figure, one panel per scenario
save_flux_network_deltas(flux_nets, "Results/figures/flux_network_change_all_scenarios.png")

# introductory figure: how the model is wired (every flux, no values), with
# dashed lines from each animal to the flux its trait-mediated indirect
# effect acts on
save_flux_network_wiring("Results/figures/flux_network_wiring.png")
write_csv(animal_eq_effect, "Results/animal_eq_effect.csv")

#------------------------------------------------------------------------#
# summary table: total and trait-mediated indirect effect, one row per scenario ----
#------------------------------------------------------------------------#
# For each scenario the animal effect on TOTAL soil + root C is the net change
# summed over the shared pools (difference = treatment - baseline; animal
# biomass pools have no baseline so are NA and excluded). "direct_dmii" repeats
# this with the trait-mediated indirect effects (_pint / exudate slopes)
# switched off: the direct plus density-mediated indirect effects. The
# trait-mediated indirect effect (tmii) is total - direct_dmii; pct_tmii is its
# share of the total effect.
animal_effect_summary <- animal_eq_effect %>%
  filter(!is.na(difference)) %>%                       # shared soil/root pools only
  group_by(model, scenario, type) %>%
  summarise(effect = sum(difference),
            baseline = sum(baseline), .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = effect) %>%
  transmute(model, scenario,
            total_effect_gCm2       = signif(total, 2),
            tmii_effect_gCm2        = signif(tmii, 2),
            pct_tmii_of_total       = signif(100 * tmii / total, 2),
            direct_dmii_effect_gCm2 = signif(direct_dmii, 2),
            pct_total_of_baseline = signif(100 * total / baseline, 2))

write_csv(animal_effect_summary, "Results/animal_effect_summary.csv")
cat("\nAnimal effect on total soil + root C (g C m-2), by scenario:\n")
print(as.data.frame(animal_effect_summary), row.names = FALSE)
