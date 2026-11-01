#------------------------------------------------------------------------#
# TEST FOR MULTIPLE STABLE STATES (alternative equilibria) ----
#------------------------------------------------------------------------#
# TEST FOR MULTIPLE STABLE STATES (alternative equilibria) - constant forcing
# Microbial-explicit soil models can have more than one STABLE equilibrium
# (e.g. an active state vs. a microbial-collapse state). This script probes each
# scenario for that by running the constant-forcing steady-state solver from
# MANY diverse initial conditions and checking whether they settle on distinct
# fixed points.
#
# WHY THIS WORKS: spinup_equilibrium() uses rootSolve::runsteady, which finds a
# steady state by FORWARD integration -- so it only lands on STABLE equilibria
# (unstable ones repel). Sampling initial conditions across the state space and
# clustering the converged results therefore counts the stable states and their
# basins of attraction. This does NOT prove uniqueness (a basin we never seed
# is missed), but repeated single-state convergence from wide, extreme starts is
# strong evidence of a single stable state; >=2 distinct states is proof of
# multistability.
#
# Output: Results/multistability_summary.csv (one row per scenario),
#         Results/multistability_states.csv  (every distinct state found),
#         Results/figures/multistability_states.png
# Run from the project root. Helper functions are in R/multistability.R.
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")
source("R/multistability.R")        # make_starts(), same_state(), cluster_states(), test_scenario()
source("R/millennial_model.R"); source("R/derive_millennial_parms.R")
source("R/init_millennial_state.R")
source("R/run_utils.R")
debugonce_functions()   # only active when options(sfw.debugonce = TRUE); see Scripts/0_run_all.R

set.seed(1)                                   # reproducible initial-condition draws
model     <- "millennial"
scen      <- read_scenarios("Data/scenarios.xlsx")
scenarios <- names(scen)

fig_dir <- "Results/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
res_dir <- "Results";         dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

#------------------------------------------------------------------------#
# test settings ----
#------------------------------------------------------------------------#
n_starts   <- 60        # random initial conditions per scenario (+ fixed corners)
span_dec    <- 2        # random scaling spans +/- this many orders of magnitude
rel_tol     <- 1e-2     # two equilibria are "the same" if max relative diff < this
abs_floor   <- 1e-6     # ignore pools below this when comparing (numerical dust)
use_treatment <- FALSE  # test the BASELINE (no-animal) equilibrium first

# The test helpers (initial conditions, clustering of equilibria, one-scenario
# test) are in R/multistability.R and read the settings above.

#------------------------------------------------------------------------#
# run all scenarios ----
#------------------------------------------------------------------------#
all_summary <- list(); all_states <- list()
for (sc in scenarios) {
  cat("\n=== multistability test:", sc, "===\n")
  r <- test_scenario(sc)
  print(r$summary)
  all_summary[[sc]] <- r$summary
  all_states[[sc]]  <- r$states
}
summary_tbl <- bind_rows(all_summary)
states_tbl  <- bind_rows(all_states)
write_csv(summary_tbl, file.path(res_dir, "multistability_summary.csv"))
write_csv(states_tbl,  file.path(res_dir, "multistability_states.csv"))

cat("\n================= RESULT =================\n")
print(as.data.frame(summary_tbl), row.names = FALSE)
if (any(summary_tbl$multistable)) {
  cat("\nALTERNATIVE STABLE STATES found in:",
      paste(summary_tbl$scenario[summary_tbl$multistable], collapse = ", "), "\n")
} else {
  cat("\nNo multistability detected: every scenario converged to a single stable",
      "state from all", unique(summary_tbl$n_tested), "seeds.\n")
}

#------------------------------------------------------------------------#
# FIGURE: total soil + root C of each distinct stable state per scenario ----
#------------------------------------------------------------------------#
# FIGURE: total soil + root C of each distinct stable state per scenario, sized
# by how many starts fell into it (basin share). Multiple points in a column =
# multiple stable states.
p <- ggplot(states_tbl, aes(scenario, total_soil_C, size = n_starts)) +
  scale_x_discrete(labels = pretty_scenario) +
  geom_point(alpha = 0.7, colour = "#2166ac") +
  scale_size_area(max_size = 10, name = "Starts in basin") +
  labs(title = "Distinct stable states per scenario (constant forcing)",
       subtitle = "Each point = one converged stable state; >1 point in a column = alternative stable states",
       x = NULL, y = expression("Total soil + root C at equilibrium (g C m"^-2*")")) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig_dir, "multistability_states.png"), p, width = 9, height = 6, dpi = 150)
print(p)

cat("\nWrote:\n  ", file.path(res_dir, "multistability_summary.csv"),
    "\n  ", file.path(res_dir, "multistability_states.csv"),
    "\n  ", file.path(fig_dir, "multistability_states.png"), "\n")
