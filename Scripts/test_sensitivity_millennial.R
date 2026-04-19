# A script to conduct a sensitivity analysis on the models:

if (!"pacman" %in% rownames(installed.packages())) {
  install.packages("pacman", dependencies = TRUE)
}

library(pacman)
p_load(deSolve, FME, tidyverse, yaml)
verbose = F

# Load in library:
library(deSolve)
library(tidyverse)

# ---- Load models ----
source("R/tree_monomolecular.R")
source("R/herb_forcing_monomolecular.R")
source("R/millennial_model.R")

# ---- Load config utilities ----
source("R/load_config.R")
source("R/make_tree_forcing.R")
source("R/derive_millennial_parms.R")
source("R/init_millennial_state.R")
source("R/plot_ode_out.R")

# ---- Load and merge configs ----
parms  <- load_config("tree_monomolecular")
parms  <- modifyList(parms, yaml::read_yaml("config/millennial.yml"))
parms  <- derive_millennial_parms(parms)

parms$tree_forcing <- make_tree_forcing_equilibrium(parms)

millennial_eqm = rootSolve::stode(
  y     = init_millennial_state(),
  func  = millennial_model,
  parms = parms
)

# ---- Build forcing ONCE ----
parms$tree_forcing <- make_tree_forcing(parms)


# ---- Run model ----
times <- seq(0, 365 * 10, by = 1)

out_millennial <- ode(
  y = millennial_eqm$y,
  times = times, 
  func = millennial_model, 
  parms = parms)

plot_ode_output(out_millennial,
                variable_cols = names(millennial_eqm$y))


# ---- Sensitivity Analysis -----

parms  <- load_config("tree_monomolecular")
parms  <- modifyList(parms, yaml::read_yaml("config/millennial.yml"))

# ---- Look at different parameters -----

source("R/single_parameter_sensitivity_millennial.R")

single_parameter_sensitivity_millennial(
  param_name = "p_a"
)

single_parameter_sensitivity_millennial(
  param_name = "p_b"
)

single_parameter_sensitivity_millennial(
  param_name = "k_exudate"
)

