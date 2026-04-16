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

# ---- Load code ----
source("R/tree_monomolecular.R")
source("R/make_tree_forcing.R")
source("R/century_model.R")
source("R/load_config.R")
source("R/init_century_state.R")
source("R/plot_ode_out.R")

# ---- Load parameters ----
parms <- load_config("tree_monomolecular")
parms <- modifyList(parms, yaml::read_yaml("config/century.yml"))

#---- Calculate equilibrium ----
parms$tree_forcing <- make_tree_forcing_equilibrium(parms)

century_eqm = rootSolve::stode(
  y     = init_century_state(),
  func  = century_model,
  parms = parms
)


# ---- Build forcing once ----
parms$tree_forcing <- make_tree_forcing(parms)

# ---- Run ----
times <- seq(0, 365 * 10, by = 1)

out_century <- ode(
  y     = century_eqm$y,
  times = times,
  func  = century_model,
  parms = parms
)

plot_ode_output(out_century,
                variable_cols = names(century_eqm$y))


# ---- Sensitivity Analysis -----

source("R/single_parameter_sensitivity_century.R")

parms  <- load_config("tree_monomolecular")
parms  <- modifyList(parms, yaml::read_yaml("config/century.yml"))

# ---- Look at different parameters -----

single_parameter_sensitivity_century(
  param_name = "w1"
)

single_parameter_sensitivity_century(
  param_name = "c1"
)

single_parameter_sensitivity_century(
  param_name = "t3"
)
