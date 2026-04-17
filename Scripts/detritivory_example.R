# Detritivory Example:

library(pacman)
p_load(deSolve, FME, tidyverse, yaml)
verbose = F

# ---- Load models ----
source("R/tree_monomolecular.R")
source("R/millennial_model_detritivory.R")
source("R/millennial_model.R")

# ---- Load config utilities ----
source("R/make_tree_forcing.R")
source("R/derive_millennial_parms.R")
source("R/init_millennial_state.R")
source("R/plot_ode_out.R")

parms  <- yaml::read_yaml("config/common.yml")
parms  <- modifyList(parms, yaml::read_yaml("config/tree_monomolecular.yml"))
parms  <- modifyList(parms, yaml::read_yaml("config/millennial.yml"))
parms  <- derive_millennial_parms(parms)

parms$tree_forcing <- make_tree_forcing_equilibrium(parms)

# With detritivores:
millennial_eqm_det = rootSolve::stode(
  y     = init_millennial_state(T),
  func  = millennial_model_detritivory,
  parms = parms
)

millennial_eqm_det$y

# Without detritivores:
millennial_eqm = rootSolve::stode(
  y     = init_millennial_state(F),
  func  = millennial_model,
  parms = parms
)

millennial_eqm$y
