# Herbivory analysis:

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

parms  <- load_config("tree_monomolecular")
parms  <- modifyList(parms, yaml::read_yaml("config/millennial.yml"))

# All roots to mineral soil:
parms$root_to_organic = 0 






# ---- Look at different parameters -----

source("R/single_parameter_sensitivity_millennial.R")

single_parameter_sensitivity_millennial(
  param_name = "a_root_herb",
  habitat = "herb"
)

# Herbivory case study:
# 8x more stabilization of root litter than shoot litter
# 61.78 -1.21*rootherbary
# consumption is 20% of the total NPP is consumed by root herbivores


source("R/scenario_sensitivity_millennial.R")

pm = expand.grid(a_root_herb = parms$a_root_herb*seq(0.5, 1.5, length = 5),
                  TBHmax = parms$TBHmax*seq(0.5, 1.5, length = 5))

scenario_sensitivity_millennial(pm,
                                habitat_type = "field")


pm = cbind(a_root_herb = parms$a_root_herb*seq(0.5, 1.5, length = 5),
                 TBHmax = parms$TBHmax*seq(0.5, 1.5, length = 5))

scenario_sensitivity_millennial(pm,
                                habitat_type = "field")