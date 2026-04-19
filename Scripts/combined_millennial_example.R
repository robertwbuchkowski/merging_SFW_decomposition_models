# Working on the combined model to integrate all of the case studies into one usable model:

library(pacman)
p_load(deSolve, FME, tidyverse, yaml)
verbose = F

# ---- Load models ----
source("R/millennial_model.R")

# ---- Load config utilities ----
source("R/climate_forcing.R")
source("R/derive_millennial_parms.R")
source("R/init_millennial_state.R")
source("R/plot_ode_out.R")

parms  <- yaml::read_yaml("config/common.yml")
parms  <- modifyList(parms, yaml::read_yaml("config/tree_monomolecular.yml"))
parms  <- modifyList(parms, yaml::read_yaml("config/millennial.yml"))
parms  <- derive_millennial_parms(parms)

parms$climate_forcing <- make_climate_forcing_equilibrium(parms)

tibble(data.frame(t(sapply(1:365*2, parms$climate_forcing)))) %>%
  mutate(doy = 1:365*2) %>%
  pivot_longer(!doy) %>%
  ggplot(aes(x = doy, y = value)) + geom_line() + facet_wrap(.~name, scales = "free_y")


tibble(data.frame(t(sapply(1:365*2, parms$climate_forcing)))) %>%
  mutate(doy = 1:365*2) %>%
  pivot_longer(!doy) %>%
  group_by(name) %>%
  summarize(sum(value))
