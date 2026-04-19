# Herbivory analysis:

library(pacman)
p_load(deSolve, FME, tidyverse, yaml)
verbose = F

# ---- Load models ----
source("R/millennial_model_herbnem.R")
source("R/derive_millennial_parms.R")
source("R/init_millennial_state.R")
source("R/plot_ode_out.R")


parms <- yaml::read_yaml("config/common_herbivore.yml")
parms  <- modifyList(parms, yaml::read_yaml("config/millennial.yml"))
parms  <- modifyList(parms, yaml::read_yaml("config/millennial_herbnem.yml"))

# All roots to mineral soil:
parms$root_to_organic = 0 

parms <- derive_millennial_parms(parms)

parms$wood_mortality = 0

# With herbivores:
millennial_herbnem = ode(
  times = seq(1, 365*3000, by = 1),
  y     = init_millennial_state(HerbNem = T),
  func  = millennial_model_herbnem,
  parms = parms
)

plot_ode_output(millennial_herbnem,
                variable_cols = names(init_millennial_state(HerbNem = T)))

head(millennial_herbnem)
tail(millennial_herbnem)


# Stability Check (code from Janey)
timeCheck_stable <- tibble(data.frame(millennial_herbnem)) %>%
  filter(
    time %in% c(
      max(times),
      max(times)-365,
      max(times)-365*2
    )
  ) # this creates an empty tibble?

timeEND <- as.data.frame(millennial_herbnem) %>%
  slice_tail(n = 1) %>%
  select(-time) %>%
  unlist()

new_herbnem <- ode(
  y = timeEND,
  times = times,
  func = millennial_model_herbnem,
  parms = parms
)

plot_ode_output(new_herbnem,
                variable_cols = names(init_millennial_state(HerbNem = T)))
