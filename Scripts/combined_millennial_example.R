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
source("R/plot_ode_output.R")

parms  <- yaml::read_yaml("config/common.yml")
parms  <- modifyList(parms, yaml::read_yaml("config/millennial.yml"))
parms  <- derive_millennial_parms(parms)

parms$climate_forcing <- make_climate_forcing(parms)

# Plot forcing if you'd like:
if(verbose){
  tibble(data.frame(t(sapply(1:365*2, parms$climate_forcing)))) %>%
    mutate(day = 1:365*2) %>%
    pivot_longer(!day) %>%
    ggplot(aes(x = day, y = value)) + geom_line() + facet_wrap(.~name, scales = "free_y")
  
}

#Run the model:

y0 = init_millennial_state()
# y0["Detritivore"] = 0
# y0["Earthworm"] = 0

millennial_out = ode(
  times = seq(1, 365*1, by = 1),
  y     = y0,
  func  = millennial_model_wplant,
  parms = parms
)

plot_ode_output(millennial_out)
