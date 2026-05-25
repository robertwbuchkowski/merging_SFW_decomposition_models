# Working on the combined model to integrate all of the case studies into one usable model:

library(pacman)
p_load(deSolve, FME, tidyverse, yaml)
verbose = F

# ---- Load models ----
source("R/millennial_model.R")
source("R/millennial_wrapper.R")

# ---- Load config utilities ----
source("R/climate_forcing.R")
source("R/derive_millennial_parms.R")
source("R/init_millennial_state.R")
source("R/plot_ode_output.R")

parms  <- yaml::read_yaml("config/common.yml")
parms  <- modifyList(parms, yaml::read_yaml("config/millennial.yml"))
parms  <- derive_millennial_parms(parms)

# parms$climate_forcing <- make_climate_forcing(parms)

parms$climate_forcing <- make_climate_forcing_equilibrium(parms)

# Plot forcing if you'd like:
if(verbose){
  tibble(data.frame(t(sapply(1:365*2, parms$climate_forcing)))) %>%
    mutate(day = 1:365*2) %>%
    pivot_longer(!day) %>%
    ggplot(aes(x = day, y = value)) + geom_line() + facet_wrap(.~name, scales = "free_y")
  
}

# --------------------------------------------
# FULL INITIAL CONDITIONS
# --------------------------------------------
init_full <- init_millennial_state()

# ============================================
# ✅ EXAMPLE 1 — FULL MODEL
# ============================================

config_full <- list(
  herb        = TRUE,
  tree        = FALSE,
  earthworm   = FALSE,
  detritivore = FALSE,
  rootherb    = TRUE
)

parms  <- modifyList(parms, list(config = config_full))

state_full <- build_initial_state(init_full, config_full)

steady_full <- stode(
  y = state_full,
  func = millennial_wrapper,
  parms = parms
)

print("FULL MODEL STEADY STATE:")
print(steady_full$y)


# ============================================
# ✅ EXAMPLE 2 — NO ROOT HERBIVORES
# ============================================

config_no_rootherb <- list(
  herb        = TRUE,
  tree        = FALSE,
  earthworm   = FALSE,
  detritivore = FALSE,
  rootherb    = FALSE
)

parms$config <- config_no_rootherb

state_no_rootherb <- build_initial_state(init_full, config_no_rootherb)

steady_no_rootherb <- steady(
  y = state_no_rootherb,
  func = millennial_wrapper,
  parms = parms
)

print("NO ROOT HERBIVORE STEADY STATE:")
print(steady_no_rootherb$y)

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
