# Earthworm Example:

library(pacman)
p_load(deSolve, FME, tidyverse, yaml)
verbose = F

# ---- Load models ----
source("R/tree_monomolecular.R")
source("R/millennial_model_earthworm.R")
source("R/millennial_model.R")

# ---- Load config utilities ----
source("R/make_tree_forcing.R")
source("R/derive_millennial_parms.R")
source("R/init_millennial_state.R")
source("R/plot_ode_out.R")

parms  <- yaml::read_yaml("config/common.yml")
parms  <- modifyList(parms, yaml::read_yaml("config/tree_monomolecular.yml"))
parms  <- modifyList(parms, yaml::read_yaml("config/millennial.yml"))

# Consumption calculation:
# AE is the chosen site
# 2.38 gdwt/g_earthworm
# 2.06 g/m2 # Fresh biomass per area of earthworms
# cast rate: 4.9028 gdwt/m2
# 3.5% carbon in the soil
# 4.9028*0.0132*1.1 = 0.07118866 gC/m2
# Earthworms assimilate 9%: 0.1887578/(1-0.09)
# Earthworms consume: 0.2074262 gC/m2
# Mineral soil C: 2284 gC/m2
# Earthworm biomass: 0.48 gC/m2
# c = F/(Earthworm*SoilC)
# c = 0.2074262/(2284*0.48) = 0.0001892022

# Water stable aggregate improvement: based on reference Erin found, 20-35% more stable

# POM is 3x more with earthworms
# MAOM is 68% more with earthworms
# 1350mm for precipitation

parms$BD = 1390 # From AE Itatinga work
parms$MAT = 20 # From field station website
parms$T_amp = 15 # From field station website
parms$pH = 4 # Calcium chloride
parms$pct_claysilt = (100-83.8) # From George itatinga
parms$fCLAY = 11.6 # From George Itatinga
# parms$phi_por = 0.35 # Guess, fix later
parms$TBTmax = 8790*0.5*2 # Based on the Voit...

parms  <- derive_millennial_parms(parms)

parms$tree_forcing <- make_tree_forcing_equilibrium(parms)


# With detritivores:
millennial_eqm_det = rootSolve::stode(
  y     = init_millennial_state(Earthworm = T),
  func  = millennial_model_earthworm,
  parms = parms
)

# Without detritivores:
millennial_eqm = rootSolve::stode(
  y     = init_millennial_state(F),
  func  = millennial_model,
  parms = parms
)

(temp_file = cbind(EW = millennial_eqm_det$y,CT = c(millennial_eqm$y, Earthworm = 0)) %>% 
  data.frame() %>%
  rownames_to_column(var = "Pool") %>%
  tibble() %>%
  mutate(Diff = 100*(EW-CT)/CT))

temp_file %>%
  filter(Pool %in% c("P", "L", "A", "M","B")) %>%
  pivot_longer(!Pool) %>%
  filter(name != "Diff") %>%
  group_by(name) %>%
  summarize(value = sum(value)) %>%
  pivot_wider() %>%
  mutate(Diff = 100*(EW-CT)/CT)

# 3450 soil carbon from Voit.
# 780 litter layer from Voit.

parms$B0 = parms$TBTmax
parms$addherb

parms$tree_forcing <- make_tree_forcing(parms)


View(t(sapply(1:365, parms$tree_forcing)))


millennial_out_det_longrun = ode(
  times = seq(1, 365*500, by = 1),
  y     = millennial_eqm_det$y,
  func  = millennial_model_earthworm,
  parms = parms
)


millennial_out_det_longrun = ode(
  times = seq(1, 365*500, by = 1),
  y     = millennial_out_det_longrun[365*500,-1],
  func  = millennial_model_earthworm,
  parms = parms)

write_rds(millennial_out_det_longrun,"Data/millennial_out_ew_longrun.rds")

plot_ode_output(millennial_out_det_longrun,
                variable_cols = names(millennial_eqm_det$y))


millennial_out_det_longrun %>%
  data.frame() %>%
  tibble() %>%
  filter(time %in% c(
    365*500-185,
    365*500-185 - 365
  ))  %>%
  mutate(Time = c("A", "B")) %>%
  select(-time) %>%
  pivot_longer(!Time) %>%
  pivot_wider(names_from = Time) %>%
  mutate(diff = A-B) %>%
  pull(diff) %>% abs() %>% max()


write_rds(millennial_out_det_longrun[365*500,-1],"Data/millennial_out_ew_longrun.rds")


stable_state = read_rds("Data/millennial_out_ew_longrun.rds")

millennial_out_wew = ode(
  times = seq(1, 365*10, by = 1),
  y     = stable_state,
  func  = millennial_model_earthworm,
  parms = parms
)

stable_state["Earthworm"] = 0

millennial_out_noew = ode(
  times = seq(1, 365*10, by = 1),
  y     = stable_state,
  func  = millennial_model_earthworm,
  parms = parms
)

millennial_out_wew %>%
  data.frame() %>%
  tibble() %>%
  mutate(Treatment = "Earthworm") %>%
  bind_rows(
    millennial_out_noew %>%
      data.frame() %>%
      tibble() %>%
      mutate(Treatment = "No earthworm")
  ) %>%
  pivot_longer(!time & !Treatment) %>%
  ggplot(aes(x = time, y = value, color = Treatment)) + geom_line() + facet_wrap(.~name, scales = "free_y")


millennial_out_wew %>%
  data.frame() %>%
  tibble() %>%
  mutate(Treatment = "Earthworm") %>%
  bind_rows(
    millennial_out_noew %>%
      data.frame() %>%
      tibble() %>%
      mutate(Treatment = "No earthworm")
  ) %>%
  pivot_longer(!time & !Treatment) %>%
  filter(time > (365*10 - 365)) %>%
  group_by(Treatment, name) %>%
  summarize(value = mean(value)) %>%
  ggplot(aes(x = Treatment, y = value, fill = Treatment)) + geom_col() + facet_wrap(.~name, scales = "free_y")

millennial_out_wew %>%
  data.frame() %>%
  tibble() %>%
  mutate(Treatment = "Earthworm") %>%
  bind_rows(
    millennial_out_noew %>%
      data.frame() %>%
      tibble() %>%
      mutate(Treatment = "No earthworm")
  ) %>%
  pivot_longer(!time & !Treatment) %>%
  filter(time > (365*10 - 365)) %>%
  pivot_wider(names_from = Treatment) %>%
  mutate(Diff = Earthworm - `No earthworm`) %>%
  group_by(name) %>%
  summarize(Diff = mean(Diff)) %>%
  ggplot(aes(x = name, y = Diff)) + geom_col() + facet_wrap(.~name, scales = "free_y")
