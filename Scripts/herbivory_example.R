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
  times = seq(1, 365*1500, by = 1),
  y     = init_millennial_state(HerbNem = T),
  func  = millennial_model_herbnem,
  parms = parms
)

plot_ode_output(millennial_herbnem,
                variable_cols = names(init_millennial_state(HerbNem = T)))

head(millennial_herbnem)
tail(millennial_herbnem)

write.csv(millennial_herbnem, "millennial_herbnem.csv")


# End Timepoint/Equilibrium
times <- seq(1, 365*1500, by = 1)

# looking at the last few years to check for stability
timeCheck_stable <- tibble(data.frame(millennial_herbnem)) %>%
  filter(
    time %in% c(
      max(times),
      max(times)-365,
      max(times)-365*2,
      max(times)-365*3
    )
  ) 

# ending parameters for model
timeEND <- as.data.frame(millennial_herbnem) %>%
  slice_tail(n = 1) %>%
  select(-time) %>%
  unlist()

write.csv(timeEND, "timeEND_stable_herbivore.csv")
  

# new simulation
times2 <- seq(1, 365*10, by = 1)

new_herbnem <- ode(
  y = timeEND,
  times = times2,
  func = millennial_model_herbnem,
  parms = parms
)

plot_ode_output(new_herbnem,
                variable_cols = names(init_millennial_state(HerbNem = T)))

write.csv(new_herbnem, "new_herbnem.csv")


# Simulate the scenarios (with and without nematodes)
timeEND2 <- as.data.frame(new_herbnem) %>%
  slice_tail(n = 1) %>%
  select(-time) %>%
  unlist()

write.csv(timeEND, "timeEND2_scenarios_herbivore.csv")


fullModel <- timeEND2
noHerbivore    <- replace(timeEND2, "RootHerb", 0) #setting herbivore = 0 

# Run 10 year scenarios
out_full <- ode(y = fullModel,  
                times = times2, 
                func = millennial_model_herbnem, 
                parms = parms
                )

write.csv(out_full, "out_full.csv")

out_noHerbivore <- ode(y = noHerbivore,
                       times = times2, 
                       func = millennial_model_herbnem, 
                       parms = parms
                       )

write.csv(out_noHerbivore, "out_noHerbivore.csv")


# Run the three scenarios - 100 years
times3 <- seq(1, 365*100, by = 1)

out_full100 <- ode(y = fullModel,  
                   times = times3, 
                   func = millennial_model_herbnem, 
                   parms = parms
                   )

write.csv(out_full100, "out_full100.csv")

out_noHerbivore100 <- ode(y = noHerbivore,
                          times = times3,
                          func = millennial_model_herbnem, 
                          parms = parms
                          )

write.csv(out_noHerbivore100, "out_noHerbivore100.csv")



### FROM JANEY ##
vars <- colnames(out_full)[2:14]


# Helper to reshape an ode output the same way the function does internally
to_long <- function(out, scenario) {
  df <- as.data.frame(out)
  names(df)[1] <- "time"
  df |>
    select(time, all_of(vars)) |>
    pivot_longer(-time, names_to = "state", values_to = "value") |>
    mutate(scenario = scenario)
}

# Base plot from function (full food web)
p <- plot_ode_output(out_full, variable_cols = vars)

# other two scenarios
p <- p +
  geom_line(data = to_long(out_full, "Full food web"), 
            aes(time, value, color = scenario)) +
  geom_line(data = to_long(out_noHerbivore, "No predator"),
            aes(time, value, color = scenario)) +
  scale_color_manual(values = c("Full food web" = "black",
                                "No Herbivores"   = "blue")) +
  labs(color = "Scenario")
p
ggsave("outputs/simulations.png",  plot = p,  width = 9.5, height = 6, dpi = 300)

#plot the 100 year simulation scenarios, but 
# Keep only one day per year (default: day 1)
keep_yearly <- function(out, day_of_year = 1) {
  out[out[, "time"] %% 365 == day_of_year, ]
}

p2 <- plot_ode_output(keep_yearly(out_full2), variable_cols = vars)

p2 <- p2 +
  geom_line(data = to_long(keep_yearly(out_full2),      "Full food web"), aes(time, value, color = scenario)) +
  geom_line(data = to_long(keep_yearly(out_noPred2),    "No predator"),   aes(time, value, color = scenario)) +
  geom_line(data = to_long(keep_yearly(out_noPredDet2), "No pred + det"), aes(time, value, color = scenario)) +
  scale_color_manual(values = c("Full food web" = "black",
                                "No predator"   = "blue",
                                "No pred + det" = "lightgreen")) +
  labs(color = "Scenario")
p2

#ggsave("outputs/simulations_100yrs.png",  plot = p2,  width = 9.5, height = 6, dpi = 300)


# Combine the three scenario outputs into one long data frame
to_long <- function(out, scenario) {
  as.data.frame(out) |>
    pivot_longer(-time, names_to = "Pool", values_to = "Value") |>
    mutate(Scenario = scenario)
}

scenario_long <- bind_rows(
  to_long(out_full,      "Full model"),
  to_long(out_noPred,    "No predators"),
  to_long(out_noPredDet, "No predators/detritivores")
)

scenario_colors <- c("No predators"              = "blue",
                     "No predators/detritivores" = "lightgreen")
