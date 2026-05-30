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
# ROOT HERBIVORE MODEL
# ============================================

config_rootherb <- list(
  herb        = TRUE,
  tree        = FALSE,
  earthworm   = FALSE,
  detritivore = FALSE,
  detpredator = FALSE,
  rootherb    = TRUE
)

parms  <- modifyList(parms, list(config = config_rootherb))

state_rootherb <- build_initial_state(init_full, config_rootherb)

steady_rootherb <- stode(
  y = state_rootherb,
  func = millennial_wrapper,
  parms = parms,
)

print("FULL MODEL STEADY STATE:")
print(steady_rootherb$y)

# ============================================
# NO ROOT HERBIVORE MODEL
# ============================================

config_no_rootherb <- list(
  herb        = TRUE,
  tree        = FALSE,
  earthworm   = FALSE,
  detritivore = FALSE,
  detpredator = FALSE,
  rootherb    = FALSE
)

parms$config <- config_no_rootherb

state_no_rootherb <- build_initial_state(init_full, config_no_rootherb)

steady_no_rootherb <- stode(
  y = state_no_rootherb,
  func = millennial_wrapper,
  parms = parms
)

print("NO ROOT HERBIVORE STEADY STATE:")
print(steady_no_rootherb$y)

comp_eqm(v1 = steady_rootherb$y,
         v2 = steady_no_rootherb$y)


#Run the model:

parms$climate_forcing <- make_climate_forcing(parms)

millennial_out = ode(
  times = seq(1, 365*500, by = 1),
  y     = read_rds("Results/stable_y.rds"), #steady_rootherb$y,
  func  = millennial_wrapper,
  parms = parms
)

plot_ode_output(millennial_out)

# looking at the last few years to check for stability
tibble(data.frame(millennial_out)) %>%
  filter(
    time %in% c(
      max(time),
      max(time)-365,
      max(time)-365*2
    )
  ) %>%
  mutate(time = paste0("t", time)) %>%
  pivot_longer(!time) %>%
  pivot_wider(names_from = time) %>%
  print(n = Inf)


# ending parameters for model
stable_y <- as.data.frame(millennial_out) %>%
  slice_tail(n = 1) %>%
  select(-time) %>%
  unlist()

write_rds(stable_y, "Results/stable_y.rds")

fullModel <- stable_y
noHerbivore    <- replace(stable_y, "RootHerb", 0) #setting herbivore = 0 

# Run 10 year scenarios
times2 <- seq(1, 365*10, by = 1)

out_full <- ode(y = fullModel,  
                times = times2, 
                func = millennial_wrapper, 
                parms = parms
)

write.csv(out_full, "Results/out_full.csv")

out_noHerbivore <- ode(y = noHerbivore,
                       times = times2, 
                       func = millennial_wrapper, 
                       parms = parms
)

write.csv(out_noHerbivore, "Results/out_noHerbivore.csv")


# Run the three scenarios - 100 years
times3 <- seq(1, 365*100, by = 1)

out_full100 <- ode(y = fullModel,  
                   times = times3, 
                   func = millennial_wrapper, 
                   parms = parms
)

write.csv(out_full100, "out_full100.csv")

out_noHerbivore100 <- ode(y = noHerbivore,
                          times = times3,
                          func = millennial_wrapper, 
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
    pivot_longer(-time, names_to = "state", values_to = "value") |>
    mutate(scenario = scenario)
}

# Base plot from function (full food web)
p <- plot_ode_output(out_full)

# other two scenarios
p <- p +
  geom_line(data = to_long(out_full, "Root Herbivory"), 
            aes(time, value, color = scenario)) +
  geom_line(data = to_long(out_noHerbivore, "No Herbivores"),
            aes(time, value, color = scenario)) +
  scale_color_manual(values = c("Root Herbivory" = "black",
                                "No Herbivores"   = "blue")) +
  labs(color = "Scenario")
p
ggsave("Plots/simulations.png",  plot = p,  width = 9.5, height = 6, dpi = 300)

#plot the 100 year simulation scenarios, but 
# Keep only one day per year (default: day 1)
keep_yearly <- function(out, day_of_year = 1) {
  out[out[, "time"] %% 365 == day_of_year, ]
}

p2 <- plot_ode_output(keep_yearly(out_full100))

p2 <- p2 +
  geom_line(data = to_long(keep_yearly(out_full100), "Root Herbivory"), 
            aes(time, value, color = scenario)) +
  geom_line(data = to_long(keep_yearly(out_noHerbivore100), "No Herbivores"),
            aes(time, value, color = scenario)) +
  scale_color_manual(values = c("Root Herbivory" = "black",
                                "No Herbivores"   = "blue")) +
  labs(color = "Scenario")
p2

ggsave("Plots/simulations_100yrs.png",  plot = p2,  width = 9.5, height = 6, dpi = 300)


# Combine the three scenario outputs into one long data frame
to_long <- function(out, scenario) {
  as.data.frame(out) |>
    pivot_longer(-time, names_to = "Pool", values_to = "Value") |>
    mutate(Scenario = scenario)
}

scenario_long <- bind_rows(
  to_long(out_full,      "Root Herbivory"),
  to_long(out_noHerbivore,    "No Herbivores")
)

scenario_colors <- c("Root Herbivory" = "blue",
                     "No Herbivores" = "lightgreen")


# ---- Fig 1: endpoint (max time) difference from Full model ----
full_endpoint <- scenario_long %>%
  filter(Scenario == "Root Herbivory", time == max(time)) %>%
  select(Pool, full_val = Value)

diff_end <- scenario_long %>%
  filter(Scenario != "Root Herbivory", time == max(time)) %>%
  left_join(full_endpoint, by = "Pool") %>%
  mutate(Diff = Value - full_val)

fig_end <- ggplot(diff_end, aes(x = Scenario, y = Diff, fill = Scenario)) +
  geom_col(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  facet_wrap(~ Pool, scales = "free_y") +
  scale_fill_manual(values = scenario_colors) +
  labs(y = expression(Delta ~ "endpoint pool size"),
       x = NULL, fill = NULL) +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

fig_end
ggsave("Plots/fig_end.png",  plot = fig_end,  width = 8, height = 6, dpi = 300)

# ---- Fig 2: whole-trajectory mean difference from Full model ----
full_mean <- scenario_long %>%
  filter(Scenario == "Root Herbivory") %>%
  group_by(Pool) %>%
  summarise(full_mean = mean(Value), .groups = "drop")

diff_mean <- scenario_long %>%
  filter(Scenario != "Root Herbivory") %>%
  group_by(Scenario, Pool) %>%
  summarise(mean_value = mean(Value), .groups = "drop") %>%
  left_join(full_mean, by = "Pool") %>%
  mutate(Diff = mean_value - full_mean)

fig_mean <- ggplot(diff_mean, aes(x = Scenario, y = Diff, fill = Scenario)) +
  geom_col(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  facet_wrap(~ Pool, scales = "free_y") +
  scale_fill_manual(values = scenario_colors) +
  labs(y = expression(Delta ~ "mean pool size"),
       x = NULL, fill = NULL) +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

fig_mean
ggsave("fig_mean.png", plot = fig_mean, width = 8, height = 6, dpi = 300)


#effect sizes

# Per-scenario, per-pool summary (endpoint + time-mean)
summary_df <- scenario_long %>%
  group_by(Scenario, Pool) %>%
  summarise(
    endpoint  = Value[time == max(time)],
    time_mean = mean(Value),
    .groups = "drop"
  )

# Baseline values (Full model)
baseline <- summary_df %>%
  filter(Scenario == "Root Herbivory") %>%
  select(Pool, baseline_end = endpoint, baseline_mean = time_mean)

# Effect sizes for the removal scenarios
effect_sizes <- summary_df %>%
  filter(Scenario != "Root Herbivory") %>%
  left_join(baseline, by = "Pool") %>%
  mutate(
    # Endpoint effects
    delta_end       = endpoint - baseline_end,
    pct_end         = 100 * delta_end / baseline_end,
    LRR_end         = log(endpoint / baseline_end),
    # Time-averaged effects
    delta_mean      = time_mean - baseline_mean,
    pct_mean        = 100 * delta_mean / baseline_mean,
    LRR_mean        = log(time_mean / baseline_mean)
  ) %>%
  select(Scenario, Pool,
         baseline_end, endpoint, delta_end, pct_end, LRR_end,
         baseline_mean, time_mean, delta_mean, pct_mean, LRR_mean)

effect_sizes
write.csv(effect_sizes, "effect_sizes.csv", row.names = FALSE)

effect_sizes <- effect_sizes %>%
  filter(!(Scenario == "No Herbivores" & Pool == "RootHerb"),
         !(Scenario == "No Herbivores" & Pool == "CWD"))

library(scales)

percent_change <- ggplot(effect_sizes, aes(x = Scenario, y = Pool, fill = pct_mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%+.1f%%", pct_mean)), size = 3) +
  scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#D7191C",
    midpoint = 0,
    limits = c(-50, 50),     # cap the color scale here
    oob = squish,            # anything beyond limits gets the extreme color
    name = "% change\n(time-mean)"
  ) +
  labs(x = NULL, y = NULL,
       title = "Effect of soil animals on soil C pools") +
  theme_minimal() +
  theme(panel.grid = element_blank())

percent_change

ggsave("Plots/percent_change.png", plot = percent_change, width = 4.5, height = 7, dpi = 300)