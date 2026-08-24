# ============================================================
# SENSITIVITY OF ANIMAL EFFECTS TO MODEL PARAMETERS (at equilibrium)
# ------------------------------------------------------------
# Like Scripts/spinup_dynamic.R but WITHOUT the long dynamic spin-up: for each
# scenario it only computes the equilibrium animal effect, and it repeats that
# calculation while sweeping chosen model parameters, one at a time, over a
# grid. For every parameter value it records:
#   TOTAL  effect = treatment (animals, all pathways) - baseline (no animals)
#   DIRECT effect = treatment with indirect (_pint / exudate) slopes zeroed
#                   - baseline
# (both as % change vs the no-animal baseline equilibrium, per pool), then plots
# how those effects respond to each parameter.
#
# Nothing is spun up dynamically and nothing is saved to Data/spinup -- this is
# a fast, equilibrium-only sensitivity sweep. Run from the project root.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R"); source("R/derive_millennial_parms.R");
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")
source("R/scenario_uncertainty.R")

model   <- "millennial"
scen    <- read_scenarios("Data/scenarios.xlsx")
scen$MitePredator <- NULL

use_fitted_params <- TRUE
fitted_params <- if (use_fitted_params && file.exists("Results/fitted_animal_params.csv"))
  load_fitted_params("Results/fitted_animal_params.csv") else NULL

fig_dir <- "Results/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
res_dir <- "Results";         dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# animal biomass pools (excluded from the "effect on pools" output)
animal_pools <- c("Earthworm", "Detritivore", "DetPredator", "RootHerb")
derive_fn    <- match.fun(model_table[[model]]$derive)

# ------------------------------------------------------------
# set_param(): change ONE parameter on a setup object and rebuild everything
# that depends on it (derived params + climate forcing), so a sweep over a
# parameter that feeds phi_por / the forcing uses consistent values.
# ------------------------------------------------------------
set_param <- function(obj, param, value) {
  obj$parms[[param]]        <- value
  obj$parms                 <- derive_fn(obj$parms)
  obj$parms$climate_forcing <- make_climate_forcing(obj$parms)
  obj
}

# ------------------------------------------------------------
# animal_effect_eq(): the core calculation (no dynamic spin-up). Given a
# treatment/baseline pair, returns a long data frame of the TOTAL and DIRECT
# animal effect on every shared pool, at equilibrium.
# ------------------------------------------------------------
animal_effect_eq <- function(pair) {
  pair$baseline  <- spinup_equilibrium(pair$baseline, verbose = FALSE)
  pair$treatment <- spinup_equilibrium(pair$treatment,
                                       warm_start = pair$baseline$init_state_spin,
                                       verbose = FALSE)
  eq_b <- pair$baseline$init_state_spin
  eq_t <- pair$treatment$init_state_spin

  # direct-only: zero the indirect (_pint + exudate) slopes, re-equilibrate
  pair_d <- zero_indirect_effects(pair, verbose = FALSE)
  eq_d   <- spinup_equilibrium(pair_d$treatment,
                               warm_start = eq_b, verbose = FALSE)$init_state_spin

  tot <- compare_vectors(eq_t, eq_b); tot$type <- "total"
  dir <- compare_vectors(eq_d, eq_b); dir$type <- "direct"
  out <- rbind(tot, dir)
  out[!is.na(out$baseline) & !(out$name %in% animal_pools), ]     # shared, non-animal pools
}

# ------------------------------------------------------------
# sweep_param(): rebuild the scenario pair, apply fitted params, then vary ONE
# parameter over `values` and recompute the equilibrium animal effect each time.
# ------------------------------------------------------------
sweep_param <- function(scenario, param, values) {
  base_pair <- setup_scenario_pair(model, scen, scenario)
  if (!is.null(fitted_params))
    base_pair$treatment <- apply_fitted_params(base_pair$treatment, fitted_params,
                                               model, scenario, verbose = FALSE)
  default <- base_pair$treatment$parms[[param]]

  rows <- list()
  for (v in values) {
    pair <- base_pair
    pair$treatment <- set_param(pair$treatment, param, v)
    pair$baseline  <- set_param(pair$baseline,  param, v)
    eff <- tryCatch(animal_effect_eq(pair), error = function(e) {
      message("  ", scenario, " ", param, "=", signif(v, 4), " failed: ", conditionMessage(e)); NULL })
    if (is.null(eff)) next
    eff$param    <- param
    eff$value    <- v
    eff$default  <- default
    eff$scenario <- scenario
    rows[[length(rows) + 1]] <- eff
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# ============================================================
# CHOOSE what to sweep.
#   sweep_params  the model parameters to test (must exist in parms)
#   n_points      grid points per parameter
#   buffer        +/- fractional range around each parameter's default (linear)
#   scenarios     which animal scenarios to run
# ============================================================
sweep_params <- c("k_frag_litter", "k_frag_organic", "k_l_o", "k_l",
                  "k_b", "k_pa", "k_ma", "pct_claysilt",
                  "k_MICd", "k_bd", "root_to_organic", "a_root_herb",
                  "MAT", "MAtheta", "NPP_herb", "NPP_tree", "k_exudate_tree", "k_exudate_herb", "k_frag_CWD", "K_ol","alpha_ol", "alpha_ob", "K_ob", "BD", "pH", "p_c", "rho_p","psi_matric", "lambda_mat", "k_a_min", "alpha_pl", "alpha_lb", "K_pl", "K_lb", "p1", "p2", "K_ld", "CUE_T", "p_a", "p_b")
n_points  <- 7
buffer    <- 0.5
scenarios <- names(scen)

# ------------------------------------------------------------
# Pretty display names for the swept parameters (code name -> figure label),
# used to relabel the axes/legends in PLOT 1 (heat map) and PLOT 2 (overview).
# pretty_param() falls back to the raw code for anything not listed.
# ------------------------------------------------------------
param_labels <- c(
  k_frag_litter   = "Litter fragmentation",
  k_frag_organic  = "Organic fragmentation",
  k_frag_CWD      = "CWD fragmentation",
  k_l_o           = "DOM -> POC turnover",
  k_l             = "LMWC leaching",
  k_b             = "Aggregate breakdown",
  k_pa            = "POC -> aggregate",
  k_ma            = "MAOC -> aggregate",
  k_MICd          = "Microbial turnover (organic)",
  k_bd            = "Microbial turnover (mineral)",
  k_a_min         = "Min. aeration factor",
  k_exudate_tree  = "Tree root exudation",
  k_exudate_herb  = "Herb root exudation",
  pct_claysilt    = "Clay + silt (%)",
  root_to_organic = "Root death -> Organic",
  a_root_herb     = "Herb root allocation",
  MAT             = "Mean annual temperature",
  MAtheta         = "Mean annual moisture",
  NPP_herb        = "Herbaceous NPP",
  NPP_tree        = "Tree NPP",
  BD              = "Bulk density",
  rho_p           = "Particle density",
  pH              = "Soil pH",
  psi_matric      = "Matric potential",
  lambda_mat      = "Matric sensitivity",
  CUE_T           = "CUE temperature slope",
  K_ol            = "Organic->DOM half-saturation",
  K_ob            = "DOM->microbe half-saturation",
  K_pl            = "POC->LMWC half-saturation",
  K_lb            = "LMWC->microbe half-saturation",
  K_ld            = "LMWC->MAOC max sorption",
  alpha_ol        = "Organic->DOM pre-exponential",
  alpha_ob        = "DOM->microbe pre-exponential",
  alpha_pl        = "POC->LMWC pre-exponential",
  alpha_lb        = "LMWC->microbe pre-exponential",
  p1              = "pH sorption coef. 1",
  p2              = "pH sorption coef. 2",
  p_a             = "Aggregate->POC partition",
  p_b             = "Necromass->MAOC partition",
  p_c             = "Clay+silt protection coef.")

pretty_param <- function(x) ifelse(x %in% names(param_labels), param_labels[x], x)

# ------------------------------------------------------------
# PHYSICAL BOUNDS for the swept grid: parameters constrained to a fixed range
# must not be evaluated outside it. a_*/p_* (assimilation / production
# efficiencies) and the named fractions are in [0, 1]; pct_claysilt in [0, 100].
# clamp_grid() trims the +/-buffer grid to these bounds (dropping duplicates).
# ------------------------------------------------------------
param_bounds <- function(param) {
  # a_*/p_* ASSIMILATION / PRODUCTION efficiencies are in [0, 1]. The soil
  # partition/scaling coefficients p_a, p_b, p_c are NOT efficiencies and are
  # left unbounded here.
  if (grepl("^(a_|p_)", param) && !param %in% c("p_a", "p_b", "p_c"))
    return(c(0, 1))
  switch(param,
         pct_claysilt = c(0, 100),
         p_a = , p_b = c(0, 1),                # partition fractions, still [0,1]
         LigFrac = , a_root_herb = , root_to_organic = ,
         prop_feaces_earthworm_LMWC = c(0, 1),
         c(-Inf, Inf))
}
clamp_grid <- function(vals, param) {
  b <- param_bounds(param)
  v <- vals[vals >= b[1] & vals <= b[2]]
  if (!length(v)) v <- pmin(pmax(vals, b[1]), b[2])   # fallback: clip, keep points
  unique(v)
}

# run the full sweep: scenario x parameter
all_rows <- list()
for (scenario in scenarios) {
  pair0 <- setup_scenario_pair(model, scen, scenario)
  avail <- intersect(sweep_params, names(pair0$treatment$parms))
  for (param in avail) {
    d0   <- pair0$treatment$parms[[param]]
    vals <- fit_param_grid(d0, buffer = buffer, n = n_points, scale = "linear")
    vals <- clamp_grid(vals, param)
    cat("sweeping", scenario, "/", param, "around", signif(d0, 4), "\n")
    s <- sweep_param(scenario, param, vals)
    if (!is.null(s)) all_rows[[paste(scenario, param)]] <- s
  }
}
sens <- bind_rows(all_rows)
write_csv(sens, file.path(res_dir, "animal_effect_sensitivity.csv"))

temp1 = read_csv("Results/animal_effect_sensitivity.csv") %>%
  group_by(type, param, rel, scenario) %>%
  summarize(difference = sum(difference)) %>%
  filter(rel %in% c("0.5", "1.5", "1")) %>%
  mutate(rel = ifelse(rel == "0.5", "Half", 
                      ifelse(rel == "1", "Base", "Double"))) %>%
  pivot_wider(names_from = rel, values_from = difference) %>%
  filter(type == "total")

temp1 %>%
  ggplot(aes(y = param, x = Half, xend = Double, yend = param)) + geom_segment() + facet_wrap(.~scenario, scales = "free") +
  geom_vline(aes(xintercept = Base), linetype = 2) + theme_minimal(base_size = 11) +
  geom_point(aes(x = Half)) + geom_point(aes(x = Double))

temp1 %>%
  mutate(range = abs(Double-Half)) %>%
  group_by(scenario) %>%
  slice_max(range, n = 3) %>%
  ungroup() %>%
  select(param) %>% table()
  
# ------------------------------------------------------------
# SCENARIO UNCERTAINTY OVERLAY: where does each parameter's reported
# uncertainty (SD, else Min/Max, else CV = 2) sit relative to the +/-50% sweep?
# For parameters present in scenarios.xlsx we compute the reported lo/hi as a
# FRACTION of the default, so it can be drawn on the same "relative to default"
# axis as the sweep (1.0 = default; 0.5 / 1.5 = the sweep ends).
# ------------------------------------------------------------
unc <- tryCatch(read_param_uncertainty("Data/scenarios.xlsx"),
                error = function(e) { message("uncertainty read failed: ",
                                              conditionMessage(e)); NULL })
if (!is.null(unc)) {
  defaults_by <- sens %>% distinct(scenario, param, default)
  unc_rel <- unc %>%
    rename(param = parameter) %>%
    inner_join(defaults_by, by = c("scenario", "param")) %>%
    mutate(rel_lo = lo / default, rel_hi = hi / default,
           rel_value = value / default) %>%
    select(scenario, param, unc_source, rel_lo, rel_hi, rel_value)
  write_csv(unc_rel, file.path(res_dir, "animal_effect_sensitivity_uncertainty.csv"))
}

# ------------------------------------------------------------
# SUMMARY: average magnitude of the animal effect at the +/-50% ends, per
# scenario x parameter x effect type. `effect_at` = mean |effect on total C|
# summed across pools, at the low and high sweep ends and at the default.
# ------------------------------------------------------------
per_value <- sens %>%
  group_by(scenario, param, type, value, default) %>%
  summarise(effect_totalC = sum(abs(difference), na.rm = TRUE), .groups = "drop") %>%
  mutate(rel = value / default)

summary_tbl <- per_value %>%
  group_by(scenario, param, type) %>%
  summarise(
    effect_lo50  = effect_totalC[which.min(abs(rel - 0.5))],
    effect_hi50  = effect_totalC[which.min(abs(rel - 1.5))],
    effect_at_default = effect_totalC[which.min(abs(rel - 1))],
    mean_effect_pm50  = mean(c(effect_totalC[which.min(abs(rel - 0.5))],
                               effect_totalC[which.min(abs(rel - 1.5))])),
    effect_range = max(effect_totalC) - min(effect_totalC),
    .groups = "drop")
write_csv(summary_tbl, file.path(res_dir, "animal_effect_sensitivity_summary.csv"))

# ------------------------------------------------------------
# PLOT 1 - HEAT MAP highlighting the most impactful parameters: parameters
# (rows, ordered by impact) x scenario (columns), fill = mean |animal effect on
# total C| at the +/-50% sweep ends (total effect). The hottest rows are the
# parameters the animal effect is most sensitive to.
#
# Each tile is annotated with a SYMBOL showing which uncertainty source
# scenarios.xlsx provides for that parameter IN THAT SCENARIO (the source can
# differ across scenarios). Axis labels flag parameters with no reported
# uncertainty at all. Symbols (explain in the caption):
#     *   SD reported            (range = value +/- SD)
#     []  Min/Max reported       (range = [Min, Max])
#     ~   CV = 2 placeholder     (no SD/Min-Max; range = value +/- 2|value|)
#   (blank) parameter absent from scenarios.xlsx -> swept on +/-50% only; its
#           axis label is marked with a leading dot.
# ------------------------------------------------------------
# source per scenario x parameter (unc read earlier; may be NULL)
src_tbl <- if (!is.null(unc))
  unc %>% transmute(scenario, param = parameter, unc_source) else
  data.frame(scenario = character(), param = character(), unc_source = character())
src_sym <- c(SD = "*", MinMax = "[]", CV2 = "~")

heat <- summary_tbl %>% filter(type == "total") %>%
  left_join(src_tbl, by = c("scenario", "param")) %>%
  mutate(sym = ifelse(is.na(unc_source), "", src_sym[unc_source]))

# parameters that appear NOWHERE in scenarios.xlsx (no uncertainty reported)
params_no_unc <- setdiff(unique(heat$param), unique(src_tbl$param))
axis_lab <- function(p) ifelse(p %in% params_no_unc,
                               paste0("\u00b7 ", pretty_param(p)),  # leading middle dot
                               pretty_param(p))

heat <- heat %>%
  mutate(param_pretty = axis_lab(param),
         param_pretty = fct_reorder(param_pretty, mean_effect_pm50, .fun = max, .desc = FALSE))

p_heat <- ggplot(heat, aes(scenario, param_pretty, fill = mean_effect_pm50)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = signif(mean_effect_pm50, 2)), size = 2.4,
            colour = "grey15", nudge_y = 0.14) +
  geom_text(aes(label = sym), size = 3.0, colour = "grey15",
            fontface = "bold", nudge_y = -0.2) +
  scale_fill_viridis_c(option = "magma", direction = -1, trans = "sqrt") +
  labs(title = "Sensitivity of the animal effect to model parameters",
       subtitle = paste0("Fill = mean |effect on total C| at +/-50% of default (total effect).\n",
                         "Uncertainty source per cell:  * = SD,  [] = Min/Max,  ~ = CV=2 placeholder;",
                         "  leading dot = not in scenarios.xlsx."),
       x = NULL, y = NULL, fill = expression("|effect| (g C m"^-2*")")) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())
ggsave(file.path(fig_dir, "animal_effect_sensitivity_heatmap.png"),
       p_heat, width = 9, height = 11, dpi = 150)
print(p_heat)

# ------------------------------------------------------------
# PLOT 2 - overview curves WITH the reported-uncertainty range shaded. x =
# parameter relative to default; y = |effect on total C|; solid = total,
# dashed = direct; the grey band marks where scenarios.xlsx says the parameter
# actually ranges (SD / Min-Max / CV2), so a reader sees whether the swept
# range is realistic. Facet by scenario; too many params to colour, so this is
# faceted per parameter in the detail plot below -- here we keep the top movers.
# ------------------------------------------------------------
top_params <- summary_tbl %>% filter(type == "total") %>%
  group_by(param) %>% summarise(m = max(mean_effect_pm50), .groups = "drop") %>%
  slice_max(m, n = 12) %>% pull(param)

overview <- per_value %>% filter(param %in% top_params) %>%
  rename(rel_param = rel) %>%
  mutate(param_pretty = pretty_param(param))

unc_layer <- if (!is.null(unc))
  geom_rect(data = unc_rel %>% filter(param %in% top_params) %>%
                     mutate(param_pretty = pretty_param(param)),
            aes(xmin = pmax(rel_lo, 0.4), xmax = pmin(rel_hi, 1.6),
                ymin = -Inf, ymax = Inf, fill = param_pretty),
            inherit.aes = FALSE, alpha = 0.12) else NULL

p_overview <- ggplot(overview, aes(rel_param, effect_totalC, colour = param_pretty)) +
  unc_layer +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey70") +
  geom_line(linewidth = 0.7) + geom_point(size = 0.7) +
  facet_wrap(~scenario, scales = "free_y") +
  scale_colour_viridis_d(guide = guide_legend(ncol = 2)) +
  scale_fill_viridis_d(guide = "none") +
  labs(title = "Animal effect vs parameter (top movers), with reported range shaded",
       subtitle = "Shaded band = reported uncertainty from scenarios.xlsx (SD / Min-Max / CV2)",
       x = "Parameter value (relative to default)",
       y = expression("Animal effect on total C (g C m"^-2*")"),
       colour = "Parameter") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "animal_effect_sensitivity_overview.png"),
       p_overview, width = 12, height = 8, dpi = 150)
print(p_overview)

# ------------------------------------------------------------
# PLOT 2 - detail for ONE scenario: the animal effect on EACH pool as a
# parameter varies (one panel per parameter, one line per pool, total effect).
# Change `focus_scenario` to inspect others.
# ------------------------------------------------------------
for(iii in 1:4){
  focus_scenario <- scenarios[iii]
  detail <- sens %>% filter(scenario == focus_scenario, type == "total") %>%
    mutate(rel_param = value / default,
           param_pretty = pretty_param(param),
           pool_pretty  = relabel_pools(name))

  p_detail <- ggplot(detail, aes(rel_param, percent_change, colour = pool_pretty)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_line(linewidth = 0.6) +
    facet_wrap(.~param_pretty, scales = "free") +
    labs(title = paste0("Animal effect per pool vs parameter -- ", focus_scenario, " (total effect)"),
         x = "Parameter value (relative to default)", y = "Effect on pool (%)",
         colour = "Pool") +
    theme_minimal(base_size = 10) + theme(legend.position = "right")
  ggsave(file.path(fig_dir, paste0("animal_effect_sensitivity_detail_", focus_scenario, ".png")),
         p_detail, width = 12, height = 8, dpi = 150)
  print(p_detail)
}

cat("\nWrote:\n  ", file.path(res_dir, "animal_effect_sensitivity.csv"),
    "\n  ", file.path(res_dir, "animal_effect_sensitivity_summary.csv"),
    "\n  ", file.path(res_dir, "animal_effect_sensitivity_uncertainty.csv"),
    "\n  ", file.path(fig_dir, "animal_effect_sensitivity_heatmap.png"),
    "\n  ", file.path(fig_dir, "animal_effect_sensitivity_overview.png"),
    "\n  ", file.path(fig_dir, "animal_effect_sensitivity_detail_<scenario>.png"), "\n")
