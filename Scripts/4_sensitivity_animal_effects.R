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
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl, ggrepel, ggpubr)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R"); source("R/derive_millennial_parms.R");
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")
source("R/scenario_uncertainty.R")
source("R/morris_sensitivity.R")

model   <- "millennial"
scen    <- read_scenarios("Data/scenarios.xlsx")

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
                  "k_b", "k_pa", "k_ma","k_litterfall_ann", "k_litterfall_herb_ann",
                  "k_MICd", "k_bd", "root_to_organic", "a_root_herb", "k_exudate_tree", "k_exudate_herb", "k_frag_CWD", "K_ol","alpha_ol", "alpha_ob", "K_ob", "p_c", "rho_p","psi_matric", "lambda_mat", "k_a_min", "alpha_pl", "alpha_lb", "K_pl", "K_lb", "p1", "p2", "K_ld", "CUE_T", "p_a", "p_b", "k_mort_root_tree", "k_mort_root_herb", "root_to_organic")
n_points  <- 3
buffer    <- 0.1
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
  p_c             = "Clay+silt protection coef.",
  k_litterfall_ann= "Tree annual litterfall prop.",
  k_litterfall_herb_ann = "Herb. annual litterfall prop.",
  k_mort_root_tree= "Tree root mortality",
  k_mort_root_herb= "Herb. root mortality")

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
         winter_root_act_prop = c(0,1),
         root_to_organic = c(0,1),
         pct_claysilt = c(0, 100),
         k_litterfall_ann = c(0,1),
         k_litterfall_herb_ann = c(0,1),
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


sens %>%
  tibble() %>%
  mutate(rel = signif((value-default)/default),1) %>%
  filter(type == "total") %>%
  select(name, difference, param, scenario, rel) %>%
  group_by(param, scenario, rel) %>%
  summarize(difference = sum(difference)) %>%
  left_join(
    tibble(rel = c(-0.1, 0, 0.1),
           rel2 = c("L", "B", "U"))
  ) %>%
  select(-rel) %>%
  pivot_wider(names_from = rel2, values_from = difference) %>%
  mutate(Lchange = (L-B)/B,
         Uchange = (U-B)/B) %>%
  mutate(Uchange = replace_na(Uchange)) %>%
  mutate(change = (abs(Lchange) + abs(Uchange))/2) %>%
  filter(change >= 0.1) %>%
  select(param, scenario, change) %>%
  mutate(change = round(100*change)) %>%
  arrange(-change) %>%
  mutate(param = pretty_param(param)) %>%
  pivot_wider(names_from = scenario, values_from = change) %>%
  write_csv(file.path(res_dir, "animal_effect_sensitivity_MStable.csv"))

# ============================================================
# MORRIS ELEMENTARY-EFFECTS SCREENING (parameters varied in COMBINATION)
# ------------------------------------------------------------
# The one-at-a-time sweep above moves each parameter alone. Morris instead moves
# the parameters together along random trajectories and reports, per parameter:
#   mu_star = mean |elementary effect|  -> TOTAL sensitivity (rank on this)
#   sigma   = sd of elementary effect   -> interactions / non-linearity
# Ranges are the same +/- buffer around each default used by the sweep, clamped
# to the physical guard rails. Output: Results/animal_effect_morris.csv.
# ============================================================
morris_r      <- 20     # number of trajectories (raise for stable estimates)
morris_levels <- 4L     # grid levels p in the Morris design

# scalar output: total animal effect on total soil + root C at equilibrium
effect_scalar <- function(base_pair, scenario, pv) {
  pair <- base_pair
  for (nm in names(pv)) {
    pair$treatment <- set_param(pair$treatment, nm, pv[[nm]])
    pair$baseline  <- set_param(pair$baseline,  nm, pv[[nm]])
  }
  eff <- animal_effect_eq(pair)                 # long: shared soil pools, total+direct
  tot <- eff[eff$type == "total", ]
  sum(tot$difference, na.rm = TRUE)             # net effect on total soil C (g C m-2)
}

set.seed(20260826)
morris_rows <- list()
for (scenario in scenarios) {
  base_pair <- setup_scenario_pair(model, scen, scenario)
  if (!is.null(fitted_params))
    base_pair$treatment <- apply_fitted_params(base_pair$treatment, fitted_params,
                                               model, scenario, verbose = FALSE)
  pnames <- intersect(sweep_params, names(base_pair$treatment$parms))

  # per-parameter bounds = default +/- buffer, clamped to physical guard rails
  lo <- hi <- setNames(numeric(length(pnames)), pnames)
  for (p in pnames) {
    d0 <- base_pair$treatment$parms[[p]]
    b  <- param_bounds(p)
    lo[p] <- max(d0 * (1 - buffer), b[1])
    hi[p] <- min(d0 * (1 + buffer), b[2])
    if (!is.finite(lo[p]) || !is.finite(hi[p]) || hi[p] <= lo[p]) { lo[p] <- hi[p] <- NA }
  }
  keep <- names(lo)[is.finite(lo) & is.finite(hi)]
  lo <- lo[keep]; hi <- hi[keep]

  cat("Morris:", scenario, "-", length(keep), "parameters,", morris_r, "trajectories\n")
  trajs <- morris_trajectories(length(keep), morris_r, levels = morris_levels)
  ee    <- morris_run(trajs, lo, hi,
                      eval_fun = function(pv) effect_scalar(base_pair, scenario, pv),
                      verbose = FALSE)
  sm    <- morris_summary(ee)
  sm$scenario <- scenario
  morris_rows[[scenario]] <- sm
}

morris_tbl <- bind_rows(morris_rows) %>%
  mutate(parameter_label = pretty_param(parameter)) %>%
  group_by(scenario) %>% arrange(scenario, desc(mu_star)) %>%
  mutate(rank = row_number()) %>% ungroup() %>%
  select(scenario, parameter, parameter_label, mu_star, sigma, n_ee, rank)
write_csv(morris_tbl, file.path(res_dir, "animal_effect_morris.csv"))

# mu_star vs sigma diagnostic plot (interaction map), faceted by scenario
p_morris <- ggplot(morris_tbl, aes(mu_star, sigma)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(colour = "#2166ac", alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = parameter)) +
  facet_wrap(~scenario, scales = "free") +
  labs(title = "Morris screening of model parameters (varied in combination)",
       subtitle = "mu* = total sensitivity (mean |EE|); sigma = interactions / non-linearity",
       x = expression(mu*"* (mean |elementary effect|, g C m"^-2*")"),
       y = expression(sigma*" (sd of elementary effect)")) +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "animal_effect_morris.png"), p_morris,
       width = 12, height = 8, dpi = 150)
print(p_morris)

cat("\nWrote:\n  ", file.path(res_dir, "animal_effect_morris.csv"),
    "\n  ", file.path(fig_dir, "animal_effect_morris.png"), "\n")
