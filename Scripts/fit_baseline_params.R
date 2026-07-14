# ============================================================
# BASELINE (no-animal) PARAMETER FITTING + SENSITIVITY
# ------------------------------------------------------------
# Everything here runs on the BASELINE arm (animals off), so it only touches
# the plant + soil parameters -- the ones NOT associated with the animals.
# (In the baseline the animal pools are zero, so animal parameters have no
# effect at all; they are excluded explicitly, see animal_params below.)
#
#   PART 1  SCAN  - sweep each parameter over a grid, re-equilibrate the
#                   baseline, and record EVERY state variable. One clean
#                   faceted plot of all state variables over those parameters.
#   PART 2  FIT   - tune one parameter to hit a target value for one pool
#                   (e.g. total SOC), by root-finding on the equilibrium.
#
# Run from project root.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")

model    <- "millennial"
scenario <- "Earthworm"          # any scenario: the BASELINE arm has no animals
scen     <- read_scenarios("Data/scenarios.xlsx")

scen$MitePredator = NULL

# ------------------------------------------------------------
# Baseline setup (animals OFF) -- the object every run below starts from.
# ------------------------------------------------------------

baseline_eq = list()
for(scenario in names(scen)){
  base <- setup_scenario(model, scen, scenario, animals = FALSE)
  base <- spinup_equilibrium(base, verbose = FALSE)
  eq0  <- base$init_state_spin                     # default-parameter equilibrium
  baseline_eq[[length(baseline_eq) + 1]] = cbind(compare_vectors(base$init_state_spin, base$init_state), scenario = scenario)
}

do.call("rbind", baseline_eq) %>% filter(!is.na(treatment)) %>%
  tibble() %>%
  mutate(pretty_name = name_lookup[name]) %>%
  mutate(pretty_name = factor(pretty_name, levels = plot_order)) %>%
  ggplot(aes(x = scenario, y = treatment, fill = pretty_name)) + geom_col()

do.call("rbind", baseline_eq) %>% filter(!is.na(treatment)) %>%
  tibble() %>%
  mutate(pretty_name = name_lookup[name]) %>%
  mutate(pretty_name = factor(pretty_name, levels = plot_order)) %>%
  ggplot(aes(x = scenario, y = percent_change, fill = pretty_name)) + geom_col(position = "dodge") + geom_hline(yintercept = c(-100,100), linetype = 2)

# ------------------------------------------------------------
# WHICH PARAMETERS ARE "NOT ASSOCIATED WITH THE ANIMALS"
# Animal parameters are listed explicitly (safer than pattern-matching: note
# a_root_herb / k_mort_root_herb / k_exudate_herb are PLANT parameters for the
# herbaceous roots, while a_rootherb / E_root_herb belong to the root
# HERBIVORE). Everything else that is a single finite number is fittable.
# ------------------------------------------------------------
animal_params <- c(
  "d_detritivores","a_detritivores","p_detritivores","E_detritivores","c_detritivores",
  "d_detpredator","a_detpredator","p_detpredator","E_detpredator","c_detpredator",
  "d_earthworm","a_earthworm","a_earthworm_soil","p_earthworm","E_earthworm",
  "c_earthworm_litter","c_earthworm_soil","c_earthworm_om","prop_feaces_earthworm_LMWC",
  "c_rootherb","d_rootherb","a_rootherb","p_rootherb","E_root_herb",
  "adj_earthworm","adj_detritivores","adj_detpredator","adj_rootherb",
  indirect_effect_params)                        # the animal-effect slopes (setup.R)

# structural / non-numeric parameters that must not be swept
non_tunable <- c("climate_forcing", "kinetics", "year_length_d", "N_theta_peaks",
                 "litter_peak_doy", "litter_width_d", "phi_por")   # phi_por is DERIVED

# CONSTRAINED: tree allocation must sum to 1 (derive_millennial_parms() stops
# otherwise), so these cannot be swept one at a time -- vary them only as a set
# that still sums to 1.
constrained <- c("a_leaf_tree", "a_wood_tree", "a_root_tree")

is_scalar_num <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
fittable_params <- setdiff(
  names(base$parms)[vapply(base$parms, is_scalar_num, logical(1))],
  c(animal_params, non_tunable, constrained))
cat("fittable (non-animal) parameters:", length(fittable_params), "\n")
print(fittable_params)

# ------------------------------------------------------------
# set_param(): change ONE parameter and rebuild everything that depends on it.
# Critical: derived parameters (phi_por <- 1 - BD/rho_p) and the climate forcing
# are functions of the parameter list, so they must be recomputed -- otherwise a
# sweep over e.g. BD or MAT silently uses stale derived values.
# ------------------------------------------------------------
derive_fn <- match.fun(model_table[[model]]$derive)
set_param <- function(obj, param, value) {
  obj$parms[[param]]        <- value
  obj$parms                 <- derive_fn(obj$parms)          # re-derive (phi_por, ...)
  obj$parms$climate_forcing <- make_climate_forcing(obj$parms)
  obj
}

# ------------------------------------------------------------
# PART 1 - SCAN: sweep a parameter, equilibrate, record ALL state variables.
# Warm-starts each step from the previous equilibrium (faster + more robust).
# ------------------------------------------------------------
scan_baseline_param <- function(obj, param, values, warm = TRUE, verbose = TRUE) {
  out <- list(); ws <- obj$init_state_spin
  for (v in values) {
    o <- tryCatch(spinup_equilibrium(set_param(obj, param, v),
                                     warm_start = if (warm) ws else NULL,
                                     verbose = FALSE),
                  error = function(e) NULL)
    if (is.null(o) || is.null(o$init_state_spin)) {
      if (verbose) message("  ", param, " = ", signif(v, 4), " : FAILED")
      next
    }
    eq   <- o$init_state_spin
    conv <- isTRUE(o$spin_info$converged) && all(is.finite(eq))
    if (conv && warm) ws <- eq                     # warm-start only from good states
    out[[length(out) + 1]] <- data.frame(
      param = param, value = v, state = names(eq), amount = as.numeric(eq),
      converged = conv, stringsAsFactors = FALSE)
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# --- parameters to sweep (edit freely; must be in `fittable_params`) ---
sweep_params <- c("NPP_tree", "k_frag_litter", "k_mort_wood_tree")
sweep_params <- intersect(sweep_params, fittable_params)

n_points <- 9        # values per parameter
buffer   <- 0.3      # +/- 30% of the default (linear); log grid if you prefer

scan_all <- list()
for (p in sweep_params) {
  d0   <- base$parms[[p]]
  vals <- fit_param_grid(d0, buffer = buffer, n = n_points, scale = "linear")
  cat("scanning", p, "around", signif(d0, 4), "\n")
  s <- scan_baseline_param(base, p, vals)
  if (is.null(s)) { message("  no converged runs for ", p); next }
  s$default <- d0
  scan_all[[p]] <- s
}
scan <- bind_rows(scan_all)
write_csv(scan, file.path(res_dir, "baseline_param_scan.csv"))

# ------------------------------------------------------------
# PLOT - all state variables over all swept parameters, in one clean figure.
# x is the parameter RELATIVE to its default (so parameters on wildly different
# scales share one axis); y is the pool relative to its default equilibrium
# (so pools spanning orders of magnitude share one axis). One panel per state
# variable, one coloured line per parameter.
# ------------------------------------------------------------
eq0_df <- tibble(state = names(eq0), eq_default = as.numeric(eq0))

plot_df <- scan %>%
  filter(converged, state != "mass_balance_check") %>%
  left_join(eq0_df, by = "state") %>%
  mutate(rel_param = value / default,
         rel_state = ifelse(abs(eq_default) > 1e-8, amount / eq_default, NA_real_),
         state = factor(state, levels = names(eq0)))

p_all <- ggplot(plot_df, aes(rel_param, rel_state, colour = param)) +
  geom_hline(yintercept = 1, linewidth = 0.3, colour = "grey70") +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey70") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 0.9) +
  facet_wrap(~state, scales = "free_y") +
  scale_colour_brewer(palette = "Dark2") +
  labs(title = "Millennial baseline (no animals): equilibrium sensitivity to non-animal parameters",
       subtitle = "Each panel = one state variable. x = parameter / its default, y = pool / its default equilibrium",
       x = "Parameter value (relative to default)",
       y = "Equilibrium pool (relative to default)",
       colour = "Parameter") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "baseline_param_scan.png"), p_all,
       width = 12, height = 8, dpi = 150)
print(p_all)

# Absolute-units version, one page per parameter (uncomment to also write these)
# for (p in unique(scan$param)) {
#   pp <- ggplot(filter(scan, param == p, converged, state != "mass_balance_check"),
#                aes(value, amount)) +
#     geom_line() + geom_point(size = 0.9) +
#     facet_wrap(~state, scales = "free_y") +
#     labs(title = paste("Baseline equilibrium vs", p), x = p, y = "Equilibrium pool") +
#     theme_minimal()
#   ggsave(file.path(fig_dir, paste0("baseline_scan_", p, ".png")), pp,
#          width = 11, height = 7, dpi = 150)
# }

# ------------------------------------------------------------
# PART 2 - FIT one non-animal parameter to a target.
#   target_pool  a state variable name, or "SOC" for the summed soil pools
#   target_value the equilibrium value you want it to take
# Root-finds on the equilibrium; returns the fitted value and the achieved target.
# ------------------------------------------------------------
soil_pools <- setdiff(names(eq0),
                      c("C_leaf_herb","C_root_herb","C_leaf_tree","C_wood_tree","C_root_tree"))

fit_baseline_param <- function(obj, param, target_pool, target_value,
                               buffer = 1, scale = c("log", "linear"), verbose = TRUE) {
  scale <- match.arg(scale)
  if (!param %in% fittable_params)
    stop("'", param, "' is not a non-animal fittable parameter.")
  x0 <- obj$parms[[param]]
  ws <- obj$init_state_spin

  value_of <- function(eq) if (identical(target_pool, "SOC"))
    sum(eq[soil_pools], na.rm = TRUE) else unname(eq[target_pool])

  f <- function(v) {
    o <- tryCatch(spinup_equilibrium(set_param(obj, param, v),
                                     warm_start = ws, verbose = FALSE),
                  error = function(e) NULL)
    if (is.null(o) || !isTRUE(o$spin_info$converged)) return(NA_real_)
    value_of(o$init_state_spin) - target_value
  }

  step <- if (scale == "log") abs(x0) * 0.25 else max(abs(x0) * 0.25, 1e-8)
  br   <- .bracket_root(f, x0, step)
  if (is.null(br))
    stop("Could not bracket a root for ", param, " -> ", target_pool, " = ", target_value,
         ". Widen the range, or check the target is achievable (see the scan).")

  r   <- uniroot(f, interval = br, tol = 1e-8)
  fit <- spinup_equilibrium(set_param(obj, param, r$root), warm_start = ws, verbose = FALSE)
  ach <- value_of(fit$init_state_spin)
  if (verbose)
    cat(sprintf("fit %s: %.6g -> %.6g   (%s: target %.4g, achieved %.4g)\n",
                param, x0, r$root, target_pool, target_value, ach))
  list(param = param, default = x0, fitted = r$root,
       target_pool = target_pool, target = target_value, achieved = ach,
       obj = fit)
}

# --- example: tune litter fragmentation so total SOC hits a measured value ---
# (target here = 10% above the default SOC, just as a runnable demonstration)
soc0 <- sum(eq0[soil_pools])
cat("\ndefault baseline SOC =", signif(soc0, 5), "\n")

fit1 <- fit_baseline_param(base, "k_frag_litter", "SOC", soc0 * 1.10)

fit_tbl <- tibble(param = fit1$param, default = fit1$default, fitted = fit1$fitted,
                  target_pool = fit1$target_pool, target = fit1$target,
                  achieved = fit1$achieved)
write_csv(fit_tbl, file.path(res_dir, "baseline_param_fit.csv"))
print(fit_tbl)

cat("\nWrote:\n  ", file.path(res_dir, "baseline_param_scan.csv"),
    "\n  ", file.path(res_dir, "baseline_param_fit.csv"),
    "\n  ", file.path(fig_dir, "baseline_param_scan.png"), "\n")
