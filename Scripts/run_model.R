# ============================================================
# RUN ONE MODEL — manual, linear, transparent
# ------------------------------------------------------------
# Pick a model in the SELECT block below, then step through the script.
# Switch plants/animals OFF by setting their initial values to 0 (a pool
# with zero biomass produces zero flux, so it stays off). No config DSL.
#
# Run from the project ROOT (so the relative paths work).
# ============================================================

library(pacman)
p_load(deSolve, rootSolve, tidyverse, yaml)

# ---- shared utilities ----
source("R/climate_forcing.R")
source("R/plot_ode_output.R")
source("R/spinup.R")
source("R/make_model_wrapper.R")

# ============================================================
# SELECT A MODEL  — change these 5 lines, nothing else.
# ============================================================
## --- Millennial ---
source("R/millennial_model.R"); source("R/derive_millennial_parms.R"); source("R/init_millennial_state.R")
model_fn   <- millennial_model_wplant
model_yaml <- "config/millennial.yml"
derive_fn  <- derive_millennial_parms
init_state <- init_millennial_state()

## --- Century ---
source("R/century_model.R"); source("R/derive_century_parms.R"); source("R/init_century_state.R")
model_fn   <- century_model
model_yaml <- "config/century.yml"
derive_fn  <- derive_century_parms
init_state <- init_century_state()

## --- MIMICS ---
source("R/MIMICS_model.R"); source("R/derive_MIMICS_parms.R"); source("R/init_MIMICS_state.R")
model_fn   <- MIMICS_model
model_yaml <- "config/MIMICS.yml"
derive_fn  <- derive_MIMICS_parms
init_state <- init_MIMICS_state()

# ============================================================
# PARAMETERS  (common + model-specific, then derive)
# ============================================================
parms <- yaml::read_yaml("config/common.yml")
parms <- modifyList(parms, yaml::read_yaml(model_yaml))
parms <- derive_fn(parms)

# ============================================================
# TURN GROUPS ON/OFF  — just zero the initial values you don't want.
# (Keep at least one plant group on so there is litter input.)
# Examples — uncomment what you need:
# ============================================================
# init_state[c("C_leaf_tree","C_wood_tree","C_root_tree")] <- 0   # no trees
# init_state[c("C_leaf_herb","C_root_herb")]               <- 0   # no herbs
init_state["Earthworm"]   <- 0                                  # no earthworms
init_state["Detritivore"] <- 0                                  # no detritivores
init_state["DetPredator"] <- 0                                  # no predators
init_state["RootHerb"]    <- 0                                  # no root herbivores

# ============================================================
# Create the wrapped model with the correct state variables
# ============================================================
wrapped_model <- make_model_wrapper(
  model_fun   = model_fn,
  full_names  = names(init_state),
  state_groups = names(init_state)[which(init_state > 0)]
)

# Keep only correct variables:
working_state = init_state[which(init_state > 0)]

print(working_state)

# ============================================================
# QUICK MASS-BALANCE CHECK at t = 0 (should be ~0)
# ============================================================
parms$climate_forcing <- make_climate_forcing(parms)
mb0 <- wrapped_model(0, working_state, parms)[[2]]
cat(sprintf("mass_balance_check at t=0: %.3e\n", mb0))

# ============================================================
# Fast warm-start to constant-forcing steady state.
# ============================================================

cfss <- stode(y = working_state,
              func = wrapped_model,
              parms = parms)

if(attr(cfss, "steady")){
  cat("Result is steady")
  init_state_spin = cfss[[1]]
  
}else{
  cat("Result did not reach stability. Using input state.")
  init_state_spin = working_state
}

# ============================================================
# SPIN-UP  — ONE continuous run. Increase n_years until stable.
#   by = 1   -> daily output (see within-year dynamics)
#   by = 365 -> yearly snapshots (lighter output for very long runs)
# ============================================================

spinup_until_stable <- function(init_state, parms,
                                model_fn = wrapped_model,
                                n_years = 50,
                                by = 1,
                                max_iter = 10,
                                tol = 1e-4,
                                verbose = TRUE) {
  
  state <- init_state
  parms$climate_forcing <- make_climate_forcing(parms)
  
  for (i in seq_len(max_iter)) {
    
    if (verbose) cat("\n--- Spin-up iteration", i, "---\n")
    
    times <- seq(0, 365 * n_years, by = by)
    out <- deSolve::ode(y = state, times = times, func = model_fn, parms = parms)
    
    # check stability
    stab <- check_stability(out)
    if (verbose) print(stab); plot_ode_output(out)
    
    # here assuming check_stability returns numeric drifts per pool
    max_drift <- max(abs(stab$rel_drift), na.rm = TRUE)
    
    if (verbose) cat("Max drift:", max_drift, "\n")
    
    if (max_drift < tol) {
      if (verbose) cat("System stabilized.\n")
      return(list(out = out,
                  final_state = final_state(out),
                  converged = TRUE,
                  iterations = i))
    }
    
    # update state for next loop
    state <- final_state(out)
  }
  
  if (verbose) cat("Reached max iterations without full stability.\n")
  
  return(list(out = out,
              final_state = final_state(out),
              converged = FALSE,
              iterations = max_iter))
}

spinup_until_stable(init_state, parms)

parms$climate_forcing <- make_climate_forcing(parms)
n_years <- 100
times <- seq(0, 365 * n_years, by = 1)
out <- deSolve::ode(y = init_state_spin, times = times, func = wrapped_model, parms = parms)

# ---- inspect ----
plot_ode_output(out)                       # trajectories (incl. mass_balance_check)
print(check_stability(out))                # per-pool drift, last year vs previous
y_spun <- final_state(out)                 # spun-up state for further runs
print(y_spun)

# Not stable yet? Just increase n_years and re-run the SPIN-UP block,
# optionally starting from y_spun:
#   out <- spinup_run(y_spun, model_fn, parms, n_years = 200, by = 1)
