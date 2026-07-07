# ============================================================
# dynamic_spinup.R - seasonal spin-up + save/restore + follow-up experiments
# ------------------------------------------------------------
# WORKFLOW
#   1. Spin up to the constant-forcing equilibrium (fast; spinup_equilibrium).
#   2. From there, run the SEASONAL dynamic spin-up to the annual limit cycle
#      (slow; dynamic_spinup). Save the result to Data/ with save_spinup().
#   3. Later, load_spinup() the saved state and run a short (~100 yr) follow-up
#      where you ADD animals to a baseline or REMOVE animals from a treatment
#      (followup_add_animals / followup_remove_animals + run_followup).
#
# Requires setup.R, spinup.R, climate_forcing.R sourced.
# ============================================================

# ------------------------------------------------------------
# dynamic_spinup(): SEASONAL integration to the limit cycle, started from an
# equilibrium state (or the object's stored equilibrium, or its raw init).
# This is the slow step -- save the result and reuse it.
#   obj      a setup object (setup_scenario / fitted treatment)
#   from     starting state (default: obj$init_state_spin, else working_state)
#   n_years  length per attempt; lengthened automatically if still drifting
# Returns list(out, final_state, converged, iterations).
# ------------------------------------------------------------
dynamic_spinup <- function(obj, from = NULL, n_years = 600, by = 30,
                           max_iter = 6, tol = 1e-3, abs_floor = 1e-3, verbose = TRUE) {
  start <- if (!is.null(from)) from
           else if (!is.null(obj$init_state_spin)) obj$init_state_spin
           else obj$working_state
  start <- start[names(obj$working_state)]            # align to active pools
  spinup_until_stable(start, obj$parms, obj$wrapped_model,
                      n_years = n_years, by = by, max_iter = max_iter,
                      tol = tol, abs_floor = abs_floor, verbose = verbose)
}

# ------------------------------------------------------------
# save_spinup() / load_spinup(): persist a spun-up state to Data/.
# Stores the state vector PLUS metadata (model, scenario, arm, the parameter
# list, and the active-pool names) so a follow-up can be rebuilt exactly.
# ------------------------------------------------------------
save_spinup <- function(obj, state, scenario, arm,
                        dir = "Data/spinup", tag = NULL) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  tag  <- if (is.null(tag)) sprintf("%s_%s_%s", obj$model, scenario, arm) else tag
  file <- file.path(dir, paste0(tag, ".rds"))
  saveRDS(list(state = state, model = obj$model, scenario = scenario, arm = arm,
               parms = obj$parms, active = names(obj$working_state),
               saved = Sys.time()), file)
  message("saved spun-up state -> ", file)
  invisible(file)
}

load_spinup <- function(file) readRDS(file)

# ------------------------------------------------------------
# run_followup(): shorter SEASONAL simulation from a given starting state,
# using a target setup object's model/parms/wrapper. Returns the deSolve out.
# ------------------------------------------------------------
run_followup <- function(start_state, target, n_years = 100, by = 30,
                         verbose = TRUE) {
  parms <- target$parms
  parms$climate_forcing <- make_climate_forcing(parms)
  y0    <- start_state[names(target$working_state)]
  if (any(is.na(y0)))
    stop("start_state is missing pools needed by the target setup: ",
         paste(names(target$working_state)[is.na(y0)], collapse = ", "))
  if (verbose) message(sprintf("follow-up: %d yr seasonal run on %d pools (%s)",
                               n_years, length(y0), target$model))
  deSolve::ode(y = y0, times = seq(0, 365 * n_years, by = by),
               func = target$wrapped_model, parms = parms)
}

# ------------------------------------------------------------
# followup_add_animals(): take a saved BASELINE (no-animal) spun-up state and
# introduce animals. The target is the treatment setup (animals active, with
# the calibrated parameters); shared plant/soil pools start from the baseline
# limit cycle, the animal pools are seeded with `seed` (default: their input
# values from the treatment setup).
# ------------------------------------------------------------
followup_add_animals <- function(baseline_saved, treatment_setup,
                                 seed = NULL, seed_mod = 0.01, n_years = 100, by = 30, verbose = TRUE) {
  ws     <- treatment_setup$working_state
  shared <- intersect(names(ws), names(baseline_saved$state))
  ws[shared] <- baseline_saved$state[shared]               # carry over spun-up pools
  animals <- intersect(c("Earthworm","Detritivore","DetPredator","RootHerb"), names(ws))
  if (is.null(seed)) seed <- treatment_setup$working_state[animals]*seed_mod  # input biomass
  for (a in animals) if (a %in% names(seed)) ws[a] <- seed[a]
  if (verbose) message("adding animals: ", paste(animals, collapse = ", "),
                       " (seed = ", paste(signif(ws[animals],3), collapse=", "), ")")
  list(start = ws,
       out = run_followup(ws, treatment_setup, n_years = n_years, by = by, verbose = verbose))
}

# ------------------------------------------------------------
# followup_remove_animals(): take a saved TREATMENT (with-animal) spun-up
# state and remove the animals. The target is the baseline setup (animals
# inactive); shared plant/soil pools start from the treatment limit cycle, and
# the animal pools are simply dropped (instantaneous removal).
# ------------------------------------------------------------
followup_remove_animals <- function(treatment_saved, baseline_setup,
                                     n_years = 100, by = 30, verbose = TRUE) {
  ws     <- baseline_setup$working_state          # no animal pools here
  shared <- intersect(names(ws), names(treatment_saved$state))
  ws[shared] <- treatment_saved$state[shared]
  if (verbose) message("removing animals; carrying over ", length(shared), " shared pools")
  list(start = ws,
       out = run_followup(ws, baseline_setup, n_years = n_years, by = by, verbose = verbose))
}

# ------------------------------------------------------------
# followup_continue_baseline(): CONTROL run for followup_add_animals() -- the
# saved baseline (no-animal) limit cycle, continued forward on the SAME
# baseline setup (animals stay off) for the same n_years/by as the add-animals
# run. This gives a time-matched "what if animals were never added" trajectory
# to plot against the add-animals result (see plot_followup_comparison()).
# ------------------------------------------------------------
followup_continue_baseline <- function(baseline_saved, baseline_setup,
                                       n_years = 100, by = 30, verbose = TRUE) {
  ws     <- baseline_setup$working_state
  shared <- intersect(names(ws), names(baseline_saved$state))
  ws[shared] <- baseline_saved$state[shared]
  if (verbose) message("continuing baseline (no animals) for ", n_years, " yr")
  list(start = ws,
       out = run_followup(ws, baseline_setup, n_years = n_years, by = by, verbose = verbose))
}

# ------------------------------------------------------------
# save_followup() / load_followup(): persist a follow-up run (add / remove /
# continue_baseline) keyed by model x scenario x kind, so the two loops in
# Scripts/followup_analysis.R (add+control, remove) can run independently of
# each other and of the plotting step.
# ------------------------------------------------------------
save_followup <- function(model, scenario, kind, result,
                          dir = "Data/followup", tag = NULL) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  tag  <- if (is.null(tag)) sprintf("%s_%s_%s", model, scenario, kind) else tag
  file <- file.path(dir, paste0(tag, ".rds"))
  saveRDS(list(start = result$start, out = result$out, model = model,
               scenario = scenario, kind = kind, saved = Sys.time()), file)
  message("saved follow-up (", kind, ") -> ", file)
  invisible(file)
}

load_followup <- function(model, scenario, kind, dir = "Data/followup") {
  file <- file.path(dir, sprintf("%s_%s_%s.rds", model, scenario, kind))
  if (!file.exists(file)) stop("Follow-up file not found: ", file)
  readRDS(file)
}
