# ============================================================
# setup.R — easy switching between models and scenarios
# ------------------------------------------------------------
# Two functions:
#   make_model_wrapper()  wraps a model so deSolve integrates only the
#                         ACTIVE (non-zero) pools; inactive pools are held
#                         at 0. Passes mass_balance_check through as output.
#   setup_model()         pick a model + a scenario (which groups to turn
#                         off) in ONE call; returns everything ready to run.
#
# Assumes the working directory is the project root and that the shared
# helpers (R/climate_forcing.R, R/spinup.R, R/plot_ode_output.R) are sourced
# in your run script. setup_model() sources the chosen model's own files.
# ============================================================

# ---- model registry: one entry per model (add a model = add a row) ----
model_table <- list(
  millennial = list(
    src    = c("R/millennial_model.R", "R/derive_millennial_parms.R", "R/init_millennial_state.R"),
    fn     = "millennial_model_wplant",
    yaml   = "config/millennial.yml",
    derive = "derive_millennial_parms",
    init   = "init_millennial_state"
  ),
  century = list(
    src    = c("R/century_model.R", "R/derive_century_parms.R", "R/init_century_state.R"),
    fn     = "century_model",
    yaml   = "config/century.yml",
    derive = "derive_century_parms",
    init   = "init_century_state"
  ),
  MIMICS = list(
    src    = c("R/MIMICS_model.R", "R/derive_MIMICS_parms.R", "R/init_MIMICS_state.R"),
    fn     = "MIMICS_model",
    yaml   = "config/MIMICS.yml",
    derive = "derive_MIMICS_parms",
    init   = "init_MIMICS_state"
  )
)

# ------------------------------------------------------------
# make_model_wrapper(): run the model on the reduced (active) state.
#   model_fun    the model function, e.g. MIMICS_model
#   full_names   names of ALL state variables (= names(init_state))
#   state_groups names of the ACTIVE states to integrate
#                (e.g. names(init_state)[init_state > 0])
# The wrapper rebuilds the full state (inactive = 0), calls the model, and
# returns derivatives for the active states only (in the order deSolve
# supplies them). mass_balance_check rides along as an output column.
# ------------------------------------------------------------
make_model_wrapper <- function(model_fun, full_names, state_groups) {
  active <- state_groups
  function(time, y, PARMS) {
    full <- setNames(numeric(length(full_names)), full_names)
    full[names(y)] <- y                       # fill active; rest stay 0

    d    <- model_fun(time, full, PARMS)
    dvec <- d[[1]]
    names(dvec) <- full_names

    if (length(d) >= 2) {
      list(dvec[names(y)], mass_balance_check = unname(d[[2]]))
    } else {
      list(dvec[names(y)])
    }
  }
}

# ------------------------------------------------------------
# setup_model(): pick a model and a scenario in one call.
#   model  one of names(model_table): "millennial", "century", "MIMICS"
#   off    character vector of groups/pools to turn OFF (zeroed). Examples:
#            off = c("C_leaf_tree","C_wood_tree","C_root_tree")  # no trees
#            off = "Earthworm"                                    # no earthworms
#            off = c("Detritivore","DetPredator")                # no det/predators
#   source_files  source the model's R files (TRUE) or assume already sourced
#
# Returns a list with: parms, full init_state, working_state (active only),
# wrapped_model (ready for deSolve), plus model_fn / full_names / active.
# ------------------------------------------------------------
setup_model <- function(model, off = character(0), source_files = TRUE) {

  m <- model_table[[model]]
  if (is.null(m))
    stop("Unknown model '", model, "'. Choose one of: ",
         paste(names(model_table), collapse = ", "))

  if (source_files) invisible(lapply(m$src, source))   # sources into globalenv

  model_fn  <- match.fun(m$fn)
  derive_fn <- match.fun(m$derive)
  init_fn   <- match.fun(m$init)

  # ---- parameters: common + model-specific, then derive ----
  parms <- yaml::read_yaml("config/common.yml")
  parms <- modifyList(parms, yaml::read_yaml(m$yaml))
  parms <- derive_fn(parms)

  # ---- initial state, then turn requested groups off ----
  init_state <- init_fn()
  bad <- setdiff(off, names(init_state))
  if (length(bad))
    stop("These 'off' names are not state variables of ", model, ": ",
         paste(bad, collapse = ", "),
         "\n  Valid names: ", paste(names(init_state), collapse = ", "))
  init_state[off] <- 0

  # ---- reduce to the active (non-zero) states ----
  full_names    <- names(init_state)
  active        <- names(init_state)[init_state > 0]
  working_state <- init_state[active]

  if (!any(c("C_leaf_herb", "C_root_herb", "C_leaf_tree", "C_wood_tree",
             "C_root_tree") %in% active))
    warning("No plant pools are active -> no litter input; soil will decay to ~0.")

  wrapped_model <- make_model_wrapper(model_fn, full_names, active)

  list(
    model         = model,
    model_fn      = model_fn,
    parms         = parms,
    init_state    = init_state,     # full vector (with zeros)
    working_state = working_state,  # active states only -> pass this to ode/spinup
    full_names    = full_names,
    active        = active,
    wrapped_model = wrapped_model   # pass this as `func` to ode/spinup
  )
}
