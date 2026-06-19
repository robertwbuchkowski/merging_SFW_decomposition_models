# ============================================================
# setup.R - switch models, scenarios, and site parameters easily
# ------------------------------------------------------------
# Functions:
#   make_model_wrapper()   integrate only the ACTIVE (non-zero) pools
#   setup_model()          pick a model + turn groups off + override params
#   read_scenarios()       read a scenario CSV (animals on/off + site params)
#   setup_scenario()       build ONE arm of a scenario (treatment or baseline)
#   setup_scenario_pair()  build treatment + matched no-animal baseline
#   compare_to_baseline()  treatment - baseline on shared pools
#
# Scenario CSV layout (rows = parameters, columns = scenarios):
#   * Flag rows  (Tree, Herb, earthworm, RootHerb, Detritivore, DetPredator)
#     use 1/0 for on/off. There is a single Detritivore pool; different
#     detritivore *scenarios* (e.g. isopods vs mites) turn that one pool on
#     and supply different detritivore PARAMETERS via override rows below.
#   * All other rows are PARAMETER overrides (site/climate AND per-scenario
#     animal parameters), e.g. MAT, T_amp, MAtheta, theta_amp, pct_claysilt,
#     fCLAY, LigFrac, GPPmax_*, pH, BD, and detritivore parameters such as
#     d_detritivores / c_detritivores (these let isopod vs mite scenarios use
#     the one Detritivore pool with different values). Names are matched to the
#     parameter list case-insensitively, so minor case differences still
#     resolve (e.g. 'fClay' -> 'fCLAY'). Parameters used by only one model
#     (e.g. Millennial's *_det_k_frag_* or pH/BD) are harmless extras elsewhere.
#
# Each scenario is compared to its OWN baseline: same plants + same climate,
# but ALL animals off.
#
# Assumes working dir = project root and that R/climate_forcing.R, R/spinup.R,
# R/plot_ode_output.R are sourced in your run script.
# ============================================================

# ---- model registry: one entry per model ----
model_table <- list(
  millennial = list(
    src    = c("R/millennial_model.R", "R/derive_millennial_parms.R", "R/init_millennial_state.R"),
    fn = "millennial_model_wplant", yaml = "config/millennial.yml",
    derive = "derive_millennial_parms", init = "init_millennial_state"),
  century = list(
    src    = c("R/century_model.R", "R/derive_century_parms.R", "R/init_century_state.R"),
    fn = "century_model", yaml = "config/century.yml",
    derive = "derive_century_parms", init = "init_century_state"),
  MIMICS = list(
    src    = c("R/MIMICS_model.R", "R/derive_MIMICS_parms.R", "R/init_MIMICS_state.R"),
    fn = "MIMICS_model", yaml = "config/MIMICS.yml",
    derive = "derive_MIMICS_parms", init = "init_MIMICS_state")
)

# which state pools each flag controls
flag_pools <- list(
  Tree        = c("C_leaf_tree", "C_wood_tree", "C_root_tree", "CWD"),
  Herb        = c("C_leaf_herb", "C_root_herb"),
  earthworm   = "Earthworm",
  RootHerb    = "RootHerb",
  Detritivore = "Detritivore",
  DetPredator = "DetPredator"
)
animal_flags <- c("earthworm", "RootHerb", "Detritivore", "DetPredator")

# ------------------------------------------------------------
# make_model_wrapper(): run the model on the reduced (active) state.
# ------------------------------------------------------------
make_model_wrapper <- function(model_fun, full_names, state_groups) {
  function(time, y, parms) {
    full <- setNames(numeric(length(full_names)), full_names)
    full[names(y)] <- y                        # active states; rest stay 0
    d    <- model_fun(time, full, parms)
    dvec <- d[[1]]; names(dvec) <- full_names
    if (length(d) >= 2)
      list(dvec[names(y)], mass_balance_check = unname(d[[2]]))
    else
      list(dvec[names(y)])
  }
}

# resolve override names to the canonical parameter keys (case-sensitive)
.apply_overrides <- function(parms, overrides) {
  if (!length(overrides)) return(parms)
  
  for (nm in names(overrides)) {
    if (nm %in% names(parms)) {
      key <- nm
    } else {
      warning(sprintf("Parameter '%s' not found in 'parms'; adding as new key.", nm))
      key <- nm
    }
    
    parms[[key]] <- overrides[[nm]]
  }
  
  parms
}

# ------------------------------------------------------------
# setup_model(): pick a model, turn groups off, override parameters.
#   off              character vector of state pools to zero
#   param_overrides  named list/vector of parameter values (applied BEFORE
#                    derive and before climate forcing is built)
# ------------------------------------------------------------
setup_model <- function(model, off = character(0), param_overrides = list(),
                        source_files = TRUE) {

  m <- model_table[[model]]
  if (is.null(m))
    stop("Unknown model '", model, "'. Choose: ",
         paste(names(model_table), collapse = ", "))
  if (source_files) invisible(lapply(m$src, source))

  model_fn  <- match.fun(m$fn)
  derive_fn <- match.fun(m$derive)
  init_fn   <- match.fun(m$init)

  # parameters: common + model yaml + overrides, THEN derive
  parms <- yaml::read_yaml("config/common.yml")
  parms <- modifyList(parms, yaml::read_yaml(m$yaml))
  parms <- .apply_overrides(parms, param_overrides)
  parms <- derive_fn(parms)
  
  # Add in the climate forcing:
  parms$climate_forcing <- make_climate_forcing(parms)
  
  # initial state, then zero requested groups
  init_state <- init_fn()
  bad <- setdiff(off, names(init_state))
  if (length(bad))
    stop("These 'off' names are not state variables of ", model, ": ",
         paste(bad, collapse = ", "))
  init_state[off] <- 0

  full_names    <- names(init_state)
  active        <- names(init_state)[init_state > 0]
  working_state <- init_state[active]

  if (!any(c("C_leaf_herb","C_root_herb","C_leaf_tree","C_wood_tree","C_root_tree") %in% active))
    warning("No plant pools active -> no litter input; soil will decay to ~0.")

  list(model = model, model_fn = model_fn, parms = parms,
       init_state = init_state, working_state = working_state,
       full_names = full_names, active = active,
       wrapped_model = make_model_wrapper(model_fn, full_names, active))
}

# ------------------------------------------------------------
# read_scenarios(): parse the scenario CSV.
# Returns a named list (one per scenario column); each element is
#   list(flags = c(Tree=,Herb=,earthworm=,RootHerb=,Detritivore=,DetPredator=),
#        params = named numeric of site/climate overrides)
# ------------------------------------------------------------
read_scenarios <- function(csv_path = "config/scenarios.csv") {
  raw  <- utils::read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  pname <- raw[[1]]
  scen  <- names(raw)[-1]
  flag_names <- names(flag_pools)
  is_flag <- tolower(pname) %in% tolower(flag_names)

  dup <- unique(pname[duplicated(pname)])
  if (length(dup))
    warning("Duplicate parameter rows in CSV (last value used): ",
            paste(dup, collapse = ", "))

  out <- lapply(scen, function(s) {
    v  <- suppressWarnings(as.numeric(raw[[s]]))
    # flags, in canonical order; 0 if a flag row is absent
    fi    <- match(tolower(flag_names), tolower(pname))
    flags <- setNames(ifelse(is.na(fi), 0L, as.integer(v[fi])), flag_names)
    # parameter overrides = every non-flag row (last value wins if repeated)
    pn   <- pname[!is_flag]; pv <- v[!is_flag]
    keep <- !duplicated(pn, fromLast = TRUE)
    params <- setNames(pv[keep], pn[keep])
    list(flags = flags, params = params)
  })
  setNames(out, scen)
}

# ------------------------------------------------------------
# setup_scenario(): build one arm of a scenario.
#   animals = TRUE  -> treatment (animals per the CSV flags)
#   animals = FALSE -> baseline  (all animals off; same plants + climate)
# ------------------------------------------------------------
setup_scenario <- function(model, scenarios, scenario, animals = TRUE,
                           source_files = TRUE) {
  sc <- scenarios[[scenario]]
  if (is.null(sc))
    stop("Unknown scenario '", scenario, "'. Available: ",
         paste(names(scenarios), collapse = ", "))

  flags <- sc$flags
  if (!animals) flags[animal_flags] <- 0L          # baseline: remove all fauna

  if (animals && flags["DetPredator"] == 1L && flags["Detritivore"] == 0L)
    warning("Scenario '", scenario, "': DetPredator on but Detritivore off (no prey).")

  off <- unlist(flag_pools[names(flags)[flags == 0L]], use.names = FALSE)

  s <- setup_model(model, off = off, param_overrides = sc$params,
                   source_files = source_files)
  s$scenario <- scenario
  s$arm      <- if (animals) "treatment" else "baseline"
  s
}

# ------------------------------------------------------------
# setup_scenario_pair(): treatment + matched no-animal baseline.
# Same plants and same climate/site parameters; baseline has no animals.
# ------------------------------------------------------------
setup_scenario_pair <- function(model, scenarios, scenario, source_files = TRUE) {
  list(
    treatment = setup_scenario(model, scenarios, scenario, animals = TRUE,
                               source_files = source_files),
    baseline  = setup_scenario(model, scenarios, scenario, animals = FALSE,
                               source_files = FALSE)
  )
}

# ------------------------------------------------------------
# compare_to_baseline(): treatment - baseline on the pools they share
# (plants + soil; animal pools exist only in the treatment).
# ------------------------------------------------------------
compare_to_baseline <- function(out_treatment, out_baseline) {
  ft <- final_state(out_treatment); fb <- final_state(out_baseline)
  shared <- intersect(names(ft), names(fb))
  data.frame(
    pool      = shared,
    treatment = as.numeric(ft[shared]),
    baseline  = as.numeric(fb[shared]),
    diff      = as.numeric(ft[shared] - fb[shared]),
    pct       = 100 * as.numeric((ft[shared] - fb[shared]) / pmax(abs(fb[shared]), 1e-8)),
    row.names = NULL
  )
}

# ============================================================
# USAGE
# ============================================================
# library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml)
# source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
# source("R/setup.R")
#
# scen <- read_scenarios("config/scenarios.csv")
# names(scen)            # "RootHerbivore" "Earthworm" "Isopod" "Mite" "MitePredator"
#
# # one scenario, both arms (same model + same climate; baseline has no animals)
# pair  <- setup_scenario_pair("MIMICS", scen, "Mite")
# out_t <- spinup_run(pair$treatment$working_state, pair$treatment$wrapped_model,
#                     pair$treatment$parms, n_years = 200, by = 1)
# out_b <- spinup_run(pair$baseline$working_state,  pair$baseline$wrapped_model,
#                     pair$baseline$parms,  n_years = 200, by = 1)
# compare_to_baseline(out_t, out_b)
#
# # loop over every scenario for a model:
# for (nm in names(scen)) {
#   p  <- setup_scenario_pair("century", scen, nm)
#   ot <- spinup_run(p$treatment$working_state, p$treatment$wrapped_model, p$treatment$parms, n_years = 200)
#   ob <- spinup_run(p$baseline$working_state,  p$baseline$wrapped_model,  p$baseline$parms,  n_years = 200)
#   cat("\n==", nm, "==\n"); print(compare_to_baseline(ot, ob))
# }
