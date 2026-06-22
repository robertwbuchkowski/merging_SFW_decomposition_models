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
#     fCLAY, LigFrac, NPP_herb, NPP_tree, pH, BD, and detritivore parameters such as
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
.apply_overrides <- function(parms, overrides, addkey = F) {
  if (!length(overrides)) return(parms)
  
  for (nm in names(overrides)) {
    if (nm %in% names(parms)) {
      key <- nm
    } else {
      if(addkey){
        warning(sprintf("Parameter '%s' not found in 'parms'; adding as new key.", nm))
        key <- nm
      }
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
  if (length(bad)){
    warning("These 'off' names are not state variables of ", model, "  and will be dropped: ",
            paste(bad, collapse = ", "))
    
    off = off[off %in% names(init_state)]
  }
    
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
# read_scenarios(): parse a LONG/TIDY scenario CSV.
# Expected columns: Parameter, Scenario, Value [, SD].
# One row per (Parameter, Scenario); the SD column (if present) is ignored.
#
# Returns a named list (one element per scenario) of:
#   list(flags  = c(Tree=, Herb=, earthworm=, RootHerb=, Detritivore=, DetPredator=),
#        params = named numeric of all non-flag parameter overrides)
# Depends on `flag_pools` (defined elsewhere in setup.R) for the canonical
# flag names and order.
# ------------------------------------------------------------
read_scenarios <- function(csv_path = "Data/scenarios.csv") {
  raw <- utils::read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  needed <- c("Parameter", "Scenario", "Value")
  if (!all(needed %in% names(raw)))
    stop("scenarios CSV must be long form with columns: ",
         paste(needed, collapse = ", "),
         "\n  (found: ", paste(names(raw), collapse = ", "), ")")
  
  raw$Parameter <- trimws(as.character(raw$Parameter))
  raw$Scenario  <- trimws(as.character(raw$Scenario))
  raw$Value     <- suppressWarnings(as.numeric(raw$Value))
  
  flag_names <- names(flag_pools)
  scen       <- unique(raw$Scenario)          # first-seen order
  
  out <- lapply(scen, function(s) {
    sub <- raw[raw$Scenario == s, , drop = FALSE]
    
    dup <- unique(sub$Parameter[duplicated(sub$Parameter)])
    if (length(dup))
      warning("Scenario '", s, "': duplicate parameter rows (last value used): ",
              paste(dup, collapse = ", "))
    
    # flags in canonical order; 0 if a flag row is absent for this scenario
    fi    <- match(tolower(flag_names), tolower(sub$Parameter))
    flags <- setNames(ifelse(is.na(fi), 0L, as.integer(sub$Value[fi])), flag_names)
    
    # parameter overrides = every non-flag row (last value wins if repeated)
    psub   <- sub[!(tolower(sub$Parameter) %in% tolower(flag_names)), , drop = FALSE]
    keep   <- !duplicated(psub$Parameter, fromLast = TRUE)
    params <- setNames(psub$Value[keep], psub$Parameter[keep])
    
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


# ------------------------------------------------------------
# run_scenario(): run setup, checks, and steady-state spinup
# for a single scenario name.
#
#   nm      scenario name (must match names(scen))
#   scen    scenario table/list (from read_scenarios)
#   model   model identifier ("millennial", "century", "MIMICS")
#
# Returns:
#   A list with:
#     pair        full baseline/treatment objects after spinup
#     eq_compare  comparison of equilibrium states between
#                 treatment and baseline
# ------------------------------------------------------------
run_scenario <- function(nm, scen, model) {
  message("\n==================  ", nm, "  ==================")
  
  pair <- setup_scenario_pair(model, scen, nm)
  
  # -----------------------------
  # Mass balance check function
  # -----------------------------
  check_mb <- function(obj) {
    mb0 <- obj$wrapped_model(0, obj$working_state, obj$parms)[[2]]
    cat(sprintf("mass_balance_check at t=0: %.3e\n", mb0))
  }
  
  # -----------------------------
  # Spin-up function
  # -----------------------------
  spinup_once <- function(obj) {
    obj$parms$climate_forcing <- make_climate_forcing_equilibrium(obj$parms)
    
    cfss <- stode(
      y     = obj$working_state,
      func  = obj$wrapped_model,
      parms = obj$parms
    )
    
    obj$parms$climate_forcing <- make_climate_forcing(obj$parms)
    
    if (attr(cfss, "steady")) {
      cat("Result is steady\n")
      obj$init_state_spin <- cfss[[1]]
    } else {
      cat("Result did not reach stability. Using input state.\n")
      obj$init_state_spin <- obj$working_state
    }
    
    obj
  }
  
  # -----------------------------
  # Treatment
  # -----------------------------
  cat("Treatment\n")
  check_mb(pair$treatment)
  check_mb(pair$baseline)
  
  pair$treatment <- spinup_once(pair$treatment)
  
  # -----------------------------
  # Baseline
  # -----------------------------
  cat("Baseline\n")
  check_mb(pair$baseline)
  check_mb(pair$baseline)  # keeping original duplicate check
  
  pair$baseline <- spinup_once(pair$baseline)
  
  # -----------------------------
  # Output
  # -----------------------------
  list(
    pair = pair,
    eq_compare = compare_vectors(
      pair$treatment$init_state_spin,
      pair$baseline$init_state_spin
    )
  )
}