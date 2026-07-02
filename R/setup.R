# ============================================================
# setup.R - switch models, scenarios, and site parameters easily
# ------------------------------------------------------------
# Functions:
#   make_model_wrapper()   integrate only the ACTIVE (non-zero) pools
#   setup_model()          pick a model + turn groups off + override params
#   read_scenarios()       read the scenarios workbook (Excel .xlsx, 3 sheets)
#   setup_scenario()       build ONE arm of a scenario (treatment or baseline)
#   setup_scenario_pair()  build treatment + matched no-animal baseline
#
# Scenarios live in Data/scenarios.xlsx (see read_scenarios() for the sheet
# layout): a "scenarios" sheet of per-scenario PARAMETER overrides, a
# "state_variable_include" sheet of on/off FLAGS (Tree, Herb, earthworm,
# RootHerb, Detritivore, DetPredator), and a "state_variables_value" sheet of
# initial pool values. There is a single Detritivore pool; isopod vs mite
# scenarios turn it on and supply different detritivore parameters. Parameter
# and pool names are matched case-insensitively; parameters or pools used by
# only one model are harmless extras elsewhere.
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
                        init_overrides = list(), source_files = TRUE) {

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
  
  # initial state from the model's init_*_state(); then override values supplied
  # from the workbook (matched by pool name; blanks/unknowns ignored), then zero
  # the requested groups.
  init_state <- init_fn()
  if (length(init_overrides)) {
    iv <- init_overrides[names(init_overrides) %in% names(init_state)]
    iv <- iv[is.finite(iv)]
    if (length(iv)) init_state[names(iv)] <- iv
  }
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
# read_scenarios(): read the scenarios workbook (Excel .xlsx) with 3 sheets.
#
#   "scenarios"              long-form PARAMETER overrides:
#                            columns Parameter, Scenario, Value (extra metadata
#                            columns such as SD/Min/Max/Reference are ignored).
#   "state_variable_include" inclusion FLAGS (1/0) per StateVariable x Scenario,
#                            StateVariable being a flag group (Tree, Herb,
#                            earthworm, RootHerb, Detritivore, DetPredator).
#   "state_variables_value"  INITIAL pool values: Model, Scenario, StateVariable,
#                            InitialEq. Model is "All" (plants/animals) or a
#                            model name (its soil pools); Scenario is "All" or a
#                            specific one. A blank InitialEq keeps the
#                            init_*_state.R default for that pool.
#
# Scenario names are matched across sheets ignoring case and punctuation, so
# "Earthworm - Temperate" and "Earthworm-Temperate" refer to the same scenario.
#
# Returns a named list (one per scenario, in "scenarios"-sheet order) of:
#   list(flags  = c(Tree=, Herb=, earthworm=, RootHerb=, Detritivore=, DetPredator=),
#        params = named numeric parameter overrides,
#        init   = named numeric initial values overriding the model defaults)
#
# A legacy long-form .csv path still works (flags + params from the one sheet).
# ------------------------------------------------------------
read_scenarios <- function(path = "Data/scenarios.xlsx",
                           sheet_params  = "scenarios",
                           sheet_include = "state_variable_include",
                           sheet_init    = "state_variables_value") {

  if (!requireNamespace("readxl", quietly = TRUE))
    stop("read_scenarios() needs the 'readxl' package to read ", path,
         "  (install.packages('readxl')).")

  rd  <- function(sh) as.data.frame(readxl::read_excel(path, sheet = sh),
                                    check.names = FALSE, stringsAsFactors = FALSE)
  nrm <- function(x) tolower(gsub("[^[:alnum:]]+", "", trimws(as.character(x))))

  # ---- parameters ----
  pr <- rd(sheet_params)
  need <- c("Parameter", "Scenario", "Value")
  if (!all(need %in% names(pr)))
    stop("Sheet '", sheet_params, "' needs columns: ", paste(need, collapse = ", "))
  pr$Parameter <- trimws(as.character(pr$Parameter))
  pr$Scenario  <- trimws(as.character(pr$Scenario))
  pr$Value     <- suppressWarnings(as.numeric(as.character(pr$Value)))
  pr <- pr[nzchar(pr$Parameter) & !is.na(pr$Parameter) &
           nzchar(pr$Scenario)  & !is.na(pr$Scenario), , drop = FALSE]

  # ---- inclusion flags ----
  inc <- rd(sheet_include)
  need <- c("StateVariable", "Scenario", "Inclusion")
  if (!all(need %in% names(inc)))
    stop("Sheet '", sheet_include, "' needs columns: ", paste(need, collapse = ", "))
  inc$sv  <- trimws(as.character(inc$StateVariable))
  inc$key <- nrm(inc$Scenario)

  # ---- initial values ----
  iv <- rd(sheet_init)
  need <- c("Model", "Scenario", "StateVariable", "InitialEq")
  if (!all(need %in% names(iv)))
    stop("Sheet '", sheet_init, "' needs columns: ", paste(need, collapse = ", "))
  iv$sv        <- trimws(as.character(iv$StateVariable))
  iv$key       <- nrm(iv$Scenario)
  iv$InitialEq <- suppressWarnings(as.numeric(as.character(iv$InitialEq)))

  flag_names <- names(flag_pools)
  scen       <- unique(pr$Scenario)                 # canonical names + order

  out <- lapply(scen, function(s) {
    skey <- nrm(s)
    sub  <- pr[pr$Scenario == s, , drop = FALSE]

    dup <- unique(sub$Parameter[duplicated(sub$Parameter)])
    if (length(dup))
      warning("Scenario '", s, "': duplicate parameter rows (last value used): ",
              paste(dup, collapse = ", "))
    keep   <- !duplicated(sub$Parameter, fromLast = TRUE)
    params <- setNames(sub$Value[keep], sub$Parameter[keep])
    params <- params[!is.na(params)]

    # flags (absent -> 0)
    fsub  <- inc[inc$key == skey, , drop = FALSE]
    fi    <- match(tolower(flag_names), tolower(fsub$sv))
    flags <- setNames(ifelse(is.na(fi), 0L, as.integer(fsub$Inclusion[fi])), flag_names)

    # initial values: "All" scenario first, this scenario overrides; drop blanks
    isub  <- iv[iv$key %in% c(nrm("All"), skey), , drop = FALSE]
    isub  <- isub[order(isub$key == skey), , drop = FALSE]   # specific rows last
    isub  <- isub[!is.na(isub$InitialEq), , drop = FALSE]
    init  <- setNames(isub$InitialEq, isub$sv)
    init  <- init[!duplicated(names(init), fromLast = TRUE)]

    list(flags = flags, params = params, init = init)
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
                   init_overrides = if (is.null(sc$init)) list() else sc$init,
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
