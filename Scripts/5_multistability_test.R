# ============================================================
# TEST FOR MULTIPLE STABLE STATES (alternative equilibria) - constant forcing
# ------------------------------------------------------------
# Microbial-explicit soil models can have more than one STABLE equilibrium
# (e.g. an active state vs. a microbial-collapse state). This script probes each
# scenario for that by running the constant-forcing steady-state solver from
# MANY diverse initial conditions and checking whether they settle on distinct
# fixed points.
#
# WHY THIS WORKS: spinup_equilibrium() uses rootSolve::runsteady, which finds a
# steady state by FORWARD integration -- so it only lands on STABLE equilibria
# (unstable ones repel). Sampling initial conditions across the state space and
# clustering the converged results therefore counts the stable states and their
# basins of attraction. This does NOT prove uniqueness (a basin we never seed
# is missed), but repeated single-state convergence from wide, extreme starts is
# strong evidence of a single stable state; >=2 distinct states is proof of
# multistability.
#
# Output: Results/multistability_summary.csv (one row per scenario),
#         Results/multistability_states.csv  (every distinct state found),
#         Results/figures/multistability_states.png
# Run from the project root.
# ============================================================
library(pacman); p_load(deSolve, rootSolve, tidyverse, yaml, readxl)
source("R/climate_forcing.R"); source("R/spinup.R"); source("R/plot_ode_output.R")
source("R/setup.R");           source("R/compare_functions.R")
source("R/fit_animals.R");     source("R/dynamic_spinup.R")

set.seed(1)                                   # reproducible initial-condition draws
model     <- "millennial"
scen      <- read_scenarios("Data/scenarios.xlsx")
scen$MitePredator <- NULL
scenarios <- names(scen)

fig_dir <- "Results/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
res_dir <- "Results";         dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# ---- test settings -------------------------------------------------------
n_starts   <- 60        # random initial conditions per scenario (+ fixed corners)
span_dec    <- 2        # random scaling spans +/- this many orders of magnitude
rel_tol     <- 1e-2     # two equilibria are "the same" if max relative diff < this
abs_floor   <- 1e-6     # ignore pools below this when comparing (numerical dust)
use_treatment <- FALSE  # test the BASELINE (no-animal) equilibrium first

# ------------------------------------------------------------
# make_starts(): a matrix of initial states (rows = pools, cols = starts).
# Columns are: the default; an all-low and an all-high corner; a microbial-
# collapse corner (tiny MIC/B); a microbial-bloom corner; then random draws
# where every pool is multiplied by an independent log-uniform factor. This
# spreads seeds across basins, including the ones microbial models flip between.
# ------------------------------------------------------------
make_starts <- function(y0, n_random, span = 2) {
  pools <- names(y0)
  cols  <- list(default = y0)

  cols$all_low  <- y0 * 10^(-span)
  cols$all_high <- y0 * 10^(+span)
  micdead <- y0; micdead[intersect(c("MIC","B"), pools)] <- abs_floor * 10
  cols$microbial_collapse <- micdead
  micbloom <- y0; micbloom[intersect(c("MIC","B"), pools)] <-
    y0[intersect(c("MIC","B"), pools)] * 10^span
  cols$microbial_bloom <- micbloom

  for (i in seq_len(n_random)) {
    fac <- 10^runif(length(y0), -span, span)
    cols[[paste0("rand", i)]] <- y0 * fac
  }
  M <- do.call(cbind, cols)
  rownames(M) <- pools
  M
}

# ------------------------------------------------------------
# distinct_states(): cluster converged equilibria into unique states by the
# max RELATIVE pool difference (pools below abs_floor ignored). Greedy: each
# state joins the first cluster it matches, else starts a new one.
# ------------------------------------------------------------
same_state <- function(a, b) {
  keep <- (abs(a) > abs_floor) | (abs(b) > abs_floor)
  if (!any(keep)) return(TRUE)
  denom <- pmax(abs(a[keep]), abs(b[keep]), abs_floor)
  max(abs(a[keep] - b[keep]) / denom) < rel_tol
}
cluster_states <- function(states) {          # list of named numeric vectors
  reps <- list(); members <- integer(0)
  for (i in seq_along(states)) {
    hit <- NA_integer_
    for (k in seq_along(reps))
      if (same_state(states[[i]], reps[[k]])) { hit <- k; break }
    if (is.na(hit)) { reps[[length(reps)+1]] <- states[[i]]; members <- c(members, 1L) }
    else members[hit] <- members[hit] + 1L
  }
  list(reps = reps, counts = members)
}

# ------------------------------------------------------------
# run one scenario: seed many starts, keep converged states, cluster them.
# ------------------------------------------------------------
soil_pools_of <- function(nm) setdiff(nm, c("Earthworm","Detritivore","DetPredator","RootHerb"))

test_scenario <- function(scenario) {
  obj <- setup_scenario(model, scen, scenario, animals = use_treatment)
  y0  <- obj$working_state
  S   <- make_starts(y0, n_starts, span = span_dec)

  eqs <- list(); conv <- logical(ncol(S))
  for (j in seq_len(ncol(S))) {
    start <- S[, j]; names(start) <- rownames(S)
    o <- tryCatch(spinup_equilibrium(obj, warm_start = start, verbose = FALSE),
                  error = function(e) NULL)
    ok <- !is.null(o) && isTRUE(o$spin_info$converged) && all(is.finite(o$init_state_spin))
    conv[j] <- ok
    if (ok) eqs[[length(eqs)+1]] <- o$init_state_spin
  }

  cl <- cluster_states(eqs)
  n_states <- length(cl$reps)

  # per-state total soil + root C, to describe how different the states are
  state_tab <- map_dfr(seq_len(n_states), function(k) {
    s  <- cl$reps[[k]]; sp <- soil_pools_of(names(s))
    hdr <- tibble(scenario = scenario, state = k, n_starts = cl$counts[k],
                  total_soil_C = sum(s[sp]))
    pool_cols <- as_tibble(as.list(round(s, 4)))       # one column per pool
    bind_cols(hdr, pool_cols)
  })

  list(
    summary = tibble(
      scenario         = scenario,
      n_converged      = sum(conv),
      n_tested         = ncol(S),
      n_stable_states  = n_states,
      multistable      = n_states > 1,
      min_total_C      = if (n_states) min(state_tab$total_soil_C) else NA_real_,
      max_total_C      = if (n_states) max(state_tab$total_soil_C) else NA_real_,
      spread_pct       = if (n_states > 1)
        round(100 * (max(state_tab$total_soil_C) - min(state_tab$total_soil_C)) /
                max(abs(state_tab$total_soil_C)), 1) else 0),
    states = state_tab)
}

# ---- run all scenarios ---------------------------------------------------
all_summary <- list(); all_states <- list()
for (sc in scenarios) {
  cat("\n=== multistability test:", sc, "===\n")
  r <- test_scenario(sc)
  print(r$summary)
  all_summary[[sc]] <- r$summary
  all_states[[sc]]  <- r$states
}
summary_tbl <- bind_rows(all_summary)
states_tbl  <- bind_rows(all_states)
write_csv(summary_tbl, file.path(res_dir, "multistability_summary.csv"))
write_csv(states_tbl,  file.path(res_dir, "multistability_states.csv"))

cat("\n================= RESULT =================\n")
print(as.data.frame(summary_tbl), row.names = FALSE)
if (any(summary_tbl$multistable))
  cat("\nALTERNATIVE STABLE STATES found in:",
      paste(summary_tbl$scenario[summary_tbl$multistable], collapse = ", "), "\n")
else
  cat("\nNo multistability detected: every scenario converged to a single stable",
      "state from all", unique(summary_tbl$n_tested), "seeds.\n")

# ------------------------------------------------------------
# FIGURE: total soil + root C of each distinct stable state per scenario, sized
# by how many starts fell into it (basin share). Multiple points in a column =
# multiple stable states.
# ------------------------------------------------------------
p <- ggplot(states_tbl, aes(scenario, total_soil_C, size = n_starts)) +
  geom_point(alpha = 0.7, colour = "#2166ac") +
  scale_size_area(max_size = 10, name = "Starts in basin") +
  labs(title = "Distinct stable states per scenario (constant forcing)",
       subtitle = "Each point = one converged stable state; >1 point in a column = alternative stable states",
       x = NULL, y = expression("Total soil + root C at equilibrium (g C m"^-2*")")) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig_dir, "multistability_states.png"), p, width = 9, height = 6, dpi = 150)
print(p)

cat("\nWrote:\n  ", file.path(res_dir, "multistability_summary.csv"),
    "\n  ", file.path(res_dir, "multistability_states.csv"),
    "\n  ", file.path(fig_dir, "multistability_states.png"), "\n")
