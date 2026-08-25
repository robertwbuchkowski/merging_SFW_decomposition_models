# ============================================================
# SCENARIO PARAMETER UNCERTAINTY
# ------------------------------------------------------------
# Reads the per-scenario parameter uncertainty (Value, SD, Min, Max, Category)
# from the "scenarios" sheet of Data/scenarios.xlsx and turns it into a usable
# +/- range for each scenario x parameter, following this precedence:
#
#   1. SD present (and > 0)     -> range = Value +/- SD          (source "SD")
#   2. else Min AND Max present -> range = [Min, Max]            (source "MinMax")
#   3. else                     -> CV = 2 assumed: SD = 2*|Value|,
#                                  range = Value +/- 2*|Value|,  (source "CV2")
#                                  with a WARNING (very wide, placeholder only).
#
# USAGE
#   source("R/scenario_uncertainty.R")
#   u <- read_param_uncertainty()                 # tidy data frame
#   u <- read_param_uncertainty(category = "Animal")   # animal params only
# Columns: scenario, parameter, units, category, value, sd, min, max,
#          lo, hi, unc_source (SD / MinMax / CV2), bounded (TRUE if the range
#          was truncated to a physical bound, e.g. a_*/p_* to [0,1]).
# ============================================================

read_param_uncertainty <- function(path = "Data/scenarios.xlsx",
                                    sheet = "scenarios",
                                    category = NULL,
                                    cv_fallback = 2,
                                    warn_cv = TRUE) {
  if (!requireNamespace("readxl", quietly = TRUE))
    stop("read_param_uncertainty() needs the 'readxl' package.")

  d <- as.data.frame(readxl::read_excel(path, sheet = sheet),
                     check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("Parameter", "Scenario", "Value")
  if (!all(need %in% names(d)))
    stop("Sheet '", sheet, "' needs columns: ", paste(need, collapse = ", "))

  num <- function(x) suppressWarnings(as.numeric(as.character(x)))
  get <- function(nm) if (nm %in% names(d)) d[[nm]] else rep(NA, nrow(d))

  out <- data.frame(
    scenario  = trimws(as.character(d$Scenario)),
    parameter = trimws(as.character(d$Parameter)),
    units     = as.character(get("Units")),
    category  = as.character(get("Category")),
    value     = num(d$Value),
    sd        = num(get("SD")),
    min       = num(get("Min")),
    max       = num(get("Max")),
    stringsAsFactors = FALSE)

  out <- out[nzchar(out$parameter) & !is.na(out$value), , drop = FALSE]
  if (!is.null(category))
    out <- out[!is.na(out$category) & out$category %in% category, , drop = FALSE]

  # resolve a lo/hi range per row via the precedence above
  has_sd <- is.finite(out$sd) & out$sd > 0
  has_mm <- is.finite(out$min) & is.finite(out$max) & (out$max > out$min)

  out$unc_source <- ifelse(has_sd, "SD", ifelse(has_mm, "MinMax", "CV2"))
  out$lo <- NA_real_; out$hi <- NA_real_

  out$lo[has_sd] <- out$value[has_sd] - out$sd[has_sd]
  out$hi[has_sd] <- out$value[has_sd] + out$sd[has_sd]

  mm <- has_mm & !has_sd
  out$lo[mm] <- out$min[mm]
  out$hi[mm] <- out$max[mm]

  cv <- out$unc_source == "CV2"
  if (any(cv)) {
    out$sd[cv] <- cv_fallback * abs(out$value[cv])
    out$lo[cv] <- out$value[cv] - cv_fallback * abs(out$value[cv])
    out$hi[cv] <- out$value[cv] + cv_fallback * abs(out$value[cv])
    if (warn_cv)
      warning("No SD or Min/Max for ", sum(cv), " parameter row(s); assumed ",
              "CV = ", cv_fallback, " (range = value +/- ", cv_fallback,
              "*|value|). These are placeholders: ",
              paste(unique(paste0(out$scenario[cv], ":", out$parameter[cv])),
                    collapse = ", "))
  }

  # ----------------------------------------------------------
  # PHYSICAL BOUNDS. Some parameters are constrained to a fixed range and their
  # uncertainty range must not exceed it:
  #   a_* (assimilation eff.), p_* (production eff.)  in [0, 1]
  #   LigFrac, a_root_herb, root_to_organic,
  #     prop_feaces_earthworm_LMWC                    in [0, 1]
  #   pct_claysilt                                    in [0, 100]
  # Any other strictly-positive rate/pool parameter gets a lower floor > 0
  # (negative rate constants are non-physical), while genuinely signed params
  # (e.g. k_b_slope_pint < 0) are left alone.
  # ----------------------------------------------------------
  # exact-name bounds
  bounds_named <- list(
    LigFrac                    = c(0, 1),
    a_root_herb                = c(0, 1),
    root_to_organic            = c(0, 1),
    prop_feaces_earthworm_LMWC = c(0, 1),
    p_a                        = c(0, 1),
    p_b                        = c(0, 1),
    pH                         = c(1, 14),
    MAtheta                    = c(0, 1),
    pct_claysilt               = c(0, 100))
  # prefix-family bounds: assimilation (a_*) and production (p_*) efficiencies.
  # The soil partition/scaling coefficients p_a, p_b (fractions, [0,1]) and p_c
  # (an unbounded scaling coefficient) are handled by name below, not here.
  in_unit_family <- grepl("^(a_|p_)", out$parameter) &
                    !out$parameter %in% c("p_a", "p_b", "p_c")

  lo_bound <- rep(-Inf, nrow(out)); hi_bound <- rep(Inf, nrow(out))
  lo_bound[in_unit_family] <- 0;  hi_bound[in_unit_family] <- 1
  for (nm in names(bounds_named)) {
    hit <- out$parameter == nm
    lo_bound[hit] <- bounds_named[[nm]][1]
    hi_bound[hit] <- bounds_named[[nm]][2]
  }

  # clamp the resolved range to the physical bounds
  out$lo <- pmax(out$lo, lo_bound)
  out$hi <- pmin(out$hi, hi_bound)
  # flag rows whose reported range was truncated by a bound (for transparency)
  out$bounded <- (is.finite(lo_bound) & lo_bound > -Inf) |
                 (is.finite(hi_bound) & hi_bound <  Inf)

  # remaining strictly-positive rate/pool params: keep a small positive floor
  # when a bound has not already been applied.
  positive_only <- grepl("^(k_|K_|alpha_|d_|c_|NPP|BD|depth|rho_|CUE|theta_opt)",
                         out$parameter) & is.finite(out$value) & out$value > 0
  needs_floor <- positive_only & !out$bounded & is.finite(out$lo) & out$lo <= 0
  out$lo[needs_floor] <- pmax(out$value[needs_floor] * 1e-3, .Machine$double.eps)

  # guard: never let clamping invert the interval
  bad_int <- is.finite(out$lo) & is.finite(out$hi) & out$lo > out$hi
  out$lo[bad_int] <- pmin(out$value[bad_int], out$hi[bad_int])

  rownames(out) <- NULL
  out
}

# convenience: named list value/lo/hi for one scenario x parameter, or NULL
param_range <- function(unc, scenario, parameter) {
  r <- unc[unc$scenario == scenario & unc$parameter == parameter, , drop = FALSE]
  if (!nrow(r)) return(NULL)
  list(value = r$value[1], lo = r$lo[1], hi = r$hi[1], source = r$unc_source[1])
}
