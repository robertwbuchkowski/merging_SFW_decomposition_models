#------------------------------------------------------------------------#
# SCENARIO PARAMETER UNCERTAINTY ----
#------------------------------------------------------------------------#
# Reads the per-scenario parameter values and uncertainty (Value, SD, Min, Max,
# Category) from the "scenarios" sheet of Data/scenarios.xlsx and gives every
# row a standard deviation:
#
#   1. SD reported (> 0)  -> that SD                                 ("SD")
#   2. otherwise          -> CV = cv_fallback (default from cv_default, else 1) ("CV-default")
#   3. ... but if Min AND Max are reported and imply a smaller spread,
#      the CV is reduced so that Min..Max spans +/- 2 SD:
#      SD = (Max - Min) / 4                                         ("CV-MinMax")
#   lo / hi = Value +/- 2 SD, cut to the reported Min / Max when present.
#
# Physical bounds (proportions in [0, 1], sign, pH ...) are applied later by
# param_bounds() in R/sensitivity_helpers.R, for spreadsheet and model
# parameters alike.
#
# USAGE
#   source("R/scenario_uncertainty.R")
#   u <- read_param_uncertainty()                      # tidy data frame
#   u <- read_param_uncertainty(category = "Animal")   # animal params only
# Columns: scenario, parameter, units, category, value, sd_reported, sd (used),
#          cv (used), min, max, lo, hi, unc_source (SD / CV-MinMax / CV-default).

read_param_uncertainty <- function(path = "Data/scenarios.xlsx",
                                    sheet = "scenarios",
                                    category = NULL,
                                    cv_fallback = if (exists("cv_default")) cv_default else 1,
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
    scenario  = if (exists("rename_scenarios", mode = "function"))
                  rename_scenarios(d$Scenario) else trimws(as.character(d$Scenario)),
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

  # ---- standard deviation for every row (SD rule) ----
  #   SD reported (> 0)        -> sd = SD                          source "SD"
  #   otherwise                -> sd = cv_fallback * |value|       source "CV-default"
  #     ... unless Min/Max is reported and implies a smaller spread:
  #         sd = (Max - Min) / 4, i.e. Min..Max spans +/- 2 SD     source "CV-MinMax"
  # lo / hi = value +/- 2 sd, cut to the reported Min / Max when present.
  # Physical bounds (proportions in [0, 1], sign constraints, ...) are NOT
  # applied here; param_bounds() in R/sensitivity_helpers.R does that, so the
  # same rules apply to parameters that are not in the spreadsheet.
  out$sd_reported <- ifelse(is.finite(out$sd) & out$sd > 0, out$sd, NA_real_)
  has_mm <- is.finite(out$min) & is.finite(out$max) & (out$max > out$min)
  sd_cv  <- cv_fallback * abs(out$value)
  sd_mm  <- ifelse(has_mm, (out$max - out$min) / 4, NA_real_)
  use_mm <- is.na(out$sd_reported) & has_mm & sd_mm < sd_cv

  out$unc_source <- ifelse(!is.na(out$sd_reported), "SD",
                           ifelse(use_mm, "CV-MinMax", "CV-default"))
  out$sd <- ifelse(out$unc_source == "SD", out$sd_reported,
                   ifelse(use_mm, sd_mm, sd_cv))
  out$cv <- ifelse(out$value != 0, out$sd / abs(out$value), NA_real_)
  out$lo <- out$value - 2 * out$sd
  out$hi <- out$value + 2 * out$sd
  out$lo[has_mm] <- pmax(out$lo[has_mm], out$min[has_mm])
  out$hi[has_mm] <- pmin(out$hi[has_mm], out$max[has_mm])

  cv <- out$unc_source == "CV-default"
  if (warn_cv && any(cv))
    warning("No SD and no narrowing Min/Max for ", sum(cv), " parameter row(s); assumed ",
            "CV = ", cv_fallback, ": ",
            paste(unique(paste0(out$scenario[cv], ":", out$parameter[cv])), collapse = ", "))

  rownames(out) <- NULL
  out
}

# convenience: named list value/lo/hi for one scenario x parameter, or NULL
param_range <- function(unc, scenario, parameter) {
  r <- unc[unc$scenario == scenario & unc$parameter == parameter, , drop = FALSE]
  if (!nrow(r)) return(NULL)
  list(value = r$value[1], lo = r$lo[1], hi = r$hi[1], source = r$unc_source[1])
}
