compare_vectors <- function(treatment, baseline) {
  # All names
  all_names <- union(names(treatment), names(baseline))
  
  # Align vectors
  treatmenta <- treatment[all_names]
  baselinea <- baseline[all_names]
  
  # Identify common elements
  common <- intersect(names(treatment), names(baseline))
  
  # Calculations
  diff <- treatmenta - baselinea
  pct <- (diff / treatmenta) * 100
  
  # Set non-common to NA
  diff[!all_names %in% common] <- NA
  pct[!all_names %in% common] <- NA
  
  # Return result
  data.frame(
    name = all_names,
    treatment = treatmenta,
    baseline = baselinea,
    difference = diff,
    percent_change = pct,
    row.names = NULL
  )
}