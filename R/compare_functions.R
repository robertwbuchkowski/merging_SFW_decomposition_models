compare_vectors <- function(v1, v2) {
  # All names
  all_names <- union(names(v1), names(v2))
  
  # Align vectors
  v1a <- v1[all_names]
  v2a <- v2[all_names]
  
  # Identify common elements
  common <- intersect(names(v1), names(v2))
  
  # Calculations
  diff <- v2a - v1a
  pct <- (diff / v1a) * 100
  
  # Set non-common to NA
  diff[!all_names %in% common] <- NA
  pct[!all_names %in% common] <- NA
  
  # Return result
  data.frame(
    name = all_names,
    v1 = v1a,
    v2 = v2a,
    difference = diff,
    percent_change = pct,
    row.names = NULL
  )
}