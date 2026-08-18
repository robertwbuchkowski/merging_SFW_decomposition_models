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
  pct <- (diff / baselinea) * 100
  
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

# ============================================================
# POOL NAMES  -  the single source of truth for pool labels and order.
# Every script and plot renames / orders / selects pools from THIS one vector.
# Names are raw pool codes; values are the display labels. The vector order is
# also the display order (used to order and to select pools for plotting).
# ============================================================
pool_names <- c(
  C_root_herb = "Herbaceous Root C",
  C_root_tree = "Tree Root C",
  Earthworm   = "Earthworms",
  Litter      = "Litter",
  CWD         = "Coarse Woody Debris (Organic horizon)",
  Organic     = "Particulate Organic C (Organic horizon)",
  DOM         = "Dissolved Organic Matter (Organic horizon)",
  MIC         = "Microbial Biomass (Organic horizon)",
  P           = "Particulate Organic C",
  L           = "Low Molecular Weight C",
  A           = "Aggregate C",
  M           = "Mineral Associated Organic C",
  B           = "Microbial Biomass",
  Detritivore = "Detritivores",
  RootHerb    = "Root Herbivores",
  TotalC      = "Total C"
)

# animal biomass pools (the fauna themselves, not the soil C they act on)
animal_pool_names <- c("Earthworm", "Detritivore", "RootHerb")

# relabel a vector of raw pool codes to display labels (unknown codes unchanged)
relabel_pools <- function(x)
  ifelse(x %in% names(pool_names), pool_names[x], x)

# display order = the order of pool_names (as a factor-levels vector of labels)
pool_order <- unname(pool_names)

# soil/plant pools to show in effect plots: everything except the animals
plot_pools <- unname(pool_names[!names(pool_names) %in% animal_pool_names])