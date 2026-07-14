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

# Name look ups:
name_lookup <- c(
  C_leaf_herb = "Herbaceous Leaf C",
  C_root_herb = "Herbaceous Root C",
  C_leaf_tree = "Tree Leaf C",
  C_wood_tree = "Tree Wood C",
  C_root_tree = "Tree Root C",
  Earthworm = "Earthworms",
  Litter = "Litter",
  CWD = "Coarse Woody Debris",
  Organic = "Organic Matter",
  DOM = "Dissolved Organic Matter",
  MIC =   "Microbial Biomass (Organic horizon)",
  P = "POC",
  L = "LWMC",
  A = "Aggregate C",
  M = "MAOC",
  B = "Microbial Biomass (Mineral Soil)",
  Detritivore = "Detritivores",
  RootHerb = "Root Herbivores"
)

plot_order <- c(
  "Tree Leaf C",
  "Tree Wood C",
  "Tree Root C",
  "Herbaceous Leaf C",
  "Herbaceous Root C",
  "Litter",
  "Coarse Woody Debris",
  "Dissolved Organic Matter",
  "Organic Matter",
  "Microbial Biomass (Organic horizon)",
  "POC",
  "LWMC",
  "Aggregate C",
  "MAOC",
  "Microbial Biomass (Mineral Soil)",
  "Earthworms",
  "Detritivores",
  "Root Herbivores"
)
