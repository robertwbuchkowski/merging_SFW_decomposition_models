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
  C_root_herb = "Herbaceous Root C",
  C_root_tree = "Tree Root C",
  Earthworm = "Earthworms",
  Litter = "Litter",
  CWD = "Coarse Woody Debris",
  Organic = "Fragmented Matter (Organic horizon)",
  DOM = "Dissolved Organic Matter",
  MIC =   "Microbial Biomass (Organic horizon)",
  P = "Particulate Organic C",
  L = "Low-weight Molecular C",
  A = "Aggregate C",
  M = "Mineral Associated Organic C",
  B = "Microbial Biomass (Mineral horizon)",
  Detritivore = "Detritivores",
  RootHerb = "Root Herbivores"
)

plot_order <- c(
  "Tree Root C",
  "Herbaceous Root C",
  "Litter",
  "Coarse Woody Debris",
  "Dissolved Organic Matter",
  "Fragmented Matter (Organic horizon)",
  "Microbial Biomass (Organic horizon)",
  "Particulate Organic C",
  "Low-weight Molecular C",
  "Aggregate C",
  "Mineral Associated Organic C",
  "Microbial Biomass (Mineral horizon)",
  "Earthworms",
  "Detritivores",
  "Root Herbivores"
)

keep_plot <- c(
  "Tree Root C",
  "Herbaceous Root C",
  "Litter",
  "Coarse Woody Debris",
  "Dissolved Organic Matter",
  "Fragmented Matter (Organic horizon)",
  "Microbial Biomass (Organic horizon)",
  "Particulate Organic C",
  "Low-weight Molecular C",
  "Aggregate C",
  "Mineral Associated Organic C",
  "Microbial Biomass (Mineral horizon)"
)