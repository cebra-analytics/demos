# Load libraries/functions
devtools::install_github("cebra-analytics/bsrmap", upgrade="never")
library(bsrmap)
library(terra)

# Output template
output_template <- terra::rast("data/template_au.tif")
terra::plot(output_template, colNA = "blue", main = "Output template raster")

## Abiotic suitability (SDM) ####

# Conform global SDM layer to template
abiotic_suitability <- conform_layer(terra::rast("data/sdm_global.tif"),
                                     output_template)
terra::plot(abiotic_suitability, colNA = "blue",
            main = "Abiotic suitability (from SDM)")

## Biotic suitability (Host layers) ####

biotic_suitability_layers <- list()

# ACLUM land use
# NOTE: Download data from (as it is too large for GitHub):
# https://www.awe.gov.au/sites/default/files/documents/geotiff_clum_50m1220m.zip
aclum_rast <- terra::rast("data/clum_50m1220m.tif")
biotic_suitability_layers[["aclum"]] <- aggregate_categories(
  aclum_rast, output_template,
  categories = terra::cats(aclum_rast)[[1]]$ID,
  selected = c(340, 341, 342, 344, 345, 347, 348, 349, 350, 353, 365, 440,
               441, 442, 444, 445, 447, 448, 449, 450, 451, 453, 540, 541,
               542, 543, 544),
  binarize = TRUE)
terra::plot(biotic_suitability_layers[["aclum"]], colNA = "blue",
            main = "Biotic suitability: ACLUM")

# NDVI
ndvi_rast <- terra::rast("data/NDVI_Oct18_Mar19.grid")
biotic_suitability_layers[["ndvi"]] <- conform_layer(
  ndvi_rast, output_template,
  normalize = TRUE)
terra::plot(biotic_suitability_layers[["ndvi"]], colNA = "blue",
            main = "Biotic suitability: NDVI")

# Combine (multiply) biotic suitability layers
biotic_suitability <- combine_layers(
  terra::rast(biotic_suitability_layers),
  use_fun = "prod")
terra::plot(biotic_suitability, colNA = "blue",
            main = "Biotic suitability: total (product)")

## Pest suitability (abiotic x biotic layers) ####

# Combine (multiply) abiotic and biotic suitability layers
pest_suitability <- combine_layers(
  terra::rast(list(abiotic_suitability, biotic_suitability)),
  use_fun = "prod")
terra::plot(pest_suitability, colNA = "blue",
            main = "Pest suitability: abiotic x biotic")

## Pest pathways ####

# Layers

# Human population
population_density <- terra::rast("data/pop_density_1000m.tif")
population_density <- conform_layer(population_density, output_template)
terra::plot(population_density, colNA = "blue",
            main = "Human population density")

# Weighted distance from major sea ports
ports_df <- read.csv("data/ports.csv")
names(ports_df)[3:4] <- c("lon", "lat")
port_distance_weight <- distance_weight_layer(
  output_template, ports_df,
  beta = log(0.5)/20,
  weights = ports_df$Count)
terra::plot(port_distance_weight, colNA = "blue",
            main = "Pathway layer: Sea port distance weight")

# Pathways

pathway_list <- list()

# Returning residents across human population
pathway_list[["residents"]] <- pathway_likelihood(
  pathway_layers = population_density,
  leakage_rate_ci = c(3, 33),
  viability_rate_ci = c(0.001, 0.01))
terra::plot(log10(pathway_list[["residents"]] + 1e-8), colNA = "blue",
            main = "Pathway arrival probability: residents (log)")

# Vessels across layer weighted by sea port distance
pathway_list[["vessels"]] <- pathway_likelihood(
  pathway_layers = port_distance_weight,
  leakage_rate_ci = c(3, 30),
  viability_rate_ci = c(0.001, 0.01))
terra::plot(log10(pathway_list[["vessels"]] + 1e-8), colNA = "blue",
            main = "Pathway arrival probability: vessels (log)")

# Combine pathway layers
pathway_arrival_likelihood <- combine_layers(
  terra::rast(pathway_list),
  use_fun = "union",
  na.rm = TRUE)
terra::plot(log10(pathway_arrival_likelihood + 1e-6), colNA = "blue",
            main = "Pathway arrival likelihood: log total (union)")

## Establishment likelihood ####

# Combine pest suitability and pathway layers
establishment_likelihood <- combine_layers(
  terra::rast(list(pest_suitability, pathway_arrival_likelihood)),
  use_fun = "prod")
terra::plot(log10(establishment_likelihood + 1e-7), colNA = "blue",
            main = "Establishment likelihood: log total (product)")
