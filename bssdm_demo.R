# Load libraries/functions
devtools::install_github("cebra-analytics/bssdm", upgrade="never")
library(bssdm)
library(terra)

# Climate data
climate_rast <- terra::rast(sprintf("data/climate/wc2.1bio%02d.tif", 1:19))

# Occurrence data (cleaned)
occurrences_cleaned <- read.csv(file = "data/fruit_fly_cleaned.csv")

# Range bagging
sdm.model <- rangebag(climate_rast, occurrences_cleaned)
rangebag_predict <- predict(sdm.model, climate_rast, raw_output = FALSE)
terra::plot(rangebag_predict, colNA = "blue",
            main = "Climate suitability (range bagging)")

# Climatch method
sdm.model <- climatch(climate_rast, occurrences_cleaned)
climatch_predict <- predict(sdm.model, climate_rast, raw_output = FALSE)
terra::plot(climatch_predict, colNA = "blue",
            main = "Climate suitability (climatch method")

