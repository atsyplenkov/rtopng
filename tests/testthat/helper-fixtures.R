rtopng_extdata_path <- function() {
  system.file("extdata", package = "rtopng")
}

rtopng_read_sf <- function(layer) {
  sf::st_read(rtopng_extdata_path(), layer, quiet = TRUE)
}

rtopng_spatial_fixtures <- function() {
  observations <- rtopng_read_sf("observations")
  prediction_locations <- rtopng_read_sf("predictionLocations")

  observations <- as(observations, "Spatial")
  prediction_locations <- as(prediction_locations, "Spatial")

  observations <- observations[1:30, ]
  prediction_locations <- prediction_locations[1:2, ]
  observations$obs <- observations$QSUMMER_OB / observations$AREASQKM

  list(
    observations = observations,
    prediction_locations = prediction_locations,
    params = list(gDist = TRUE, cloud = FALSE, rresol = 25, hresol = 3, debug.level = -1)
  )
}

rtopng_sf_fixtures <- function() {
  observations <- rtopng_read_sf("observations")
  prediction_locations <- rtopng_read_sf("predictionLocations")
  observations$obs <- observations$QSUMMER_OB / observations$AREASQKM

  list(
    observations = observations,
    prediction_locations = prediction_locations,
    params = list(gDist = TRUE, cloud = FALSE)
  )
}
