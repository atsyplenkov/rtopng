rtopng_extdata_path <- function() {
  extdata_path <- system.file("extdata", package = "rtopng")
  if (nzchar(extdata_path)) {
    return(extdata_path)
  }

  if (requireNamespace("testthat", quietly = TRUE)) {
    extdata_path <- testthat::test_path("..", "..", "inst", "extdata")
    if (dir.exists(extdata_path)) {
      return(normalizePath(extdata_path))
    }
  }

  stop("Could not locate rtopng extdata for tests")
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
    params = list(
      gDist = TRUE,
      cloud = FALSE,
      rresol = 25,
      hresol = 3,
      debug.level = -1
    )
  )
}


rtopng_sf_subset_fixtures <- function(
  n_obs = 10,
  n_pred = 5,
  params = list(gDist = TRUE, cloud = FALSE)
) {
  observations <- rtopng_read_sf("observations")
  prediction_locations <- rtopng_read_sf("predictionLocations")
  observations$obs <- observations$QSUMMER_OB / observations$AREASQKM

  list(
    observations = observations[seq_len(n_obs), ],
    prediction_locations = prediction_locations[seq_len(n_pred), ],
    params = params
  )
}

rtopng_spatial_subset_fixtures <- function(
  n_obs = 10,
  n_pred = 2,
  params = list(
    gDist = TRUE,
    cloud = FALSE,
    rresol = 25,
    hresol = 3,
    debug.level = -1
  )
) {
  fixtures <- rtopng_spatial_fixtures()

  list(
    observations = fixtures$observations[seq_len(n_obs), ],
    prediction_locations = fixtures$prediction_locations[seq_len(n_pred), ],
    params = params
  )
}
