test_that("spacetime STSDF objects work end-to-end", {
  skip_if_not_installed("spacetime")

  fixtures <- utop_spacetime_fixtures(n_obs = 8, n_pred = 4, n_time = 3)

  rtop_obj <- createRtopObject(
    fixtures$observations,
    fixtures$prediction_locations,
    formulaString = "obs ~ 1",
    params = list(
      rresol = 4,
      rstype = "regular",
      debug.level = -1,
      nugget = FALSE
    )
  )

  expect_s3_class(rtop_obj, "rtop")
  expect_equal(unname(dim(rtop_obj$observations))[1], 8)
  expect_equal(unname(dim(rtop_obj$observations))[2], 3)

  rtop_obj <- rtopFitVariogram(rtop_obj, iprint = -1)
  expect_s3_class(rtop_obj$variogramModel, "rtopVariogramModel")

  result <- rtopKrige(rtop_obj)
  expect_s4_class(result$predictions, "STSDF")
  expect_equal(unname(dim(result$predictions))[1], 4)
  expect_equal(unname(dim(result$predictions))[2], 3)
  expect_false(anyNA(result$predictions@data$var1.pred))
})

test_that("rtopVariogram dispatches correctly for STSDF", {
  skip_if_not_installed("spacetime")

  fixtures <- utop_spacetime_fixtures(n_obs = 6, n_pred = 3, n_time = 2)
  vario <- rtopVariogram(fixtures$observations, formulaString = "obs ~ 1")

  expect_s3_class(vario, "rtopVariogram")
  expect_true("dist" %in% names(vario))
  expect_true("gamma" %in% names(vario))
})

test_that("spacetime pipeline tolerates nugget=TRUE", {
  skip_if_not_installed("spacetime")

  fixtures <- utop_spacetime_fixtures(n_obs = 6, n_pred = 3, n_time = 2)

  rtop_obj <- createRtopObject(
    fixtures$observations,
    fixtures$prediction_locations,
    formulaString = "obs ~ 1",
    params = list(
      rresol = 4,
      rstype = "regular",
      debug.level = -1,
      nugget = TRUE
    )
  )

  expect_s3_class(rtop_obj, "rtop")
  expect_true("overlapObs" %in% names(rtop_obj))
})
