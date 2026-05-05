test_that("sf kriging returns complete prediction fields", {
  fixtures <- rtopng_sf_fixtures()

  set.seed(1)
  rtop_obj <- createRtopObject(
    fixtures$observations,
    fixtures$prediction_locations,
    params = fixtures$params,
    formulaString = "obs ~1"
  )
  rtop_obj <- rtopFitVariogram(rtop_obj, iprint = -1)
  rtop_pred <- rtopKrige(rtop_obj)

  expect_s3_class(rtop_pred$predictions, "sf")
  expect_equal(nrow(rtop_pred$predictions), 235)
  expect_false(anyNA(rtop_pred$predictions$var1.pred))
  expect_false(anyNA(rtop_pred$predictions$var1.var))
})

test_that("sf kriging preserves the seeded prediction anchors", {
  fixtures <- rtopng_sf_fixtures()

  set.seed(1)
  rtop_obj <- createRtopObject(
    fixtures$observations,
    fixtures$prediction_locations,
    params = fixtures$params,
    formulaString = "obs ~1"
  )
  rtop_obj <- rtopFitVariogram(rtop_obj, iprint = -1)
  rtop_pred <- rtopKrige(rtop_obj)

  expect_lt(abs(rtop_pred$predictions$var1.pred[1] - 0.01110463), 1e-7)
  expect_lt(abs(rtop_pred$predictions$var1.pred[2] - 0.01216192), 1e-7)
})

test_that("sf cross-validation keeps the legacy correlation anchor", {
  fixtures <- rtopng_sf_fixtures()

  set.seed(1)
  rtop_obj <- createRtopObject(
    fixtures$observations,
    fixtures$prediction_locations,
    params = fixtures$params,
    formulaString = "obs ~1"
  )
  rtop_obj <- rtopFitVariogram(rtop_obj, iprint = -1)
  rtop_cv <- rtopKrige(rtop_obj, cv = TRUE)

  expect_equal(
    cor(rtop_cv$predictions$observed, rtop_cv$predictions$var1.pred),
    0.7483928,
    tolerance = 1e-6
  )
})
