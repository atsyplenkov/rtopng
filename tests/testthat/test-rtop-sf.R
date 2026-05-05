test_that("sf workflow matches legacy outputs", {
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
  rtop_cv <- rtopKrige(rtop_obj, cv = TRUE)

  expect_s3_class(rtop_pred$predictions, "sf")
  expect_equal(nrow(rtop_pred$predictions), 235)
  expect_lt(abs(rtop_pred$predictions$var1.pred[1] - 0.01110463), 1e-7)
  expect_lt(abs(rtop_pred$predictions$var1.pred[2] - 0.01216192), 1e-7)
  expect_false(anyNA(rtop_pred$predictions$var1.pred))
  expect_false(anyNA(rtop_pred$predictions$var1.var))

  rtop_cor <- cor(rtop_cv$predictions$observed, rtop_cv$predictions$var1.pred)
  expect_lt(abs(rtop_cor - 0.7483928), 1e-6)
})
