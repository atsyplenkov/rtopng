test_that("spatial workflow matches legacy outputs", {
  fixtures <- rtopng_spatial_fixtures()

  set.seed(1501)
  rtop_obj <- createRtopObject(
    fixtures$observations,
    fixtures$prediction_locations,
    params = fixtures$params,
    formulaString = "obs ~ 1"
  )
  rtop_obj <- rtopFitVariogram(rtop_obj, iprint = -1)

  expect_s3_class(rtop_obj$variogramModel, "rtopVariogramModel")

  rtop_cv <- rtopKrige(rtop_obj, cv = TRUE)
  rtop_pred <- rtopKrige(rtop_obj)
  varmat <- varMat(
    fixtures$observations,
    fixtures$prediction_locations,
    variogramModel = rtop_obj$variogramModel,
    gDistEst = TRUE,
    gDistPred = TRUE,
    rresol = 25,
    hresol = 3
  )
  rtop_reuse <- rtopKrige(rtop_cv)

  expect_s4_class(rtop_cv$predictions, "SpatialPolygonsDataFrame")
  expect_s4_class(rtop_pred$predictions, "SpatialPolygonsDataFrame")
  expect_s4_class(rtop_reuse$predictions, "SpatialPolygonsDataFrame")
  expect_equal(nrow(rtop_cv$predictions), 30)
  expect_equal(nrow(rtop_pred$predictions), 2)
  expect_equal(nrow(rtop_reuse$predictions), 2)
  expect_true(all(
    c("observed", "var1.pred", "var1.var") %in% names(rtop_cv$predictions)
  ))
  expect_true(all(c("var1.pred", "var1.var") %in% names(rtop_pred$predictions)))
  expect_true(isTRUE(all.equal(varmat$varMatObs, rtop_cv$varMatObs)))
  expect_true(isTRUE(all.equal(rtop_reuse$predictions, rtop_pred$predictions)))
  expect_lt(
    abs(
      cor(rtop_cv$predictions$observed, rtop_cv$predictions$var1.pred) -
        0.1678744283
    ),
    1e-7
  )

  rtop_updated <- varMat(rtop_reuse)
  rtop_updated <- updateRtopVariogram(rtop_updated, exp = 1.5, action = "mult")
  rtop_updated_mat <- varMat(rtop_updated)

  expect_false(isTRUE(all.equal(
    rtop_updated$varMatObs,
    rtop_updated_mat$varMatObs
  )))

  set.seed(1501)
  rtop_obj <- createRtopObject(
    fixtures$observations,
    fixtures$prediction_locations,
    params = fixtures$params,
    formulaString = "obs~1"
  )
  rtop_obj <- rtopFitVariogram(rtop_obj, iprint = -1)

  rtop_sim_5 <- rtopSim(rtop_obj, nsim = 5, logdist = TRUE, debug.level = -1)
  rtop_sim_input <- rtop_obj
  rtop_sim_input$predictionLocations <- rtop_sim_input$observations
  rtop_sim_input$observations$unc <- var(rtop_sim_5$observations$obs) *
    min(rtop_sim_5$observations$area) /
    rtop_sim_5$observations$area
  rtop_sim_input$predictionLocations$replaceNumber <- seq_len(nrow(
    rtop_sim_input$predictionLocations
  ))
  rtop_sim_10 <- rtopSim(
    rtop_sim_input,
    nsim = 10,
    replace = TRUE,
    debug.level = -1
  )

  expect_lt(abs(rtop_sim_5$simulations@data$sim1[1] - 0.01161453402), 1e-7)
  expect_lt(abs(rtop_sim_5$simulations@data$sim2[1] - 0.01064349883), 1e-7)
  expect_lt(abs(rtop_sim_5$simulations@data$sim4[2] - 0.01376182189), 1e-7)
  expect_lt(abs(rtop_sim_10$simulations@data$sim1[1] - 0.0139201299), 1e-7)
  expect_lt(abs(rtop_sim_10$simulations@data$sim6[1] - 0.0142161127), 1e-7)
  expect_lt(abs(rtop_sim_10$simulations@data$sim7[14] - 0.02216472038), 1e-7)
})

test_that("intamap branch matches the direct kriging output", {
  skip_if_not_installed("intamap")
  probe <- try(useRtopWithIntamap(), silent = TRUE)
  if (inherits(probe, "try-error")) {
    skip("intamap does not expose estimateParameters() on this system")
  }

  fixtures <- rtopng_spatial_fixtures()

  set.seed(1501)
  rtop_obj <- createRtopObject(
    fixtures$observations,
    fixtures$prediction_locations,
    params = fixtures$params,
    formulaString = "obs ~ 1"
  )
  rtop_obj <- rtopFitVariogram(rtop_obj, iprint = -1)
  rtop_pred <- rtopKrige(rtop_obj)

  output <- interpolate(
    fixtures$observations,
    fixtures$prediction_locations,
    optList = list(
      formulaString = obs ~ 1,
      gDist = TRUE,
      cloud = FALSE,
      nmax = 10,
      rresol = 25,
      hresol = 3
    ),
    methodName = "rtop",
    iprint = -1
  )

  expect_true(isTRUE(all.equal(
    rtop_pred$predictions@data$var1.pred,
    output$predictions@data$var1.pred
  )))
  expect_true(isTRUE(all.equal(
    rtop_pred$predictions@data$var1.var,
    output$predictions@data$var1.var
  )))
})
