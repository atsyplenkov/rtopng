# demo_ok.R
# Ordinary-kriging MAF baseline using automap.
#
# Literature reference from the hydro-kriging skill:
# - ch04 ordinary-kriging streamflow literature kriges log unit streamflow,
#   z = ln(Q / A), then back-transforms and rescales by target drainage area.
# - ch01 recommends automap::autoKrige() / autoKrige.cv() as a fast baseline
#   for ordinary kriging with automatic variogram fitting and CV.
#
# This file reads demo.gpkg created by demo_prep.R. Predicted/CV MAF values are
# kept in memory and are not written back to the GeoPackage.

library(sf)
library(automap)
library(yardstick)
library(tidyhydro)
library(dplyr)

nselog_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  nse_vec(log(truth), log(estimate), na_rm = na_rm, ...)
}

nselog <- function(data, ...) UseMethod("nselog")
nselog <- new_numeric_metric(nselog, direction = "maximize")

nselog.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  numeric_metric_summarizer(
    name = "nselog",
    fn = nselog_vec,
    data = data,
    truth = !!rlang::enquo(truth),
    estimate = !!rlang::enquo(estimate),
    na_rm = na_rm
  )
}

input_gpkg <- "demo.gpkg"

GaugedCatchments <- st_read(input_gpkg, layer = "gauged_catchments", quiet = TRUE)
UngaugedCatchments <- st_read(input_gpkg, layer = "ungauged_catchments", quiet = TRUE)
GaugedStations <- st_read(input_gpkg, layer = "gauged_stations", quiet = TRUE)

# Re-compute areas with sf. The OK literature baseline works with unit
# streamflow, so the kriged variable is log(MAF / drainage area).
GaugedCatchments$Area_km2 <- as.numeric(st_area(GaugedCatchments)) / (1000 * 1000)
UngaugedCatchments$Are_km2 <- as.numeric(st_area(UngaugedCatchments)) / (1000 * 1000)

GaugedStations <- GaugedStations %>%
  select(-any_of("MAF")) %>%
  left_join(
    st_drop_geometry(GaugedCatchments) %>% select(Cod, Area_km2, MAF),
    by = "Cod"
  )
GaugedStations$log_unit_maf <- log(GaugedStations$MAF / GaugedStations$Area_km2)

# Use the target outlet coordinates stored in the ungauged catchment attributes.
UngaugedPoints <- st_as_sf(
  st_drop_geometry(UngaugedCatchments),
  coords = c("X_3035", "Y_3035"),
  crs = st_crs(UngaugedCatchments)
)

# automap uses Spatial* objects internally.
gauged_sp <- as(GaugedStations, "Spatial")
ungauged_sp <- as(UngaugedPoints, "Spatial")

# Ordinary kriging: z ~ 1.
ok_model <- autoKrige(log_unit_maf ~ 1, gauged_sp, ungauged_sp, verbose = FALSE)

ungauged_maf <- st_drop_geometry(UngaugedPoints) %>%
  mutate(
    log_unit_maf_pred = ok_model$krige_output@data$var1.pred,
    log_unit_maf_var = ok_model$krige_output@data$var1.var,
    MAF_pred = exp(log_unit_maf_pred) * Are_km2
  ) %>%
  select(Locatin, MAF_pred, log_unit_maf_pred, log_unit_maf_var)

print(ungauged_maf)

# Leave-one-out CV, matching the validation pattern used in demo_rtop.R.
ok_cv <- autoKrige.cv(
  log_unit_maf ~ 1,
  gauged_sp,
  nfold = nrow(GaugedStations),
  verbose = c(FALSE, FALSE)
)

gauged_cv <- st_drop_geometry(GaugedStations) %>%
  transmute(
    Cod = as.character(Cod),
    MAF_obs = MAF,
    MAF_pred = exp(ok_cv$krige.cv_output@data$var1.pred) * Area_km2,
    MAF_resid = MAF_pred - MAF_obs
  ) %>%
  filter(is.finite(MAF_obs), is.finite(MAF_pred), MAF_obs > 0, MAF_pred > 0)

maf_metrics <- metric_set(kge2012, pbias, rmse, nse, nselog)
cv_metrics <- maf_metrics(gauged_cv, truth = MAF_obs, estimate = MAF_pred)

print(cv_metrics)
