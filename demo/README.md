# Demo material

The current demo workflow is based on the TNDTK example by Simone Persiano and collaborators:

- GitHub repository: <https://github.com/SimonePersiano/TNDTK>
- Original tutorial <>

The demo uses the Tyrol / South Tyrol subset distributed with the TNDTK tutorial: 27 gauged catchments, 27 gauging stations, 3 ungauged target catchments, and daily streamflow. Daily streamflow were used to estimate MAF. The `demo.gpkg` is a complete copy of the shapefiles used by Persiano et al. with the only difference of hardcoded MAF values. 

### `demo_rtop.R`

Runs MAF top-kriging using [`rtop`](https://cran.r-project.org/package=rtop). This follows the MAF branch of TNDTK Tutorial:

1. Fit the log-linear scaling relationship between observed MAF and catchment area.
2. Krige `obs = MAF / Area^c2` with catchment polygons as supports.
3. Rescale predictions by target catchment area.
4. Run leave-one-out CV for gauged catchments.
5. Print metrics with `yardstick::metric_set()` and `tidyhydro` metrics.

### `demo_ok.R`

Runs an ordinary-kriging baseline using [`automap`](https://cran.r-project.org/package=automap) on gauging station **points**.

The hydrological variable follows Farmer (2016):

```r
z = log(MAF / drainage_area)
```

Predictions are back-transformed and multiplied by target drainage area.

### `demo_uk.R`

Runs universal kriging using [`automap`](https://cran.r-project.org/package=automap) on gauging station **points**, with altitude as a covariate/drift term:

```r
log_unit_maf ~ Altitud
```

The `ungauged_catchments` layer has no altitude field. For this standalone demo, the script assigns each target the nearest gauged-station altitude. This is a placeholder. For production use, replace it with outlet altitude, station altitude, or catchment-mean altitude from a DEM or trusted metadata source.

## Leave-one-out CV comparison on gauged MAF

| Method | KGE2012 | PBIAS | RMSE | NSE | NSElog |
|---|---:|---:|---:|---:|---:|
| Top-Kriging | 0.963 | 2.45 | 0.697 | 0.978 | 0.945 |
| Ordinary Kriging | 0.960 | 0.162 | 0.916 | 0.961 | 0.929 |
| Universal Kriging | 0.978 | -1.60 | 0.816 | 0.969 | 0.955 |

## Ungauged MAF predictions (approximate as the ungauged stations altitude is unkown :-\ )

| Target | Top-Kriging | Ordinary Kriging | Universal Kriging |
|---|---:|---:|---:|
| Ahr_3 | 10.436 | 11.454 | 10.766 |
| Gader_1 | 6.178 | 5.789 | 5.512 |
| Isel_4 | 2.226 | 2.153 | 2.245 |


## References

- Castellarin, A., Persiano, S., Pugliese, A., Aloe, A., Skøien, J. O., and Pistocchi, A. (2018). Prediction of streamflow regimes over large geographical areas: interpolated flow-duration curves for the Danube Region. *Hydrological Sciences Journal*, 63(6), 845-861. <https://doi.org/10.1080/02626667.2018.1445855>
- Persiano, S., Pugliese, A., Aloe, A., Skøien, J. O., Castellarin, A., and Pistocchi, A. (2022). Streamflow data availability in Europe: a detailed dataset of interpolated flow-duration curves. *Earth System Science Data Discussions*. <https://doi.org/10.5194/essd-2021-469>
- Pugliese, A., Castellarin, A., and Brath, A. (2014). Geostatistical prediction of flow-duration curves in an index-flow framework. *Hydrology and Earth System Sciences*, 18, 3801-3816. <https://doi.org/10.5194/hess-18-3801-2014>
- Pugliese, A., Farmer, W. H., Castellarin, A., Archfield, S. A., and Vogel, R. M. (2016). Regional flow duration curves: geostatistical techniques versus multivariate regression. *Advances in Water Resources*, 96, 11-22. <https://doi.org/10.1016/j.advwatres.2016.06.008>
- Pugliese, A., Persiano, S., Bagli, S., Mazzoli, P., Parajka, J., Arheimer, B., Capell, R., Montanari, A., Blöschl, G., and Castellarin, A. (2018). A geostatistical data-assimilation technique for enhancing macro-scale rainfall-runoff simulations. *Hydrology and Earth System Sciences*, 22, 4633-4648. <https://doi.org/10.5194/hess-22-4633-2018>
- Skøien, J. O., Merz, R., and Blöschl, G. (2006). Top-kriging: geostatistics on stream networks. *Hydrology and Earth System Sciences*, 10, 277-287. <https://doi.org/10.5194/hess-10-277-2006>
- Farmer, W. H. (2016). Ordinary kriging as a tool to estimate historical daily streamflow records. *Hydrology and Earth System Sciences*, 20, 2721-2735. <https://doi.org/10.5194/hess-20-2721-2016>
- Laaha, G., Skøien, J. O., Nobilis, F., and Blöschl, G. (2013). Spatial prediction of stream temperatures using top-kriging with an external drift. *Environmental Modeling & Assessment*, 18, 671-683. <https://doi.org/10.1007/s10666-013-9373-3>
