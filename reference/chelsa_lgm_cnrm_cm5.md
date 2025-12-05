# SpatRaster Representing LGM Conditions (GCM: CNRM-CM5)

Raster layer containing bioclimatic variables representing Last Glacial
Maximum (LGM) climatic conditions based on the CNRM-CM5 General
Circulation Model. The variables were resampled to 10arc-minutes and
masked using the `m` provided in the package. Data sourced from CHELSA:
<https://chelsa-climate.org/last-glacial-maximum-climate/>

## Format

A `SpatRaster` object.

## Value

No return value. Used with function
[`rast`](https://rspatial.github.io/terra/reference/rast.html) to bring
raster variables to analysis.

## Examples

``` r
chelsa_lgm_cnrm_cm5 <- terra::rast(system.file("extdata",
                                               "CHELSA_LGM_CNRM-CM5.tif",
                                               package = "kuenm2"))
terra::plot(chelsa_lgm_cnrm_cm5)
```
