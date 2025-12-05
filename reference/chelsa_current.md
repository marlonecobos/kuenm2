# SpatRaster Representing present-day Conditions (CHELSA)

Raster layer containing bioclimatic variables representing present-day
climatic conditions. The variables were resampled to a 10 arc-minute
resolution and masked using the `m` region provided in the package. Data
sourced from CHELSA: <https://chelsa-climate.org/>

## Format

A `SpatRaster` object.

## Value

No return value. Used with function
[`rast`](https://rspatial.github.io/terra/reference/rast.html) to bring
raster variables to analysis.

## Examples

``` r
chelsa_current <- terra::rast(system.file("extdata",
                                           "Current_CHELSA.tif",
                                           package = "kuenm2"))
terra::plot(chelsa_current)
```
