# SpatRaster Representing present-day Conditions (WorldClim)

Raster layer containing bioclimatic variables representing present-day
climatic conditions. The variables were obtained at a 10 arc-minute
resolution and masked using the `m` region provided in the package. Data
sourced from WorldClim: <https://worldclim.org/data/worldclim21.html>

## Format

A `SpatRaster` object.

## Value

No return value. Used with function
[`rast`](https://rspatial.github.io/terra/reference/rast.html) to bring
raster variables to analysis.

## Examples

``` r
var <- terra::rast(system.file("extdata",
                               "Current_variables.tif",
                                package = "kuenm2"))
terra::plot(var)
```
