# Example Bias File

A `SpatRaster` object representing a bias layer used for extracting
background points with the
[`prepare_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)
function.

## Format

A `SpatRaster` object.

## Value

No return value. Used with function
[`rast`](https://rspatial.github.io/terra/reference/rast.html) to bring
raster variables to analysis.

## Examples

``` r
bias <- terra::rast(system.file("extdata", "bias_file.tif",
                                package = "kuenm2"))

terra::plot(bias)
```
