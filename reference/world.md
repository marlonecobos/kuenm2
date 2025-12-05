# World country polygons from Natural Earth

A spatial vector of the world countries. This is a simplified version of
the `countries110` from rnaturalearth R package.

## Format

A `Spatvector` object.

## Value

No return value. Used with function
[`vect`](https://rspatial.github.io/terra/reference/vect.html) to bring
vector variables to analysis.

## Examples

``` r
m <- terra::vect(system.file("extdata",
                             "world.gpkg",
                              package = "kuenm2"))
terra::plot(m)
```
