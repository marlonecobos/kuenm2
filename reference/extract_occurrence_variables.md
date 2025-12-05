# Extracts Environmental Variables for Occurrences

This function extracts values from environmental or predictor variables
(`SpatRaster`) for georeferenced occurrence points. It also adds a
column indicating that these are presence points(pr_bg = 1).

## Usage

``` r
extract_occurrence_variables(occ, x, y, raster_variables)
```

## Arguments

- occ:

  A data.frame containing occurrence data. It must include columns with
  longitude (x) and latitude (y) coordinates.

- x:

  (character) a string specifying the name of the column in occ that
  contains the longitude values.

- y:

  (character) a string specifying the name of the column in occ that
  contains the latitude values.

- raster_variables:

  (SpatRaster) predictor variables used to calibrate the models.

## Value

A data.frame containing the original x and y coordinates of the
occurrence points (`x` and `y`), the values of the variables extracted
from `raster_variables`, and a new column `pr_bg` with a value of 1 for
all occurrences.

## Examples

``` r
# Import occurrences
data(occ_data, package = "kuenm2")

# Import raster layers
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))

# Extracts environmental variables for occurrences
occ_var <- extract_occurrence_variables(occ = occ_data, x = "x", y = "y",
                                        raster_variables = var)
```
