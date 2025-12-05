# Advanced occurrence data cleaning

Advanced processes of data cleaning involving duplicate removal and
movement of records.

## Usage

``` r
advanced_cleaning(data, x, y, raster_layer, cell_duplicates = TRUE,
                  move_points_inside = FALSE, move_limit_distance = NULL,
                  verbose = TRUE)

remove_cell_duplicates(data, x, y,
                       raster_layer)

move_2closest_cell(data, x, y, raster_layer,
                   move_limit_distance, verbose = TRUE)
```

## Arguments

- data:

  data.frame with occurrence records. Rows with NA values will be
  omitted.

- x:

  (character) name of the column in `data` containing longitude values.

- y:

  (character) name of the column in `data` containing latitude values.

- raster_layer:

  a raster layer (object of class
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)).

- cell_duplicates:

  (logical) whether to remove duplicate coordinates considering raster
  cells. Default = TRUE.

- move_points_inside:

  (logical) whether to move records outside of raster cells with valid
  values to the closest cell with values. Default = FALSE.

- move_limit_distance:

  maximum distance to move records outside cells with valid values.
  Default = NULL. Must be defined if `move_points_inside` = TRUE.

- verbose:

  (logical) whether to print messages of progress. Default = TRUE.

## Value

A data.frame with occurrence records resulting from advanced cleaning
procedures. Other columns will be added to describe changes made in the
original data.

## Details

Data used in this functions should have gone through initial processes
of cleaning and filtering.

## See also

[`initial_cleaning()`](https://marlonecobos.github.io/kuenm2/reference/initial_cleaning.md)

## Examples

``` r
# Import occurrences
data(occ_data_noclean, package = "kuenm2")

# Import raster layers
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                   package = "kuenm2"))

# Keep only one layer
var <- var$bio_1

# all basic cleaning steps
clean_init <- initial_cleaning(data = occ_data_noclean, species = "species",
                               x = "x", y = "y", remove_na = TRUE,
                               remove_empty = TRUE, remove_duplicates = TRUE,
                               by_decimal_precision = TRUE,
                               decimal_precision = 2)

# Advanced cleaning steps
# exclude duplicates based on raster cell (pixel)
celldup <- remove_cell_duplicates(data = clean_init, x = "x", y = "y",
                                  raster_layer = var)

# move records to valid pixels
moved <- move_2closest_cell(data = celldup, x = "x", y = "y",
                            raster_layer = var, move_limit_distance = 10)
#> Moving occurrences to closest pixels...

# the steps at a time
clean_data <- advanced_cleaning(data = clean_init, x = "x", y = "y",
                                raster_layer = var, cell_duplicates = TRUE,
                                move_points_inside = TRUE,
                                move_limit_distance = 10)
#> Moving occurrences to closest pixels...
```
