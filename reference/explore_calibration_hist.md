# Explore variable distribution for occurrence and background points

This function prepares data to generate overlaid histograms to visualize
the distribution of predictor variables for occurrence (presence) and
background points.

## Usage

``` r
explore_calibration_hist(data, include_m = FALSE, raster_variables = NULL,
                         magnify_occurrences = 2, breaks = 15)
```

## Arguments

- data:

  an object of class `prepared_data` returned by the
  [`prepare_data`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)()
  function.

- include_m:

  (logical) whether to include data for plotting the histogram of the
  entire area from which background points were sampled. Default is
  FALSE, meaning only background and presence information will be
  plotted.

- raster_variables:

  (SpatRaster) predictor variables used to prepared the data with
  `prepared_data`. Only applicable if `include_m` is TRUE.

- magnify_occurrences:

  (numeric) factor by which the frequency of occurrences is magnified
  for better visualization. Default is 2, meaning occurrence frequencies
  in the plot will be doubled.

- breaks:

  (numeric) a single number giving the desired number of intervals in
  the histogram.

## Value

A list of with information to plot informative histograms to explore
data to be used in the modeling process. Histogram plots can be plotted
with the function
[`plot_calibration_hist()`](https://marlonecobos.github.io/kuenm2/reference/plot_calibration_hist.md).

## See also

[`plot_calibration_hist()`](https://marlonecobos.github.io/kuenm2/reference/plot_calibration_hist.md)

## Examples

``` r
# Import raster layers
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))

# Import occurrences
data(sp_swd, package = "kuenm2")

# Explore calibration data
calib_hist <- explore_calibration_hist(data = sp_swd,
                                       raster_variables = var,
                                       include_m = TRUE)

# To visualize results use the function plot_calibration_hist()
```
