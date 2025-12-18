# Project selected models to multiple sets of new data (scenarios)

This function performs predictions of selected models on multiple
scenarios, as specified in a `projection_data` object created with the
[`prepare_projection()`](https://marlonecobos.github.io/kuenm2/reference/prepare_projection.md)
function. In addition to generating predictions for each replicate, the
function calculates consensus measures (e.g., mean, median) across
replicates and models.

## Usage

``` r
project_selected(models, projection_data, out_dir, mask = NULL,
                 consensus_per_model = TRUE, consensus_general = TRUE,
                 consensus = c("median", "range", "mean", "stdev"),
                 write_replicates = FALSE, extrapolation_type = "E",
                 var_to_clamp = NULL, type = NULL, overwrite = FALSE,
                 parallel = FALSE, ncores = NULL,
                 progress_bar = TRUE, verbose = TRUE)
```

## Arguments

- models:

  an object of class `fitted_models` returned by the
  [`fit_selected()`](https://marlonecobos.github.io/kuenm2/reference/fit_selected.md)
  function.

- projection_data:

  an object of class `projection_data` returned by the
  [`prepare_projection()`](https://marlonecobos.github.io/kuenm2/reference/prepare_projection.md)
  function. This file contains the paths to the rasters representing
  each scenario.

- out_dir:

  (character) a path to a root directory for saving the raster file of
  each projection.

- mask:

  (SpatRaster, SpatVector, or SpatExtent) spatial object used to mask
  the variables before predict. Default is NULL.

- consensus_per_model:

  (logical) whether to calculate consensus across replicates when there
  are more than one replicate per model. Default is TRUE.

- consensus_general:

  (logical) whether to calculate consensus across models when there are
  more than one selected model. Default is TRUE.

- consensus:

  (character) consensus measures to calculate. Options available are
  'median', 'range', 'mean' and 'stdev' (standard deviation). Default is
  c("median", "range", "mean", "stdev").

- write_replicates:

  (logical) whether to write the projections for each replicate. Default
  is FALSE.

- extrapolation_type:

  (character) extrapolation type of model. Models can be transferred
  with three options: free extrapolation ('E'), extrapolation with
  clamping ('EC'), and no extrapolation ('NE'). Default = 'E'. See
  details.

- var_to_clamp:

  (character) vector specifying which variables to clamp. Only
  applicable if extrapolation_type is "EC" or "NE". Default is `NULL`,
  meaning all variables will be clamped or not extrapolated.

- type:

  (character) the format of prediction values. For `maxnet` models,
  valid options are `"raw"`, `"cumulative"`, `"logistic"`, and
  `"cloglog"`. For `glm` models, valid options are `"response"` and
  `"raw"`. If `NULL` (default), the function uses `"cloglog"` for
  `maxnet` models and `"response"` for `glm` models.

- overwrite:

  (logical) whether to overwrite SpatRaster if they already exists. Only
  applicable if `write_files` is set to TRUE. Default is FALSE.

- parallel:

  (logical) whether to fit the candidate models in parallel. Default is
  FALSE.

- ncores:

  (numeric) number of cores to use for parallel processing. Default is
  NULL and uses available cores - 1. This is only applicable if
  `parallel = TRUE`.

- progress_bar:

  (logical) whether to display a progress bar during processing. Default
  is TRUE.

- verbose:

  (logical) whether to display messages during processing. Default is
  TRUE.

## Value

A `model_projections` object that provides the paths to the raster files
with the projection results and the corresponding thresholds used to
binarize the predictions.

## See also

[`organize_future_worldclim()`](https://marlonecobos.github.io/kuenm2/reference/organize_future_worldclim.md),
[`prepare_projection()`](https://marlonecobos.github.io/kuenm2/reference/prepare_projection.md)

## Examples

``` r
# Step 1: Organize variables for current projection
## Import current variables (used to fit models)
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))

## Create a folder in a temporary directory to copy the variables
out_dir_current <- file.path(tempdir(), "Current_raw_wc")
dir.create(out_dir_current, recursive = TRUE)

## Save current variables in temporary directory
terra::writeRaster(var, file.path(out_dir_current, "Variables.tif"))


# Step 2: Organize future climate variables (example with WorldClim)
## Directory containing the downloaded future climate variables (example)
in_dir <- system.file("extdata", package = "kuenm2")

## Create a folder in a temporary directory to copy the future variables
out_dir_future <- file.path(tempdir(), "Future_raw_wc")

## Organize and rename the future climate data (structured by year and GCM)
### 'SoilType' will be appended as a static variable in each scenario
organize_future_worldclim(input_dir = in_dir, output_dir = out_dir_future,
                          name_format = "bio_", fixed_variables = var$SoilType)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%
#> 
#> Variables successfully organized in directory:
#> /tmp/Rtmp3qLrFS/Future_raw_wc

# Step 3: Prepare data to run multiple projections
## An example with maxnet models
## Import example of fitted_models (output of fit_selected())
data(fitted_model_maxnet, package = "kuenm2")

## Prepare projection data using fitted models to check variables
pr <- prepare_projection(models = fitted_model_maxnet,
                         present_dir = out_dir_current,
                         future_dir = out_dir_future,
                         future_period = c("2081-2100"),
                         future_pscen = c("ssp126", "ssp585"),
                         future_gcm = c("ACCESS-CM2", "MIROC6"),
                         raster_pattern = ".tif*")

# Step 4: Run multiple model projections
## A folder to save projection results
out_dir <- file.path(tempdir(), "Projection_results/maxnet_projections")
dir.create(out_dir, recursive = TRUE)

## Project selected models to multiple scenarios
p <- project_selected(models = fitted_model_maxnet, projection_data = pr,
                      out_dir = out_dir)
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
```
