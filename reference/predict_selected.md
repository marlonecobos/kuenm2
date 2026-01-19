# Predict selected models for a single scenario

This function predicts selected models for a single set of new data
using either `maxnet` or `glm` It provides options to save the output
and compute consensus results (mean, median, etc.) across replicates and
models.

## Usage

``` r
predict_selected(models, new_variables, mask = NULL, write_files = FALSE,
                 write_replicates = FALSE, out_dir = NULL,
                 consensus_per_model = TRUE, consensus_general = TRUE,
                 consensus = c("median", "range", "mean", "stdev"),
                 extrapolation_type = "E", var_to_restrict = NULL,
                 type = "cloglog", overwrite = FALSE, progress_bar = TRUE)
```

## Arguments

- models:

  an object of class `fitted_models` returned by the
  [`fit_selected`](https://marlonecobos.github.io/kuenm2/reference/fit_selected.md)()
  function.

- new_variables:

  a SpatRaster or data.frame of predictor variables. The names of these
  variables must match those used to calibrate the models or those used
  to run PCA if `do_pca = TRUE` in the
  [`prepare_data`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)()
  function.

- mask:

  (SpatRaster, SpatVector, or SpatExtent) spatial object used to mask
  the variables before predict. Default is NULL.

- write_files:

  (logical) whether to save the predictions (SpatRasters or data.frame)
  to disk. Default is FALSE.

- write_replicates:

  (logical) whether to save the predictions for each replicates to disk.
  Only applicable if `write_files` is TRUE. Default is FALSE.

- out_dir:

  (character) directory path where predictions will be saved. Only
  relevant if `write_files = TRUE`.

- consensus_per_model:

  (logical) whether to compute consensus (mean, median, etc.) for each
  model across its replicates. Default is TRUE.

- consensus_general:

  (logical) whether to compute a general consensus across all models.
  Default is TRUE.

- consensus:

  (character) vector specifying the types of consensus to calculate
  across replicates and models. Available options are `"median"`,
  `"range"`, `"mean"`, and `"stdev"` (standard deviation). Default is
  `c("median", "range", "mean", "stdev")`.

- extrapolation_type:

  (character) extrapolation type of model. Models can be transferred
  with three options: free extrapolation ('E'), extrapolation with
  clamping ('EC'), and no extrapolation ('NE'). Default = 'E'. See
  details.

- var_to_restrict:

  (character) vector specifying which variables to clamp or not to
  extrapolate for. Only applicable if extrapolation_type is "EC" or
  "NE". Default is `NULL`, clamping and no extrapolation will be done
  for all variables.

- type:

  (character) the format of prediction values. For `maxnet` models,
  valid options are `"raw"`, `"cumulative"`, `"logistic"`, and
  `"cloglog"`. For `glm` models, valid options are `"cloglog"`,
  `"response"`, `"raw"`, `"cumulative"` and `"link"`. Default is
  "cloglog.

- overwrite:

  (logical) whether to overwrite SpatRasters if they already exist. Only
  applicable if `write_files = TRUE`. Default is FALSE.

- progress_bar:

  (logical) whether to display a progress bar during processing. Default
  is TRUE.

## Value

A list containing SpatRaster or data.frames predictions for each
replicate, long with the consensus results for each model and the
overall general consensus.

## Details

When predicting to areas where the variables are beyond the lower or
upper limits of the calibration data, users can choose to free
extrapolate the predictions (`extrapolation_type = "E"`), extrapolate
with clamping (extrapolation_type = "EC"), or not extrapolate
(extrapolation_type = "NE"). When clamping, the variables are set to
minimum and maximum values established for the maximum and minimum
values within calibration data. In the no extrapolation approach, any
cell with at least one variable listed in `var_to_restrict` falling
outside the calibration range is assigned a suitability value of 0.

## Examples

``` r
# Import variables to predict on
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))

# Example with maxnet
# Import example of fitted_models (output of fit_selected())
data("fitted_model_maxnet", package = "kuenm2")

# Predict to single scenario
p <- predict_selected(models = fitted_model_maxnet, new_variables = var)
#>   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%

# Example with GLMs
# Import example of fitted_models (output of fit_selected()) without replicates
data("fitted_model_glm", package = "kuenm2")

# Predict to single scenario
p_glm <- predict_selected(models = fitted_model_glm, new_variables = var)
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%

# Plot predictions
terra::plot(c(p$General_consensus$median, p_glm$General_consensus),
            col = rev(terrain.colors(240)), main = c("MAXNET", "GLM"),
            zlim = c(0, 1))
```
