# Evaluate models with independent data

This function evaluates the selected models using independent data
(i.e., data not used during model calibration). The function computes
omission rate and pROC, and optionally assesses whether the
environmental conditions in the independent data are analogous (i.e.,
within the range) to those in the calibration data.

## Usage

``` r
independent_evaluation(fitted_models, new_data,
                              consensus = c("mean", "median"),
                              type = "cloglog", extrapolation_type = "E",
                              var_to_clamp = NULL, perform_mop = TRUE,
                              mop_type = "detailed",
                              calculate_distance = TRUE,
                              where_distance = "all",
                              return_predictions = TRUE,
                              return_binary = TRUE,
                              progress_bar = FALSE, ...)
```

## Arguments

- fitted_models:

  an object of class `fitted_models` returned by the
  [`fit_selected`](https://marlonecobos.github.io/kuenm2/reference/fit_selected.md)()
  function.

- new_data:

  a `data.frame` containing environmental variables for independent test
  records. The column names must correspond exactly to the environmental
  variables used to fit the selected models, and each row to an
  individual test record.

- consensus:

  (character) vector specifying the types of consensus to use. Available
  options are `"median"` and `"mean"`. Default is `c("median", "mean")`.

- type:

  (character) the format of prediction values. For `maxnet` models,
  valid options are `"raw"`, `"cumulative"`, `"logistic"`, and
  `"cloglog"`. For `glm` models, valid options are `"response"` and
  `"cloglog"`. Default is `"cloglog"`.

- extrapolation_type:

  (character) extrapolation type of model. Models can be transferred
  with three options: free extrapolation ('E'), extrapolation with
  clamping ('EC'), and no extrapolation ('NE'). Default = 'E'. See
  details.

- var_to_clamp:

  (character) vector specifying which variables to clamp or not
  extrapolate. Only applicable if extrapolation_type is "EC" or "NE".
  Default is `NULL`, meaning all variables will be clamped or not
  extrapolated.

- perform_mop:

  (logical) whether to execute a Mobility-Oriented Parity (MOP)
  analysis. This analysis assesses if the environmental conditions in
  the `new_data` are analogous (within ranges) to those in the
  calibration data. Defaults to `TRUE`.

- mop_type:

  (character) type of MOP analysis to be performed. Options available
  are "basic", "simple" and "detailed". Default is 'simples'. See
  [`projection_mop`](https://marlonecobos.github.io/kuenm2/reference/projection_mop.md)()
  for more details.

- calculate_distance:

  (logical) whether to calculate distances (dissimilarities) between
  new_data and calibration data. Default is TRUE.

- where_distance:

  (character) specifies which values in `new_data` should be used to
  calculate distances. Options are: "in_range" (only conditions within
  the calibration range), "out_range" (only conditions outside the
  calibration range), and "all" (all conditions). Default is "all".

- return_predictions:

  (logical) whether to return continuous predictions at the locations of
  independent records in `new_data`. Default is TRUE.

- return_binary:

  (logical) whether to return binary predictions at the locations of
  independent records in `new_data`. The predictions are binarized using
  the respective thresholds stores in `fitted_models`. Default is TRUE.

- progress_bar:

  (logical) whether to display a progress bar during mop processing.
  Default is FALSE.

- ...:

  additional arguments passed to
  [`mop()`](https://rdrr.io/pkg/mop/man/mop.html).

## Value

A list containing the following elements:

- **evaluation**: A `data.frame` with omission rate and pROC values for
  each selected model and for the consensus.

- **mop_results**: (Only if `perform_mop = TRUE`) An object of class
  `mop_results`, with metrics of environmental similarity between
  calibration and independent data.

- **predictions**: (Only if `return_predictions = TRUE`) A `list` of
  `data.frames` containing continuous and binary predictions at the
  independent record locations, along with MOP distances, an indicator
  of whether environmental conditions at each location fall within the
  calibration range, and the identity of the variables that have lower
  and higher values than the calibration range. If the `fitted_models`
  object includes categorical variables, the returned `data.frame` will
  also contain columns indicating which values in `new_data` were not
  present in the calibration data.

## Examples

``` r
# Example with maxnet
# Import example of fitted_models (output of fit_selected())
data("fitted_model_maxnet", package = "kuenm2")

# Import independent records to evaluate the models
data("new_occ", package = "kuenm2")

# Import raster layers
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))

#Extract variables to occurrences
new_data <- extract_occurrence_variables(occ = new_occ, x = "x", y = "y",
                                         raster_variables = var)

#Add some fake data beyond the limits of calibration ranges
fake_data <- data.frame("pr_bg" = c(1, 1, 1),
                        "x" = c(NA, NA, NA),
                        "y" = c(NA, NA, NA),
                        "bio_1" = c(10, 15, 23),
                        "bio_7" = c(12, 16, 20),
                        "bio_12" = c(2300, 2000, 1000),
                        "bio_15" = c(30, 40, 50),
                        "SoilType" = c(1, 1, 1))
new_data <- rbind(new_data, fake_data)


# Evaluate models with independent data
res_ind <- independent_evaluation(fitted_models = fitted_model_maxnet,
                                  new_data = new_data)
```
