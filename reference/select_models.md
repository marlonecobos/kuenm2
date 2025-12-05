# Select models that perform the best among candidates

This function selects the best models according to user-defined
criteria, evaluating statistical significance (partial ROC), predictive
ability (omission rates), and model complexity (AIC).

## Usage

``` r
select_models(calibration_results = NULL, candidate_models = NULL, data = NULL,
              algorithm = NULL, compute_proc = FALSE,
              addsamplestobackground = TRUE, weights = NULL,
              remove_concave = FALSE, omission_rate = NULL,
              allow_tolerance = TRUE, tolerance = 0.01,
              significance = 0.05, delta_aic = 2, parallel = FALSE,
              ncores = NULL, progress_bar = FALSE,verbose = TRUE)
```

## Arguments

- calibration_results:

  an object of class `calibration_results` returned by the
  [`calibration()`](https://marlonecobos.github.io/kuenm2/reference/calibration.md)
  function. Default is NULL.

- candidate_models:

  (data.frame) a summary of the evaluation metrics for each candidate
  model. Required only if `calibration_results` is NULL. In the output
  of the
  [`calibration()`](https://marlonecobos.github.io/kuenm2/reference/calibration.md),
  this data.frame is located in `$calibration_results$Summary`. Default
  is NULL.

- data:

  an object of class `prepared_data` returned by the
  [`prepare_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)
  function. Required only if `calibration_results` is NULL and
  `compute_proc` is TRUE.

- algorithm:

  (character) model algorithm, either "glm" or "maxnet". The default,
  NULL, uses the one defined as part of `calibration_results`, or
  `data`. If those arguments are not used, `algorithm` must be defined.

- compute_proc:

  (logical) whether to compute partial ROC tests for the selected
  models. This is required when partial ROC is not calculated for all
  candidate models during calibration. Default is FALSE.

- addsamplestobackground:

  (logical) whether to add to the background any presence sample that is
  not already there. Required only if `compute_proc` is TRUE and
  `calibration_results` is NULL.Default is TRUE.

- weights:

  (numeric) a numeric vector specifying weights for the occurrence
  records. Required only if `compute_proc` is TRUE and
  `calibration_results` is NULL. Default is NULL.

- remove_concave:

  (logical) whether to remove candidate models presenting concave
  curves. Default is FALSE.

- omission_rate:

  (numeric) the maximum omission rate a candidate model can have to be
  considered as a potentially selected model. The default, NULL, uses
  the value provided as part of `calibration_results`. For purposes of
  selection in existing results of evaluation, this value must match one
  of the values used in omission tests, and must be manually defined.

- allow_tolerance:

  (logical) whether to allow selection of models with minimum values of
  omission rates even if their omission rate surpasses the
  `omission_rate`. This is only applicable if all candidate models have
  omission rates higher than the `omission_rate`. Default is TRUE.

- tolerance:

  (numeric) The value added to the minimum omission rate if it exceeds
  the `omission_rate`. If `allow_tolerance = TRUE`, selected models will
  have an omission rate equal to or less than the minimum rate plus this
  tolerance. Default is 0.01.

- significance:

  (numeric) the significance level to select models based on the partial
  ROC (pROC). Default is 0.05. See Details.

- delta_aic:

  (numeric) the value of delta AIC used as a threshold to select models.
  Default is 2.

- parallel:

  (logical) whether to calculate the PROC of the candidate models in
  parallel. Default is FALSE.

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

If calibration_results is provided, it returns a new calibration_results
with the new selected models and summary. If calibration_results is
NULL, it returns a list containing the following elements:

- selected_models: data frame with the ID and the summary of evaluation
  metrics for the selected models.

- summary: A list containing the delta AIC values for model selection,
  and the ID values of models that failed to fit, had concave curves,
  non-significant pROC values, omission rates above the threshold, delta
  AIC values above the threshold, and the selected models.

## Details

Partial ROC is calculated following Peterson et al. (2008).

## Examples

``` r
# Import example of calibration results (output of calibration function)
## GLM
data(calib_results_glm, package = "kuenm2")

#Select new best models based on another value of omission rate
new_best_model <- select_models(calibration_results = calib_results_glm,
                                algorithm = "glm", compute_proc = TRUE,
                                omission_rate = 10)  # Omission error of 10
#> Selecting best among 122 models.
#> Calculating pROC...
#> 
#> Filtering 122 models.
#> Removing 0 model(s) because they failed to fit.
#> 21 model(s) were selected with omission rate below 10%.
#> Selecting 1 final model(s) with delta AIC <2.
#> Validating pROC of selected models...
#> 
#> All selected models have significant pROC values.

# Compare with best models selected previously
calib_results_glm$summary$Selected  # Model 86 selected
#> [1] 85
new_best_model$summary$Selected  # Models 64, 73 and 86 selected
#> [1] 85
```
