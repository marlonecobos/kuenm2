# Fitting and evaluation of models, and selection of the best ones

This function fits and evaluates candidate models using the data and
grid of formulas prepared with
[`prepare_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md).
It supports both algorithms `glm` and `maxnet`. The function then
selects the best models based on unimodality (optional), partial ROC,
omission rate, and AIC values.

## Usage

``` r
calibration(data, error_considered, remove_concave = FALSE,
            proc_for_all = FALSE, omission_rate = NULL, delta_aic = 2,
            allow_tolerance = TRUE, tolerance = 0.01,
            addsamplestobackground = TRUE, use_weights = NULL,
            write_summary = FALSE, output_directory = NULL,
            skip_existing_models = FALSE, return_all_results = TRUE,
            parallel = FALSE, ncores = NULL, progress_bar = TRUE,
            verbose = TRUE)
```

## Arguments

- data:

  an object of class `prepared_data` returned by the
  [`prepare_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)
  function. It contains the calibration data, formulas grid, kfolds, and
  model type.

- error_considered:

  (numeric) values from 0 to 100 representing the percentage of
  potential error due to any source of uncertainty in your data. This
  value is used to calculate omission rates and partial ROC. See
  details.

- remove_concave:

  (logical) whether to remove candidate models presenting concave
  curves. Default is FALSE.

- proc_for_all:

  (logical) whether to apply partial ROC tests to all candidate models
  or only to the selected models. Default is FALSE, meaning that tests
  are applied only to the selected models.

- omission_rate:

  (numeric) values from 0 - 100, the maximum omission rate a candidate
  model can have to be considered as a potentially selected model. The
  default, NULL, uses the value in `error_considered`. If more that one
  value is used in `error_considered`, `omission_rate` must be defined.

- delta_aic:

  (numeric) the value of delta AIC used as a threshold to select models.
  Default is 2.

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

- addsamplestobackground:

  (logical) whether to add to the background any presence sample that is
  not already there. Default is TRUE.

- use_weights:

  (logical) whether to apply the weights present in the data. The
  default, NULL, uses weights provided in `data`. If they are not
  present in `data`, NULL weights are 1 for presences and 100 for
  background. If turned to FALSE, it uses NULL weights even if present
  in `data`.

- write_summary:

  (logical) whether to save the evaluation results for each candidate
  model to disk. Default is FALSE.

- output_directory:

  (character) the file name, with or without a path, for saving the
  evaluation results for each candidate model. This is only applicable
  if `write_summary = TRUE`.

- skip_existing_models:

  (logical) whether to check for and skip candidate models that have
  already been fitted and saved in `output_directory`. This is only
  applicable if `write_summary = TRUE`. Default is FALSE.

- return_all_results:

  (logical) whether to return the evaluation results for each replicate.
  Default is TRUE, meaning evaluation results for each replicate will be
  returned.

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

An object of class 'calibration_results' containing the following
elements:

- species: a character string with the name of the species.

- calibration data: a data.frame containing a column (`pr_bg`) that
  identifies occurrence points (1) and background points (0), along with
  the corresponding values of predictor variables for each point.

- formula_grid: data frame containing the calibration grid with possible
  formulas and parameters.

- kfolds: a list of vectors with row indices corresponding to each fold.

- data_xy: a data.frame with occurrence and background coordinates.

- continuous_variables: a character indicating the continuous variables.

- categorical_variables: a character, categorical variable names (if
  used).

- weights: a numeric vector specifying weights for data_xy (if used).

- pca: if a principal component analysis was performed with variables, a
  list of class "prcomp". See
  [`prcomp`](https://rdrr.io/r/stats/prcomp.html)() for details.

- algorithm: the model type (glm or maxnet)

- calibration_results: a list containing a data frame with all
  evaluation metrics for all partitions (if `return_all_results = TRUE`)
  and a summary of the evaluation metrics for each candidate model.

- omission_rate: The omission rate used to select models.

- addsamplestobackground: a logical value indicating whether any
  presence sample not already in the background was added.

- selected_models: data frame with the ID and the summary of evaluation
  metrics for the selected models.

- summary: A list containing the delta AIC values for model selection,
  and the ID values of models that failed to fit, had concave curves,
  non-significant pROC values, omission rates above the threshold, delta
  AIC values above the threshold, and the selected models.

## Details

Partial ROC is calculated using the values defined in `error_considered`
following Peterson et al. (2008).

Omission rates are calculated using separate testing data subsets. Users
can specify multiple values of `error_considered` to calculate this
metric (e.g., c(5, 10)), but only one can be used as the omission rate
for model selection.

Model fitting and complexity (AICc) is assessed using models generated
with the complete set of occurrences. AICc values are computed as
proposed by Warren and Seifert (2011).

## References

Ninomiya, Yoshiyuki, and Shuichi Kawano. "AIC for the Lasso in
generalized linear models." (2016): 2537-2560.

Warren, D. L., & Seifert, S. N. (2011). Ecological niche modeling in
Maxent: the importance of model complexity and the performance of model
selection criteria. Ecological applications, 21(2), 335-342.

## Examples

``` r
# Import occurrences
data(occ_data, package = "kuenm2")

# Import raster layers
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))
# Use only continuous variables
var <- var[[c("bio_1", "bio_7", "bio_12", "bio_15")]]

# Prepare data for maxnet model
sp_swd <- prepare_data(algorithm = "maxnet", occ = occ_data,
                       x = "x", y = "y",
                       raster_variables = var,
                       species = occ_data[1, 1],
                       n_background = 100,
                       features = c("l", "lq"),
                       r_multiplier = 1,
                       partition_method = "kfolds")
# Model calibration (maxnet)
m <- calibration(data = sp_swd, error_considered = 10)
#> Task 1/1: fitting and evaluating models:
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   5%  |                                                                              |======                                                                |   9%  |                                                                              |==========                                                            |  14%  |                                                                              |=============                                                         |  18%  |                                                                              |================                                                      |  23%  |                                                                              |===================                                                   |  27%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |=============================                                         |  41%  |                                                                              |================================                                      |  45%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%  |                                                                              |=========================================                             |  59%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |===================================================                   |  73%  |                                                                              |======================================================                |  77%  |                                                                              |=========================================================             |  82%  |                                                                              |============================================================          |  86%  |                                                                              |================================================================      |  91%  |                                                                              |===================================================================   |  95%  |                                                                              |======================================================================| 100%
#> 
#> 
#> Model selection step:
#> Selecting best among 22 models.
#> Calculating pROC...
#> 
#> Filtering 22 models.
#> Removing 0 model(s) because they failed to fit.
#> 7 model(s) were selected with omission rate below 10%.
#> Selecting 1 final model(s) with delta AIC <2.
#> Validating pROC of selected models...
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#> 
#> All selected models have significant pROC values.
m
#> calibration_results object summary (maxnet)
#> =============================================================
#> Species: Myrcia hatschbachii 
#> Number of candidate models: 22 
#>   - Models removed because they failed to fit: 0 
#>   - Models identified with concave curves: 0 
#>   - Model with concave curves not removed 
#>   - Models removed with non-significant values of pROC: 0 
#>   - Models removed with omission error > 10%: 15 
#>   - Models removed with delta AIC > 2: 6 
#> Selected models: 1 
#>   - Up to 5 printed here:
#>   ID                   Formulas Features R_multiplier pval_pROC_at_10.mean
#> 8  8 ~bio_1 + bio_7 + bio_15 -1        l            1                    0
#>   Omission_rate_at_10.mean dAIC Parameters
#> 8                   0.0978    0          3

# Prepare data for glm model
sp_swd_glm <- prepare_data(algorithm = "glm", occ = occ_data,
                           x = "x", y = "y",
                           raster_variables = var,
                           species = occ_data[1, 1],
                           n_background = 100,
                           features = c("l", "lq"),
                           partition_method = "kfolds")
m_glm <- calibration(data = sp_swd_glm, error_considered = 10)
#> Task 1/1: fitting and evaluating models:
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   5%  |                                                                              |======                                                                |   9%  |                                                                              |==========                                                            |  14%  |                                                                              |=============                                                         |  18%  |                                                                              |================                                                      |  23%  |                                                                              |===================                                                   |  27%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |=============================                                         |  41%  |                                                                              |================================                                      |  45%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%  |                                                                              |=========================================                             |  59%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |===================================================                   |  73%  |                                                                              |======================================================                |  77%  |                                                                              |=========================================================             |  82%  |                                                                              |============================================================          |  86%  |                                                                              |================================================================      |  91%  |                                                                              |===================================================================   |  95%  |                                                                              |======================================================================| 100%
#> 
#> 
#> Model selection step:
#> Selecting best among 22 models.
#> Calculating pROC...
#> 
#> Filtering 22 models.
#> Removing 0 model(s) because they failed to fit.
#> 13 model(s) were selected with omission rate below 10%.
#> Selecting 2 final model(s) with delta AIC <2.
#> Validating pROC of selected models...
#>   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%
#> 
#> All selected models have significant pROC values.
m_glm
#> calibration_results object summary (glm)
#> =============================================================
#> Species: Myrcia hatschbachii 
#> Number of candidate models: 22 
#>   - Models removed because they failed to fit: 0 
#>   - Models identified with concave curves: 1 
#>   - Model with concave curves not removed 
#>   - Models removed with non-significant values of pROC: 0 
#>   - Models removed with omission error > 10%: 9 
#>   - Models removed with delta AIC > 2: 11 
#> Selected models: 2 
#>   - Up to 5 printed here:
#>    ID                                                        Formulas Features
#> 12 12                        ~bio_1 + bio_7 + I(bio_1^2) + I(bio_7^2)       lq
#> 18 18 ~bio_1 + bio_7 + bio_12 + I(bio_1^2) + I(bio_7^2) + I(bio_12^2)       lq
#>    pval_pROC_at_10.mean Omission_rate_at_10.mean       dAIC Parameters
#> 12                    0                   0.0962 0.02361168          4
#> 18                    0                   0.0769 0.00000000          6
```
