# Fit models selected after calibration

This function fits models selected during model
[`calibration()`](https://marlonecobos.github.io/kuenm2/reference/calibration.md).

## Usage

``` r
fit_selected(calibration_results, replicate_method = "kfolds",
             n_replicates = 1, sample_proportion = 0.7, type = "cloglog",
             write_models = FALSE,
             file_name = NULL, parallel = FALSE, ncores = NULL,
             progress_bar = TRUE, verbose = TRUE, seed = 1)
```

## Arguments

- calibration_results:

  an object of class `calibration_results` returned by the
  [`calibration()`](https://marlonecobos.github.io/kuenm2/reference/calibration.md)
  function.

- replicate_method:

  (character) method used for producing replicates. Available options
  are `"kfolds"`, `"subsample"`, and `"bootstrap"`. See **Details** for
  more information.

- n_replicates:

  (numeric) number of replicates or folds to generate. If
  `replicate_method` is `"subsample"` or `"bootstrap"`, this defines the
  number of replicates. If `"kfolds"`, it specifies the number of folds.
  Default is 4.

- sample_proportion:

  (numeric) proportion of occurrence and background points to be used to
  fit model replicates. Only applicable when `replicate_method` is
  `"subsample"` or `"bootstrap"`. Default is 0.7 (i.e., 70% of the
  data).

- type:

  (character) the format of prediction values for computing thresholds.
  For maxnet models, valid options are "raw", "cumulative", "logistic",
  and "cloglog". For glm models, valid options are "cloglog", "response"
  and "raw". Default is "cloglog".

- write_models:

  (logical) whether to save the final fitted models to disk. Default is
  FALSE.

- file_name:

  (character) the file name, with or without a path, for saving the
  final models. This is only applicable if `write_models = TRUE`.

- parallel:

  (logical) whether to fit the final models in parallel. Default is
  FALSE.

- ncores:

  (numeric) number of cores to use for parallel processing. Default is
  NULL and uses available cores - 1. This is only applicable if
  `parallel = TRUE`.

- progress_bar:

  (logical) whether to display a progress bar during processing. Default
  is TRUE.

- verbose:

  (logical) whether to display detailed messages during processing.
  Default is TRUE.

- seed:

  (numeric) integer value used to specify an initial seed to split the
  data. Default is 1.

## Value

An object of class 'fitted_models' containing the following elements:

- species:

  a character string with the name of the species.

- Models:

  a list of fitted models, including replicates (fitted with part of the
  data) and full models (fitted with all data).

- calibration_data:

  a data.frame containing a column (`pr_bg`) that identifies occurrence
  points (1) and background points (0), along with the corresponding
  values of predictor variables for each point.

- selected_models:

  a data frame with the ID and summary of evaluation metrics for the
  selected models.

- weights:

  a numeric vector specifying weights for the predictor variables (if
  used).

- pca:

  a list of class
  [`prcomp`](https://rspatial.github.io/terra/reference/prcomp.html)
  representing the result of principal component analysis (if
  performed).

- addsamplestobackground:

  a logical value indicating whether any presence sample not already in
  the background was added.

- omission_rate:

  the omission rate determined during the calibration step.

- thresholds:

  the thresholds to binarize each replicate and the consensus (mean and
  median), calculated based on the omission rate set in
  [`calibration()`](https://marlonecobos.github.io/kuenm2/reference/calibration.md).

## Details

This function also computes model consensus (mean and median), the
thresholds to binarize model predictions based on the omission rate set
during model calibration to select models.

## Examples

``` r
# An example with maxnet models
data(calib_results_maxnet, package = "kuenm2")

# Fit models using calibration results
fm <- fit_selected(calibration_results = calib_results_maxnet,
                   n_replicates = 4)
#> Fitting replicates...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%
#> 
#> Fitting full models...
#>   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%

# Output the fitted models
fm
#> fitted_models object summary
#> ============================
#> Species: Myrcia hatschbachii 
#> Algortihm: maxnet 
#> Number of fitted models: 2 
#> Models fitted with 4 replicates

# An example with GLMs
data(calib_results_glm, package = "kuenm2")

# Fit models using calibration results
fm_glm <- fit_selected(calibration_results = calib_results_glm,
                       replicate_method = "subsample",
                       n_replicates = 5)
#> Fitting replicates...
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#> 
#> Fitting full models...
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%

# Output the fitted models
fm_glm
#> fitted_models object summary
#> ============================
#> Species: Myrcia hatschbachii 
#> Algortihm: glm 
#> Number of fitted models: 1 
#> Models fitted with 5 replicates
```
