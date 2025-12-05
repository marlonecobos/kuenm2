# Partial ROC calculation for multiple candidate models

Computes partial ROC tests for multiple candidate models.

## Usage

``` r
partial_roc(formula_grid, data, omission_rate = 10,
            addsamplestobackground = TRUE, weights = NULL,
            algorithm = "maxnet", parallel = FALSE, ncores = NULL,
            progress_bar = TRUE)
```

## Arguments

- formula_grid:

  a data.frame with the grid of formulas defining the candidate models
  to test.

- data:

  an object of class `prepared_data` returned by the
  [`prepare_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)
  function or an object of class calibration_results returned by the
  [`calibration()`](https://marlonecobos.github.io/kuenm2/reference/calibration.md)
  function. It contains the calibration data and k-folds.

- omission_rate:

  (numeric) values from 0 to 100 representing the percentage of
  potential error due to any source of uncertainty. This value is used
  to calculate the omission rate. Default is 10. See details.

- addsamplestobackground:

  (logical) whether to add to the background any presence sample that is
  not already there. Default is TRUE.

- weights:

  (numeric) a numeric vector specifying weights for the occurrence
  records. Default is NULL.

- algorithm:

  (character) type algorithm, either "glm" or "maxnet". Default is
  "maxnet".

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

## Value

A data frame with summary statistics of the and AUC ratios and
significance calculated from the replicates of each candidate model.
Specifically, it includes the mean and standard deviation of these
metrics for each model.

## Details

Partial ROC is calculated following Peterson et al. (2008)
<doi:10.1016/j.ecolmodel.2007.11.008>.

## Examples

``` r
# Import prepared data to get model formulas
data(sp_swd, package = "kuenm2")

# Calculate proc for the first 5 candidate models
res_proc <- partial_roc(formula_grid = sp_swd$formula_grid[1:2,],
                        data = sp_swd, omission_rate = 10,
                        algorithm = "maxnet")
#>   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%
```
