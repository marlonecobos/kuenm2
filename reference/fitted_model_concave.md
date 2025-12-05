# Fitted model with concave curves

A maxnet `fitted_models` object resulting from
[`fit_selected()`](https://marlonecobos.github.io/kuenm2/reference/fit_selected.md)
with a model presenting concave curves.

## Usage

``` r
data("fitted_model_concave")
```

## Format

A `fitted_models` with the following elements:

- species:

  Species names

- Models:

  A `list` with the fitted maxnet models (replicates and full models)

- calibration_data:

  A `data.frame` containing the variables extracted for presence and
  background points

- continuous_variables:

  A `character` indicating the names of the continuous variables

- categorical_variables:

  A `character` indicating the names of the categorical variables

- selected_models:

  A `data.frame` with formulas and evaluation metrics for each selected
  model

- weights:

  A `numeric` vector specifying weights for the occurrence records.
  `NULL` if no weights were set.

- pca:

  A `prcomp` object containing PCA results. `NULL` if PCA was not
  performed.

- addsamplestobackground:

  A `logical` value indicating whether to add any presence point not
  already included to the background.

- omission_rate:

  A `numeric` value indicating the omission rate used to evaluate
  models.

- thresholds:

  A `numeric` vector with thresholds used to binarize each replicate and
  the consensus (mean and median), calculated based on the omission rate
  defined in
  [`calibration()`](https://marlonecobos.github.io/kuenm2/reference/calibration.md).

- algorithm:

  A `character` string indicating the algorithm used (maxnet).

- partition_method:

  A `character` string indicating the partitioning method used.

- n_replicates:

  A `numeric` value indicating the number of replicates or folds.
