# Calibration Results (Maxnet)

A `calibration_results` object resulted from
[`calibration()`](https://marlonecobos.github.io/kuenm2/reference/calibration.md)
using maxnet algorithm

## Usage

``` r
data("calib_results_maxnet")
```

## Format

A `calibration_results` with the following elements:

- species:

  Species names

- calibration_data:

  A `data.frame` with the variables extracted to presence and background
  points

- formula_grid:

  A `data.frame` with the ID, formulas, and regularization multipliers
  of each candidate model

- part_data:

  A `list` with the partition data, where each element corresponds to a
  replicate and contains the **indices of the test points** for that
  replicate

- partition_method:

  A `character` indicating the partition method

- n_replicates:

  A `numeric` value indicating the number of replicates or k-folds

- train_proportion:

  A `numeric` value indicating the proportion of occurrences used as
  train points when the partition method is 'subsample' or 'boostrap'

- data_xy:

  A `data.frame` with the coordinates of the occurrence and bakground
  points

- continuous_variables:

  A `character` indicating the names of the continuous variables

- categorical_variables:

  A `character` indicating the names of the categorical variables

- weights:

  A `numeric` value specifying weights for the occurrence records. It's
  NULL, meaning it was not set weights.

- pca:

  A `prcomp` object storing PCA information. Is NULL, meaning PCA was
  not performed

- algorithm:

  A `character` indicanting the algorithm (maxnet)

- calibration_results:

  A `list` containing the evaluation metrics for each candidate model

- omission_rate:

  A `numeric` value indicating the omission rate used to evaluate the
  models (10%)

- addsamplestobackground:

  A `logical` value indicating whether to add to the background any
  presence sample that is not already there.

- selected_models:

  A `data.frame` with the formulas and evaluation metrics for each
  selected model

- summary:

  A `list` with the number and the ID of the models removed and selected
  during selection procedure
