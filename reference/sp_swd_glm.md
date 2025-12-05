# Prepared Data for glm models

A `prepared_data` object resulted from
[`prepare_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)
to calibrate models using 'glm' algorithm.

## Usage

``` r
data("sp_swd_glm")
```

## Format

A `prepared_data` object with the following elements:

- species:

  Species names

- calibration_data:

  A `data.frame` containing the variables extracted for presence and
  background points

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

  A `character` indicanting the algorithm (glm)
