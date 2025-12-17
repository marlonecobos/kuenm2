# Prepare data for model calibration with user-prepared calibration data

This function prepares data for model calibration using user-prepared
calibration data. It includes optional PCA, training/testing
partitioning, and the creation of a grid parameter combinations,
including distinct regularization multiplier values, various feature
classes, and different sets of environmental variables.

## Usage

``` r
prepare_user_data(algorithm, user_data, pr_bg, species = NULL, x = NULL,
                  y = NULL, features = c("lq", "lqp"),
                  r_multiplier = c(0.1, 0.5, 1, 2, 3),
                  user_formulas = NULL,
                  partition_method = "kfolds", n_partitions = 4,
                  train_proportion = 0.7, user_part = NULL,
                  categorical_variables = NULL,
                  do_pca = FALSE, center = TRUE, scale = TRUE,
                  exclude_from_pca = NULL, variance_explained = 95,
                  min_explained = 5, min_number = 2, min_continuous = NULL,
                  weights = NULL, include_xy = TRUE, write_pca = FALSE,
                  pca_directory = NULL, write_file = FALSE, file_name = NULL,
                  seed = 1)
```

## Arguments

- algorithm:

  (character) modeling algorithm, either "glm" or "maxnet".

- user_data:

  (data frame) A data.frame with a column with presence (1) and
  background (0) records, together with variable values (one variable
  per column). See an example with
  `data("user_data", package = "kuenm2")`.

- pr_bg:

  (character) the name of the column in `user_data` that contains the
  presence/background records.

- species:

  (character) string specifying the species name (optional). Default is
  NULL.

- x:

  (character) a string specifying the name of the column in `user_data`
  that contains the longitude values. Default is NULL. Must be defined
  if present in `user_data` otherwise it will be considered as another
  predictor variable.

- y:

  (character) a string specifying the name of the column in `user_data`
  that contains the latitude values. Default is NULL. Must be defined if
  present in `user_data` otherwise it will be considered as another
  predictor variable.

- features:

  (character) a vector of feature classes. Default is c("q", "lq", "lp",
  "qp", "lqp").

- r_multiplier:

  (numeric) a vector of regularization parameters for maxnet. Default is
  c(0.1, 1, 2, 3, 5).

- user_formulas:

  (character) Optional character vector with custom formulas provided by
  the user. See Details. Default is NULL.

- partition_method:

  (character) method used for data partitioning. Available options are
  `"kfolds"`, `"subsample"`, and `"bootstrap"`. See **Details** for more
  information. Default = "kfolds".

- n_partitions:

  (numeric) number of partitions to generate. If `partition_method` is
  `"subsample"` or `"bootstrap"`, this defines the number of training
  testing replicates. If `"kfolds"`, it specifies the number of folds.
  Must be \> 1; default = 4.

- train_proportion:

  (numeric) proportion of occurrence and background points to be used
  for model training in each partition. Only applicable when
  `partition_method` is `"subsample"` or `"bootstrap"`. Default is 0.7
  (i.e., 70% for training and 30% for testing).

- user_part:

  a user provided list with partitions or folds for cross-validation to
  be used in model calibration. Each element of the list should contain
  a vector of indices indicating the test points, which will be used to
  split `user_data` into training and testing sets. Useful in
  experiments that require exactly the same partition sets.

- categorical_variables:

  (character) names of the variables that are categorical. Default is
  NULL.

- do_pca:

  (logical) whether to perform a principal component analysis (PCA) with
  the set of variables. Default is FALSE.

- center:

  (logical) whether the variables should be zero-centered. Default is
  TRUE.

- scale:

  (logical) whether the variables should be scaled to have unit variance
  before the analysis takes place. Default is FALSE.

- exclude_from_pca:

  (character) variable names within raster_variables that should not be
  included in the PCA transformation. Instead, these variables will be
  added directly to the final set of output variables without being
  modified. The default is NULL, meaning all variables will be used
  unless specified otherwise.

- variance_explained:

  (numeric) the cumulative percentage of total variance that must be
  explained by the selected principal components. Default is 95.

- min_explained:

  (numeric) the minimum percentage of total variance that a principal
  component must explain to be retained. Default is 5.

- min_number:

  (numeric) the minimum number of variables to be included in the model
  formulas to be generated.

- min_continuous:

  (numeric) the minimum number of continuous variables required in a
  combination. Default is NULL.

- weights:

  (numeric) a numeric vector specifying weights for the occurrence
  records. Default is NULL.

- include_xy:

  (logical) whether to include the coordinates (longitude and latitude)
  in the results from preparing data. Default is TRUE.

- write_pca:

  (logical) whether to save the PCA-derived raster layers (principal
  components) to disk. Default is FALSE.

- pca_directory:

  (character) the path or name of the folder where the PC raster layers
  will be saved. This is only applicable if `write_pca = TRUE`. Default
  is NULL.

- write_file:

  (logical) whether to write the resulting prepared_data list in a local
  directory. Default is FALSE.

- file_name:

  (character) the path or name of the folder where the resulting list
  will be saved. This is only applicable if `write_file = TRUE`. Default
  is NULL.

- seed:

  (numeric) integer value to specify an initial seed to split the data.
  Default is 1.

## Value

An object of class `prepared_data` containing all elements necessary to
perform further explorations of data and run a model calibration
routine.

## Details

Training and testing are performed multiple times (i.e., the number set
in `n_partitions`), and model selection is based on the average
performance of models after running this routine. A description of the
available data partitioning methods is below:

- **"kfolds"**: Splits the dataset into *K* subsets (folds) of
  approximately equal size, keeping proportion of 0 and 1 stable
  compared to the full set. In each training/test run, one fold is used
  as the test set, while the remaining folds are combined to form the
  training set.

- **"bootstrap"**: Creates the training dataset by sampling observations
  from the original dataset *with replacement* (i.e., the same
  observation can be selected multiple times). The test set consists of
  the observations that were not selected in that specific sampling.

- **"subsample"**: Similar to bootstrap, but the training set is created
  by sampling *without replacement* (i.e., each observation is selected
  at most once).

`user_formulas` must be a character vector of model formulas. Supported
terms include linear effects, quadratic terms (e.g., `I(bio_7^2)`),
products (e.g., `bio_1:bio_7`), hinge (e.g., `hinge(bio_1)`), threshold
(e.g., `thresholds(bio_2)`), and categorical predictors (e.g.,
`categorical(SoilType)`). Example of a valid formula:
`~ bio_1 + bio_7 + I(bio_7^2) + bio_1:bio_7 + hinge(bio_1) + thresholds(bio_2) + categorical(SoilType)`.
All variables appearing in the formulas must exist in the data.frame
supplied as `user_data`.

## See also

[`calibration()`](https://marlonecobos.github.io/kuenm2/reference/calibration.md),
[`explore_calibration_hist()`](https://marlonecobos.github.io/kuenm2/reference/explore_calibration_hist.md),
[`explore_partition_env()`](https://marlonecobos.github.io/kuenm2/reference/explore_partition_env.md),
[`explore_partition_geo()`](https://marlonecobos.github.io/kuenm2/reference/explore_partition_geo.md),
[`explore_partition_extrapolation()`](https://marlonecobos.github.io/kuenm2/reference/explore_partition_extrapolation.md),
[`plot_calibration_hist()`](https://marlonecobos.github.io/kuenm2/reference/plot_calibration_hist.md),
[`plot_explore_partition()`](https://marlonecobos.github.io/kuenm2/reference/plot_explore_partition.md)

## Examples

``` r
# Import user-prepared data
data("user_data", package = "kuenm2")

# Prepare data for maxnet model
maxnet_swd_user <- prepare_user_data(algorithm = "maxnet",
                                     user_data = user_data, pr_bg = "pr_bg",
                                     species = "Myrcia hatschbachii",
                                     categorical_variables = "SoilType",
                                     features = c("l", "q", "p", "lq", "lqp"),
                                     r_multiplier = c(0.1, 1, 2, 3, 5))
maxnet_swd_user
#> prepared_data object summary
#> ============================
#> Species: Myrcia hatschbachii 
#> Number of Records: 527 
#>   - Presence: 51 
#>   - Background: 476 
#> Partition Method: kfolds 
#>   - Number of kfolds: 4 
#> Continuous Variables:
#>   - bio_1, bio_7, bio_12, bio_15 
#> Categorical Variables:
#>   - SoilType 
#> PCA Information: PCA not performed
#> Weights: No weights provided
#> Calibration Parameters:
#>   - Algorithm: maxnet 
#>   - Number of candidate models: 610 
#>   - Features classes (responses): l, q, p, lq, lqp 
#>   - Regularization multipliers: 0.1, 1, 2, 3, 5 

# Prepare data for glm model
glm_swd_user <- prepare_user_data(algorithm = "glm",
                                  user_data = user_data, pr_bg = "pr_bg",
                                  species = "Myrcia hatschbachii",
                                  categorical_variables = "SoilType",
                                  features = c("l", "q", "p", "lq", "lqp"))
glm_swd_user
#> prepared_data object summary
#> ============================
#> Species: Myrcia hatschbachii 
#> Number of Records: 527 
#>   - Presence: 51 
#>   - Background: 476 
#> Partition Method: kfolds 
#>   - Number of kfolds: 4 
#> Continuous Variables:
#>   - bio_1, bio_7, bio_12, bio_15 
#> Categorical Variables:
#>   - SoilType 
#> PCA Information: PCA not performed
#> Weights: No weights provided
#> Calibration Parameters:
#>   - Algorithm: glm 
#>   - Number of candidate models: 122 
#>   - Features classes (responses): l, q, p, lq, lqp 
```
