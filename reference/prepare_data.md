# Prepare data for model calibration

This function prepares data for model calibration, including optional
PCA, background point generation, training/testing partitioning, and the
creation of a grid of parameter combinations, including regularization
multiplier values, feature classes, and sets of environmental variables.

## Usage

``` r
prepare_data(algorithm, occ, x, y, raster_variables, species = NULL,
             n_background = 1000, features = c("lq", "lqp"),
             r_multiplier = c(0.1, 0.5, 1, 2, 3),
             user_formulas = NULL,
             partition_method = "kfolds",
             n_partitions = 4, train_proportion = 0.7,
             categorical_variables = NULL,
             do_pca = FALSE, center = TRUE, scale = TRUE,
             exclude_from_pca = NULL, variance_explained = 95,
             min_explained = 5, min_number = 2, min_continuous = NULL,
             bias_file = NULL, bias_effect = NULL, weights = NULL,
             include_xy = TRUE, write_pca = FALSE, pca_directory = NULL,
             write_file = FALSE, file_name = NULL, seed = 1)
```

## Arguments

- algorithm:

  (character) modeling algorithm, either "glm" or "maxnet".

- occ:

  (data frame) a data.frame containing the coordinates (longitude and
  latitude) of the occurrence records.

- x:

  (character) a string specifying the name of the column in `occ` that
  contains the longitude values.

- y:

  (character) a string specifying the name of the column in `occ` that
  contains the latitude values.

- raster_variables:

  (SpatRaster) predictor variables from which environmental values will
  be extracted using `occ` and a background will be sampled. Must
  correspond geographically with the area where model is calibrated.

- species:

  (character) string specifying the species name (optional). Default is
  NULL.

- n_background:

  (numeric) number of points to represent the background for the model.
  Default is 1000.

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

  (numeric) the minimum number of variables to be included in model
  formulas to be generated. Default = 2.

- min_continuous:

  (numeric) the minimum number of continuous variables required in a
  combination. Default is NULL.

- bias_file:

  (SpatRaster) a raster containing bias values (probability weights)
  that influence the selection of background points. It must have the
  same extent, resolution, and number of cells as the raster variables.
  Default is NULL.

- bias_effect:

  (character) a string specifying how the values in the `bias_file`
  should be interpreted. Options are "direct" or "inverse". If "direct",
  higher values in bias file increase the likelihood of selecting
  background points. If "inverse", higher values decrease the
  likelihood. Default = NULL. Must be defined if `bias_file` is
  provided.

- weights:

  (numeric) a numeric vector specifying weights for the occurrence
  records. The default, NULL, uses 1 for presence and 100 for
  background.

- include_xy:

  (logical) whether to include the coordinates (longitude and latitude)
  in the results from preparing data. Columns containing coordinates
  will be renamed as "x" and "y". Default is TRUE.

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

  (character) name of file (no extension needed) to write resulting
  object in a local directory. Only needed if `write_file = TRUE`.
  Default is NULL.

- seed:

  (numeric) integer value to specify an initial seed to split the data
  and extract background. Default is 1.

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
All variables appearing in the formulas must exist in the raster
supplied as `raster_variables`.

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
# Import occurrences
data(occ_data, package = "kuenm2")

# Import raster layers
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))

# Import a bias file
bias <- terra::rast(system.file("extdata", "bias_file.tif",
                                package = "kuenm2"))

# Prepare data for maxnet model
sp_swd <- prepare_data(algorithm = "maxnet", occ = occ_data,
                       x = "x", y = "y",
                       raster_variables = var,
                       species = occ_data[1, 1],
                       categorical_variables = "SoilType",
                       n_background = 500, bias_file = bias,
                       bias_effect = "direct",
                       features = c("l", "q", "p", "lq", "lqp"),
                       r_multiplier = c(0.1, 1, 2, 3, 5))
#> Warning: 27 rows were excluded from database because NAs were found.
print(sp_swd)
#> prepared_data object summary
#> ============================
#> Species: Myrcia hatschbachii 
#> Number of Records: 524 
#>   - Presence: 51 
#>   - Background: 473 
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
sp_swd_glm <- prepare_data(algorithm = "glm", occ = occ_data,
                           x = "x", y = "y",
                           raster_variables = var,
                           species = occ_data[1, 1],
                           categorical_variables = "SoilType",
                           n_background = 500, bias_file = bias,
                           bias_effect = "direct",
                           features = c("l", "q", "p", "lq", "lqp"))
#> Warning: 27 rows were excluded from database because NAs were found.
print(sp_swd_glm)
#> prepared_data object summary
#> ============================
#> Species: Myrcia hatschbachii 
#> Number of Records: 524 
#>   - Presence: 51 
#>   - Background: 473 
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
