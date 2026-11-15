# Analysis of extrapolation risks in partitions using the MOP metric

This function calculates environmental dissimilarities and identifies
non-analogous conditions by comparing the training data against the test
data for each partition, using the MOP (Mobility-Oriented Parity)
metric.

## Usage

``` r
explore_partition_extrapolation(data, include_train_background = TRUE,
                                       include_test_background = FALSE,
                                       variables = NULL,
                                       mop_type = "detailed",
                                       calculate_distance = TRUE,
                                       where_distance = "all",
                                       progress_bar = FALSE,
                                       progress_type = "timer", ...)
```

## Arguments

- data:

  an object of class `prepared_data` returned by the
  [`prepare_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)
  function.

- include_train_background:

  (logical) whether to include the background points used in training to
  define the environmental range of the training data. If set to FALSE,
  only the environmental conditions of the training presence records
  will be considered. Default is TRUE, meaning both presence and
  background points are used.

- include_test_background:

  (logical) whether to compute MOP for both the test presence records
  and the background points not used during training. Default is FALSE,
  meaning MOP will be calculated only for the test presences.

- variables:

  (character) names of the variables to be used in the MOP calculation.
  Default is NULL, meaning all variables in `data` will be used.

- mop_type:

  (character) type of MOP analysis to be performed. Options available
  are `"basic"`, `"simple"`, and `"detailed"`. Default is `"detailed"`.
  See
  [`projection_mop`](https://marlonecobos.github.io/kuenm2/reference/projection_mop.md)()
  for more details.

- calculate_distance:

  (logical) whether to calculate distances (dissimilarities) between
  train and test data. Default is TRUE.

- where_distance:

  (character) specifies which values in train data should be used to
  calculate distances. Options are: `"in_range"` (only conditions within
  the train range), `"out_range"` (only conditions outside the train
  range), and `"all"` (all conditions). Default is `"all"`.

- progress_bar:

  (logical) whether to display a progress bar during processing. Default
  is FALSE.

- progress_type:

  (character) type of the progress bar. See
  [`pbapply::pbtypes()`](https://peter.solymos.org/pbapply/reference/pboptions.html)
  for available progress bar types depending on operating system.
  Default is `"timer"`.

- ...:

  additional arguments passed to
  [`mop()`](https://rdrr.io/pkg/mop/man/mop.html).

## Value

An object of class `explore_partition` containing:

- `Mop_results`: a `data.frame` with MOP distances (if
  `calculate_distance = TRUE`), indicators of non-analogous conditions,
  number and names of variables outside the training range, and novel
  categories (if categorical variables are present) for each partition;

- `calibration_data`: a `data.frame` with the calibration data ordered
  by partition;

- `categorical_variables`: character vector with categorical variable
  names (if present);

- `pca`: PCA object (if `do_pca = TRUE` was used in
  [`prepare_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)).

## Examples

``` r
# Prepare data
# Import occurrences
data(occ_data, package = "kuenm2")

# Import raster layers
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))

# Prepare data for maxnet model
sp_swd <- prepare_data(algorithm = "maxnet", occ = occ_data,
                       x = "x", y = "y",
                       raster_variables = var,
                       species = occ_data[1, 1],
                       n_background = 100,
                       categorical_variables = "SoilType",
                       features = c("l", "lq"),
                       r_multiplier = 1,
                       partition_method = "kfolds")
#> Warning: 3 rows were excluded from database because NAs were found.

# Analysis of extrapolation risks in partitions
res <- explore_partition_extrapolation(data = sp_swd)
```
