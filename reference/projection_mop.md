# Analysis of extrapolation risks in projections using the MOP metric

Calculates the mobility-oriented parity metric and other sub-products to
represent dissimilarities and non-analogous conditions when comparing a
set of reference conditions (M) against model projection conditions (G).

## Usage

``` r
projection_mop(data, projection_data, out_dir,
               subset_variables = FALSE, mask = NULL, type = "basic",
               na_in_range = TRUE, calculate_distance = FALSE,
               where_distance = "in_range", distance = "euclidean",
               scale = FALSE, center = FALSE, fix_NA = TRUE, percentage = 1,
               comp_each = 2000, tol = NULL, rescale_distance = FALSE,
               parallel = FALSE, ncores = NULL, progress_bar = TRUE,
               overwrite = FALSE)
```

## Arguments

- data:

  an object of class `fitted_models` returned by the
  [`enmpa::fit_selected()`](https://rdrr.io/pkg/enmpa/man/fit_glms.html)
  function or an object of class `prepared_data` returned by the
  [`prepare_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)
  function.

- projection_data:

  an object of class `projection_data` returned by the
  [`prepare_projection()`](https://marlonecobos.github.io/kuenm2/reference/prepare_projection.md)function.
  This file contains the paths to the rasters representing each
  scenario.

- out_dir:

  (character) a path to a root directory for saving the raster file of
  each projection.

- subset_variables:

  (logical) whether to include in the analysis only the variables
  present in the selected models. Default is FALSE

- mask:

  (SpatRaster, SpatVector, or SpatExtent) spatial object used to mask
  the variables (optional). Default is NULL.

- type:

  (character) type of MOP analysis to be performed. Options available
  are "basic", "simple" and "detailed". See Details for further
  information.

- na_in_range:

  (logical) whether to assign `NA` to regions within the projected
  area (G) where environmental conditions fall within the range of the
  calibration data (M). If `TRUE` (default), these regions are assigned
  `NA`. If `FALSE`, they are assigned `0` in the simple and basic MOP
  outputs, and `"within ranges"` in the detailed MOP output.

- calculate_distance:

  (logical) whether to calculate distances (dissimilarities) between m
  and g. The default, FALSE, runs rapidly and does not assess
  dissimilarity levels.

- where_distance:

  (character) where to calculate distances, considering how conditions
  in g are positioned in comparison to the range of conditions in m.
  Options available are "in_range", "out_range" and "all". Default is
  "in_range".

- distance:

  (character) which distances are calculated, euclidean or mahalanobis.
  Only applicable if calculate_distance = TRUE.

- scale:

  (logical or numeric) whether to scale as in
  [`scale`](https://rdrr.io/r/base/scale.html). Default is FALSE.

- center:

  (logical or numeric) whether to center as in
  [`scale`](https://rdrr.io/r/base/scale.html). Default is FALSE.

- fix_NA:

  (logical) whether to fix layers so cells with NA values are the same
  in all layers. Setting to FALSE may save time if the rasters are big
  and have no NA matching problems. Default is TRUE.

- percentage:

  (numeric) percentage of `m` closest conditions used to derive mean
  environmental distances to each combination of conditions in `g`.

- comp_each:

  (numeric) number of combinations in `g` to be used for distance
  calculations at a time. Increasing this number requires more RAM

- tol:

  (numeric) tolerance to detect linear dependencies when calculating
  Mahalanobis distances. The default, NULL, uses `.Machine$double.eps`.

- rescale_distance:

  (logical) whether to re-scale distances 0-1. Re-scaling prevents
  comparisons of dissimilarity values obtained from runs with different
  values of `percentage`.

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

- overwrite:

  (logical) whether to overwrite SpatRaster if they already exists. Only
  applicable if `write_files` is set to TRUE. Default is FALSE.

## Value

An object of class `mop_projections`, with the root directory and the
dataframe containing the file paths where the results were stored for
each scenario. The paths contain the following files:

- **summary** - a data.frame with details of the data used in the
  analysis:

  - *variables* - names of variables considered.

  - *type* - type of MOP analysis performed.

  - *scale* - value according to the argument `scale`.

  - *center* - value according to the argument `center`.

  - *calculate_distance* - value according to the argument
    `calculate_distance`.

  - *distance* - option regarding distance used.

  - *percentage* - percentage of `m` used as reference for distance
    calculation.

  - *rescale_distance* - value according to the argument
    `rescale_distance`.

  - *fix_NA* - value according to the argument `fix_NA`.

  - *N_m* - total number of elements (cells with values or valid rows)
    in `m`.

  - *N_g* - total number of elements (cells with values or valid rows)
    in `g`.

  - *m_min* - minimum values (lower limit) of the variables in reference
    conditions (`m`).

  - *m_max* - maximum values (upper limit) of the variables in reference
    conditions (`m`).

- **mop_distances** - if `calculate_distance` = TRUE, a SpatRaster or
  vector with distance values for the set of interest (`g`). Higher
  values represent greater dissimilarity compared to the set of
  reference (`m`).

- **mop_basic** - a SpatRaster or vector, for the set of interest,
  representing conditions in which at least one of the variables is
  non-analogous to the set of reference. Values should be: 1 for
  non-analogous conditions, and NA for conditions inside the ranges of
  the reference set.

- **mop_simple** - a SpatRaster or vector, for the set of interest,
  representing how many variables in the set of interest are
  non-analogous to those in the reference set. NA is used for conditions
  inside the ranges of the reference set.

- **mop_detailed** - a list containing:

  - *interpretation_combined* - a data.frame to help identify
    combinations of variables in *towards_low_combined* and
    *towards_high_combined* that are non-analogous to `m`.

  - *towards_low_end* - a SpatRaster or matrix for all variables
    representing where non-analogous conditions were found towards low
    values of each variable.

  - *towards_high_end* - a SpatRaster or matrix for all variables
    representing where non-analogous conditions were found towards high
    values of each variable.

  - *towards_low_combined* - a SpatRaster or vector with values
    representing the identity of the variables found to have
    non-analogous conditions towards low values. If vector,
    interpretation requires the use of the data.frame
    *interpretation_combined*.

  - *towards_high_combined* - a SpatRaster or vector with values
    representing the identity of the variables found to have
    non-analogous conditions towards high values. If vector,
    interpretation requires the use of the data.frame
    *interpretation_combined*.

## Details

`type` options return results that differ in the detail of how
non-analogous conditions are identified.

- **basic** - makes calculation as proposed by Owens et al. (2013)
  <doi:10.1016/j.ecolmodel.2013.04.011>.

- **simple** - calculates how many variables in the set of interest are
  non-analogous to those in the reference set.

- **detailed** - calculates five additional extrapolation metrics. See
  `mop_detailed` under `Value` below for full details.

`where_distance` options determine what values should be used to
calculate dissimilarity

- **in_range** - only conditions inside `m` ranges

- **out_range** - only conditions outside `m` ranges

- **all** - all conditions

When the variables used to represent conditions have different units,
scaling and centering are recommended. This step is only valid when
Euclidean distances are used.

## See also

[`organize_future_worldclim()`](https://marlonecobos.github.io/kuenm2/reference/organize_future_worldclim.md),
[`prepare_projection()`](https://marlonecobos.github.io/kuenm2/reference/prepare_projection.md)

## Examples

``` r
# Step 1: Organize variables for current projection
## Import current variables (used to fit models)
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))

## Create a folder in a temporary directory to copy the variables
out_dir_current <- file.path(tempdir(), "Current_raw4")
dir.create(out_dir_current, recursive = TRUE)

## Save current variables in temporary directory
terra::writeRaster(var, file.path(out_dir_current, "Variables.tif"))


# Step 2: Organize future climate variables (example with WorldClim)
## Directory containing the downloaded future climate variables (example)
in_dir <- system.file("extdata", package = "kuenm2")

## Create a folder in a temporary directory to copy the future variables
out_dir_future <- file.path(tempdir(), "Future_raw4")

## Organize and rename the future climate data (structured by year and GCM)
### 'SoilType' will be appended as a static variable in each scenario
organize_future_worldclim(input_dir = in_dir, output_dir = out_dir_future,
                          name_format = "bio_", fixed_variables = var$SoilType)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%
#> 
#> Variables successfully organized in directory:
#> /tmp/RtmphsTnoz/Future_raw4

# Step 3: Prepare data to run multiple projections
## An example with maxnet models
## Import example of fitted_models (output of fit_selected())
data(fitted_model_maxnet, package = "kuenm2")

## Prepare projection data using fitted models to check variables
pr <- prepare_projection(models = fitted_model_maxnet,
                         present_dir = out_dir_current,
                         future_dir = out_dir_future,
                         future_period = c("2041-2060", "2081-2100"),
                         future_pscen = c("ssp126", "ssp585"),
                         future_gcm = c("ACCESS-CM2", "MIROC6"),
                         raster_pattern = ".tif*")

# Step 4: Perform MOP for all projection scenarios
## Create a folder to save MOP results
out_dir <- file.path(tempdir(), "MOPresults")
dir.create(out_dir, recursive = TRUE)

## Run MOP
kmop <- projection_mop(data = fitted_model_maxnet, projection_data = pr,
                       out_dir = out_dir, type = "detailed")
#>   |                                                                              |                                                                      |   0%  |                                                                              |========                                                              |  11%  |                                                                              |================                                                      |  22%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================                                       |  44%  |                                                                              |=======================================                               |  56%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================                |  78%  |                                                                              |==============================================================        |  89%  |                                                                              |======================================================================| 100%
```
