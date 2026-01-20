# Analysis of extrapolation risks using the MOP metric (for single scenario)

Calculates the mobility-oriented parity metric and other sub-products to
represent dissimilarities and non-analogous conditions when comparing a
set of reference conditions (M) against another set of scenario
conditions (G).

## Usage

``` r
single_mop(
  data,
  new_variables,
  subset_variables = FALSE,
  mask = NULL,
  type = "basic",
  na_in_range = TRUE,
  calculate_distance = FALSE,
  where_distance = "in_range",
  distance = "euclidean",
  scale = FALSE,
  center = FALSE,
  fix_NA = TRUE,
  percentage = 1,
  comp_each = 2000,
  tol = NULL,
  rescale_distance = FALSE,
  parallel = FALSE,
  ncores = NULL,
  progress_bar = TRUE,
  write_files = FALSE,
  out_dir = NULL,
  overwrite = FALSE
)
```

## Arguments

- data:

  an object of class `fitted_models` returned by the
  [`fit_selected()`](https://marlonecobos.github.io/kuenm2/reference/fit_selected.md)
  function, an object of class `prepared_data` returned by the
  [`prepare_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)
  function, or an object of class `calibration_results` returned by the
  [`calibration()`](https://marlonecobos.github.io/kuenm2/reference/calibration.md)
  function.

- new_variables:

  a SpatRaster or data.frame of predictor variables. The names of these
  variables must match those used to prepare the date or calibrate the
  models provided in `data`.

- subset_variables:

  (logical) whether to include in the analysis only the variables
  present in the selected models. Only applicable if `data` is a
  `fitted_models` object. Default is FALSE

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

  (logical) whether to compute MOP in parallel. Default is FALSE.

- ncores:

  (numeric) number of cores to use for parallel processing. Default is
  NULL and uses available cores - 1. This is only applicable if
  `parallel = TRUE`.

- progress_bar:

  (logical) whether to display a progress bar during distance
  calculation. Only applicable if calculate_distance is TRUE. Default is
  TRUE.

- write_files:

  (logical) whether to save the MOP results (SpatRasters and
  data.frames) to disk. Default is FALSE.

- out_dir:

  (character) directory path where results will be saved. Only relevant
  if `write_files = TRUE`.

- overwrite:

  (logical) whether to overwrite SpatRasters if they already exist. Only
  applicable if `write_files = TRUE`. Default is FALSE.

## Value

An object of class `mop_results` containing:

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

  - *m_ranges* - the range (minimum and maximum values) of the variable
    in reference conditions (`m`)

- **mop_distances** - if `calculate_distance` = TRUE, a SpatRaster with
  distance values for the set of interest (`g`). Higher values represent
  greater dissimilarity compared to the set of reference (`m`).

- **mop_basic** - a SpatRaster for the set of interest representing
  conditions in which at least one of the variables is non-analogous to
  the set of reference. Values should be: 1 for non-analogous
  conditions, and NA for conditions inside the ranges of the reference
  set.

- **mop_simple** - a SpatRaster, for the set of interest, representing
  how many variables in the set of interest are non-analogous to those
  in the reference set. NA is used for conditions inside the ranges of
  the reference set.

- **mop_detailed** - a list containing:

  - *interpretation_combined* - a data.frame to help identify
    combinations of variables in *towards_low_combined* and
    *towards_high_combined* that are non-analogous to `m`.

  - *towards_low_end* - a SpatRaster for all variables representing
    where non-analogous conditions were found towards low values of each
    variable.

  - *towards_high_end* - a SpatRaster for all variables representing
    where non-analogous conditions were found towards high values of
    each variable.

  - *towards_low_combined* - a SpatRaster with values representing the
    identity of the variables found to have non-analogous conditions
    towards low values.

  - *towards_high_combined* - a SpatRaster with values representing the
    identity of the variables found to have non-analogous conditions
    towards high values.

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

[`projection_mop()`](https://marlonecobos.github.io/kuenm2/reference/projection_mop.md)

## Examples

``` r
# Import an example of fitted models (output of fit_selected())
data("fitted_model_maxnet", package = "kuenm2")

# Import variables under a new set of conditions
# Here, future climate data
future_scenario <- terra::rast(system.file("extdata",
                                          "wc2.1_10m_bioc_ACCESS-CM2_ssp585_2081-2100.tif",
                                          package = "kuenm2"))



# Rename variables to match the variable names in the fitted models
names(future_scenario) <- sub("bio0", "bio", names(future_scenario))
names(future_scenario) <- sub("bio", "bio_", names(future_scenario))

# Run MOP
sm <- single_mop(data = fitted_model_maxnet, new_variables = future_scenario,
                 type = "detailed")
```
