# Preparation of data for model projections

This function prepared data for model projections to multiple scenarios,
storing the paths to the rasters representing each scenario.

## Usage

``` r
prepare_projection(models = NULL, variable_names = NULL, present_dir = NULL,
                   past_dir = NULL, past_period = NULL, past_gcm = NULL,
                   future_dir = NULL, future_period = NULL,
                   future_pscen = NULL, future_gcm = NULL,
                   write_file = FALSE, filename = NULL,
                   raster_pattern = ".tif*")
```

## Arguments

- models:

  an object of class `fitted_models` returned by the
  [`fit_selected`](https://marlonecobos.github.io/kuenm2/reference/fit_selected.md)()
  function. Default is NULL.

- variable_names:

  (character) names of the variables used to fit the model or do the PCA
  in the
  [`prepare_data`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)()
  function. Only applicable if `models` argument is not provided.
  Default is NULL.

- present_dir:

  (character) path to the folder containing variables that represent the
  current scenario for projection. Default is NULL.

- past_dir:

  (character) path to the folder containing subfolders with v ariables
  representing past scenarios for projection. Default is NULL.

- past_period:

  (character) names of the subfolders within `past_dir`, representing
  specific time periods (e.g., 'LGM' or 'MID').

- past_gcm:

  (character) names of the subfolders within `past_period` folders,
  representing specific General Circulation Models (GCMs).

- future_dir:

  (character) path to the folder containing subfolders with variables
  representing future scenarios for projection. Default is NULL.

- future_period:

  (character) names of the subfolders within `future_dir`, representing
  specific time periods (e.g., '2041-2060' or '2081-2100'). Default is
  NULL.

- future_pscen:

  (character) names of the subfolders within `future_period`,
  representing specific emission scenarios (e.g., 'ssp126' or 'ssp585').
  Default is NULL.

- future_gcm:

  (character) names of the subfolders within `future_pscen` folders,
  representing specific General Circulation Models (GCMs). Default is
  NULL.

- write_file:

  (logical) whether to write the object containing the paths to the
  structured folders. This object is required for projecting models
  across multiple scenarios using the
  [`project_selected`](https://marlonecobos.github.io/kuenm2/reference/project_selected.md)()
  function. Default is FALSE.

- filename:

  (character) the path or name of the folder where the object will be
  saved. This is only applicable if `write_file = TRUE`. Default is
  NULL.

- raster_pattern:

  (character) pattern used to identify the format of raster files within
  the folders. Default is ".tif\*".

## Value

An object of class `prepared_projection` containing the following
elements:

- Present, Past, and Future: paths to the variables structured in
  subfolders.

- Raster_pattern: the pattern used to identify the format of raster
  files within the folders.

- PCA: if a principal component analysis (PCA) was performed on the set
  of variables with
  [`prepare_data`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)(),
  a list with class "prcomp" will be returned. See `?stats::prcomp()`
  for details.

- variables: names of the raw predictos variables used to project.

## See also

[`organize_future_worldclim()`](https://marlonecobos.github.io/kuenm2/reference/organize_future_worldclim.md)

## Examples

``` r
# Import example of fitted_models (output of fit_selected())
data("fitted_model_maxnet", package = "kuenm2")

# Organize and structure future climate variables from WorldClim
# Import the current variables used to fit the model.
# In this case, SoilType will be treated as a static variable (constant
# across future scenarios).
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))

# Create a "Current_raw" folder in a temporary directory and copy the raw
# variables there.
out_dir_current <- file.path(tempdir(), "Current_raw")
dir.create(out_dir_current, recursive = TRUE)

# Save current variables in temporary directory
terra::writeRaster(var, file.path(out_dir_current, "Variables.tif"))

# Set the input directory containing the raw future climate variables.
# For this example, the data is located in the "inst/extdata" folder.
in_dir <- system.file("extdata", package = "kuenm2")

# Create a "Future_raw" folder in a temporary directory and copy the raw
# variables there.
out_dir_future <- file.path(tempdir(), "Future_raw")

# Organize and rename the future climate data, structuring it by year and GCM.
# The 'SoilType' variable will be appended as a static variable in each scenario.
# The files will be renamed following the "bio_" format
organize_future_worldclim(input_dir = in_dir,
                          output_dir = out_dir_future,
                          name_format = "bio_", variables = NULL,
                          fixed_variables = var$SoilType, mask = NULL,
                          overwrite = TRUE)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%
#> 
#> Variables successfully organized in directory:
#> /tmp/RtmpuJ37XZ/Future_raw

# Prepare projections using fitted models to check variables
pr <- prepare_projection(models = fitted_model_maxnet,
                         present_dir = out_dir_current,
                         past_dir = NULL,
                         past_period = NULL,
                         past_gcm = NULL,
                         future_dir = out_dir_future,
                         future_period = c("2041-2060", "2081-2100"),
                         future_pscen = c("ssp126", "ssp585"),
                         future_gcm = c("ACCESS-CM2", "MIROC6"),
                         write_file = FALSE,
                         filename = NULL,
                         raster_pattern = ".tif*")
pr
#> projection_data object summary
#> =============================
#> Variables prepared to project models for Present and Future 
#> Future projections contain the following periods, scenarios and GCMs:
#>   - Periods: 2041-2060 | 2081-2100 
#>   - Scenarios: ssp126 | ssp585 
#>   - GCMs: ACCESS-CM2 | MIROC6 
#> All variables are located in the following root directory:
#> /tmp/RtmpuJ37XZ

# Prepare projections using variables names
pr_b <- prepare_projection(models = NULL,
                           variable_names = c("bio_1", "bio_7", "bio_12"),
                           present_dir = out_dir_current,
                           past_dir = NULL,
                           past_period = NULL,
                           past_gcm = NULL,
                           future_dir = out_dir_future,
                           future_period = c("2041-2060", "2081-2100"),
                           future_pscen = c("ssp126", "ssp585"),
                           future_gcm = c("ACCESS-CM2", "MIROC6"),
                           write_file = FALSE,
                           filename = NULL,
                           raster_pattern = ".tif*")
pr_b
#> projection_data object summary
#> =============================
#> Variables prepared to project models for Present and Future 
#> Future projections contain the following periods, scenarios and GCMs:
#>   - Periods: 2041-2060 | 2081-2100 
#>   - Scenarios: ssp126 | ssp585 
#>   - GCMs: ACCESS-CM2 | MIROC6 
#> All variables are located in the following root directory:
#> /tmp/RtmpuJ37XZ
```
