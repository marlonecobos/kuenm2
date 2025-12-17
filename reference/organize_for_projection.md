# Organize and structure variables for past and future projections

This function helpts to organize climate variable files from past and
future scenarios into folders categorized by time period ("Past" or
"Future"), specific period (e.g., "LGM" or "2081–2100"), emission
scenario (e.g., "ssp585"), and GCMs. This structure simplifies the
preparation of climate data and ensures compatibility with the
[`prepare_projection()`](https://marlonecobos.github.io/kuenm2/reference/prepare_projection.md)
function, making the variables properly organized for modeling
projections. See **Details** for more information.

## Usage

``` r
organize_for_projection(output_dir, models = NULL,
                               variable_names = NULL,
                               categorical_variables = NULL,
                               present_file = NULL,
                               past_files = NULL, past_period = NULL,
                               past_gcm = NULL, future_files = NULL,
                               future_period = NULL, future_pscen = NULL,
                               future_gcm = NULL, fixed_variables = NULL,
                               check_extent = TRUE,
                               resample_to_present = TRUE, mask = NULL,
                               overwrite = FALSE)
```

## Arguments

- output_dir:

  (character) path to the folder where the organized data will be saved.

- models:

  an object of class fitted_models returned by the
  [`fit_selected()`](https://marlonecobos.github.io/kuenm2/reference/fit_selected.md)
  function. Default is NULL.

- variable_names:

  (character) names of the variables used to fit the model or do the PCA
  in the
  [`prepare_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)
  function. Only applicable if 'models' argument is not provided.
  Default is NULL.

- categorical_variables:

  (character) names of the variables that are categorical. Default is
  NULL.

- present_file:

  (character) **full paths** to the variables from the present scenario.
  Default is NULL.

- past_files:

  (character) **full paths** to the variables from the past scenario(s).
  Default is NULL.

- past_period:

  (character) names of the subfolders within 'past_files', representing
  specific time periods (e.g., 'LGM' or 'MID'). Only applicable if
  'past_files' is provided. Default is NULL.

- past_gcm:

  (character) names of the subfolders within 'past_files', representing
  specific General Circulation Models (GCMs). Only applicable if
  'past_files' is provided. Default is NULL.

- future_files:

  (character) **full paths** to the variables from the future
  scenario(s). Default is NULL.

- future_period:

  (character) names of the subfolders within 'future_files',
  representing specific time periods (e.g., '2041-2060' or '2081-2100').
  Only applicable if 'future_files' is provided. Default is NULL.

- future_pscen:

  (character) names of the subfolders within 'future_files',
  representing specific emission scenarios (e.g., 'ssp126' or 'ssp585').
  Only applicable if 'future_files' is provided. Default is NULL.

- future_gcm:

  (character) names of the subfolders within 'future_files',
  representing specific General Circulation Models (GCMs). Only
  applicable if 'future_files' is provided. Default is NULL.

- fixed_variables:

  (SpatRaster) optional static variables (i.e., soil type) used in the
  model, which will remain unchanged in past or future scenarios. This
  variable will be included with each scenario. Default is NULL.

- check_extent:

  (logical) whether to ensure that the 'fixed_variables' have the same
  spatial extent as the bioclimatic variables. Applicable only if
  'fixed_variables' is provided. Default is TRUE.

- resample_to_present:

  (logical) whether to resample past or future variables so they match
  the extent of the present variables. Only used when 'present_file' is
  provided. Default is TRUE.

- mask:

  (SpatRaster, SpatVector, or SpatExtent) spatial object used to mask
  the variables (optional). Default is NULL.

- overwrite:

  whether to overwrite existing files in the output directory. Default
  is FALSE.

## Value

A message indicating that the variables were successfully organized in
the 'output_dir' directory.

## Details

The listed input rasters must be stored as `.tif` files, with one file
per scenario. Filenames should include identifiable patterns for time
period, GCM, and (for future scenarios) the emission scenario (SSP).

For example:

- A file representing "Past" conditions for the "LGM" period using the
  "MIROC6" GCM should be named: `"Past_LGM_MIROC6.tif"`

- A file representing "Future" conditions for the period "2081–2100"
  under the emission scenario "ssp585" and the GCM "ACCESS-CM2" should
  be named: `"Future_2081-2100_ssp585_ACCESS-CM2.tif"`

All scenario files must contain the same variable names (e.g., `bio1`,
`bio2`, etc.) and units as those used for model calibration with
present-day data. Tip: When listing the files, use
`list.files(path, full.names = TRUE)` to obtain the full file paths
required by the function.

## See also

[prepare_projection](https://marlonecobos.github.io/kuenm2/reference/prepare_projection.md)
[organize_future_worldclim](https://marlonecobos.github.io/kuenm2/reference/organize_future_worldclim.md)

## Examples

``` r
# Set the input directory containing the climate variables.
# In this example, we use present and LGM variables from CHELSA
# located in the "inst/extdata" folder of the package.
present_lgm_dir <- system.file("extdata", package = "kuenm2")

# Define an output directory (here, using a temporary folder)
# Replace with your own working directory if needed.
out_dir <- file.path(tempdir(), "Projection_variables")

# List files for present-day conditions
present_list <- list.files(path = present_lgm_dir,
                           pattern = "Current_CHELSA", # Select only CHELSA present-day files
                           full.names = TRUE)

# List files for LGM conditions
lgm_list <- list.files(path = present_lgm_dir,
                       pattern = "LGM", # Select only LGM files
                       full.names = TRUE)

# Organize variables for projection
organize_for_projection(output_dir = out_dir,
                        variable_names = c("bio1", "bio7", "bio12", "bio15"),
                        present_file = present_list,
                        past_files = lgm_list,
                        past_period = "LGM",
                        past_gcm = c("CCSM4", "CNRM-CM5", "FGOALS-g2",
                                     "IPSL-CM5A-LR", "MIROC-ESM", "MPI-ESM-P",
                                     "MRI-CGCM3"),
                        resample_to_present = TRUE,
                        overwrite = TRUE)
#> 
#> Variables successfully organized in directory:
#> /tmp/RtmppvPMA9/Projection_variables
```
