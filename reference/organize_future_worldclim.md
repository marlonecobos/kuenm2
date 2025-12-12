# Organize and structure future climate variables from WorldClim

This function imports future climate variables downloaded from
WorldClim, renames the files, and organizes them into folders
categorized by year, emission scenario (SSP) and General Circulation
Model (GCM). It simplifies the preparation of climate data, making it
compatible with the
[`prepare_projection()`](https://marlonecobos.github.io/kuenm2/reference/prepare_projection.md)
function, ensuring that all required variables are properly structured
for modeling projections.

## Usage

``` r
organize_future_worldclim(input_dir, output_dir, name_format = "bio_",
                          variables = NULL, fixed_variables = NULL,
                          check_extent = TRUE, mask = NULL,
                          progress_bar = TRUE, overwrite = FALSE)
```

## Arguments

- input_dir:

  (character) path to the folder containing the future climate variables
  downloaded from WorldClim.

- output_dir:

  (character) path to the folder where the organized data will be saved.

- name_format:

  (character) the format for renaming variable. Options are "bio\_",
  "Bio\_", "bio_0", and "Bio_0". See details for more information.
  Default is "bio\_".

- variables:

  (character) the names of the variables to retain. Default is NULL,
  meaning all variables will be kept.

- fixed_variables:

  (SpatRaster) optional static variables (i.e., soil type) used in the
  model, which will remain unchanged in future scenarios. This variable
  will be included with each future scenario. Default is NULL.

- check_extent:

  (logical) whether to ensure that the `fixed_variables` have the same
  spatial extent as the bioclimatic variables. Applicable only if
  `fixed_variables` is provided. Default is TRUE.

- mask:

  (SpatRaster, SpatVector, or SpatExtent) spatial object used to mask
  the variables (optional). Default is NULL.

- progress_bar:

  (logical) whether to display a progress bar during processing. Default
  is TRUE.

- overwrite:

  whether to overwrite existing files in the output directory. Default
  is FALSE.

## Value

A list of paths to the folders where the organized climate data has been
saved.

## Details

The raw variables downloaded from WorldClim are named as "Bio01",
"Bio02", "Bio03", "Bio10", etc. The `name_format` parameter controls how
these variables will be renamed:

- "bio\_": the variables will be renamed to bio_1, bio_2, bio_3, bio_10,
  etc.

- "bio_0": the variables will be renamed to bio_01, bio_02, bio_03,
  bio_10, etc

- "Bio\_": the variables will be renamed to Bio_1, Bio_2, Bio_3, Bio_10,
  etc.

- "Bio_0": the variables will be renamed to Bio_01, Bio_02, Bio_03,
  Bio_10, etc.

## See also

[`prepare_projection()`](https://marlonecobos.github.io/kuenm2/reference/prepare_projection.md)

## Examples

``` r
# Import the current variables used to fit the model.
# In this case, SoilType will be treated as a static variable (constant
# across future scenarios).
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))

# Set the input directory containing the raw future climate variables.
# For this example, the data is located in the "inst/extdata" folder.
in_dir <- system.file("extdata", package = "kuenm2")

# Create a "Future_raw" folder in a temporary directory and copy the raw
# variables there.
out_dir <- file.path(tempdir(), "Future_raw")

# Organize and rename the future climate data, structuring it by year and GCM.
# The 'SoilType' variable will be appended as a static variable in each scenario.
# The files will be renamed following the "bio_" format
organize_future_worldclim(input_dir = in_dir, output_dir = out_dir,
                          name_format = "bio_",
                          fixed_variables = var$SoilType)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%
#> 
#> Variables successfully organized in directory:
#> /tmp/RtmpOvuyJ6/Future_raw

# Check files organized
dir(out_dir, recursive = TRUE)
#> [1] "2041-2060/ssp126/ACCESS-CM2/Variables.tif"
#> [2] "2041-2060/ssp126/MIROC6/Variables.tif"    
#> [3] "2041-2060/ssp585/ACCESS-CM2/Variables.tif"
#> [4] "2041-2060/ssp585/MIROC6/Variables.tif"    
#> [5] "2081-2100/ssp126/ACCESS-CM2/Variables.tif"
#> [6] "2081-2100/ssp126/MIROC6/Variables.tif"    
#> [7] "2081-2100/ssp585/ACCESS-CM2/Variables.tif"
#> [8] "2081-2100/ssp585/MIROC6/Variables.tif"    
```
