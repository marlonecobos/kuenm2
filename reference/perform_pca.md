# Principal Component Analysis for raster layers

This function performs principal component analysis (PCA) with a set of
raster variables.

## Usage

``` r
perform_pca(raster_variables, exclude_from_pca = NULL, project = FALSE,
            projection_data = NULL, out_dir = NULL, overwrite = FALSE,
            progress_bar = FALSE, center = TRUE, scale = FALSE,
            variance_explained = 95, min_explained = 5)
```

## Arguments

- raster_variables:

  (SpatRaster) set of predictor variables that the function will
  summarize into a set of orthogonal, uncorrelated components based on
  PCA.

- exclude_from_pca:

  (character) variable names within raster_variables that should not be
  included in the PCA transformation. Instead, these variables will be
  added directly to the final set of output variables without being
  modified. The default is NULL, meaning all variables will be used
  unless specified otherwise.

- project:

  (logical) whether the function should project new data from different
  scenarios (e.g. future variables) onto the PCA coordinates generated
  by the initial analysis. If TRUE, the argument projection_data needs
  to be defined. Default is FALSE.

- projection_data:

  an object of class `prepared_projection` returned by the
  [`prepare_projection()`](https://marlonecobos.github.io/kuenm2/reference/prepare_projection.md)
  function. This file contains the paths to the raster files
  representing each scenario. Only applicable if `project = TRUE`.
  Default is NULL.

- out_dir:

  (character) a path to a root directory for saving the raster files of
  each projection. Default = NULL.

- overwrite:

  (logical) whether to overwrite SpatRaster if they already exists when
  projecting. Only applicable if `write_files` is set to TRUE. Default
  is FALSE.

- progress_bar:

  (logical) whether to display a progress bar during processing
  projections. Only applicable if `project = TRUE`. Default is FALSE

- center:

  (logical) whether the variables should be zero-centered. Default is
  TRUE.

- scale:

  (logical) whether the variables should be scaled to have unit variance
  before the analysis takes place. Default is FALSE.

- variance_explained:

  (numeric) the cumulative percentage of total variance that must be
  explained by the selected principal components. Default is 95.

- min_explained:

  (numeric) the minimum percentage of total variance that a principal
  component must explain to be retained. Default is 5.

## Value

A list containing the following elements:

- env: A SpatRaster object that contains the orthogonal components
  derived from the PCA. PCs correspond to the variables used to perform
  the analysis.

- pca: an object of class prcomp, containing the details of the PCA
  analysis. See [`prcomp`](https://rdrr.io/r/stats/prcomp.html)().

- variance_explained_cum_sum: The cumulative percentage of total
  variance explained by each of the selected principal components. This
  value indicates how much of the data’s original variability is
  captured by the PCA transformation.

- projection_directory: the root directory where projection files were
  saved. Not NULL only if `project` was set to TRUE. This directory
  contains the projected raster files for each scenario.

## Examples

``` r
# PCA with current variables
# Import raster layers
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))

# PCA
pca_var <- perform_pca(raster_variables = var, exclude_from_pca = "SoilType",
                       center = TRUE, scale = TRUE)

pca_var
#> $env
#> class       : SpatRaster 
#> size        : 52, 40, 5  (nrow, ncol, nlyr)
#> resolution  : 0.1666667, 0.1666667  (x, y)
#> extent      : -53.5, -46.83333, -30.83333, -22.16667  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> sources     : memory  (4 layers) 
#>               Current_variables.tif  
#> varnames    : Current_variables 
#>               Current_variables 
#> names       :       PC1,       PC2,       PC3,       PC4, SoilType 
#> min values  : -3.621362, -2.041276, -3.923471, -1.730859,        1 
#> max values  :  2.929786,  3.029667,  1.752452,  1.601162,       23 
#> 
#> $pca
#> Standard deviations (1, .., p=4):
#> [1] 1.4574175 0.9290330 0.8020786 0.6078666
#> 
#> Rotation (n x k) = (4 x 4):
#>               PC1        PC2         PC3         PC4
#> bio_1  -0.5327531 -0.1500983 -0.66580466  0.50034864
#> bio_7   0.3338164 -0.9359166 -0.09775690 -0.05541088
#> bio_12  0.5032916  0.2775767 -0.73525490 -0.35923383
#> bio_15 -0.5928223 -0.1564666 -0.08091963 -0.78583200
#> 
#> $variance_explained_cumsum
#>       PC1       PC2       PC3       PC4 
#>  53.10165  74.67920  90.76246 100.00000 
#> 
#> $projection_directory
#> NULL
#> 

# Project PCA for new scenarios (future)
# First, organize and prepare future variables
# Set the input directory containing the raw future climate variables
# For this example, the data is located in the "inst/extdata" folder.
in_dir <- system.file("extdata", package = "kuenm2")

# Create a "Future_raw" folder in a temporary directory and copy the variables.
out_dir_future <- file.path(tempdir(), "Future_raw1")

# Organize and rename the future climate data, structuring it by year and GCM.
# The 'SoilType' variable will be appended as a static variable in each scenario.
# The files will be renamed following the "bio_" format
organize_future_worldclim(input_dir = in_dir, output_dir = out_dir_future,
                          name_format = "bio_", fixed_variables = var$SoilType)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%
#> 
#> Variables successfully organized in directory:
#> /tmp/RtmpWvlY6q/Future_raw1

# Prepare projections
pr <- prepare_projection(variable_names = c("bio_1", "bio_7", "bio_12",
                                            "bio_15", "SoilType"),
                         future_dir = out_dir_future,
                         future_period = c("2041-2060", "2081-2100"),
                         future_pscen = c("ssp126", "ssp585"),
                         future_gcm = c("ACCESS-CM2", "MIROC6"),
                         raster_pattern = ".tif*")

# Create folder to save projection results
out_dir <- file.path(tempdir(), "PCA_projections")
dir.create(out_dir, recursive = TRUE)

# Perform and project PCA for new scenarios (future)
proj_pca <- perform_pca(raster_variables = var, exclude_from_pca = "SoilType",
                        project = TRUE, projection_data = pr,
                        out_dir = out_dir, center = TRUE, scale = TRUE)

proj_pca$projection_directory  # Directory with projected PCA-variables
#> [1] "/tmp/RtmpWvlY6q/PCA_projections"
```
