% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project_selected.R
\name{project_selected}
\alias{project_selected}
\title{Project selected models for multiple scenarios}
\usage{
project_selected(models, projection_data,
                               out_dir,
                               mask = NULL,
                               consensus_per_model = TRUE,
                               consensus_general = TRUE,
                               consensus = c("median", "range", "mean", "stdev"),
                               write_replicates = FALSE, clamping = FALSE,
                               var_to_clamp = NULL, type = "cloglog",
                               overwrite = FALSE, parallel = FALSE,
                               ncores = 1, parallel_option = "doSNOW",
                               progress_bar = TRUE, verbose = TRUE)
}
\arguments{
\item{models}{an object of class \code{fitted_models} returned by the
\code{\link{fit_selected}}() function.}

\item{projection_data}{an object of class \code{prepared_proj} returned by the \code{\link{prepare_proj()}}() function. This file contains the paths to the rasters representing each scenario.}

\item{out_dir}{(character) a path to a root directory for saving the raster file of each projection.}

\item{mask}{(SpatRaster, SpatVector, or SpatExtent) spatial object used to
mask the variables before predict. Default is NULL.}

\item{consensus_per_model}{(logical) whether to calculate consensus across replicates when there are more than one replicate per model. Default is TRUE.}

\item{consensus_general}{(logical) whether to calculate consensus across models when there are more than one selected model. Default is TRUE.}

\item{consensus}{(character) consensus measures to calculate. Options available are 'median', 'range', 'mean' and 'stdev' (standard deviation). Default is c("median", "range", "mean", "stdev").}

\item{write_replicates}{(logical) whether to write the projections for each replicate. Default is FALSE.}

\item{clamping}{(logical) whether to restricts variable values to the range of the calibration data to avoid extrapolation. Default is \code{TRUE} (free extrapolation).}

\item{var_to_clamp}{(character) vector specifying which variables to clamp. Only applicable if \code{clamping = TRUE}. Default is \code{NULL}, meaning all variables will be clamped.}

\item{type}{(character) the format of the prediction values. Available options are \code{"raw"}, \code{"cumulative"}, \code{"logistic"}, and \code{"cloglog"}. Default is \code{"cloglog"}.}

\item{overwrite}{(logical) whether to overwrite SpatRaster if they already exists. Only applicable if \code{write_files} is set to TRUE. Default is FALSE.}

\item{parallel}{(logical) whether to fit the candidate models in parallel.
Default is FALSE.}

\item{ncores}{(numeric) number of cores to use for parallel processing.
Default is 1. This is only applicable if \code{parallel = TRUE}.}

\item{parallel_option}{(character) the package to use for parallel processing:
"doParallel" or "doSNOW". Default is "doSNOW". This is only applicable if
\code{parallel = TRUE}.}

\item{progress_bar}{(logical) whether to display a progress bar during processing. Default is TRUE.}

\item{verbose}{(logical) whether to display messages during processing. Default is TRUE.}
}
\value{
A \code{model_projections} object that provides the paths to the raster files with the projection results and the corresponding thresholds used to binarize the predictions.
}
\description{
This function predicts selected models across multiple scenarios, as specified in a \code{prepared_proj} object created with the \code{\link{prepare_proj()}}() function. In addition to generating predictions for each replicate, the function calculates consensus measures (e.g., mean, median) across replicates and models.
}
\examples{
# Organize and structure future climate variables from WorldClim
# Import the current variables used to fit the model.
# In this case, SoilType will be treated as a static variable (constant across future scenarios).
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))
# Create a "Current_raw" folder in a temporary directory and copy the raw variables there.
out_dir_current <- file.path(tempdir(), "Current_raw")
dir.create(out_dir_current, recursive = TRUE)
# Save current variables in temporary directory
writeRaster(var, file.path(out_dir_current, "Variables.tif"))

# Set the input directory containing the raw future climate variables.
# For this example, the data is located in the "inst/extdata" folder.
in_dir <- system.file("extdata", package = "kuenm2")
# Create a "Future_raw" folder in a temporary directory and copy the raw variables there.
out_dir_future <- file.path(tempdir(), "Future_raw")
# Organize and rename the future climate data, structuring it by year and GCM.
# The 'SoilType' variable will be appended as a static variable in each scenario.
# The files will be renamed following the "bio_" format
organize_future_worldclim(input_dir = in_dir,
                          output_dir = out_dir_future,
                          name_format = "bio_", variables = NULL,
                          fixed_variables = var$SoilType, mask = NULL,
                          overwrite = TRUE)

#Example with GLMNET
# Import example of fitted_models (output of fit_selected())
data("fitted_model_maxnet", package = "kuenm2")

# Prepare projections using fitted models to check variables
pr <- prepare_proj(models = fitted_model_maxnet,
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

#Create folder to save projection results
out_dir <- file.path(tempdir(), "Projection_results/glmnet")
dir.create(out_dir, recursive = TRUE)

# Project selected models for multiple scenarios
p <- project_selected(models = fitted_model_maxnet,
                      projection_data = pr,
                      out_dir = out_dir,
                      consensus_per_model = TRUE,
                      consensus_general = TRUE,
                      consensus = c("median", "range", "mean", "stdev"),
                      write_replicates = FALSE,
                      clamping = FALSE,
                      var_to_clamp = NULL,
                      type = "cloglog",
                      overwrite = TRUE,
                      parallel = FALSE,
                      ncores = 1,
                      parallel_option = "doSNOW",
                      progress_bar = TRUE,
                      verbose = TRUE)

#Example with GLMN
# Import example of fitted_models (output of fit_selected())
data("fitted_model_glm", package = "kuenm2")

# Prepare projections using fitted models to check variables
pr <- prepare_proj(models = fitted_model_glm,
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

#Create folder to save projection results
out_dir <- file.path(tempdir(), "Projection_results/glm")
dir.create(out_dir, recursive = TRUE)

# Project selected models for multiple scenarios
p_glm <- project_selected(models = fitted_model_glm,
                          projection_data = pr,
                          out_dir = out_dir,
                          consensus_per_model = TRUE,
                          consensus_general = TRUE,
                          consensus = c("median", "range", "mean", "stdev"),
                          write_replicates = FALSE,
                          clamping = FALSE,
                          var_to_clamp = NULL,
                          type = "cloglog",
                          overwrite = TRUE,
                          parallel = FALSE,
                          ncores = 1,
                          parallel_option = "doSNOW",
                          progress_bar = TRUE,
                          verbose = TRUE)

}
