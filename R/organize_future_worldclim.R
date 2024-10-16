#' Organize and structure future climate variables from WorldClim
#'
#' @description
#' This function imports future climate variables downloaded from WorldClim, renames the files, and organizes them into folders categorized by year and General Circulation Model (GCM). It simplifies the preparation of climate data, making it compatible with the \code{\link{prepare_proj}}() function, ensuring that all required variables are properly structured for modeling projections.
#'
#' @param input_dir (character) path to the folder containing the future climate variables downloaded from WorldClim.
#' @param output_dir (character) path to the folder where the organized data will be saved.
#' @param name_format (character) the format for renaming variable. Options are "bio_", "Bio_", "bio_0", and "Bio_0". See details for more information. Default is "bio_".
#' @param variables (character) the names of the variables to retain. Default is NULL, meaning all variables will be kept.
#' @param fixed_variables (SpatRaster) optional static variables (i.e., soil type) used in the model, which will remain unchanged in future scenarios. This variable will be included with each future scenario. Default is NULL.
#' @param check_extent (logical) whether to ensure that the `fixed_variables` have the same spatial extent as the bioclimatic variables. Applicable only if `fixed_variables` is provided. Default is TRUE.
#' @param mask (SpatRaster, SpatVector, or SpatExtent) spatial object used to mask the variables (optional). Default is NULL.
#' @param overwrite whether to overwrite existing files in the output directory. Default is FALSE.
#'
#' @return
#' A list of paths to the folders where the organized climate data has been saved.
#'
#' @importFrom terra rast crop ext writeRaster
#'
#' @details
#' The raw variables downloaded from WorldClim are named as "Bio01", "Bio02", "Bio03", "Bio10", etc. The `name_format` parameter controls how these variables will be renamed:
#' - "bio_": the variables will be renamed to bio_1, bio_2, bio_3, bio_10, etc.
#' - "bio_0": the variables will be renamed to bio_01, bio_02, bio_03, bio_10, etc
#' - "Bio_": the variables will be renamed to Bio_1, Bio_2, Bio_3, Bio_10, etc.
#' - "Bio_0": the variables will be renamed to Bio_01, Bio_02, Bio_03, Bio_10, etc.
#'
#' @export
#'
#' @usage organize_future_worldclim(input_dir, output_dir, name_format = "bio_",
#'                                variables = NULL, fixed_variables = NULL,
#'                                check_extent = TRUE, mask = NULL,
#'                                overwrite = FALSE)
#'
#' @examples
#' # Import the current variables used to fit the model.
#' # In this case, SoilType will be treated as a static variable (constant across future scenarios).
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#' # Set the input directory containing the raw future climate variables.
#' # For this example, the data is located in the "inst/extdata" folder.
#' in_dir <- "inst/extdata/"
#' # Create a "Future_raw" folder in a temporary directory and copy the raw variables there.
#' out_dir <- file.path(tempdir(), "Future_raw")
#' # Organize and rename the future climate data, structuring it by year and GCM.
#' # The 'SoilType' variable will be appended as a static variable in each scenario.
#' # The files will be renamed following the "bio_" format
#' organize_future_worldclim(input_dir = in_dir,
#'                                output_dir = out_dir,
#'                                name_format = "bio_", variables = NULL,
#'                                fixed_variables = var$SoilType, mask = NULL,
#'                                overwrite = TRUE)
#'
organize_future_worldclim <- function(input_dir, output_dir,
                                  name_format = "bio_",
                                  variables = NULL,
                                  fixed_variables = NULL,
                                  check_extent = TRUE,
                                  mask = NULL, overwrite = FALSE){
  #Check data
  if (!inherits(input_dir, "character")) {
    stop(paste0("Argument input_dir must be a character, not ",
                class(input_dir)))
  }
  if (!inherits(output_dir, "character")) {
    stop(paste0("Argument output_dir must be a character, not ",
                class(output_dir)))
  }
  if(length(name_format) != 1){
    stop("'name_format' must be a single value: 'bio_', 'Bio_', 'bio_0', or 'Bio_0'")
  }
  if (!name_format %in% c("bio_", "Bio_", "bio_0", "Bio_0")) {
    stop("'name_format' must be 'bio_', 'Bio_', 'bio_0', or 'Bio_0'")
  }

  if(!is.null(variables) & !inherits(variables, "character")){
      stop(paste0("Argument variables must be NULL or a character, not ",
                  class(variables)))
  }

  if(!is.null(fixed_variables)){
    if (!inherits(fixed_variables, "SpatRaster")) {
      stop(paste0("Argument fixed_variables must be NULL or a SpatRaster, not ",
                  class(variables)))}
    if(!inherits(check_extent, "logical")){
      stop(paste0("Argument check_extent must be logical, not ",
                  class(check_extent)))
    }}

  if (!is.null(mask) & !inherits(mask, c("SpatVector", "SpatVector", "SpatExtent"))) {
    stop(paste0("Argument mask must be a SpatVector, SpatVector, or SpatExtent, not ",
                class(mask)))
  }

  if (!is.logical(overwrite)) {
    stop(paste0("Argument overwrite must be logical, not ",
                class(overwrite)))
  }

  #Create folder
  if(!file.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }

  #List files
  lf <- list.files(input_dir, pattern = ".tif*")
  #Subset only patterns of worldclim future
  lf <- lf[grepl("wc2.1_", lf)]

  if(length(lf) == 0)
    stop("The input_dir does not contain any variable representing future scenarios from WorldClim 2.1")
  #For each file, rename, extract gcm and year and cut if necessary
  p <- lapply(lf, function(x){
    #Read file
    var_x <- terra::rast(file.path(input_dir, x))
    #Rename variables
    vnumber <- as.numeric(gsub("[^0-9]", "", names(var_x)))
    if(name_format %in% c("bio_0", "bio0", "Bio_0", "Bio0")){
      vnumber[vnumber < 10] <- paste0(0, vnumber[vnumber < 10])
    }
    names(var_x) <- paste0(name_format, vnumber)

    #If there is a mask, crop
    if(!is.null(mask)){
      var_x <- terra::crop(var_x, mask, mask = TRUE)
    }
    #Extract gcm and year
    gcm <- gsub(".*bioc_(.*?)_ssp.*", "\\1", x)
    year <- gsub(".*_(.*?)\\.tif", "\\1", x)
    ssp <-  gsub(".*_([^_]+)_.*", "\\1", x)

    #If variables is not null, subset variables
    if(!is.null(variables)){
      #Check if there are variables that do not match
      var_out <- setdiff(variables, names(var_x))
      if(length(var_out) > 0){
        stop("The variable names do not match those in the future climate data")
      }
      var_x <- var_x[[variables]]
    }

    #Append fixed variables
    if(!is.null(fixed_variables)){
      #If there is a mask, crop
      if(!is.null(mask)){
        fixed_variables <- terra::crop(fixed_variables, mask, mask = TRUE)
      }
      if(check_extent){
        if(terra::ext(fixed_variables) != terra::ext(var_x)){
          terra::ext(fixed_variables) <- terra::ext(var_x)
        }
      }
      #Append
      var_x <- c(var_x, fixed_variables)
    }

    #Create folders to save results
    res_save <- file.path(output_dir, year, ssp, gcm)
    dir.create(res_save, recursive = T, showWarnings = FALSE)
    #Save
    terra::writeRaster(var_x, file.path(res_save, "Variables.tiff"),
                overwrite = overwrite)
    return(res_save)
    })
  #Return folders
  return(p)
}

#Test function
# #Get folder with future variables (here, in the inst/extdata folder)
# future_dir <- system.file("extdata", package = "kuenm2")
# input_dir <- future_dir
# output_dir <- file.path(tempdir(), "Future_raw")
# name_format = "bio_1"
# variables = NULL
# fixed_variables = var$SoilType
# mask = NULL
# overwrite = FALSE

