#' Organize and structure future climate variables from WorldClim
#'
#' @description
#' This function imports future climate variables downloaded from WorldClim,
#' renames the files, and organizes them into folders categorized by year and
#' General Circulation Model (GCM). It simplifies the preparation of climate
#' data, making it compatible with the [prepare_projection()] function,
#' ensuring that all required variables are properly structured for modeling
#' projections.
#'
#' @usage
#' organize_future_worldclim(input_dir, output_dir, name_format = "bio_",
#'                           variables = NULL, fixed_variables = NULL,
#'                           check_extent = TRUE, mask = NULL,
#'                           progress_bar = TRUE, overwrite = FALSE)
#'
#' @param input_dir (character) path to the folder containing the future climate
#' variables downloaded from WorldClim.
#' @param output_dir (character) path to the folder where the organized data
#' will be saved.
#' @param name_format (character) the format for renaming variable. Options are
#' "bio_", "Bio_", "bio_0", and "Bio_0". See details for more information.
#' Default is "bio_".
#' @param variables (character) the names of the variables to retain. Default
#' is NULL, meaning all variables will be kept.
#' @param fixed_variables (SpatRaster) optional static variables (i.e., soil
#' type) used in the model, which will remain unchanged in future scenarios.
#' This variable will be included with each future scenario. Default is NULL.
#' @param check_extent (logical) whether to ensure that the `fixed_variables`
#' have the same spatial extent as the bioclimatic variables. Applicable only
#' if `fixed_variables` is provided. Default is TRUE.
#' @param mask (SpatRaster, SpatVector, or SpatExtent) spatial object used to
#' mask the variables (optional). Default is NULL.
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is TRUE.
#' @param overwrite whether to overwrite existing files in the output directory.
#' Default is FALSE.
#'
#' @return
#' A list of paths to the folders where the organized climate data has been
#' saved.
#'
#' @importFrom terra rast crop ext writeRaster
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @details
#' The raw variables downloaded from WorldClim are named as "Bio01", "Bio02",
#' "Bio03", "Bio10", etc. The `name_format` parameter controls how these
#' variables will be renamed:
#' - "bio_": the variables will be renamed to bio_1, bio_2, bio_3, bio_10, etc.
#' - "bio_0": the variables will be renamed to bio_01, bio_02, bio_03, bio_10, etc
#' - "Bio_": the variables will be renamed to Bio_1, Bio_2, Bio_3, Bio_10, etc.
#' - "Bio_0": the variables will be renamed to Bio_01, Bio_02, Bio_03, Bio_10, etc.
#'
#' @export
#'
#' @seealso
#' [prepare_projection()]
#'
#' @examples
#' # Import the current variables used to fit the model.
#' # In this case, SoilType will be treated as a static variable (constant
#' # across future scenarios).
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' # Set the input directory containing the raw future climate variables.
#' # For this example, the data is located in the "inst/extdata" folder.
#' in_dir <- system.file("extdata", package = "kuenm2")
#'
#' # Create a "Future_raw" folder in a temporary directory and copy the raw
#' # variables there.
#' out_dir <- file.path(tempdir(), "Future_raw")
#'
#' # Organize and rename the future climate data, structuring it by year and GCM.
#' # The 'SoilType' variable will be appended as a static variable in each scenario.
#' # The files will be renamed following the "bio_" format
#' organize_future_worldclim(input_dir = in_dir, output_dir = out_dir,
#'                           name_format = "bio_",
#'                           fixed_variables = var$SoilType)
#'
#' # Check files organized
#' dir(out_dir, recursive = TRUE)

organize_future_worldclim <- function(input_dir,
                                      output_dir,
                                      name_format = "bio_",
                                      variables = NULL,
                                      fixed_variables = NULL,
                                      check_extent = TRUE,
                                      mask = NULL,
                                      progress_bar = TRUE,
                                      overwrite = FALSE) {
  #Check data
  if (missing(input_dir)) {
    stop("Argument 'input_dir' must be defined.")
  }
  if (missing(output_dir)) {
    stop("Argument 'output_dir' must be defined.")
  }
  if (!inherits(input_dir, "character")) {
    stop("Argument input_dir must be a character.")
  }
  if (!inherits(output_dir, "character")) {
    stop("Argument output_dir must be a character.")
  }
  if (length(name_format) != 1) {
    stop("'name_format' must be a single value: 'bio_', 'Bio_', 'bio_0', or 'Bio_0'.")
  }
  if (!name_format %in% c("bio_", "Bio_", "bio_0", "Bio_0", "bio", "Bio")) {
    stop("'name_format' must be 'bio_', 'Bio_', 'bio_0', ,'Bio_0', 'bio', or 'Bio'.")
  }

  if (!is.null(variables) & !inherits(variables, "character")) {
    stop("Argument 'variables' must be NULL or a 'character'.")
  }

  if (!is.null(fixed_variables)) {
    if (!inherits(fixed_variables, "SpatRaster")) {
      stop("Argument 'fixed_variables' must be NULL or a 'SpatRaster'.")
    }
  }

  if (!is.null(mask) & !inherits(mask, c("SpatVector", "SpatVector", "SpatExtent"))) {
    stop("Argument 'mask' must be a 'SpatVector', 'SpatVector', or 'SpatExtent'.")
  }

  #Create folder
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  #List files
  lf <- list.files(input_dir, pattern = ".tif*")
  #Subset only patterns of worldclim future
  lf <- lf[grepl("wc2.1_", lf)]

  if (length(lf) == 0)
    stop("'input_dir' does not contain any variable from WorldClim 2.1 future scenarios.")

  if (progress_bar) {
    pb <- utils::txtProgressBar(0, length(lf), style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n) }

  #For each file, rename, extract gcm and year and cut if necessary
  p <- list()

  for(x in 1:length(lf)) {
    #Read file
    var_x <- terra::rast(file.path(input_dir, lf[x]))
    #Rename variables
    #vnumber <- as.numeric(gsub("[^0-9]", "", names(var_x)))
    #or
    vnumber <- as.integer(sub(".*?(\\d+)$", "\\1", names(var_x))) #Work with geodata
    if (name_format %in% c("bio_0", "bio0", "Bio_0", "Bio0")) {
      vnumber[vnumber < 10] <- paste0(0, vnumber[vnumber < 10])
    }
    names(var_x) <- paste0(name_format, vnumber)

    #If there is a mask, crop
    if (!is.null(mask)) {
      var_x <- terra::crop(var_x, mask, mask = TRUE)
    }
    #Extract gcm and year
    gcm <- gsub(".*bioc_(.*?)_ssp.*", "\\1", lf[x])
    year <- gsub(".*_(.*?)\\.tif", "\\1", lf[x])
    ssp <-  gsub(".*_([^_]+)_.*", "\\1", lf[x])

    #If variables is not null, subset variables
    if (!is.null(variables)) {
      #Check if there are variables that do not match
      var_out <- setdiff(variables, names(var_x))
      if (length(var_out) > 0) {
        stop("Variable names do not match those in future climate data.")
      }
      var_x <- var_x[[variables]]
    }

    #Append fixed variables
    if (!is.null(fixed_variables)) {
      #If there is a mask, crop
      if (!is.null(mask)) {
        fixed_variables <- terra::crop(fixed_variables, mask, mask = TRUE)
      }
      if (check_extent) {
        if (terra::ext(fixed_variables) != terra::ext(var_x)) {
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

    # Sets the progress bar to the current state
    if (progress_bar) {
      utils::setTxtProgressBar(pb, x) }

    p[[x]] <- res_save
  }

  #Return folders
  return(cat("\nVariables successfully organized in the root directory:\n",
             output_dir))
}


