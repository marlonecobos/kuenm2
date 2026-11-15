#' Organize and structure future climate variables from WorldClim
#'
#' @description
#' This function imports future climate variables downloaded from WorldClim,
#' renames the files, and organizes them into folders categorized by year,
#' emission scenario (SSP) and General Circulation Model (GCM). It simplifies
#' the preparation of climate data, making it compatible with the
#' [prepare_projection()] function, ensuring that all required variables are
#' properly structured for modeling projections.
#'
#' @usage
#' organize_future_worldclim(input_dir, output_dir, name_format = "bio_",
#'                           variables = NULL, static_variables = NULL,
#'                           check_extent = TRUE, mask = NULL,
#'                           progress_bar = TRUE, progress_type = "timer",
#'                           overwrite = FALSE)
#'
#' @param input_dir (character) path to the folder containing the future climate
#' variables downloaded from WorldClim.
#' @param output_dir (character) path to the folder where the organized data
#' will be saved.
#' @param name_format (character) the format for renaming variables. Options are
#' `"bio_"`, `"Bio_"`, `"bio_0"`, `"Bio_0"`, `"bio"`, and `"Bio"`. See details
#' for more information. Default is `"bio_"`.
#' @param variables (character) the names of the variables to retain. Default
#' is NULL, meaning all variables will be kept.
#' @param static_variables (SpatRaster) optional static variables (i.e., soil
#' type) used in the model, which will remain unchanged in future scenarios.
#' This variable will be included with each future scenario. Default is NULL.
#' @param check_extent (logical) whether to ensure that the `static_variables`
#' have the same spatial extent as the bioclimatic variables. Applicable only
#' if `static_variables` is provided. Default is TRUE.
#' @param mask (SpatRaster, SpatVector, or SpatExtent) spatial object used to
#' mask the variables (optional). Default is NULL.
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is TRUE.
#' @param progress_type (character) type of the progress bar. See
#' `pbapply::pbtypes()` for available progress bar types depending on operating
#' system. Default is `"timer"`.
#' @param overwrite (logical) whether to overwrite existing files in the output
#' directory. Default is FALSE.
#'
#' @return
#' Invisibly returns a character vector of paths to the folders where the
#' organized climate data has been saved.
#'
#' @importFrom terra rast crop ext writeRaster
#' @importFrom pbapply pboptions pblapply
#'
#' @details
#' The raw variables downloaded from WorldClim are named as "Bio01", "Bio02",
#' "Bio03", "Bio10", etc. The `name_format` parameter controls how these
#' variables will be renamed:
#' - `"bio_"`: the variables will be renamed to bio_1, bio_2, bio_3, bio_10, etc.
#' - `"bio_0"`: the variables will be renamed to bio_01, bio_02, bio_03, bio_10, etc.
#' - `"Bio_"`: the variables will be renamed to Bio_1, Bio_2, Bio_3, Bio_10, etc.
#' - `"Bio_0"`: the variables will be renamed to Bio_01, Bio_02, Bio_03, Bio_10, etc.
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
#'                           static_variables = var$SoilType)
#'
#' # Check files organized
#' dir(out_dir, recursive = TRUE)
organize_future_worldclim <- function(input_dir,
                                      output_dir,
                                      name_format = "bio_",
                                      variables = NULL,
                                      static_variables = NULL,
                                      check_extent = TRUE,
                                      mask = NULL,
                                      progress_bar = TRUE,
                                      progress_type = "timer",
                                      overwrite = FALSE) {
  # Check data
  if (missing(input_dir)) {
    stop("Argument 'input_dir' must be defined.")
  }
  if (missing(output_dir)) {
    stop("Argument 'output_dir' must be defined.")
  }
  if (!inherits(input_dir, "character")) {
    stop("Argument 'input_dir' must be a character.")
  }
  if (!inherits(output_dir, "character")) {
    stop("Argument 'output_dir' must be a character.")
  }
  if (length(name_format) != 1L) {
    stop("'name_format' must be a single value: 'bio_', 'Bio_', 'bio_0', 'Bio_0', 'bio', or 'Bio'.")
  }
  if (!name_format %in% c("bio_", "Bio_", "bio_0", "Bio_0", "bio", "Bio")) {
    stop("'name_format' must be 'bio_', 'Bio_', 'bio_0', 'Bio_0', 'bio', or 'Bio'.")
  }
  if (!is.null(variables) && !inherits(variables, "character")) {
    stop("Argument 'variables' must be NULL or a 'character'.")
  }
  if (!is.null(static_variables) && !inherits(static_variables, "SpatRaster")) {
    stop("Argument 'static_variables' must be NULL or a 'SpatRaster'.")
  }
  if (!is.null(mask) && !inherits(mask, c("SpatRaster", "SpatVector", "SpatExtent"))) {
    stop("Argument 'mask' must be a 'SpatRaster', 'SpatVector', or 'SpatExtent'.")
  }

  # Create folder
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # List files
  lf <- list.files(input_dir, pattern = "\\.tif{1,2}$", ignore.case = TRUE)
  lf <- lf[grepl("wc2\\.1_", lf)]

  if (length(lf) == 0L) {
    stop("'input_dir' does not contain any variable from WorldClim 2.1 future scenarios.")
  }

  # Pre-crop static_variables once outside the loop
  if (!is.null(static_variables) && !is.null(mask)) {
    static_variables <- terra::crop(static_variables, mask, mask = TRUE)
  }

  # Clean prefix when zero-padded format is requested (avoids "bio_001")
  use_zero_pad <- name_format %in% c("bio_0", "bio0", "Bio_0", "Bio0")
  base_prefix <- sub("0$", "", name_format)

  # Setup progress bar options safely for CRAN
  old_pbopt <- pbapply::pboptions(
    type = if (progress_bar) progress_type else "none"
  )
  on.exit(pbapply::pboptions(old_pbopt), add = TRUE)

  # Process each file with pblapply
  saved_paths <- pbapply::pblapply(
    X = lf,
    FUN = function(file_x) {
      # Read file
      var_x <- terra::rast(file.path(input_dir, file_x))

      # Rename variables
      vnumber <- as.integer(sub(".*?(\\d+)$", "\\1", names(var_x)))
      if (use_zero_pad) {
        vnumber <- sprintf("%02d", vnumber)
      }
      names(var_x) <- paste0(base_prefix, vnumber)

      # Subset variables before cropping to avoid unnecessary spatial operations
      if (!is.null(variables)) {
        var_out <- setdiff(variables, names(var_x))
        if (length(var_out) > 0L) {
          stop("Variable names do not match those in future climate data.")
        }
        var_x <- var_x[[variables]]
      }

      # Crop bioclimatic variables if mask is provided
      if (!is.null(mask)) {
        var_x <- terra::crop(var_x, mask, mask = TRUE)
      }

      # Extract GCM, SSP, and year safely
      gcm <- sub(".*bioc_(.*?)_ssp.*", "\\1", file_x)
      ssp <- sub(".*_(ssp[^_]+)_.*", "\\1", file_x)
      year <- sub(".*_([^_]+)\\.tif{1,2}$", "\\1", file_x, ignore.case = TRUE)

      # Append static variables
      if (!is.null(static_variables)) {
        sv_x <- static_variables
        if (check_extent && terra::ext(sv_x) != terra::ext(var_x)) {
          terra::ext(sv_x) <- terra::ext(var_x)
        }
        var_x <- c(var_x, sv_x)
      }

      # Create folders and save results
      res_save <- file.path(output_dir, year, ssp, gcm)
      if (!dir.exists(res_save)) {
        dir.create(res_save, recursive = TRUE, showWarnings = FALSE)
      }
      terra::writeRaster(var_x, file.path(res_save, "Variables.tif"),
                         overwrite = overwrite)

      return(res_save)
    }
  )

  message("\nVariables successfully organized in directory:\n", output_dir)
  return(invisible(unlist(saved_paths, use.names = FALSE)))
}
