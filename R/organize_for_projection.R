#' Organize and structure variables for past and future projections
#'
#' @description
#' This function helpts to organize climate variable files from past and future
#' scenarios into folders categorized by time period ("Past" or "Future"),
#' specific period (e.g., "LGM" or "2081–2100"), emission scenario (e.g.,
#' "ssp585"), and GCMs. This structure simplifies the preparation of climate
#' data and ensures compatibility with the `prepare_projection()` function,
#' making the variables properly organized for modeling projections.
#' See **Details** for more information.
#'
#' @usage organize_for_projection(output_dir, models = NULL,
#'                                variable_names = NULL, present_file = NULL,
#'                                past_files = NULL, past_period = NULL,
#'                                past_gcm = NULL, future_files = NULL,
#'                                future_period = NULL, future_pscen = NULL,
#'                                future_gcm = NULL, fixed_variables = NULL,
#'                                check_extent = TRUE,
#'                                resample_to_present = TRUE, mask = NULL,
#'                                overwrite = FALSE)
#'
#' @param output_dir (character) path to the folder where the organized data
#' will be saved.
#' @param models an object of class fitted_models returned by the
#' `fit_selected()` function. Default is NULL.
#' @param variable_names (character) names of the variables used to fit the
#' model or do the PCA in the `prepare_data()` function. Only applicable if
#' 'models' argument is not provided. Default is NULL.
#' @param present_file (character) **full paths** to the variables from the
#' present scenario. Default is NULL.
#' @param past_files (character) **full paths** to the variables from the past
#'  scenario(s). Default is NULL.
#' @param past_period (character) names of the subfolders within 'past_files',
#'  representing specific time periods (e.g., 'LGM' or 'MID'). Only applicable
#'  if 'past_files' is provided. Default is NULL.
#' @param past_gcm (character) names of the subfolders within 'past_files',
#'  representing specific General Circulation Models (GCMs). Only applicable if
#'  'past_files' is provided. Default is NULL.
#' @param future_files (character) **full paths** to the variables from the
#'  future scenario(s). Default is NULL.
#' @param future_period (character) names of the subfolders within
#' 'future_files', representing specific time periods (e.g., '2041-2060' or
#' '2081-2100'). Only applicable if 'future_files' is provided. Default is NULL.
#' @param future_pscen (character) names of the subfolders within 'future_files',
#' representing specific emission scenarios (e.g., 'ssp126' or 'ssp585'). Only
#' applicable if 'future_files' is provided. Default is NULL.
#' @param future_gcm (character) names of the subfolders within 'future_files',
#' representing specific General Circulation Models (GCMs). Only applicable if
#' 'future_files' is provided. Default is NULL.
#' @param fixed_variables (SpatRaster) optional static variables (i.e., soil
#' type) used in the model, which will remain unchanged in past or future
#' scenarios. This variable will be included with each scenario. Default is NULL.
#' @param check_extent (logical) whether to ensure that the 'fixed_variables'
#' have the same spatial extent as the bioclimatic variables. Applicable only if
#' 'fixed_variables' is provided. Default is TRUE.
#' @param resample_to_present (logical) whether to resample past or future
#' variables so they match the extent of the present variables. Only used when
#' 'present_file' is provided. Default is TRUE.
#' @param mask (SpatRaster, SpatVector, or SpatExtent) spatial object used to
#' mask the variables (optional). Default is NULL.
#' @param overwrite whether to overwrite existing files in the output directory.
#' Default is FALSE.
#'
#' @details
#' The listed input rasters must be stored as `.tif` files, with one file per
#' scenario. Filenames should include identifiable patterns for time period,
#' GCM, and (for future scenarios) the emission scenario (SSP).
#'
#' For example:
#' - A file representing "Past" conditions for the "LGM" period using the
#' "MIROC6" GCM should be named: `"Past_LGM_MIROC6.tif"`
#' - A file representing "Future" conditions for the period "2081–2100" under
#' the emission scenario "ssp585" and the GCM "ACCESS-CM2" should be named:
#' `"Future_2081-2100_ssp585_ACCESS-CM2.tif"`
#'
#' All scenario files must contain the same variable names (e.g., `bio1`,
#' `bio2`, etc.) and units as those used for model calibration with present-day
#' data.
#' Tip: When listing the files, use `list.files(path, full.names = TRUE)` to
#' obtain the full file paths required by the function.
#'
#' @return
#' A message indicating that the variables were successfully organized in the
#' 'output_dir' directory.
#'
#' @export
#'
#' @importFrom terra rast writeRaster crop res ext resample
#'
#' @seealso [prepare_projection] [organize_future_worldclim]
#'
#' @examples
#' # Set the input directory containing the climate variables.
#' # In this example, we use present and LGM variables from CHELSA
#' # located in the "inst/extdata" folder of the package.
#' present_lgm_dir <- system.file("extdata", package = "kuenm2")
#'
#' # Define an output directory (here, using a temporary folder)
#' # Replace with your own working directory if needed.
#' out_dir <- file.path(tempdir(), "Projection_variables")
#'
#' # List files for present-day conditions
#' present_list <- list.files(path = present_lgm_dir,
#'                            pattern = "Current_CHELSA", # Select only CHELSA present-day files
#'                            full.names = TRUE)
#'
#' # List files for LGM conditions
#' lgm_list <- list.files(path = present_lgm_dir,
#'                        pattern = "LGM", # Select only LGM files
#'                        full.names = TRUE)
#'
#' # Organize variables for projection
#' organize_for_projection(output_dir = out_dir,
#'                         variable_names = c("bio1", "bio7", "bio12", "bio15"),
#'                         present_file = present_list,
#'                         past_files = lgm_list,
#'                         past_period = "LGM",
#'                         past_gcm = c("CCSM4", "CNRM-CM5", "FGOALS-g2",
#'                                      "IPSL-CM5A-LR", "MIROC-ESM", "MPI-ESM-P",
#'                                      "MRI-CGCM3"),
#'                         resample_to_present = TRUE,
#'                         overwrite = TRUE)
#'
organize_for_projection <- function(output_dir,
                                    models = NULL,
                                    variable_names = NULL,
                                    present_file = NULL,
                                    past_files = NULL,
                                    past_period = NULL,
                                    past_gcm = NULL,
                                    future_files = NULL,
                                    future_period = NULL,
                                    future_pscen = NULL,
                                    future_gcm = NULL,
                                    fixed_variables = NULL,
                                    check_extent = TRUE,
                                    resample_to_present = TRUE,
                                    mask = NULL,
                                    overwrite = FALSE){
  #Check data
  if(!inherits(output_dir, "character")){
    stop("Argument 'output_dir' must be a 'character'.")
  }

  if(!is.null(models) && !inherits(models, "fitted_models")){
    stop("Argument 'models' must be a 'fitted_models' object or NULL.")
  }

  if (!is.null(variable_names) && !inherits(variable_names, "character")) {
    stop("Argument 'variable_names' must be NULL or a 'character'.")
  }

  if(!is.null(variable_names) && !is.null(models)){
    stop("You must provide 'variable_names' or 'models'")
  }

  #Extract variable_names from models
  if(!is.null(models)){
    variable_names <- c(models$continuous_variables,
                        models$categorical_variables)
  }

  if(!is.null(present_file) && !inherits(present_file, "character")){
    stop("Argument 'present_file' must be NULL or a 'character'.")
  }

  if(!is.null(past_files)){
    if(!inherits(past_files, "character")){
      stop("Argument 'past_files' must be NULL or a 'character'.")}
    if(!inherits(past_period, "character")){
      stop("Argument 'past_period' must be NULL or a 'character'.")}
    if(!inherits(past_gcm, "character")){
      stop("Argument 'past_gcm' must be NULL or a 'character'.")}

    #Check if files exists
    existe <- !file.exists(past_files)
    if(any(existe)){
      stop("Some of the files listed in 'past_files' do not exist.
Please ensure you used list.files(full.names = TRUE) to provide the correct file paths.")
    }

    #Check if periods exist
    period_exist <- !sapply(past_period, function(i) any(grepl(i, past_files)))
    if(any(period_exist)){
      stop("Some values in 'past_period' were not found in the file names listed in 'past_files'.
Please ensure the past_period match the file names correctly.")
    }

    #Check if gcms exist
    gcm_exist <- !sapply(past_gcm, function(i) any(grepl(i, past_files)))

    if(any(gcm_exist)){
      stop("Some values in 'past_gcm' were not found in the file names listed in 'past_files'.
Please ensure the past_gcm match the file names correctly.")
    }

  }

  if(!is.null(future_files)){
    if(!inherits(future_files, "character")){
      stop("Argument 'future_files' must be NULL or a 'character'.")}
    if(!inherits(future_period, "character")){
      stop("Argument 'future_period' must be NULL or a 'character'.")}
    if(!inherits(future_pscen, "character")){
      stop("Argument 'future_pscen' must be NULL or a 'character'.")}
    if(!inherits(future_gcm, "character")){
      stop("Argument 'future_gcm' must be NULL or a 'character'.")}

    #Check if files exists
    existe <- !file.exists(future_files)
    if(any(existe)){
      stop("Some of the files listed in 'future_files' do not exist.
Please ensure you used list.files(full.names = TRUE) to provide the correct file paths.")
    }

    #Check if periods exist
    pscen_exist <- !sapply(future_pscen, function(i) any(grepl(i, future_files)))
    if(any(period_exist)){
      stop("Some values in 'future_pscen' were not found in the file names listed in 'future_files'.
Please ensure the future_pscen match the file names correctly.")
    }

    #Check if gcms exist
    gcm_exist <- !sapply(future_gcm, function(i) any(grepl(i, future_files)))
    if(any(gcm_exist)){
      stop("Some values in 'future_gcm' were not found in the file names listed in 'future_files'.
Please ensure the future_gcm match the file names correctly.")
    }

    }

  if (!is.null(fixed_variables)) {
    if (!inherits(fixed_variables, "SpatRaster")) {
      stop("Argument 'fixed_variables' must be NULL or a 'SpatRaster'.")
    }
  }
  if (!is.null(mask) & !inherits(mask, c("SpatVector", "SpatRaster", "SpatExtent"))) {
    stop("Argument 'mask' must be a 'SpatVector', 'SpatRaster', or 'SpatExtent'.")
  }

  if(!inherits(check_extent, "logical")){
    stop("Argument 'check_extent' must be logical")
  }

  if(!inherits(overwrite, "logical")){
    stop("Argument 'overwrite' must be logical")
  }

  if(resample_to_present && is.null(present_file)){
    stop("If 'resample_to_present = TRUE', you must provide a 'present_file'")
  }

  #### Start function ####
  #Create folder
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  #### Present ####
  #Create grid of scenarios
  if(!is.null(present_file)){
    #Check if file exists
    if(!file.exists(present_file)){
      stop("The file in 'present_file' does not exist.
Please ensure you used list.files(full.names = TRUE) to provide the correct file path.")
    }

    #Read file
    r_present <- terra::rast(present_file)
    #Check names, mask and append fixed_variables
    r_present <- helper_organize_proj(r = r_present, mask, variable_names,
                                      fixed_variables, check_extent,
                                      resample_to_present = FALSE,
                                      r_present = NULL,
                                      file_name = "present_files")
    #Write results
    dir.create(file.path(output_dir, "Present"))
    terra::writeRaster(r_present,
                       file.path(output_dir, "Present/Variables.tif"),
                       overwrite = overwrite)
  } else {
    r_present <- NULL
  }

  #### Past ####
  if(!is.null(past_files)){
    #Check if file exists
    if(any(!file.exists(past_files))){
      stop("One or more files in 'past_files' do not exist.
Please ensure you used list.files(full.names = TRUE) to provide the correct file path.")
    }

    #Read file
    r_past <- lapply(past_files, terra::rast)
    names(r_past) <- gsub("\\.tif$", "", basename(past_files))

    #Get grids of scenarios in past
    g_past <- expand.grid("period" = past_period,
                          "gcm" = past_gcm)
    #For each scenario, organize and save
    for(i in 1:nrow(g_past)){
      period_i <- g_past$period[i]
      gcm_i <- g_past$gcm[i]
      #Get files
      r_past_i <- r_past[grepl(period_i, names(r_past)) &
                            grepl(gcm_i, names(r_past))][[1]]
      r_past_i <- helper_organize_proj(r = r_past_i, mask, variable_names,
                                       fixed_variables, check_extent,
                                       resample_to_present = resample_to_present,
                                       r_present = r_present,
                                       file_name = "past_files")
      #Create recursive folder to save
      dir_past_i <- file.path(output_dir, "Past", period_i, gcm_i)
      dir.create(dir_past_i, showWarnings = FALSE, recursive = TRUE)
      #Save
      terra::writeRaster(r_past_i,
                         filename = file.path(dir_past_i, "Variables.tif"),
                         overwrite = overwrite)
    }
  }

  #### Future ####
  if(!is.null(future_files)){
    #Check if file exists
    if(any(!file.exists(future_files))){
      stop("One or more files in 'future_files' do not exist.
Please ensure you used list.files(full.names = TRUE) to provide the correct file path.")
    }

    #Read file
    r_future <- lapply(future_files, terra::rast)
    names(r_future) <- gsub("\\.tif$", "", basename(future_files))

    #Get grids of scenarios in future
    g_future <- expand.grid("period" = future_period,
                            "ssp" = future_pscen,
                            "gcm" = future_gcm)
    #For each scenario, organize and save
    for(i in 1:nrow(g_future)){
      period_i <- g_future$period[i]
      ssp_i <- g_future$ssp[i]
      gcm_i <- g_future$gcm[i]
      #Get files
      r_future_i <- r_future[grepl(period_i, names(r_future)) &
                           grepl(ssp_i, names(r_future)) &
                             grepl(gcm_i, names(r_future))][[1]]
      r_future_i <- helper_organize_proj(r = r_future_i, mask, variable_names,
                                       fixed_variables, check_extent,
                                       resample_to_present = resample_to_present,
                                       r_present = r_present,
                                       file_name = "future_files")
      #Create recursive folder to save
      dir_future_i <- file.path(output_dir, "future", period_i, ssp_i, gcm_i)
      dir.create(dir_future_i, showWarnings = FALSE, recursive = TRUE)
      #Save
      terra::writeRaster(r_future_i,
                         filename = file.path(dir_future_i, "Variables.tif"),
                         overwrite = overwrite)
    }
  }

  #Return folders
  return(message("\nVariables successfully organized in directory:\n",
                 output_dir))
}
