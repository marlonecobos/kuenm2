#' Preparation of data for model projections
#' @description
#' This function prepared data for model projections to multiple scenarios, storing the paths to the rasters representing each scenario
#'
#' @param models an object of class `fitted_models` returned by the
#' \code{\link{fit_selected_glmnetmx}}() function. Default is NULL.
#' @param variable_names (character) names of the variables used to fit the model or do the PCA in the \code{\link{prepare_data}}() function. Only applicable if `models` argument is not provided. Default is NULL.
#' @param present_dir (character) path to the folder containing variables that represent the current scenario for projection. Default is NULL.
#' @param past_dir (character) path to the folder containing subfolders with variables representing past scenarios for projection. Default is NULL.
#' @param past_period (character) names of the subfolders within `past_dir`, representing specific time periods (e.g., 'LGM' or 'MID').
#' @param past_gcm (character) names of the subfolders within `past_period` folders, representing specific General Circulation Models (GCMs).
#' @param future_dir (character) path to the folder containing subfolders with variables representing future scenarios for projection. Default is NULL.
#' @param future_period (character) names of the subfolders within `future_dir`, representing specific time periods (e.g., '2041-2060' or '2081-2100'). Default is NULL.
#' @param future_pscen (character) names of the subfolders within `future_period`, representing specific emission scenarios (e.g., 'ssp126' or 'ssp585'). Default is NULL.
#' @param future_gcm (character) names of the subfolders within `future_pscen` folders, representing specific General Circulation Models (GCMs). Default is NULL.
#' @param write_file (logical) whether to write the object containing the paths to the structured folders. This object is required for projecting models across multiple scenarios using the \code{\link{project_selected}}() function. Default is FALSE.
#' @param filename (character) the path or name of the folder where the object will be saved. This is only applicable if `write_file = TRUE`. Default is NULL.
#' @param raster_pattern (character) pattern used to identify the format of raster files within the folders. Default is ".tif*".
#'
#' @importFrom terra rast
#' @export
#' @return An object of class `prepared_proj` containing the following elements:
#' - Present, Past, and Future: paths to the variables structured in subfolders.
#' - Raster_pattern: the pattern used to identify the format of raster files within the folders.
#' - PCA: if a principal component analysis (PCA) was performed on the set of variables with \code{\link{prepare_data}}(), a list with class "prcomp" will be returned. See `?stats::prcomp()` for details.
#' - variables: names of the raw predictos variables used to project.
#'
#' @usage prepare_proj(models = NULL, variable_names = NULL, present_dir = NULL,
#'                    past_dir = NULL, past_period = NULL, past_gcm = NULL,
#'                    future_dir = NULL, future_period = NULL,
#'                    future_pscen = NULL, future_gcm = NULL,
#'                    write_file = FALSE, filename = NULL,
#'                    raster_pattern = ".tif*")
#'
#' @examples
#' # Import example of fitted_models (output of fit_selected())
#' data("fitted_model_glmnet", package = "kuenm2")
#'
#' # Organize and structure future climate variables from WorldClim
#' # Import the current variables used to fit the model.
#' # In this case, SoilType will be treated as a static variable (constant across future scenarios).
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#' # Create a "Current_raw" folder in a temporary directory and copy the raw variables there.
#' out_dir_current <- file.path(tempdir(), "Current_raw")
#' dir.create(out_dir_current, recursive = TRUE)
#' # Save current variables in temporary directory
#' writeRaster(var, file.path(out_dir_current, "Variables.tif"))
#'
#' # Set the input directory containing the raw future climate variables.
#' # For this example, the data is located in the "inst/extdata" folder.
#' in_dir <- system.file("extdata", package = "kuenm2")
#' # Create a "Future_raw" folder in a temporary directory and copy the raw variables there.
#' out_dir_future <- file.path(tempdir(), "Future_raw")
#' # Organize and rename the future climate data, structuring it by year and GCM.
#' # The 'SoilType' variable will be appended as a static variable in each scenario.
#' # The files will be renamed following the "bio_" format
#' organize_future_worldclim(input_dir = in_dir,
#'                           output_dir = out_dir_future,
#'                           name_format = "bio_", variables = NULL,
#'                           fixed_variables = var$SoilType, mask = NULL,
#'                           overwrite = TRUE)
#' # Prepare projections using fitted models to check variables
#' pr <- prepare_proj(models = fitted_model_glmnet,
#'                    present_dir = out_dir_current,
#'                    past_dir = NULL,
#'                    past_period = NULL,
#'                    past_gcm = NULL,
#'                    future_dir = out_dir_future,
#'                    future_period = c("2041-2060", "2081-2100"),
#'                    future_pscen = c("ssp126", "ssp585"),
#'                    future_gcm = c("ACCESS-CM2", "MIROC6"),
#'                    write_file = FALSE,
#'                    filename = NULL,
#'                    raster_pattern = ".tif*")
#' pr
#'
#' # Prepare projections using variables names
#' pr_b <- prepare_proj(models = NULL,
#'                      variable_names = c("bio_1", "bio_7", "bio_12"),
#'                      present_dir = out_dir_current,
#'                      past_dir = NULL,
#'                      past_period = NULL,
#'                      past_gcm = NULL,
#'                      future_dir = out_dir_future,
#'                      future_period = c("2041-2060", "2081-2100"),
#'                      future_pscen = c("ssp126", "ssp585"),
#'                      future_gcm = c("ACCESS-CM2", "MIROC6"),
#'                      write_file = FALSE,
#'                      filename = NULL,
#'                      raster_pattern = ".tif*")
#' pr_b

prepare_proj <- function(models = NULL,
                         variable_names = NULL,
                         present_dir = NULL,
                         past_dir = NULL, past_period = NULL, past_gcm = NULL,
                         future_dir = NULL, future_period = NULL,
                         future_pscen = NULL, future_gcm = NULL,
                         write_file = FALSE,
                         filename = NULL,
                         raster_pattern = ".tif*") {
  #Check data
  if(is.null(models) & is.null(variable_names)) {
    stop("You must specify models or variable_names")
  }

  if (!is.null(models) & !inherits(models, "fitted_models")) {
    stop(paste0("Argument models must be NULL or a fitted_models object, not ",
                class(models)))
  }
  if (!is.null(variable_names) & !inherits(variable_names, "character")) {
    stop(paste0("Argument variable_names must be NULL or a character, not ",
                class(variable_names)))
  }
  if (!is.null(present_dir) & !inherits(present_dir, "character")) {
    stop(paste0("Argument present_dir must be NULL or a character, not ",
                class(present_dir)))
  }
  if (!is.null(past_dir)){
    if(!inherits(past_dir, "character")) {
      stop(paste0("Argument past_dir must be NULL or a character, not ",
                  class(past_dir)))}
    if(is.null(past_period))
      stop("If past_dir is not NULL, past_period must be specified")
    if (!is.null(past_period) & !inherits(past_period, "character")) {
      stop(paste0("Argument past_period must be a character, not ",
                  class(past_period)))
    }
    if(is.null(past_gcm))
      stop("If past_dir is not NULL, past_gcm must be specified")
    if (!is.null(past_gcm) & !inherits(past_gcm, "character")) {
      stop(paste0("Argument past_gcm must be a character, not ",
                  class(past_gcm)))
    }}

  if (!is.null(future_dir)){
    if(!inherits(future_dir, "character")) {
      stop(futuree0("Argument future_dir must be NULL or a character, not ",
                  class(future_dir)))}
    if(is.null(future_period))
      stop("If future_dir is not NULL, future_period must be specified")
    if (!is.null(future_period) & !inherits(future_period, "character")) {
      stop(futuree0("Argument future_period must be a character, not ",
                  class(future_period)))
    }
    if(is.null(future_pscen))
      stop("If future_dir is not NULL, future_pscen must be specified")
    if (!is.null(future_pscen) & !inherits(future_pscen, "character")) {
      stop(futuree0("Argument future_pscen must be a character, not ",
                    class(future_pscen)))
    }
    if(is.null(future_gcm))
      stop("If future_dir is not NULL, future_gcm must be specified")
    if (!is.null(future_gcm) & !inherits(future_gcm, "character")) {
      stop(futuree0("Argument future_gcm must be a character, not ",
                  class(future_gcm)))
    }
  }

  if(!inherits(write_file, "logical")) {
    stop(paste0("Argument write_file must be logical, not ",
                class(write_file)))
  }

  if(write_file & is.null(filename)){
    stop("If write_file = TRUE, filename must be specified")
  }

  if(write_file & !inherits(filename, "character")){
    stop(paste0("Argument filename must be character, not ",
                class(filename)))
  }

  if (!inherits(raster_pattern, "character")) {
    stop(paste0("Argument raster_pattern must be a character, not ",
                class(raster_pattern)))
  }

  #Get variables used to fit models
  if(!is.null(models)){
  vars <- c(models$continuous_variables, models$categorical_variables)}

  if(is.null(models) & !is.null(variable_names)){
    vars <- variable_names
  }

  if(!is.null(models) & !is.null(variable_names)) {
    warning("You specified models and variable names. Using variables names from models")
  }

  ####Check directories with present projections####
  if(!is.null(present_dir)) {
    #Check folders
    if(!file.exists(present_dir)){
      stop(paste("present_dir", present_dir, "does not exist"))
    }

    #List internal directories or files
    internal_dir <- list.dirs(present_dir, recursive = T)[-1]

    #Check if there is a file in the directory
    fdir <- list.files(present_dir, pattern = raster_pattern)

    #Check if there are several scenarios
    pdir <- list.dirs(present_dir)[-1]

    if(length(pdir) == 0 & length(fdir) == 0 & length(internal_dir) == 0) {
      stop("present_dir ", present_dir, " has no folders or files to project")
    }

    #To project using a single file
    if(length(pdir) == 0 & length(fdir) > 0) {
      r <- terra::rast(file.path(present_dir, fdir))
      #Check absent vars
      abs_vars <- vars[!(vars %in% names(r))]
      if(length(abs_vars) > 0)
        stop("The following vars are absent from current projection folder", ":\n",
             paste(abs_vars, collapse = "\n"))
      #If everything is OK, create list with the path
      res_present <- list()
      res_present[["Present"]] <- normalizePath(present_dir)
    }

    #To project using folders
    if(length(pdir) > 0) {
      #Check if variables are in the folders
      invisible(sapply(pdir, function(x){
        r <- list.files(x, pattern = raster_pattern, full.names = TRUE)
        if(length(r) == 0) {
          stop("The directory ", x, " has no ", raster_pattern, " files")
        }
        r <- terra::rast(r)
        #Check absent vars
        abs_vars <- vars[!(vars %in% names(r))]
        if(length(abs_vars) > 0)
          stop("The following vars are absent from ", x, ":\n",
               paste(abs_vars, collapse = "\n"))
      }))

      #If everything is OK, create list with the path
      #Get scenarios
      sc <- gsub(present_dir, "", pdir, fixed = TRUE)
      res_present <- list()
      for (scenario in sc) {
        res_present[[scenario]] <- normalizePath(file.path(present_dir,
                                                           scenario))}
    }

  } else {
    res_present <- NULL }#End of check present

  #####Check directories with future projections####
  if(!is.null(future_dir)) {

    #Check folders
    if(!file.exists(future_dir)){
      stop(paste("future_dir", future_dir, "does not exist"))
    }

    #List internal directories
    internal_dir <- list.dirs(future_dir, recursive = T)[-1]

    if(length(internal_dir) == 0) {
      stop(paste("future_dir", future_dir, "is empty"))
    }

  if(sum(!is.null(future_period),
         !is.null(future_pscen),
         !is.null(future_gcm)) != 3) {
    stop("To prepare projections for future time, you must set future_period, future_pscen and future_gcm arguments")
  }

  #Check folders
  expected_folders_future <- unlist(sapply(future_period, function(i){
    lapply(future_pscen, function(x){
      file.path(future_dir, i, x, future_gcm)
    })
  }, simplify = F, USE.NAMES = FALSE))
  future_exists <- unlist(sapply(expected_folders_future, file.exists,
                                 simplify = TRUE, USE.NAMES = FALSE))

  if(any(future_exists == FALSE)){
    eff <- expected_folders_future[future_exists == FALSE]
    stop("The following folders do not exist: ", "\n", paste(eff, collapse = "\n"))
  }

  #Check if variables are in the folders
  all_dir <- expected_folders_future
  invisible(sapply(all_dir, function(x){
    r <- list.files(x, pattern = raster_pattern, full.names = TRUE)
    if(length(r) == 0) {
      stop("The directory ", x, " has no ", raster_pattern, " files")
    }
    r <- terra::rast(r)
    #Check absent vars
    abs_vars <- vars[!(vars %in% names(r))]
    if(length(abs_vars) > 0)
    stop("The following vars are absent from ", x, ":\n",
         paste(abs_vars, collapse = "\n"))
  }))

  #If everything is OK, create list with the path

  res_future <- list()
  for (year in future_period) {
    res_future[[year]] <- list()
  for (ssp in future_pscen) {
    res_future[[year]][[ssp]] <- list()
    for (gcm in future_gcm) {
      res_future[[year]][[ssp]][[gcm]] <- normalizePath(file.path(future_dir, year,
                                                                  ssp, gcm))
    }}}

  } else {
    res_future <- NULL} #End of check future

  #####Check directories with past projections####
  if(!is.null(past_dir)) {

    #Check folders
    if(!file.exists(past_dir)){
      stop(paste("past_dir", past_dir, "does not exist"))
    }

    #List internal directories
    internal_dir <- list.dirs(past_dir, recursive = T)[-1]

    if(length(internal_dir) == 0) {
      stop(paste("past_dir", past_dir, "is empty"))
    }

    if(sum(!is.null(past_period),
           !is.null(past_gcm)) != 2) {
      stop("To prepare projections for past time, you must set past_period and past_gcm arguments")
    }


  #Check folders
  expected_folders_past <- unlist(sapply(past_period, function(i){
      file.path(past_dir, i, past_gcm)
  }, simplify = F, USE.NAMES = FALSE))
  past_exists <- unlist(sapply(expected_folders_past, file.exists,
                                 simplify = TRUE, USE.NAMES = FALSE))

  if(any(past_exists == FALSE)){
    eff <- expected_folders_past[past_exists == FALSE]
    stop("The following folders do not exist: ", "\n", paste(eff,
                                                             collapse = "\n"))
  }

  #Check if variables are in the folders
  all_dir <- expected_folders_past
  invisible(sapply(all_dir, function(x){
    r <- list.files(x, pattern = raster_pattern, full.names = TRUE)
    if(length(r) == 0) {
      stop("The directory ", x, " has no ", raster_pattern, " files")
    }
    r <- terra::rast(r)
    #Check absent vars
    abs_vars <- vars[!(vars %in% names(r))]
    if(length(abs_vars) > 0)
      stop("The following vars are absent from ", x, ":\n",
           paste(abs_vars, collapse = "\n"))
  }))

  #If everything is OK, create list with the path
  res_past <- list()
  for (year in past_period) {
    res_past[[year]] <- list()
      for (gcm in past_gcm) {
        res_past[[year]][[gcm]] <- normalizePath(file.path(past_dir, year, gcm))
      }}
  } else {
    res_past <- NULL} #End of check past

  #Get pca
  if(!is.null(models)){
    if(!is.null(models$pca)){
      pca <- models$pca
    } else pca <- NULL
  } else pca <- NULL

  #Construct prepared_proj object
  res <- new_projection_data(res_present = res_present,
                           res_past = res_past,
                           res_future = res_future,
                           raster_pattern = raster_pattern,
                           variables = vars,
                           pca = pca)
  #Save results as RDS
  if(write_file){
  saveRDS(res, paste0(filename, ".RDS"))}

  return(res)
} #End of function


# #Check function
# bm <- readRDS("../test_kuenm2/Best_Models.RDS")
#
# pr <- prepare_proj(models = bm,
#              present_dir = "../test_kuenm2/Projections/Present/",
#              past_dir = "../test_kuenm2/Projections/Past/",
#              past_period = c("LGM", "MID"),
#              past_gcm = c("CCSM4", "MIROC-ESM", "MPI-ESM-P"),
#              future_dir = "../test_kuenm2/Projections/Future/",
#              future_period = c("2041-2060", "2081-2100"),
#              future_pscen = c("ssp245", "ssp585"),
#              future_gcm = c("BCC-CSM2-MR", "ACCESS-CM2", "CMCC-ESM2"),
#              filename = "../test_kuenm2/Projection_file",
#              raster_pattern = ".tif*")
#
#
#
# bm <- readRDS("../test_kuenm2/Best_Models.RDS")
# models <- bm
# future_dir <- "../test_kuenm2/Projections/Future/"
# future_period <- "2041-2060"
# future_pscen <- c("ssp245", "ssp585")
# future_gcm <- c("BCC-CSM2-MR", "ACCESS-CM2", "CMCC-ESM2")
# raster_pattern = ".tif*"
# past_dir <- "../test_kuenm2/Projections/Past/"
# past_period = c("LGM", "MID")
# past_gcm = c("CCSM4", "MIROC-ESM", "MPI-ESM-P")
# present_dir <- "../test_kuenm2/Projections/Present/"
# filename = "../test_kuenm2/Projection_file"
