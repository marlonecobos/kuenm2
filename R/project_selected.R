#' Project selected models for multiple scenarios
#'
#' @description
#' This function predicts selected models across multiple scenarios, as specified in a `prepared_proj` object created with the \code{\link{prepare_proj()}}() function. In addition to generating predictions for each replicate, the function calculates consensus measures (e.g., mean, median) across replicates and models.
#'
#' @param models an object of class `fitted_models` returned by the
#' \code{\link{fit_selected}}() function.
#' @param projection_data an object of class `prepared_proj` returned by the \code{\link{prepare_proj()}}() function. This file contains the paths to the rasters representing each scenario.
#' @param out_dir (character) a path to a root directory for saving the raster file of each projection.
#' @param mask (SpatRaster, SpatVector, or SpatExtent) spatial object used to
#'        mask the variables before predict. Default is NULL.
#' @param consensus_per_model (logical) whether to calculate consensus across replicates when there are more than one replicate per model. Default is TRUE.
#' @param consensus_general (logical) whether to calculate consensus across models when there are more than one selected model. Default is TRUE.
#' @param consensus (character) consensus measures to calculate. Options available are 'median', 'range', 'mean' and 'stdev' (standard deviation). Default is c("median", "range", "mean", "stdev").
#' @param write_replicates (logical) whether to write the projections for each replicate. Default is FALSE.
#' @param clamping (logical) whether to restricts variable values to the range of the calibration data to avoid extrapolation. Default is `TRUE` (free extrapolation).
#' @param var_to_clamp (character) vector specifying which variables to clamp. Only applicable if `clamping = TRUE`. Default is `NULL`, meaning all variables will be clamped.
#' @param type (character) the format of the prediction values. Available options are `"raw"`, `"cumulative"`, `"logistic"`, and `"cloglog"`. Default is `"cloglog"`.
#' @param overwrite (logical) whether to overwrite SpatRaster if they already exists. Only applicable if `write_files` is set to TRUE. Default is FALSE.
#' @param parallel (logical) whether to fit the candidate models in parallel.
#' Default is FALSE.
#' @param ncores (numeric) number of cores to use for parallel processing.
#' Default is 1. This is only applicable if `parallel = TRUE`.
#' @param parallelType (character) the package to use for parallel processing:
#' "doParallel" or "doSNOW". Default is "doSNOW". This is only applicable if
#' `parallel = TRUE`.
#' @param progress_bar (logical) whether to display a progress bar during processing. Default is TRUE.
#' @param verbose (logical) whether to display messages during processing. Default is TRUE.
#'
#' @export
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom terra crop wrap unwrap
#'
#' @return A `model_projections` object that provides the paths to the raster files with the projection results and the corresponding thresholds used to binarize the predictions.
#'
#' @usage project_selected(models, projection_data,
#'                                out_dir,
#'                                mask = NULL,
#'                                consensus_per_model = TRUE,
#'                                consensus_general = TRUE,
#'                                consensus = c("median", "range", "mean", "stdev"),
#'                                write_replicates = FALSE, clamping = FALSE,
#'                                var_to_clamp = NULL, type = "cloglog",
#'                                overwrite = FALSE, parallel = FALSE,
#'                                ncores = 1, parallelType = "doSNOW",
#'                                progress_bar = TRUE, verbose = TRUE)
#'
#' @examples
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
#'
#' #Example with GLMNET
#' # Import example of fitted_models (output of fit_selected())
#' data("fitted_model_glmnet", package = "kuenm2")
#'
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
#'
#' #Create folder to save projection results
#' out_dir <- file.path(tempdir(), "Projection_results/glmnet")
#' dir.create(out_dir, recursive = TRUE)
#'
#' # Project selected models for multiple scenarios
#' p <- project_selected(models = fitted_model_glmnet,
#'                       projection_data = pr,
#'                       out_dir = out_dir,
#'                       consensus_per_model = TRUE,
#'                       consensus_general = TRUE,
#'                       consensus = c("median", "range", "mean", "stdev"),
#'                       write_replicates = FALSE,
#'                       clamping = FALSE,
#'                       var_to_clamp = NULL,
#'                       type = "cloglog",
#'                       overwrite = TRUE,
#'                       parallel = FALSE,
#'                       ncores = 1,
#'                       parallelType = "doSNOW",
#'                       progress_bar = TRUE,
#'                       verbose = TRUE)
#'
#' #Example with GLMN
#' # Import example of fitted_models (output of fit_selected())
#' data("fitted_model_glm", package = "kuenm2")
#'
#' # Prepare projections using fitted models to check variables
#' pr <- prepare_proj(models = fitted_model_glm,
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
#'
#' #Create folder to save projection results
#' out_dir <- file.path(tempdir(), "Projection_results/glm")
#' dir.create(out_dir, recursive = TRUE)
#'
#' # Project selected models for multiple scenarios
#' p_glm <- project_selected(models = fitted_model_glm,
#'                           projection_data = pr,
#'                           out_dir = out_dir,
#'                           consensus_per_model = TRUE,
#'                           consensus_general = TRUE,
#'                           consensus = c("median", "range", "mean", "stdev"),
#'                           write_replicates = FALSE,
#'                           clamping = FALSE,
#'                           var_to_clamp = NULL,
#'                           type = "cloglog",
#'                           overwrite = TRUE,
#'                           parallel = FALSE,
#'                           ncores = 1,
#'                           parallelType = "doSNOW",
#'                           progress_bar = TRUE,
#'                           verbose = TRUE)
#'
project_selected <- function(models,
                             projection_data,
                             out_dir,
                             mask = NULL,
                             consensus_per_model = TRUE,
                             consensus_general = TRUE,
                             consensus = c("median", "range", "mean", "stdev"),
                             write_replicates = FALSE,
                             clamping = FALSE,
                             var_to_clamp = NULL,
                             type = "cloglog",
                             overwrite = FALSE,
                             parallel = FALSE,
                             ncores = 1,
                             parallelType = "doSNOW",
                             progress_bar = TRUE,
                             verbose = TRUE){
  #Check data
  if (!inherits(models, "fitted_models")) {
    stop(paste0("Argument models must be a fitted_models object, not ",
                class(models)))
  }
  if (!inherits(projection_data, "projection_data")) {
    stop(paste0("Argument projection_data must be a projection_data object, not ",
                class(projection_data)))
  }
  if (!inherits(out_dir, "character")) {
    stop(paste0("Argument out_dir must be a character, not ",
                class(out_dir)))
  }

  if(!is.null(mask) & !inherits(mask, c("SpatRaster", "SpatVector",
                                        "SpatExtent"))){
    stop(paste0("Argument mask must be a SpatVector, SpatExtent or SpatRaster, not ",
                class(mask)))
  }

  if (!inherits(consensus_per_model, "logical")) {
    stop(paste0("Argument consensus_per_model must be logical, not ",
                class(consensus_per_model)))
  }
  if (!inherits(consensus_general, "logical")) {
    stop(paste0("Argument consensus_general must be logical, not ",
                class(consensus_general)))
  }
  if (!inherits(consensus, "character")) {
    stop(paste0("Argument consensus must be a character, not ",
                class(consensus)))
  }
  consensus_out <- setdiff(consensus, c("median", "range", "mean", "stdev"))
  if(length(consensus_out) > 0){
    stop("Invalid consensus provided.
  The available options are: 'median', 'range', 'mean' and 'stdev'")
  }
  if(!any(c("median", "mean") %in% consensus)){
    stop("Consensus must have at least one of the options: 'median' or 'mean'")
  }
  if (!inherits(write_replicates, "logical")) {
    stop(paste0("Argument write_replicates must be logical, not ",
                class(write_replicates)))
  }
  if (!inherits(clamping, "logical")) {
    stop(paste0("Argument clamping must be logical, not ",
                class(clamping)))
  }
  if (clamping & !is.null(var_to_clamp) & !inherits(var_to_clamp, "character")) {
    stop(paste0("Argument var_to_clamp must be NULL or character, not ",
                class(var_to_clamp)))
  }
  if (!inherits(type, "character")) {
    stop(paste0("Argument type must be character, not ",
                class(type)))
  }
  if(!any(c("raw", "cumulative", "logistic", "cloglog") %in% type)){
    stop("Invalid type provided.
  The available options are: 'raw', 'cumulative', 'logistic', or 'cloglog'")
  }
  if (!inherits(overwrite, "logical")) {
    stop(paste0("Argument overwrite must be logical, not ",
                class(overwrite)))}
  if (!inherits(parallel, "logical")) {
    stop(paste0("Argument parallel must be logical, not ",
                class(parallel)))}
  if (!inherits(ncores, "numeric")) {
    stop(paste0("Argument ncores must be numeric, not ",
                class(ncores)))}
  if (!inherits(parallelType, "character")) {
    stop(paste0("Argument parallelType must be character, not ",
                class(type)))
  }
  if(!any(c("doParallel", "doSNOW") %in% parallelType)){
    stop("Invalid parallelType provided. The available options are: 'doParallel' or 'doSNOW'")
  }
  if (!inherits(progress_bar, "logical")) {
    stop(paste0("Argument progress_bar must be logical, not ",
                class(progress_bar)))}
  if (!inherits(verbose, "logical")) {
    stop(paste0("Argument verbose must be logical, not ",
                class(verbose)))}

  #Save parameters in a list to send to foreach nodes#
  par_list <- list(models = models,
                   mask = mask,
                   projection_data = projection_data,
                   consensus_per_model = consensus_per_model,
                   consensus_general = consensus_general,
                   consensus = consensus,
                   write_replicates = write_replicates,
                   clamping = clamping,
                   var_to_clamp = var_to_clamp,
                   type = type,
                   overwrite = overwrite)

  ####PREPARE DATAFRAME TO PREDICT####
  #Extract variables from best models
  vars <- names(models[["Models"]][[1]][[1]]$samplemeans)
  vars <- setdiff(vars, c("pr_bg", "fold"))

  if(!file.exists(out_dir)){
    dir.create(out_dir, showWarnings = FALSE)
  }
  #Normalize path
  out_dir <- normalizePath(out_dir)

  # #Check scenarios to predict
  # sc <- names(projection_data[sapply(projection_data, function(x) !is.null(x))])
  #
  # #Get raster pattern to read
  # raster_pattern <- projection_data$Raster_pattern
  #
  # #Get dataframe with path to predictions
  # #Present
  # if("Present" %in% sc){
  #   #Create folder
  #   present_dir <- file.path(out_dir, "Present/")
  #   present_sc <- names(projection_data[["Present"]])
  #   suppressWarnings({
  #     d_present <- data.frame(Time = "Present",
  #                             Period = "Present",
  #                             Scenario = present_sc,
  #                             input_path = unlist(projection_data[["Present"]]),
  #                             output_path = normalizePath(file.path(present_dir,
  #                                                                   present_sc)))})
  # }
  # #Past
  # if("Past" %in% sc){
  #   #Create folder
  #   past_dir <- file.path(out_dir, "Past/")
  #   #Get grid of projections
  #   df_past <- do.call(rbind, lapply(names(projection_data$Past), function(time) {
  #     time_data <- projection_data$Past[[time]]
  #     do.call(rbind, lapply(names(time_data), function(gcm) {
  #       data.frame(Time = "Past", Period = time, GCM = gcm, Path = time_data[[gcm]], stringsAsFactors = FALSE)
  #     }))
  #   }))
  #
  #   #Looping in the grid
  #
  #   #Create dataframe with path to results
  #   suppressWarnings({
  #     d_past <- data.frame(Time = "Past",
  #                          Period = df_past$Period,
  #                          GCM = df_past$GCM,
  #                          input_path = df_past$Path,
  #                          output_path = normalizePath(file.path(past_dir, df_past$Period,
  #                                                                df_past$GCM),
  #                                                      mustWork = FALSE))})
  # }
  # #Future
  # ####Project to Future scenarios####
  # if("Future" %in% sc){
  #   #Create folder
  #   future_dir <- file.path(out_dir, "Future/")
  #
  #   #Create grid of time-ssp-gcm
  #   df_future <- do.call(rbind, lapply(names(projection_data[["Future"]]), function(year_range) {
  #     year_range_data <- projection_data[["Future"]][[year_range]]
  #     do.call(rbind, lapply(names(year_range_data), function(ssp) {
  #       ssp_data <- year_range_data[[ssp]]
  #       do.call(rbind, lapply(names(ssp_data), function(gcm) {
  #         data.frame(Time = "Future", Period = year_range, ssp = ssp,
  #                    GCM = gcm, Path = ssp_data[[gcm]],
  #                    stringsAsFactors = FALSE)
  #       }))      }))    }))
  #
  #
  #   #Create dataframe with path to results
  #   suppressWarnings({
  #     d_future <- data.frame(Time = df_future$Time,
  #                            Period = df_future$Period,
  #                            ssp = df_future$ssp,
  #                            GCM = df_future$GCM,
  #                            input_path = df_future$Path,
  #                            output_path = normalizePath(file.path(future_dir, df_future$Period,
  #                                                                  df_future$ssp, df_future$GCM),
  #                                                        mustWork = FALSE))})
  # }
  #
  # #Get dataframe with path to each projection
  # if(!("Present" %in% sc)){
  #   d_present <- NULL
  # }
  # if(!("Past" %in% sc)){
  #   d_past <- NULL
  # }
  # if(!("Future" %in% sc)){
  #   d_future <- NULL
  # }
  #
  # #Return and write files with path
  # res_path <- bind_rows_projection(list(d_present, d_past, d_future))
  # #Create ID
  # res_path$id <- 1:nrow(res_path)

  #Check scenarios for predicting
  res_path <- check_pred_scenarios(projection_data = projection_data,
                                   out_dir = out_dir)
  raster_pattern <- projection_data$raster_pattern


  ####Configure parallelization####
  n_models <- nrow(res_path)
  if(n_models == 1 & isTRUE(parallel)){
    parallel <- FALSE
  } else {parallel <- TRUE}

  if(n_models < ncores & isTRUE(parallel)){
    ncores <- n_models
  } else {ncores = ncores}

  #Show progress bar?
  if (progress_bar) {
    pb <- utils::txtProgressBar(0, n_models, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n) }


  #Set parallelization
  if(parallel) {
    #Make cluster
    cl <- parallel::makeCluster(ncores)

    if (parallelType == "doParallel") {
      doParallel::registerDoParallel(cl)
      opts <- NULL
    }

    if (parallelType == "doSNOW") {
      doSNOW::registerDoSNOW(cl)
      if (progress_bar)
        opts <- list(progress = progress)
      else opts <- NULL
    }
  } else {
    opts <- NULL}
  ###################################

  #If parallel and mask, wrap variables
  if(!is.null(mask) & parallel)
    par_list$mask <- terra::wrap( par_list$mask)

  #Run predictions
  if(parallel){
    foreach(x = 1:n_models,
            .options.snow = opts
            ,
            .export = c("predict_selected", "multiple_projections")
    ) %dopar% {
      multiple_projections(i = x, res_path, raster_pattern, par_list)
    }
  } else { #Not in parallel (using %do%)
    # Loop for com barra de progresso manual
    for (x in 1:n_models) {
      # Execute a função fit_eval_models
      multiple_projections(i = x, res_path, raster_pattern, par_list)

      # Sets the progress bar to the current state
      if(progress_bar){
        utils::setTxtProgressBar(pb, x) }
    }
  }

  #Stop cluster
  if(parallel){parallel::stopCluster(cl)}

  #Append threshold to final results
  res_final <- new_model_projections(paths = res_path,
                                     thresholds = models$thresholds)

  #Save
  saveRDS(res_final, file.path(out_dir, "Projection_paths.RDS"))
  return(res_final)
} #End of function

# #Test function internally
# models = fm
# projection_data = pr
# out_dir = out_dir
# #write_path = TRUE
# consensus_per_model = TRUE
# consensus_general = TRUE
# consensus = c("median", "range", "mean", "stdev") #weighted mean
# write_replicates = FALSE
# clamping = FALSE
# var_to_clamp = NULL
# type = "cloglog"
# overwrite = TRUE
# parallel = TRUE
# ncores = 8
# parallelType = "doSNOW"
# progress_bar = TRUE
# verbose = TRUE
