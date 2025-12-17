#' Project selected models to multiple sets of new data (scenarios)
#'
#' @description
#' This function performs predictions of selected models on multiple scenarios,
#' as specified in a `projection_data` object created with the
#' [prepare_projection()] function. In addition to generating predictions
#' for each replicate, the function calculates consensus measures (e.g., mean,
#' median) across replicates and models.
#'
#' @usage
#' project_selected(models, projection_data, out_dir, mask = NULL,
#'                  consensus_per_model = TRUE, consensus_general = TRUE,
#'                  consensus = c("median", "range", "mean", "stdev"),
#'                  write_replicates = FALSE, extrapolation_type = "E",
#'                  var_to_clamp = NULL, type = NULL, overwrite = FALSE,
#'                  parallel = FALSE, ncores = NULL,
#'                  progress_bar = TRUE, verbose = TRUE)
#'
#' @param models an object of class `fitted_models` returned by the
#' [fit_selected()] function.
#' @param projection_data an object of class `projection_data` returned by the
#' [prepare_projection()] function. This file contains the paths to the
#' rasters representing each scenario.
#' @param out_dir (character) a path to a root directory for saving the raster
#' file of each projection.
#' @param mask (SpatRaster, SpatVector, or SpatExtent) spatial object used to
#'        mask the variables before predict. Default is NULL.
#' @param consensus_per_model (logical) whether to calculate consensus across
#' replicates when there are more than one replicate per model. Default is TRUE.
#' @param consensus_general (logical) whether to calculate consensus across
#' models when there are more than one selected model. Default is TRUE.
#' @param consensus (character) consensus measures to calculate. Options
#' available are 'median', 'range', 'mean' and 'stdev' (standard deviation).
#' Default is c("median", "range", "mean", "stdev").
#' @param write_replicates (logical) whether to write the projections for each
#' replicate. Default is FALSE.
#' @param extrapolation_type (character) extrapolation type of model. Models can
#' be transferred with three options: free extrapolation ('E'), extrapolation
#' with clamping ('EC'), and no extrapolation ('NE'). Default = 'E'. See details.
#' @param var_to_clamp (character) vector specifying which variables to clamp.
#' Only applicable if extrapolation_type is "EC" or "NE". Default is `NULL`, meaning all
#' variables will be clamped or not extrapolated.
#' @param type (character) the format of prediction values. For `maxnet` models,
#' valid options are `"raw"`, `"cumulative"`, `"logistic"`, and `"cloglog"`. For
#' `glm` models, valid options are `"response"` and `"raw"`. If `NULL` (default),
#' the function uses `"cloglog"` for `maxnet` models and `"response"` for `glm`
#' models.
#' @param overwrite (logical) whether to overwrite SpatRaster if they already
#' exists. Only applicable if `write_files` is set to TRUE. Default is FALSE.
#' @param parallel (logical) whether to fit the candidate models in parallel.
#' Default is FALSE.
#' @param ncores (numeric) number of cores to use for parallel processing.
#' Default is NULL and uses available cores - 1. This is only applicable if
#' `parallel = TRUE`.
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is TRUE.
#' @param verbose (logical) whether to display messages during processing.
#' Default is TRUE.
#'
#' @export
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom terra crop wrap unwrap
#'
#' @return
#' A `model_projections` object that provides the paths to the raster
#' files with the projection results and the corresponding thresholds used to
#' binarize the predictions.
#'
#' @seealso
#' [organize_future_worldclim()], [prepare_projection()]
#'
#' @examples
#' # Step 1: Organize variables for current projection
#' ## Import current variables (used to fit models)
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' ## Create a folder in a temporary directory to copy the variables
#' out_dir_current <- file.path(tempdir(), "Current_raw_wc")
#' dir.create(out_dir_current, recursive = TRUE)
#'
#' ## Save current variables in temporary directory
#' terra::writeRaster(var, file.path(out_dir_current, "Variables.tif"))
#'
#'
#' # Step 2: Organize future climate variables (example with WorldClim)
#' ## Directory containing the downloaded future climate variables (example)
#' in_dir <- system.file("extdata", package = "kuenm2")
#'
#' ## Create a folder in a temporary directory to copy the future variables
#' out_dir_future <- file.path(tempdir(), "Future_raw_wc")
#'
#' ## Organize and rename the future climate data (structured by year and GCM)
#' ### 'SoilType' will be appended as a static variable in each scenario
#' organize_future_worldclim(input_dir = in_dir, output_dir = out_dir_future,
#'                           name_format = "bio_", fixed_variables = var$SoilType)
#'
#' # Step 3: Prepare data to run multiple projections
#' ## An example with maxnet models
#' ## Import example of fitted_models (output of fit_selected())
#' data(fitted_model_maxnet, package = "kuenm2")
#'
#' ## Prepare projection data using fitted models to check variables
#' pr <- prepare_projection(models = fitted_model_maxnet,
#'                          present_dir = out_dir_current,
#'                          future_dir = out_dir_future,
#'                          future_period = c("2081-2100"),
#'                          future_pscen = c("ssp126", "ssp585"),
#'                          future_gcm = c("ACCESS-CM2", "MIROC6"),
#'                          raster_pattern = ".tif*")
#'
#' # Step 4: Run multiple model projections
#' ## A folder to save projection results
#' out_dir <- file.path(tempdir(), "Projection_results/maxnet_projections")
#' dir.create(out_dir, recursive = TRUE)
#'
#' ## Project selected models to multiple scenarios
#' p <- project_selected(models = fitted_model_maxnet, projection_data = pr,
#'                       out_dir = out_dir)

project_selected <- function(models,
                             projection_data,
                             out_dir,
                             mask = NULL,
                             consensus_per_model = TRUE,
                             consensus_general = TRUE,
                             consensus = c("median", "range", "mean", "stdev"),
                             write_replicates = FALSE,
                             extrapolation_type = "E",
                             var_to_clamp = NULL,
                             type = NULL,
                             overwrite = FALSE,
                             parallel = FALSE,
                             ncores = NULL,
                             progress_bar = TRUE,
                             verbose = TRUE) {
  #Check data
  if (missing(models)) {
    stop("Argument 'models' must be defined.")
  }
  if (missing(projection_data)) {
    stop("Argument 'projection_data' must be defined.")
  }
  if (missing(out_dir)) {
    stop("Argument 'out_dir' must be defined.")
  }

  if (!inherits(models, "fitted_models")) {
    stop("Argument 'models' must be a 'fitted_models' object.")
  }
  if (!inherits(projection_data, "projection_data")) {
    stop("Argument 'projection_data' must be a 'projection_data' object.")
  }
  if (!inherits(out_dir, "character")) {
    stop("Argument 'out_dir' must be a 'character'.")
  }

  if (!is.null(mask) & !inherits(mask, c("SpatRaster", "SpatVector",
                                         "SpatExtent"))) {
    stop("Argument 'mask' must be a 'SpatVector', 'SpatExtent' or 'SpatRaster'.")
  }

  if (!inherits(consensus, "character")) {
    stop("Argument 'consensus' must be a 'character'.")
  }
  consensus_out <- setdiff(consensus, c("median", "range", "mean", "stdev"))
  if (length(consensus_out) > 0) {
    stop("Invalid 'consensus' provided.",
         "\nAvailable options are: 'median', 'range', 'mean' and 'stdev'.")
  }
  if (!any(c("median", "mean") %in% consensus)) {
    stop("'consensus' must contain at least one of the options: 'median' or 'mean'.")
  }

  if (extrapolation_type %in% c("EC", "NE") & !is.null(var_to_clamp) &
      !inherits(var_to_clamp, "character")) {
    stop("Argument 'var_to_clamp' must be NULL or 'character'.")
  }

  if(is.null(type)){
    if(models$algorithm == "maxnet") {
      type <- "cloglog"
    } else if (models$algorithm == "glm") {
      type <- "response"
    }
  }

  if(!inherits(type, "character")){
    stop("Argument 'type' must be NULL or 'character'.")
  }

  if(models$algorithm == "maxnet"){
    if (!any(c("raw", "cumulative", "logistic", "cloglog") %in% type)) {
      stop("Invalid 'type' provided.",
           "\nAvailable options for maxnet models are: 'raw', 'cumulative',
           'logistic', or 'cloglog'.")
    }
    if(type == "raw")
      type <-  "exponential"
  }

  if(models$algorithm == "glm"){
    if (!any(c("response", "cloglog") %in% type)) {
      stop("Invalid 'type' provided.",
           "\nAvailable options for glm models are 'response' or 'cloglog'.")
    }
    if(type == "cloglog")
      type = "link"
  }

  #Save parameters in a list to send to foreach nodes#
  par_list <- list(models = models,
                   mask = mask,
                   projection_data = projection_data,
                   consensus_per_model = consensus_per_model,
                   consensus_general = consensus_general,
                   consensus = consensus,
                   write_replicates = write_replicates,
                   extrapolation_type = extrapolation_type,
                   var_to_clamp = var_to_clamp,
                   type = type,
                   overwrite = overwrite)

  ####PREPARE DATAFRAME TO PREDICT####
  #Extract variables from best models
  vars <- names(models[["Models"]][[1]][[1]]$samplemeans)
  vars <- setdiff(vars, c("pr_bg", "fold"))

  if (!file.exists(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE)  # We should show these warnings
  }
  #Normalize path
  out_dir <- normalizePath(out_dir)

  #Check scenarios for predicting
  res_path <- check_pred_scenarios(projection_data = projection_data,
                                   out_dir = out_dir)
  raster_pattern <- projection_data$raster_pattern


  ####Configure parallelization####
  n_models <- nrow(res_path)
  if (n_models == 1 & isTRUE(parallel)) {
    parallel <- FALSE
  }

  #Show progress bar?
  if (progress_bar) {
    pb <- utils::txtProgressBar(0, n_models, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
  }


  ###################################

  #If parallel and mask, wrap variables
  if (!is.null(mask) & parallel) {
    par_list$mask <- terra::wrap(par_list$mask)
  }

  #Run predictions
  if (parallel) {
    if (is.null(ncores)) {
      ncores <- max(1, parallel::detectCores() - 1)
    }
    if (n_models < ncores) {
      ncores <- n_models
    }

    #Make cluster
    cl <- parallel::makeCluster(ncores)

    doSNOW::registerDoSNOW(cl)
    if (progress_bar) {
      opts <- list(progress = progress)
    } else {
      opts <- NULL
    }

    foreach(
      x = 1:n_models, .options.snow = opts,
      .export = c("predict_selected", "multiple_projections")
    ) %dopar% {
      multiple_projections(i = x, res_path, raster_pattern, par_list)
    }
  } else { #Not in parallel (using %do%)
    # Loop for com barra de progresso manual
    for (x in 1:n_models) {
      # Execute a função fit_eval_models
      multiple_projections(i = x, res_path = res_path,
                           raster_pattern = raster_pattern,
                           par_list = par_list)

      # Sets the progress bar to the current state
      if (progress_bar) {
        utils::setTxtProgressBar(pb, x)
      }
    }
  }

  #Stop cluster
  if (parallel) {
    parallel::stopCluster(cl)
  }

  #Append threshold to final results
  res_final <- new_model_projections(paths = res_path,
                                     thresholds = models$thresholds)

  #Save
  saveRDS(res_final, file.path(out_dir, "Projection_paths.RDS"))
  return(res_final)
}
