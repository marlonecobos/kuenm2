#' Compute changes of suitable areas in other scenarios (single scenario / GCM)
#'
#' @usage
#' prediction_changes(current_predictions, new_predictions,
#'                    predicted_to = "future", fitted_models = NULL,
#'                    consensus = "mean", user_threshold = NULL,
#'                    force_resample = FALSE, gain_color = "#009E73",
#'                    loss_color = "#D55E00", stable_suitable = "#0072B2",
#'                    stable_unsuitable = "grey", write_results = FALSE,
#'                    output_dir = NULL, overwrite = FALSE,
#'                    write_bin_models = FALSE)
#'
#' @param current_predictions (SpatRaster) A `SpatRaster` object returned by
#'  `predict_selected()` with suitability predicted under current conditions.
#' @param new_predictions (SpatRaster) A `SpatRaster` object returned by
#' `predict_selected()` with suitability predicted under new conditions
#' (past or future). This SpatRaster must have the same resolution and extent
#' as `current_predictions`.
#' @param predicted_to (character) a string specifying whether `new_predictions`
#' represent "past" or "future" conditions. Default is "future".
#' @param fitted_models an object of class `fitted_models` returned by
#' [fit_selected()]
#' @param consensus (character) the consensus metric stored in `fitted_models`
#' used to binarize models. Available options are `"mean"`, `"median"`,
#' `"range"`, and `"stdev"` (standard deviation). Default is `"mean"`.
#' @param user_threshold (numeric) an optional threshold for binarizing predictions.
#' Default is `NULL`, meaning the function will apply the thresholds stored in
#' `model_projections`, which were calculated earlier using the omission rate
#' from [calibration()].
#' @param force_resample (logical) whether to force rasters to have the same
#' extent and resolution. Default is `TRUE`.
#' @param gain_color (character) color used to represent gains. Default is
#' "#009E73" (teal green).
#' @param loss_color (character) color used to represent losses. Default is
#' "#D55E00" (orange-red).
#' @param stable_suitable (character) color used for representing areas that
#' remain suitable across scenarios. Default is "#0072B2" (oxford blue).
#' @param stable_unsuitable (character) color used for representing areas that
#' remain unsuitable across scenarios. Default is "grey".
#' @param write_results (logical) whether to save the results to disk. Default
#' is FALSE.
#' @param output_dir (character) directory path where results will be saved.
#' Only relevant if `write_results = TRUE`.
#' @param overwrite (logical) whether to overwrite SpatRasters if they already
#' exist. Only applicable if `write_results = TRUE`. Default is FALSE.
#' @param write_bin_models (logical) whether to write the binarized models for
#' each scenario to the disk. Only applicable if `write_results = TRUE`. Default
#' is FALSE.
#'
#' @details
#' When projecting a niche model to different temporal scenarios (past or
#' future), species’ areas can be classified into three categories relative to
#' the current baseline: **gain**, **loss** and **stability**. The
#' interpretation of these categories depends on the temporal direction of the projection.
#' **When projecting to future scenarios**:
#' - *Gain*: Areas that are currently unsuitable become suitable in the future.
#' - *Loss*: Areas that are currently suitable become unsuitable in the future.
#' - *Stability*: Areas that retain their current classification in the future,
#' whether suitable or unsuitable.
#'
#' **When projecting to past scenarios**:
#' - *Gain*: Areas that were unsuitable in the past are now suitable in the
#' present.
#' - *Loss*: Areas that were suitable in the past are now unsuitable in the
#' present.
#' - *Stability*: Areas that retain their past classification in the present,
#' whether suitable or unsuitable.
#'
#' @returns
#' A `SpatRaster` showing the areas of **gain**, **loss** and **stability**.
#'
#' @export
#'
#' @importFrom terra ext res resample coltab writeRaster
#'
#' @examples
#' # Import an example of fitted models (output of fit_selected())
#' data("fitted_model_maxnet", package = "kuenm2")
#'
#'
#' # Import current variables for prediction
#' present_var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                        package = "kuenm2"))
#'
#'
#' # Import variables for a single future scenario for prediction
#' future_var <- terra::rast(system.file("extdata",
#'                                       "wc2.1_10m_bioc_ACCESS-CM2_ssp585_2081-2100.tif",
#'                                       package = "kuenm2"))
#'
#' # Rename variables to match the variable names used in the fitted models
#' names(future_var) <- sub("bio0", "bio", names(future_var))
#' names(future_var) <- sub("bio", "bio_", names(future_var))
#'
#' # Append the static soil variable to the future variables
#' future_var <- c(future_var, present_var$SoilType)
#'
#' # Predict under present and future conditions
#' p_present <- predict_selected(models = fitted_model_maxnet,
#'                               new_variables = present_var)
#' p_future <- predict_selected(models = fitted_model_maxnet,
#'                              new_variables = future_var)
#'
#'
#' # Compute changes between scenarios
#' p_changes <- prediction_changes(current_predictions = p_present$General_consensus$mean,
#'                                 new_predictions = p_future$General_consensus$mean,
#'                                 fitted_models = fitted_model_maxnet,
#'                                 predicted_to = "future")
#' # Plot result
#' terra::plot(p_changes)

prediction_changes <- function(current_predictions,
                               new_predictions,
                               predicted_to = "future",
                               fitted_models = NULL,
                               consensus = "mean",
                               user_threshold = NULL,
                               force_resample = FALSE,
                               gain_color = "#009E73",
                               loss_color = "#D55E00",
                               stable_suitable = "#0072B2",
                               stable_unsuitable = "grey",
                               write_results = FALSE,
                               output_dir = NULL,
                               overwrite = FALSE,
                               write_bin_models = FALSE) {
  #### Check data ####
  if (missing(current_predictions)) {
    stop("Argument 'current_predictions' must be defined.")
  }
  if (!inherits(current_predictions, "SpatRaster")) {
    stop("Argument 'current_predictions' must be a 'SpatRaster' object.")
  }

  if (missing(new_predictions)) {
    stop("Argument 'new_predictions' must be defined.")
  }
  if (!inherits(new_predictions, "SpatRaster")) {
    stop("Argument 'new_predictions' must be a 'SpatRaster' object.")
  }

  if (!inherits(predicted_to, "character") || length(predicted_to) > 1) {
    stop("Argument 'predicted_to' must be a single 'character'.")
  }
  if (!predicted_to %in% c("past", "future")){
    stop("Argument 'predicted_to' must be 'past' or 'future'")
  }

  if(!is.null(fitted_models)){
    if (!inherits(fitted_models, "fitted_models")) {
      stop("Argument 'fitted_models' must be a 'fitted_models' object.")
    }

    if (!inherits(consensus, "character")) {
      stop("Argument 'consensus' must be a 'character'.")
    }
    if (length(consensus) > 1) {
      stop("Argument 'consensus' must be a unique value.",
           "\nAvailable options are: 'median', 'range', 'mean' or 'stdev'.")
    }
    consensus_out <- setdiff(consensus, c("median", "range", "mean", "stdev"))
    if (length(consensus_out) > 0) {
      stop("Invalid 'consensus' provided.",
           "\nAvailable options are: 'median', 'range', 'mean' or 'stdev'.")
    }

  }

  if (!is.null(user_threshold) & !inherits(user_threshold, "numeric")) {
    stop("Argument 'user_threshold' must be NULL or 'numeric'.")
  }

  if(is.null(user_threshold) && is.null(fitted_models)){
    stop("You must provide a 'fitted_models' or an 'user_threshold'")
  }

  if (!inherits(force_resample, "logical")) {
    stop("Argument 'force_resample' must be 'logical' (TRUE or FALSE)")
  }

  if (!inherits(gain_color, "character") || length(gain_color) > 1) {
    stop("Argument 'gain_color' must be a single 'character' value.")
  }
  if (!inherits(loss_color, "character") || length(loss_color) > 1) {
    stop("Argument 'loss_color' must be a single 'character' value.")
  }
  if (!inherits(stable_suitable, "character") || length(stable_suitable) > 1) {
    stop("Argument 'stable_suitable' must be a single 'character' value.")
  }
  if (!inherits(stable_unsuitable, "character") || length(stable_unsuitable) > 1) {
    stop("Argument 'stable_unsuitable' must be a single 'character' value.")
  }

  if (!inherits(write_results, "logical")) {
    stop("Argument 'write_results' must be 'logical' (TRUE or FALSE)")
  }

  if(write_results){
    if(is.null(output_dir)){
      stop("If 'write_files = TRUE', 'output_dir' must be specified")
    }

    if (!inherits(output_dir, "character")) {
      stop("Argument 'output_dir' must be a 'character'.")
    }

    if(!file.exists(output_dir)){
      stop("Directory specified in 'output_dir' does not exist.")
    }

    if(!inherits(overwrite, "logical")){
      stop("Argument 'overwrite' must be 'logical'.")
    }

    if(!inherits(write_bin_models, "logical")){
      stop("Argument 'write_bin_models' must be 'logical'.")
    }
  }


  #### Extract threshold ####
  if (is.null(user_threshold)) {
    thr <- fitted_models[["thresholds"]][["consensus"]][[consensus]]
  } else {
    thr <- user_threshold
  }

  # Check resolution
  if(sum(terra::res(current_predictions) == terra::res(new_predictions)) != 2){
    stop("SpatRasters provided in 'current_predictions' and 'new_predictions' must have the same resolution")
  }

  # Check extent
  if(terra::ext(current_predictions) != terra::ext(new_predictions)){
    if(force_resample){
      new_predictions <- terra::resample(new_predictions,
                                         current_predictions,
                                         method = "average")
    } else {
      stop("The SpatRaster in 'new_predictions' has a different extent than 'current_predictions'. ",
           "Please check the raster extents or set 'force_resample = TRUE' to resample automatically.")
    }
  }



  #### Binarize models ####
  # Current predictions
  r_bin <- binarize_values(x = current_predictions, v = thr, new_value = 2)
  #Set levels
  levels(r_bin) <- data.frame(id = c(0, 2),
                              Result = c("Unsuitable", "Suitable"))
  names(r_bin) <- "Present"

  # Projected predictions
  r_change <- binarize_values(x = new_predictions, v = thr, new_value = 1)
  levels(r_change) <- data.frame(id = c(0, 1),
                                 Result = c("Unsuitable", "Suitable"))


  #### Identify changes ####
  #Table to set levels in raster
  if(predicted_to == "future"){
    #Table to set levels in raster
    cls <- data.frame(id = c(1, 2, 3, 0),
                             Result = c("Gain", "Loss", "Stable suitable",
                                        "Stable unsuitable"))
    # Tab to set colors
    df_colors <- data.frame(value = 0:3,
                            col = c(stable_unsuitable, gain_color,
                                    loss_color, stable_suitable))

  } else if(predicted_to == "past") {
    #Table to set levels in raster
    cls <- data.frame(id = c(1, 2, 3, 0),
                           Result = c("Loss", "Gain", "Stable suitable",
                                      "Stable unsuitable"))
    # Tab to set colors
    df_colors <- data.frame(value = 0:3,
                            col = c(stable_unsuitable, loss_color,
                                    gain_color, stable_suitable))
  }

  #Compute changes
  r_result <- r_bin + r_change
  #Get legend
  levels(r_result) <- cls

  # Set colors
  terra::coltab(r_result) <- df_colors


  #### Write results? ####
  if(write_results){
    terra::writeRaster(r_result,
                       paste0(output_dir, "/Changes_summary.tif"),
                       overwrite = overwrite)

    if(write_bin_models){
      terra::writeRaster(r_bin,
                         paste0(output_dir, "/Current_predictions_binarized.tif"),
                         overwrite = overwrite)
      terra::writeRaster(r_change,
                         paste0(output_dir, "/New_predictions_binarized.tif"),
                         overwrite = overwrite)
    }

  }


  # Return result
  return(r_result)

  }


# # Objects of the function
# current_predictions = p_present$General_consensus$mean
# projected_predictions = p_future$General_consensus$mean
# projected_to = "future"
# fitted_models = fm
# user_threshold = NULL
# force_resample = TRUE
# write_results = FALSE
# output_dir = NULL
# overwrite = FALSE
# write_bin_models = FALSE
