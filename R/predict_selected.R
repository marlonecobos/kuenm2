#' Predict selected models for a single scenario
#'
#' @description
#' This function predicts selected models for a single environmental scenario
#' using either `maxnet` or `glm` models. It provides options for saving the
#' output and calculating consensus measures (mean, median, etc.) across
#' replicates and models.
#'
#' @usage predict_selected(models, raster_variables, write_files = FALSE,
#'                         write_replicates = FALSE, out_dir = NULL,
#'                         consensus_per_model = TRUE, consensus_general = TRUE,
#'                         consensus = c("median", "range", "mean", "stdev"),
#'                         clamping = FALSE, var_to_clamp = NULL,
#'                         type = "cloglog", overwrite = FALSE,
#'                         progress_bar = TRUE)
#'
#' @param models an object of class `fitted_models` returned by the
#' \code{\link{fit_selected}}() function.
#' @param raster_variables a SpatRaster or data.frame of predictor variables. The names
#' of these variables must match those used to calibrate the models or those
#' used to run PCA if `do_pca = TRUE` in the \code{\link{prepare_data}}()
#' function.
#' @param mask (SpatRaster, SpatVector, or SpatExtent) spatial object used to
#'        mask the variables before predict. Default is NULL.
#' @param write_files (logical) whether to save the predictions (SpatRasters)
#' to disk. Default is FALSE.
#' @param write_replicates (logical) whether to save the predictions for each
#' replicate to disk. Only applicable if `write_files` is TRUE. Default is
#' FALSE.
#' @param out_dir (character) directory path where predictions will be saved.
#' Only relevant if `write_files = TRUE`.
#' @param consensus_per_model (logical) whether to compute consensus (mean,
#' median, etc.) for each model across its replicates. Default is TRUE.
#' @param consensus_general (logical) whether to compute a general consensus
#' across all models. Default is TRUE.
#' @param consensus (character) vector specifying the types of consensus to
#' calculate across replicates and models. Available options are `"median"`,
#' `"range"`, `"mean"`, and `"stdev"` (standard deviation). Default is
#' `c("median", "range", "mean", "stdev")`.
#' @param clamping (logical) whether to restrict variable values to the range
#' of the calibration data to avoid extrapolation. Default is FALSE (free
#' extrapolation).
#' @param var_to_clamp (character) vector specifying which variables to clamp.
#' Only applicable if `clamping = TRUE`. Default is `NULL`, meaning all
#' variables will be clamped.
#' @param type (character) the format of prediction values. For `maxnet` models,
#' available options are `"raw"`, `"cumulative"`, `"logistic"`, and `"cloglog"`,
#' with the default being `"cloglog"` for maxnet and `"response"` for glm.
#' @param overwrite (logical) whether to overwrite SpatRasters if they already
#' exist. Only applicable if `write_files = TRUE`. Default is FALSE.
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is TRUE.
#'
#' @return A list containing SpatRaster predictions for each replicate, along
#' with the consensus results for each model and the overall general consensus.
#'
#' @importFrom terra rast clamp predict median mean stdev diff writeRaster
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats glm
#' @export
#'
#' @examples
#' #Load package
#' library(terra)
#' # Import variables to predict
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#' ##Example with maxnet
#' # Import example of fitted_models (output of fit_selected())
#' data("fitted_model_maxnet", package = "kuenm2")
#' # Predict to single scenario
#' p <- predict_selected(models = fitted_model_maxnet,
#'                       raster_variables = var,
#'                       write_files = FALSE,
#'                       write_replicates = FALSE,
#'                       out_dir = NULL,
#'                       consensus_per_model = TRUE,
#'                       consensus_general = TRUE,
#'                       consensus = c("median", "range", "mean", "stdev"),
#'                       clamping = FALSE,
#'                       var_to_clamp = NULL,
#'                       type = "cloglog",
#'                       overwrite = FALSE,
#'                       progress_bar = TRUE)
#' ##Example with GLMN
#' # Import example of fitted_models (output of fit_selected())
#' data("fitted_model_glm", package = "kuenm2")
#' # Predict to single scenario
#' p_glm <- predict_selected(models = fitted_model_glm,
#'                           raster_variables = var,
#'                           write_files = FALSE,
#'                           write_replicates = FALSE,
#'                           out_dir = NULL,
#'                           consensus_per_model = TRUE,
#'                           consensus_general = TRUE,
#'                           consensus = c("median", "range", "mean", "stdev"),
#'                           clamping = FALSE,
#'                           var_to_clamp = NULL,
#'                           type = "cloglog",
#'                           overwrite = FALSE,
#'                           progress_bar = TRUE)
#' #Compare
#' plot(c(p$General_consensus$mean,
#'        p_glm$General_consensus$mean),
#'      col = rev(terrain.colors(240)), main = c("MAXNET", "GLM"),
#'      zlim = c(0, 1))
#'
#'

predict_selected <- function(models,
                             raster_variables,
                             mask = NULL,
                             write_files = FALSE,
                             write_replicates = FALSE,
                             out_dir = NULL,
                             consensus_per_model = TRUE,
                             consensus_general = TRUE,
                             consensus = c("median", "range", "mean", "stdev"),
                             clamping = FALSE,
                             var_to_clamp = NULL,
                             type = "cloglog",
                             overwrite = FALSE,
                             progress_bar = TRUE) {

  #Check data
  if (!inherits(models, "fitted_models")) {
    stop(paste0("Argument models must be a fitted_models object, not ",
                class(models)))
  }
  if (!inherits(raster_variables, "SpatRaster")) {
    stop(paste0("Argument raster_variables must be a SpatRaster, not ",
                class(raster_variables)))
  }

  if(!is.null(mask) & !inherits(mask, c("SpatRaster", "SpatVector",
                                        "SpatExtent"))){
    stop(paste0("Argument mask must be a SpatVector, SpatExtent or SpatRaster, not ",
                class(mask)))
  }

  if (!inherits(write_files, "logical")) {
    stop(paste0("Argument write_files must be logical, not ",
                class(write_files)))
  }
  if(write_files & is.null(out_dir)){
    stop("If write_files = TRUE, out_dir must be specified")
  }
  if(write_files & !inherits(out_dir, "character")) {
    stop(paste0("Argument out_dir must be a character, not ",
                class(out_dir)))
  }
  if(write_files & !inherits(write_replicates, "logical")) {
    stop(paste0("Argument write_replicates must be logical, not ",
                class(write_replicates)))
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
  if (!inherits(progress_bar, "logical")) {
    stop(paste0("Argument progress_bar must be logical, not ",
                class(progress_bar)))}

  if(!is.null(mask)){
    raster_variables <- terra::crop(raster_variables, mask, mask = TRUE)
  }

  if(!is.null(models$pca)){
    if(!("vars_out" %in% names(models$pca))) {
      raster_variables <- terra::predict(raster_variables[[models$pca$vars_in]],
                                         models$pca)
    } else {
      raster_variables <- c(terra::predict(raster_variables[[models$pca$vars_in]],
                                           models$pca),
                            raster_variables[[models$pca$vars_out]])
    }
  }

  # Extract info from fitted object
  categorical_layers <- models[["categorical_variables"]]
  models <- models[["Models"]]
  nm <- names(models)

  # Get names of the models (Replicates or full model)
  names_models <- unlist(unique(sapply(nm, function(i) {
    names(models[[i]])
  }, USE.NAMES = FALSE, simplify = FALSE)))

  # If there are replicates, remove the full model from the dataset
  if (any(grepl("Rep", names_models))) {
    models <- lapply(nm, function(i) {
      models[[i]][["Full_model"]] <- NULL
      return(models[[i]])
    })
    names(models) <- nm
    names_models <- unlist(unique(sapply(nm, function(i) {
      names(models[[i]])
    }, USE.NAMES = FALSE, simplify = FALSE)))
  }


  # Clamp variables if required
  # To do:
  #  - Add a warning if the variable is not in the calibration data

  if (clamping) {
    varmin <- models[[1]][[1]]$varmin[-1]  # Get var min
    varmax <- models[[1]][[1]]$varmax[-1]  # Get var max
    if (is.null(var_to_clamp)) {
      var_to_clamp <- setdiff(names(varmin), c("pr_bg", "fold"))
    }
    clamped_variables <- terra::rast(lapply(var_to_clamp, function(i) {
      terra::clamp(raster_variables[[i]], lower = varmin[i], upper = varmax[i], values = TRUE)
    }))
    raster_variables[[names(clamped_variables)]] <- clamped_variables
  }

  # Prepare progress bar
  n_models <- length(models)
  if (progress_bar) {
    pb <- utils::txtProgressBar(0, n_models, style = 3)
  }

  #Create empty list
  p_models <- list()

  # Fill the list with predictions
  for (i in seq_along(models)) {
    inner_list <- list()

    for (x in models[[i]]) {

      if (inherits(raster_variables, "SpatRaster")) {

        #  Convert Layers to Factors in SpatRaster
        if (!is.null(categorical_layers) && length(categorical_layers) > 0) {
          for (layer in categorical_layers) {
            if (layer %in% names(raster_variables)) {  # Check if the layer exists
              raster_variables[[layer]] <- terra::as.factor(raster_variables[[layer]])
            } else {
              warning(paste("Layer", layer, "not found in the SpatRaster."))
            }
          }
        }

        if (inherits(x, "glmnet")) {
          prediction <- terra::predict(raster_variables, x, na.rm = TRUE,
                                       type = type, fun = predict.glmnet_mx)
        } else if (inherits(x, "glm")) {
          prediction <- terra::predict(raster_variables, x, na.rm = TRUE,
                                       type = "response")
        }
      } else if (inherits(raster_variables, "data.frame")) {
        if (inherits(x, "glmnet")) {
          prediction <- as.numeric(predict.glmnet_mx(x, newdata = raster_variables,
                                                     type = type))
        } else if (inherits(x, "glm")) {
          prediction <- as.numeric(predict(x, newdata = raster_variables,
                                           type = "response"))
        }
      }
      inner_list[[length(inner_list) + 1]] <- prediction
    }

    if (inherits(raster_variables, "SpatRaster")) {
      p_models[[length(p_models) + 1]] <- terra::rast(inner_list)
    } else if (inherits(raster_variables, "data.frame")) {
      p_models[[length(p_models) + 1]] <- inner_list
    }

    # Update progress bar
    if (progress_bar) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  # Rename models and replicates
  names(p_models) <- nm
  for (i in nm) {
    names(p_models[[i]]) <- names(models[[i]])
  }

  # Create an empty list to store results
  res <- list()

  # Get consensus by model if raster_variables is a SpatRaster
  if (inherits(raster_variables, "SpatRaster")) {
    rep <- unlist(p_models)

    if (consensus_per_model) {
      if ("median" %in% consensus) {
        res$Consensus_per_model$median <- terra::rast(lapply(p_models, terra::median))
      }
      if ("mean" %in% consensus) {
        res$Consensus_per_model$mean <- terra::rast(lapply(p_models, terra::mean))
      }
      if ("stdev" %in% consensus) {
        res$Consensus_per_model$stdev <- terra::rast(lapply(p_models, terra::stdev))
      }
      if ("range" %in% consensus) {
        res$Consensus_per_model$range <- terra::rast(lapply(p_models, function(r) {
          terra::diff(range(r))
        }))
      }
    }

    gen_res <- list()
    if (consensus_general && length(p_models) == 1 && consensus_per_model) {
      if ("median" %in% consensus) {
        gen_res$median <- res$Consensus_per_model$median
      }
      if ("mean" %in% consensus) {
        gen_res$mean <- res$Consensus_per_model$mean
      }
      if ("stdev" %in% consensus) {
        gen_res$stdev <- res$Consensus_per_model$stdev
      }
      if ("range" %in% consensus) {
        gen_res$range <- res$Consensus_per_model$range
      }
    } else if (consensus_general && length(p_models) > 1) {
      all_rep <- terra::rast(p_models)

      if ("median" %in% consensus) {
        gen_res$median <- terra::median(all_rep)
      }
      if ("mean" %in% consensus) {
        gen_res$mean <- terra::mean(all_rep)
      }
      if ("stdev" %in% consensus) {
        gen_res$stdev <- terra::stdev(all_rep)
      }
      if ("range" %in% consensus) {
        gen_res$range <- terra::diff(range(all_rep))
      }
    }

    res <- lapply(1:length(nm), function(x) {
      mcs <- lapply(consensus, function(y) {
        res$Consensus_per_model[[y]][[x]]
      })
      mcs <- terra::rast(mcs)
      names(mcs) <- consensus
      list(Replicates = rep[[x]], Model_consensus = mcs)
    })

    names(res) <- nm
    res <- c(res, General_consensus = terra::rast(gen_res))
  }

  # Write results to disk if required
  if (write_files) {
    if (!file.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }

    sapply(nm, function(i) {
      if (write_replicates) {
        terra::writeRaster(res[[i]]$Replicates,
                           file.path(out_dir, paste0(i, "_replicates.tiff")),
                           overwrite = overwrite)
      }
      terra::writeRaster(res[[i]]$Model_consensus,
                         file.path(out_dir, paste0(i, "_consensus.tiff")),
                         overwrite = overwrite)
    })

    terra::writeRaster(res$General_consensus,
                       file.path(out_dir, "General_consensus.tiff"),
                       overwrite = overwrite)
  }

  return(res)
} # End of function
