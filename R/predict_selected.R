#' Predict selected models for a single scenario
#'
#' @description
#' This function predicts selected models for a single set of new data
#' using either `maxnet` or `glm` It provides options to save the
#' output and compute consensus results (mean, median, etc.) across
#' replicates and models.
#'
#' @usage
#' predict_selected(models, raster_variables, mask = NULL, write_files = FALSE,
#'                  write_replicates = FALSE, out_dir = NULL,
#'                  consensus_per_model = TRUE, consensus_general = TRUE,
#'                  consensus = c("median", "range", "mean", "stdev"),
#'                  extrapolation_type = "E", var_to_clamp = NULL,
#'                  type = "cloglog", overwrite = FALSE, progress_bar = TRUE)
#'
#' @param models an object of class `fitted_models` returned by the
#' \code{\link{fit_selected}}() function.
#' @param raster_variables a SpatRaster or data.frame of predictor variables.
#' The names of these variables must match those used to calibrate the models or
#' those used to run PCA if `do_pca = TRUE` in the \code{\link{prepare_data}}()
#' function.
#' @param mask (SpatRaster, SpatVector, or SpatExtent) spatial object used to
#' mask the variables before predict. Default is NULL.
#' @param write_files (logical) whether to save the predictions (SpatRasters or
#' data.frame) to disk. Default is FALSE.
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
#' @param extrapolation_type (character) extrapolation type of model. Models can
#' be transferred with three options: free extrapolation ('E'), extrapolation
#' with clamping ('EC'), and no extrapolation ('NE'). Default = 'E'. See details.
#' @param var_to_clamp (character) vector specifying which variables to clamp or
#' not extrapolate. Only applicable if extrapolation_type is "EC" or "NE".
#' Default is `NULL`, meaning all variables will be clamped or not extrapolated.
#' @param type (character) the format of prediction values. For `maxnet` models,
#' valid options are `"raw"`, `"cumulative"`, `"logistic"`, and `"cloglog"`.
#' For `glm` models, valid options are `"cloglog"`, `"response"`, `"raw"`,
#' `"cumulative"` and `"link"`. Default is "cloglog.
#' @param overwrite (logical) whether to overwrite SpatRasters if they already
#' exist. Only applicable if `write_files = TRUE`. Default is FALSE.
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is TRUE.
#'
#' @details
#' When predicting to areas where the variables are beyond the lower or upper
#' limits of the calibration data, users can choose to free extrapolate the
#' predictions (`extrapolation_type = "E"`), extrapolate with clamping
#' (extrapolation_type = "EC"), or not extrapolate (extrapolation_type = "NE").
#' When clamping, the variables are set to minimum and maximum values
#' established for the maximum and minimum values within calibration data. In
#' the no extrapolation approach, any cell with at least one variable listed in
#' `var_to_clamp` falling outside the calibration range is assigned a suitability
#' value of 0.
#'
#' @return
#' A list containing SpatRaster or data.frames predictions for each replicate,
#' long with the consensus results for each model and the overall general consensus.
#'
#' @importFrom terra rast clamp predict median mean stdev diff writeRaster cells
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @importFrom stats glm predict median sd
#'
#' @export
#'
#' @examples
#' # Import variables to predict on
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' # Example with maxnet
#' # Import example of fitted_models (output of fit_selected())
#' data("fitted_model_maxnet", package = "kuenm2")
#'
#' # Predict to single scenario
#' p <- predict_selected(models = fitted_model_maxnet, raster_variables = var)
#'
#' # Example with GLMs
#' # Import example of fitted_models (output of fit_selected())
#' data("fitted_model_glm", package = "kuenm2")
#'
#' # Predict to single scenario
#' p_glm <- predict_selected(models = fitted_model_glm, raster_variables = var)
#'
#' # Plot predictions
#' terra::plot(c(p$General_consensus$median, p_glm$General_consensus$median),
#'             col = rev(terrain.colors(240)), main = c("MAXNET", "GLM"),
#'             zlim = c(0, 1))

predict_selected <- function(models,
                             raster_variables,
                             mask = NULL,
                             write_files = FALSE,
                             write_replicates = FALSE,
                             out_dir = NULL,
                             consensus_per_model = TRUE,
                             consensus_general = TRUE,
                             consensus = c("median", "range", "mean", "stdev"),
                             extrapolation_type = "E",
                             var_to_clamp = NULL,
                             type = "cloglog",
                             overwrite = FALSE,
                             progress_bar = TRUE) {

  #Check data
  if (missing(models)) {
    stop("Argument 'models' must be defined.")
  }
  if (missing(raster_variables)) {
    stop("Argument 'raster_variables' must be defined.")
  }
  if (!inherits(models, "fitted_models")) {
    stop("Argument 'models' must be a 'fitted_models' object.")
  }
  if (!inherits(raster_variables, c("SpatRaster", "data.frame"))) {
    stop("Argument 'raster_variables' must be a 'SpatRaster' o 'data.frame")
  }

  if (!is.null(mask) & !inherits(mask, c("SpatRaster", "SpatVector",
                                         "SpatExtent"))) {
    stop("Argument 'mask' must be a 'SpatVector', 'SpatExtent' or 'SpatRaster'.")
  }

  if (write_files & is.null(out_dir)) {
    stop("If write_files = TRUE, out_dir must be specified")
  }
  if (write_files & !inherits(out_dir, "character")) {
    stop("Argument 'out_dir' must be a 'character'.")
  }
  if (!inherits(consensus, "character")) {
    stop("Argument 'consensus' must be a 'character'.")
  }
  consensus_out <- setdiff(consensus, c("median", "range", "mean", "stdev"))
  if (length(consensus_out) > 0) {
    stop("Invalid 'consensus' provided.",
         "\nAvailable options are: 'median', 'range', 'mean' and 'stdev'.")
  }

  if(length(extrapolation_type) > 1){
    stop("Extrapolation type accepts only one of these values: 'E', 'EC', or
         'NE'")
  }

  extrapolation_out <- setdiff(extrapolation_type, c("E", "EC", "NE"))
  if (length(extrapolation_out) > 0) {
    stop("Invalid 'extrapolation type' provided.",
         "\nAvailable options are: 'E', 'EC', and 'NE'.")
  }

  if (extrapolation_type %in% c("E", "EC") & !is.null(var_to_clamp) &
      !inherits(var_to_clamp, "character")) {
    stop("Argument 'var_to_clamp' must be NULL or 'character'.")
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
    if (!any(c("response", "cloglog", "cumulative", "link", "raw") %in% type)) {
      stop("Invalid 'type' provided.",
           "\nAvailable options for glm models are 'response', or 'cloglog', or 'cumulative', or 'link', or 'raw'.")
    }
    }

  # Analyses start here
  if (!is.null(mask) & inherits(raster_variables, "SpatRaster")) {
    raster_variables <- terra::crop(raster_variables, mask, mask = TRUE)
  }

  if (!is.null(models$pca)) {
    if(inherits(raster_variables, "SpatRaster")){
      if (!("vars_out" %in% names(models$pca))) {
      raster_variables <- terra::predict(raster_variables[[models$pca$vars_in]],
                                         models$pca)
    } else {
      raster_variables <- c(terra::predict(raster_variables[[models$pca$vars_in]],
                                           models$pca),
                            raster_variables[[models$pca$vars_out]])
    }
    }

    if(inherits(raster_variables, "data.frame")){
      if (!("vars_out" %in% names(models$pca))) {
        raster_variables <- stats::predict(raster_variables[[models$pca$vars_in]],
                                     models$pca)
      } else {
        raster_variables <- c(stats::predict(raster_variables[[models$pca$vars_in]],
                                             models$pca),
                              raster_variables[[models$pca$vars_out]])
      }
    }
  }

  # Extract info from fitted object
  categorical_layers <- models[["categorical_variables"]]
  models <- models[["Models"]]
  nm <- names(models)
  nrep <- length(models[[1]])

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

  if (extrapolation_type == "EC") {
    varmin <- models[[1]][[1]]$varmin[-1]  # Get var min
    varmax <- models[[1]][[1]]$varmax[-1]  # Get var max
    if (is.null(var_to_clamp)) {
      var_to_clamp <- setdiff(names(varmin), c("pr_bg", "fold"))
    }

    if(inherits(raster_variables, "SpatRaster")){
      clamped_variables <- terra::rast(lapply(var_to_clamp, function(i) {
        terra::clamp(raster_variables[[i]], lower = varmin[i],
                     upper = varmax[i], values = TRUE)
    }))
    raster_variables[[names(clamped_variables)]] <- clamped_variables
    } else if (inherits(raster_variables, "data.frame")){
      # Apply clamp for each column
      for (col_name in var_to_clamp) {
        # Clampar lower values
        raster_variables[[col_name]][raster_variables[[col_name]] < varmin[col_name]] <- varmin[col_name]
        # Clamp higher values
        raster_variables[[col_name]][raster_variables[[col_name]] > varmax[col_name]] <- varmax[col_name]
      }
    } #End of if raster_variables is data.frame
  } #End of if extrapolation_type == "EC"

  if(extrapolation_type == "NE"){
    varmin <- models[[1]][[1]]$varmin[-1]  # Get var min
    varmax <- models[[1]][[1]]$varmax[-1]  # Get var max
    if (is.null(var_to_clamp)) {
      var_to_clamp <- setdiff(names(varmin), c("pr_bg", "fold"))
    }

    if(inherits(raster_variables, "SpatRaster")){
      #Idenfity cells to not extrapolate
      no_extrapolate <- lapply(var_to_clamp, function(i){
        i_min <- terra::cells(x = (raster_variables[[i]] < varmin[i]) * 1, y = 1)
        i_max <- terra::cells(x = (raster_variables[[i]] > varmax[i]) * 1, y = 1)
        return(unique(i_min, i_max))
      })

      } else if(inherits(raster_variables, "data.frame")){
        no_extrapolate <- lapply(var_to_clamp, function(i){
          i_min <- which(raster_variables[,i] < varmin[i])
          i_max <- which(raster_variables[,i] > varmax[i])
          return(unique(i_min, i_max))
        })
    } #End of if raster_variables is data.frame

    #Get unique cells to no extrapolate
    no_extrapolate <- unique(unlist(no_extrapolate))
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
                                       type = type, fun = predict_glm_mx)
        }
      } else if (inherits(raster_variables, "data.frame")) {
        if (inherits(x, "glmnet")) {
          prediction <- as.numeric(predict.glmnet_mx(x, newdata = raster_variables,
                                                     type = type))
        } else if (inherits(x, "glm")) {
          prediction <- as.numeric(predict_glm_mx(x, newdata = raster_variables,
                                   type = type))
        }
      }

      #If no extrapolate, set cells beyond limits to 0
      if(extrapolation_type == "NE"){
        if(length(no_extrapolate) > 0){
          prediction[no_extrapolate] <- 0
          if(type == "cumulative") { #Recalculate cumulative values
            #For spatraster
            if(inherits(raster_variables, "SpatRaster")){
              numeric_predictions <- terra::values(prediction, na.rm = TRUE)
              prediction[!is.na(prediction)] <- cumulative_predictions(numeric_predictions) }
            #For data.frame
            if(inherits(raster_variables, "data.frame")){
              prediction <- cumulative_predictions(prediction)
            }
          } #End of type = cumulative
        } #End of length(extrapolate) > 0
        } #End of extrapolation_type == "NE"

      #Save in list
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

  #### Get consensus by model if raster_variables is a SpatRaster ####
  if (inherits(raster_variables, "SpatRaster")) {
    rep <- unlist(p_models)

    if (nrep == 1) {
      res$Consensus_per_model$Full_model <- lapply(p_models, function(x) {x})
    } else {
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
    }

    gen_res <- list()
    if (consensus_general && length(p_models) == 1 && consensus_per_model) {
      if (nrep == 1) {
        gen_res$Full_model <- res$Consensus_per_model$Full_model
      } else {
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
      if (nrep == 1) {
        mcs <- lapply("Full_model", function(y) {
          res$Consensus_per_model[[y]][[x]]
        })
        mcs <- terra::rast(mcs)
        names(mcs) <- "Full_model"
        list(Model_consensus = mcs)
      } else {
        mcs <- lapply(consensus, function(y) {
          res$Consensus_per_model[[y]][[x]]
        })
        mcs <- terra::rast(mcs)
        names(mcs) <- consensus
        list(Replicates = rep[[x]], Model_consensus = mcs)
      }
    })

    names(res) <- nm
    res <- c(res, General_consensus = terra::rast(gen_res))
  }

  #### Get consensus by model if raster_variables is a data.frame ####
  if (inherits(raster_variables, "data.frame")) {
    rep <- lapply(p_models, as.data.frame)

    if (nrep == 1) {
      res$Consensus_per_model$Full_model <- lapply(p_models, function(x) {x})
    } else {
      if (consensus_per_model) {
        if ("median" %in% consensus) {
          res$Consensus_per_model$median <- lapply(p_models, function(x){
            apply(as.data.frame(x), 1, stats::median)
          })
        }
        if ("mean" %in% consensus) {
          res$Consensus_per_model$mean <- lapply(p_models, function(x){
            apply(as.data.frame(x), 1, mean)})
        }
        if ("stdev" %in% consensus) {
          res$Consensus_per_model$stdev <- lapply(p_models, function(x){
            apply(as.data.frame(x), 1, stats::sd)})
        }
        if ("range" %in% consensus) {
          res$Consensus_per_model$range <- lapply(p_models, function(x){
            x_df <- as.data.frame(x)
            min_x <- apply(x_df, 1, min, na.rm = TRUE)
            max_x <- apply(x_df, 1, max, na.rm = TRUE)
            return(max_x - min_x)
            })
        }
      }
    }

    gen_res <- list()
    if (consensus_general && length(p_models) == 1 && consensus_per_model) {
      if (nrep == 1) {
        gen_res$Full_model <- res$Consensus_per_model$Full_model
      } else {
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
      }

    } else if (consensus_general && length(p_models) > 1) {
      all_rep <- as.data.frame(p_models)

      if ("median" %in% consensus) {
        gen_res$median <- apply(all_rep, 1, stats::median)
      }
      if ("mean" %in% consensus) {
        gen_res$mean <- apply(all_rep, 1, mean)
      }
      if ("stdev" %in% consensus) {
        gen_res$stdev <- apply(all_rep, 1, stats::sd)
      }
      if ("range" %in% consensus) {
        min_all_rep <- apply(all_rep, 1, min, na.rm = TRUE)
        max_all_rep <- apply(all_rep, 1, max, na.rm = TRUE)
        gen_res$range <- max_all_rep - min_all_rep
      }
    }

    res <- lapply(1:length(nm), function(x) {
      if (nrep == 1) {
        mcs <- lapply("Full_model", function(y) {
          res$Consensus_per_model[[y]][[x]]
        })
        names(mcs) <- "Full_model"
        list(Model_consensus = mcs)
      } else {
        mcs <- lapply(consensus, function(y) {
          res$Consensus_per_model[[y]][[x]]
        })
        names(mcs) <- consensus
        list(Replicates = rep[[x]], Model_consensus = as.data.frame(mcs))
      }
    })

    names(res) <- nm
    res$General_consensus <- as.data.frame(gen_res)
  }

  # Write results to disk if required
  if (write_files) {
    if (!file.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }

    #Save if raster_variables are spatraster
    if(inherits(raster_variables, "SpatRaster")){
      sapply(nm, function(i) {
        if (write_replicates & nrep > 1) {
          terra::writeRaster(res[[i]]$Replicates,
                           file.path(out_dir, paste0(i, "_replicates.tif")),
                           overwrite = overwrite)
          }
          terra::writeRaster(res[[i]]$Model_consensus,
                         file.path(out_dir, paste0(i, "_consensus.tif")),
                         overwrite = overwrite)
          })
      terra::writeRaster(res$General_consensus,
                       file.path(out_dir, "General_consensus.tif"),
                       overwrite = overwrite)
    }

    #Save if raster_variables are data.frame
    if(inherits(raster_variables, "data.frame")){
      sapply(nm, function(i) {
        if (write_replicates & nrep > 1) {
          utils::write.csv(res[[i]]$Replicates,
                           file.path(out_dir, paste0(i, "_replicates.csv")))
        }
        utils::write.csv(res[[i]]$Model_consensus,
                file.path(out_dir, paste0(i, "_consensus.csv")))
      })
      utils::write.csv(res$General_consensus,
                       file.path(out_dir, "General_consensus.csv"))
    }
  }

  return(res)
} # End of function
