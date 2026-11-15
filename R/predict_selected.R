#' Predict selected models for a single scenario
#'
#' @description
#' This function predicts selected models for a single set of new data
#' using either `maxnet` or `glm`. It provides options to save the
#' output and compute consensus results (mean, median, etc.) across
#' replicates and models.
#'
#' @usage
#' predict_selected(models, new_variables, mask = NULL, write_files = FALSE,
#'                  write_replicates = FALSE, out_dir = NULL,
#'                  consensus_per_model = TRUE, consensus_general = TRUE,
#'                  consensus = c("median", "range", "mean", "stdev"),
#'                  extrapolation_type = "E", var_to_restrict = NULL,
#'                  type = "cloglog", overwrite = FALSE, progress_bar = TRUE,
#'                  progress_type = "timer")
#'
#' @param models an object of class `fitted_models` returned by the
#' \code{\link{fit_selected}}() function.
#' @param new_variables a SpatRaster or data.frame of predictor variables.
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
#' @param var_to_restrict (character) vector specifying which variables to clamp
#' or not to extrapolate for. Only applicable if extrapolation_type is "EC" or "NE".
#' Default is `NULL`, clamping and no extrapolation will be done for all variables.
#' @param type (character) the format of prediction values. For `maxnet` models,
#' valid options are `"raw"`, `"cumulative"`, `"logistic"`, and `"cloglog"`.
#' For `glm` models, valid options are `"cloglog"`, `"response"`, `"raw"`,
#' `"cumulative"` and `"link"`. Default is `"cloglog"`.
#' @param overwrite (logical) whether to overwrite SpatRasters if they already
#' exist. Only applicable if `write_files = TRUE`. Default is FALSE.
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is TRUE.
#' @param progress_type (character) type of the progress bar. See
#' `pbapply::pbtypes()` for available progress bar types depending on operating
#' system. Default is "timer".
#'
#' @details
#' When predicting to areas where the variables are beyond the lower or upper
#' limits of the calibration data, users can choose to free extrapolate the
#' predictions (`extrapolation_type = "E"`), extrapolate with clamping
#' (extrapolation_type = "EC"), or not extrapolate (extrapolation_type = "NE").
#' When clamping, the variables are set to minimum and maximum values
#' established for the maximum and minimum values within calibration data. In
#' the no extrapolation approach, any cell with at least one variable listed in
#' `var_to_restrict` falling outside the calibration range is assigned a suitability
#' value of 0.
#'
#' @return
#' A list containing SpatRaster or data.frames predictions for each replicate,
#' along with the consensus results for each model and the overall general consensus.
#'
#' @importFrom terra rast clamp predict median mean stdev writeRaster cells values ncell mask
#' @importFrom utils write.csv
#' @importFrom stats glm predict median sd complete.cases
#' @importFrom pbapply pboptions pblapply
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
#' p <- predict_selected(models = fitted_model_maxnet, new_variables = var)
#'
#' # Example with GLMs
#' # Import example of fitted_models (output of fit_selected()) without replicates
#' data("fitted_model_glm", package = "kuenm2")
#'
#' # Predict to single scenario
#' p_glm <- predict_selected(models = fitted_model_glm, new_variables = var)
#'
#' # Plot predictions
#' terra::plot(c(p$General_consensus$median, p_glm$General_consensus),
#'             col = rev(terrain.colors(240)), main = c("MAXNET", "GLM"),
#'             zlim = c(0, 1))

predict_selected <- function(models,
                             new_variables,
                             mask = NULL,
                             write_files = FALSE,
                             write_replicates = FALSE,
                             out_dir = NULL,
                             consensus_per_model = TRUE,
                             consensus_general = TRUE,
                             consensus = c("median", "range", "mean", "stdev"),
                             extrapolation_type = "E",
                             var_to_restrict = NULL,
                             type = "cloglog",
                             overwrite = FALSE,
                             progress_bar = TRUE,
                             progress_type = "timer") {

  # Check data
  if (missing(models)) {
    stop("Argument 'models' must be defined.")
  }
  if (missing(new_variables)) {
    stop("Argument 'new_variables' must be defined.")
  }
  if (!inherits(models, "fitted_models")) {
    stop("Argument 'models' must be a 'fitted_models' object.")
  }
  if (!inherits(new_variables, c("SpatRaster", "data.frame"))) {
    stop("Argument 'new_variables' must be a 'SpatRaster' or 'data.frame'.")
  }

  if (!is.null(mask) && !inherits(mask, c("SpatRaster", "SpatVector",
                                          "SpatExtent"))) {
    stop("Argument 'mask' must be a 'SpatVector', 'SpatExtent' or 'SpatRaster'.")
  }

  if (write_files && is.null(out_dir)) {
    stop("If write_files = TRUE, out_dir must be specified.")
  }
  if (write_files && !inherits(out_dir, "character")) {
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

  if (length(extrapolation_type) > 1) {
    stop("Extrapolation type accepts only one of these values: 'E', 'EC', or 'NE'.")
  }

  extrapolation_out <- setdiff(extrapolation_type, c("E", "EC", "NE"))
  if (length(extrapolation_out) > 0) {
    stop("Invalid 'extrapolation_type' provided.",
         "\nAvailable options are: 'E', 'EC', and 'NE'.")
  }

  if (extrapolation_type %in% c("E", "EC", "NE") && !is.null(var_to_restrict) &&
      !inherits(var_to_restrict, "character")) {
    stop("Argument 'var_to_restrict' must be NULL or 'character'.")
  }

  if (!inherits(type, "character")) {
    stop("Argument 'type' must be 'character'.")
  }

  if (models$algorithm == "maxnet") {
    if (!any(c("raw", "cumulative", "logistic", "cloglog") %in% type)) {
      stop("Invalid 'type' provided.",
           "\nAvailable options for maxnet models are: 'raw', 'cumulative', 'logistic', or 'cloglog'.")
    }
    if (type == "raw") {
      type <- "exponential"
    }
  }

  if (models$algorithm == "glm") {
    if (!any(c("response", "cloglog", "cumulative", "link", "raw") %in% type)) {
      stop("Invalid 'type' provided.",
           "\nAvailable options for glm models are 'response', 'cloglog', 'cumulative', 'link', or 'raw'.")
    }
  }

  # Analyses start here
  is_raster <- inherits(new_variables, "SpatRaster")

  if (!is.null(mask) && is_raster) {
    new_variables <- terra::crop(new_variables, mask, mask = TRUE)
  }

  if (!is.null(models$pca)) {
    if (is_raster) {
      if (!("vars_out" %in% names(models$pca))) {
        new_variables <- terra::predict(new_variables[[models$pca$vars_in]],
                                        models$pca)
      } else {
        new_variables <- c(terra::predict(new_variables[[models$pca$vars_in]],
                                          models$pca),
                           new_variables[[models$pca$vars_out]])
      }
    } else {
      if (!("vars_out" %in% names(models$pca))) {
        new_variables <- stats::predict(new_variables[, models$pca$vars_in, drop = FALSE],
                                        models$pca)
      } else {
        new_variables <- cbind(stats::predict(new_variables[, models$pca$vars_in, drop = FALSE],
                                              models$pca),
                               new_variables[, models$pca$vars_out, drop = FALSE])
      }
      new_variables <- as.data.frame(new_variables)
    }
  }

  # Extract info from fitted object
  categorical_layers <- models[["categorical_variables"]]
  models <- models[["Models"]]
  nm <- names(models)
  nrep <- length(models[[1]])

  # Get names of the models (replicates or full model)
  names_models <- unique(unlist(lapply(models, names), use.names = FALSE))

  # If there are replicates, remove the full model from the dataset
  if (any(grepl("Replicate", names_models))) {
    models <- lapply(models, function(m) {
      m[["Full_model"]] <- NULL
      return(m)
    })
    names(models) <- nm
    nrep <- length(models[[1]])
  }

  # Convert categorical layers to factors before extracting values
  if (!is.null(categorical_layers) && length(categorical_layers) > 0) {
    for (layer in categorical_layers) {
      if (layer %in% names(new_variables)) {
        new_variables[[layer]] <- if (is_raster) {
          terra::as.factor(new_variables[[layer]])
        } else {
          as.factor(new_variables[[layer]])
        }
      } else {
        warning(paste("Layer", layer, "not found in 'new_variables'."))
      }
    }
  }


  # Prepare tabular newdata once (avoids repeated terra::predict overhead)
  template_rast <- NULL
  n_cells <- NULL
  n_rows <- NULL

  if (is_raster) {
    template_rast <- new_variables[[1]]
    n_cells <- terra::ncell(template_rast)
    vals_df <- terra::as.data.frame(new_variables, cells = TRUE, na.rm = TRUE)
    valid_cells <- vals_df$cell
    newdata_df <- vals_df[, -1, drop = FALSE]
  } else {
    valid_cells <- which(stats::complete.cases(new_variables))
    newdata_df <- new_variables[valid_cells, , drop = FALSE]
    n_rows <- nrow(new_variables)
  }

  # Handle clamping ("EC") or no-extrapolation ("NE") on tabular data
  if (extrapolation_type %in% c("EC", "NE")) {
    varmin <- models[[1]][[1]]$varmin[-1]
    varmax <- models[[1]][[1]]$varmax[-1]
    if (is.null(var_to_restrict)) {
      var_to_restrict <- setdiff(names(varmin), c("pr_bg", "fold"))
    }

    if (extrapolation_type == "EC") {
      for (col_name in var_to_restrict) {
        if (is.numeric(newdata_df[[col_name]])) {
          newdata_df[[col_name]][newdata_df[[col_name]] < varmin[col_name]] <- varmin[col_name]
          newdata_df[[col_name]][newdata_df[[col_name]] > varmax[col_name]] <- varmax[col_name]
        }
      }
    } else if (extrapolation_type == "NE") {
      no_extrapolate <- unique(unlist(lapply(var_to_restrict, function(i) {
        if (is.numeric(newdata_df[[i]])) {
          which(newdata_df[[i]] < varmin[i] | newdata_df[[i]] > varmax[i])
        } else {
          integer(0)
        }
      }), use.names = FALSE))
    }
  }

  # Setup progress bar options safely for CRAN
  old_pbopt <- pbapply::pboptions(
    type = if (progress_bar) progress_type else "none"
  )
  on.exit(pbapply::pboptions(old_pbopt), add = TRUE)

  # Fill the list with predictions (fast vectorized prediction)
  p_vecs <- pbapply::pblapply(
    X = models,
    FUN = function(model_reps) {
      lapply(model_reps, function(x) {
        if (inherits(x, "glmnet")) {
          pred <- as.numeric(predict.glmnet_mx(x, newdata = newdata_df,
                                               type = type))
        } else if (inherits(x, "glm")) {
          pred <- as.numeric(predict_glm_mx(x, newdata = newdata_df,
                                            type = type))
        }

        # If no extrapolate, set rows beyond limits to 0
        if (extrapolation_type == "NE" && length(no_extrapolate) > 0) {
          pred[no_extrapolate] <- 0
          if (type == "cumulative") {
            pred <- cumulative_predictions(pred)
          }
        }
        return(pred)
      })
    }
  )

  # Build per-model results
  res <- lapply(seq_along(nm), function(idx) {
    mod_reps <- p_vecs[[idx]]
    names(mod_reps) <- names(models[[idx]])
    mat_mod <- do.call(cbind, mod_reps)

    if (nrep == 1L) {
      mcs <- list_to_container(
        v_list = list(Full_model = mat_mod[, 1L]),
        is_raster = is_raster,
        valid_cells = valid_cells,
        template_rast = template_rast,
        n_cells = n_cells,
        n_rows = n_rows
      )
      return(list(Model_consensus = mcs))
    } else {
      rep_out <- list_to_container(
        v_list = mod_reps,
        is_raster = is_raster,
        valid_cells = valid_cells,
        template_rast = template_rast,
        n_cells = n_cells,
        n_rows = n_rows
      )
      if (consensus_per_model) {
        cons_list <- calc_mat_consensus(mat_mod, consensus)
        mcs <- list_to_container(
          v_list = cons_list,
          is_raster = is_raster,
          valid_cells = valid_cells,
          template_rast = template_rast,
          n_cells = n_cells,
          n_rows = n_rows
        )
        return(list(Replicates = rep_out, Model_consensus = mcs))
      } else {
        return(list(Replicates = rep_out))
      }
    }
  })
  names(res) <- nm

  # Build general consensus if requested
  if (consensus_general) {
    if (nrep == 1L && length(nm) == 1L) {
      if (is_raster) {
        res <- c(res, General_consensus = res[[1L]]$Model_consensus)
      } else {
        res$General_consensus <- res[[1L]]$Model_consensus
      }
    } else {
      all_mat <- do.call(cbind, unlist(p_vecs, recursive = FALSE))
      gen_list <- calc_mat_consensus(all_mat, consensus)
      gen_out <- list_to_container(
        v_list = gen_list,
        is_raster = is_raster,
        valid_cells = valid_cells,
        template_rast = template_rast,
        n_cells = n_cells,
        n_rows = n_rows
      )
      if (is_raster) {
        res <- c(res, General_consensus = gen_out)
      } else {
        res$General_consensus <- gen_out
      }
    }
  }

  # Write results to disk if required
  if (write_files) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }

    if (is_raster) {
      for (i in nm) {
        if (write_replicates && nrep > 1L) {
          terra::writeRaster(res[[i]]$Replicates,
                             file.path(out_dir, paste0(i, "_replicates.tif")),
                             overwrite = overwrite)
        }
        if (!is.null(res[[i]]$Model_consensus)) {
          terra::writeRaster(res[[i]]$Model_consensus,
                             file.path(out_dir, paste0(i, "_consensus.tif")),
                             overwrite = overwrite)
        }
      }
      if (!is.null(res$General_consensus)) {
        terra::writeRaster(res$General_consensus,
                           file.path(out_dir, "General_consensus.tif"),
                           overwrite = overwrite)
      }
    } else {
      for (i in nm) {
        if (write_replicates && nrep > 1L) {
          utils::write.csv(res[[i]]$Replicates,
                           file.path(out_dir, paste0(i, "_replicates.csv")),
                           row.names = FALSE)
        }
        if (!is.null(res[[i]]$Model_consensus)) {
          utils::write.csv(res[[i]]$Model_consensus,
                           file.path(out_dir, paste0(i, "_consensus.csv")),
                           row.names = FALSE)
        }
      }
      if (!is.null(res$General_consensus)) {
        utils::write.csv(res$General_consensus,
                         file.path(out_dir, "General_consensus.csv"),
                         row.names = FALSE)
      }
    }
  }

  return(res)
}
