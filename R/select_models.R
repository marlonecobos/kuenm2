#' Select models that perform the best among candidates
#'
#' @description
#' This function selects the best models according to user-defined criteria,
#' evaluating statistical significance (partial ROC), predictive ability
#' (omission rates), and model complexity (AIC).
#'
#' @usage
#' select_models(calibration_results = NULL, candidate_models = NULL, data = NULL,
#'               algorithm = NULL, compute_proc = FALSE,
#'               addsamplestobackground = TRUE, weights = NULL,
#'               test_concave = FALSE, omrat_threshold = 10,
#'               allow_tolerance = TRUE, tolerance = 0.01, AIC_option = "ws",
#'               significance = 0.05, delta_aic = 2, parallel = FALSE,
#'               ncores = NULL, progress_bar = FALSE,verbose = TRUE)
#'
#' @param calibration_results an object of class `calibration_results` returned
#' by the [calibration()] function. Default is NULL.
#' @param candidate_models (data.frame) a summary of the evaluation metrics for each
#' candidate model. Required only if `calibration_results` is NULL. In the
#' output of the [calibration()], this data.frame is located in
#' `$calibration_results$Summary`. Default is NULL.
#' @param data an object of class `prepared_data` returned by the
#' [prepare_data()] function. Required only if `calibration_results`
#' is NULL and `compute_proc` is TRUE.
#' @param algorithm (character) model algorithm, either "glm" or "maxnet". The
#' default, NULL, uses the one defined as part of `calibration_results`, or
#' `data`. If those arguments are not used, `algorithm` must be defined.
#' @param compute_proc (logical) whether to compute partial ROC tests for the
#' selected models. This is required when partial ROC is not calculated for all
#' candidate models during calibration. Default is FALSE.
#' @param addsamplestobackground (logical) whether to add to the background any
#' presence sample that is not already there. Required only if `compute_proc` is
#' TRUE and `calibration_results` is NULL.Default is TRUE.
#' @param weights (numeric) a numeric vector specifying weights for the
#' occurrence records. Required only if `compute_proc` is TRUE and
#' `calibration_results` is NULL. Default is NULL.
#' @param test_concave (logical) whether to remove candidate models presenting
#' concave curves. Default is FALSE.
#' @param omrat_threshold (numeric) the maximum omission rate a candidate model
#' can have to be considered a best model. Default is 10. This value must match
#' one of the values specified in `omrat` in [calibration()].
#' @param allow_tolerance (logical) whether to allow selection of models with
#' minimum values of omission rates even if their omission rate surpasses the
#' `omrat_threshold`. This is only applicable if all candidate models have
#' omission rates higher than the `omrat_threshold`. Default is TRUE.
#' @param tolerance (numeric) The value added to the minimum omission rate if it
#' exceeds the `omrat_threshold`. If `allow_tolerance = TRUE`, selected models
#' will have an omission rate equal to or less than the minimum rate plus this
#' tolerance. Default is 0.01.
#' @param AIC_option (character) the type of AIC to be calculated: "ws" for AIC
#' proposed by Warren and Seifert (2011), or "nk" for AIC proposed by Ninomiya
#' and Kawano (2016). This is only applicable if algorithm = "maxnet".
#' Default is "ws". See References for details.
#' @param significance (numeric) the significance level to select models
#' based on the partial ROC (pROC). Default is 0.05. See Details.
#' @param delta_aic (numeric) the value of delta AIC used as a threshold to
#' select models. Default is 2.
#' @param parallel (logical) whether to calculate the PROC of the candidate
#' models in parallel. Default is FALSE.
#' @param ncores (numeric) number of cores to use for parallel processing.
#' Default is NULL and uses available cores - 1. This is only applicable if
#' `parallel = TRUE`.
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is TRUE.
#' @param verbose (logical) whether to display messages during processing.
#' Default is TRUE.
#'
#' @details
#' Partial ROC is calculated following Peterson et al.
#' (2008).
#'
#' @return
#' If calibration_results is provided, it returns a new calibration_results with
#' the new selected models and summary. If calibration_results is NULL, it
#' returns a list containing the following elements:
#' - selected_models: data frame with the ID and the summary of evaluation
#' metrics for the selected models.
#' - summary: A list containing the delta AIC values for model selection, and
#' the ID values of models that failed to fit, had concave curves,
#' non-significant pROC values, omission rates above the threshold, delta AIC
#' values above the threshold, and the selected models.
#'
#' @export
#'
#' @examples
#' # Import example of calibration results (output of calibration function)
#' ## GLM
#' data(calib_results_glm, package = "kuenm2")
#'
#' #Select new best models based on another value of omission rate
#' new_best_model <- select_models(candidate_models = calib_results_glm$calibration_results$Summary,
#'                                 algorithm = "glm",
#'                                 omrat_threshold = 5)  # Omission error of 5
#'
#' # Compare with best models selected previously
#' calib_results_glm$summary$Selected  # Model 1 selected
#' new_best_model$summary$Selected  # Models 1 and 5 selected


select_models <- function(calibration_results = NULL,
                          candidate_models = NULL,
                          data = NULL,
                          algorithm = NULL,
                          compute_proc = FALSE,
                          addsamplestobackground = TRUE,
                          weights = NULL,
                          test_concave = FALSE,
                          omrat_threshold = 10,
                          allow_tolerance = TRUE,
                          tolerance = 0.01,
                          AIC_option = "ws",
                          significance = 0.05,
                          delta_aic = 2,
                          parallel = FALSE,
                          ncores = NULL,
                          progress_bar = FALSE,
                          verbose = TRUE) {
  #Check data
  if (is.null(calibration_results) & is.null(candidate_models)) {
    stop("Either 'calibration_results' or 'candidate_models' must be defined.")
  }

  if (is.null(calibration_results) & is.null(data) & compute_proc) {
    stop("If 'compute_proc' = TRUE, 'calibration_results' or 'data' must be defined.")
  }

  if (!is.null(calibration_results)) {
    candidate_models <- calibration_results$calibration_results$Summary
    algorithm <- calibration_results$algorithm
  } else {
    if (!is.null(data)) {
      algorithm <- data$algorithm
    }
    if (is.null(algorithm)) {
      stop("Argument 'algorithm' must be defined.")
    }
  }

  #Update weights and addsamplestobackground from calibration_results, if necessary
  if (!is.null(calibration_results)) {
    weights <- calibration_results$weights
    addsamplestobackground <- calibration_results$addsamplestobackground
  }

  # Adjust AIC column based on model type
  if (algorithm == "maxnet") {
    # Remove the unused AIC column in maxnet
    if (AIC_option == "nk") {
      AIC_option <- "AIC_nk"
      candidate_models$AIC_ws <- NULL
    } else if (AIC_option == "ws") {
      AIC_option <- "AIC_ws"
      candidate_models$AIC_nk <- NULL
    } else {
      stop("Unsupported AIC option. Please use 'nk' or 'ws'.")
    }
  } else if (algorithm == "glm") {
    AIC_option <- "AIC" # For glm models, we only use a single AIC column
  } else {
    stop("Unsupported algorithm. Please use 'maxnet' or 'glm'.")
  }

  # Omission rate column name
  om_thr <- paste0("Omission_rate_at_", omrat_threshold, ".mean")

  #proc-pval columns
  proc_pval <- paste0("pval_pROC_at_", omrat_threshold, ".mean")


  #Check if it's necessary calculate proc
  if (all(is.na(candidate_models[[proc_pval]]))) {
    if (!compute_proc) {
      stop("pROC values were not provided as part of the input. Set 'compute_proc' to TRUE.")
    }
  }

  # Log the number of models being filtered
  if (verbose) {
    message("Selecting best among ", nrow(candidate_models), " models.")
  }


  #### Initiate filtering if it's necessary to calculate proc ####
  if (compute_proc) {

    if (verbose) {
      message("Calculating pROC...")
    }

    any_bad <- TRUE  # To initiate looping
    id_to_remove <- 0

    while(any_bad) {
      #Remove bad models
      candidate_models <- candidate_models[!(candidate_models$ID %in% id_to_remove), ]


      # Log the number of models being filtered
      if (verbose) {
        message("\nFiltering ", nrow(candidate_models), " models.")
      }

      # Remove models with errors
      na_models <- candidate_models[is.na(candidate_models$is_concave), "ID"]

      if (verbose) {
        message("Removing ", length(na_models), " model(s) because they failed to fit.")
      }

      candidate_models <- candidate_models[!is.na(candidate_models$is_concave), ]

      # Remove concave curves if test_concave is TRUE
      if (test_concave) {
        concave_models <- candidate_models[candidate_models$is_concave, "ID"]

        if (verbose) {
          message("Removing ", length(concave_models), " model(s) with concave curves.")
        }

        candidate_models <- candidate_models[!candidate_models$is_concave, ]
      } else {
        concave_models <- 0
      }

      # Subset models by omission rate
      high_omr <- candidate_models[candidate_models[, om_thr] > omrat_threshold / 100, "ID"]
      cand_om <- candidate_models[candidate_models[, om_thr] <= omrat_threshold / 100, ]

      if (verbose) {
        message(nrow(cand_om), " models were selected with omission rate below ",
                omrat_threshold, "%.")
      }

      # Stop if no models meet the omission rate threshold and allow_tolerance is FALSE
      if (nrow(cand_om) == 0 & !allow_tolerance) {
        stop("There are no models with values of omission rate below ",
             omrat_threshold, "%. Try with 'allow_tolerance' = TRUE.")
      }

      # Apply tolerance if no models meet the omission rate threshold and allow_tolerance is TRUE
      if (nrow(cand_om) == 0 & allow_tolerance) {
        min_thr <- min(candidate_models[, om_thr])
        cand_om <- subset(candidate_models,
                          candidate_models[, om_thr] <= min_thr + tolerance)
        high_omr <- candidate_models[candidate_models[, om_thr] > min_thr + tolerance, "ID"]

        if (verbose) {
          message("Minimum value of omission rate (", round(min_thr * 100, 1),
                  "%) is above the selected threshold (", omrat_threshold,
                  "%).\nApplying tolerance and selecting ", nrow(cand_om),
                  " models with omission rate <",
                  round(min_thr * 100 + tolerance, 1), "%.")
        }
      }

      # Calculate delta AIC and select models based on delta AIC
      cand_om$dAIC <- cand_om[, AIC_option] - min(cand_om[, AIC_option],
                                                  na.rm = TRUE)
      high_aic <- cand_om[cand_om$dAIC > delta_aic | is.na(cand_om$dAIC), "ID"]
      cand_final <- cand_om[cand_om$dAIC <= delta_aic & !is.na(cand_om$dAIC), ]

      if (verbose) {
        message("Selecting ", nrow(cand_final),
                " final model(s) with delta AIC <", delta_aic, ".")
      }

      #Identify models to calculate proc
      to_calculate <- cand_final$ID[is.na(cand_final[, proc_pval])]

      #If there are models to calculate proc...
      if (length(to_calculate) > 0) {
        cand_final_with_proc <- cand_final[!(cand_final$ID %in% to_calculate), ]

        if (nrow(cand_final_with_proc) > 0) {
          cand_final <- cand_final[cand_final$ID %in% to_calculate, ]
        }

        # Validate pROC
        if (verbose) {
          message("Validating pROC of selected models...")
        }

        if (is.null(data)) {
          data <- calibration_results[1:10]
          class(data) <- "prepared_data"
        }

        proc_values <- partial_roc(formula_grid = cand_final, data = data,
                                   omission_rate = omrat_threshold,
                                   addsamplestobackground, weights,
                                   algorithm, parallel, ncores,
                                   progress_bar)

        # Create a copy of cand_final to keep the original data unchanged
        cand_final_updated <- cand_final

        # Replace NA values in cand_final_updated with corresponding values from proc_values
        ncpr <- ncol(proc_values)
        rrep <- match(proc_values$ID, cand_final_updated$ID)

        cand_final_updated[rrep, colnames(proc_values[, -ncpr])] <-
          proc_values[, -ncpr]

        #Append cand_final_with_proc
        if (nrow(cand_final_with_proc) > 0) {
          cand_final_updated <- rbind(cand_final_with_proc, cand_final_updated)
        }

        # Check if p_value is non-significative
        p_value_omr <- cand_final_updated[, paste0("pval_pROC_at_",
                                                  omrat_threshold, ".mean")]
        any_bad <- any(p_value_omr > significance | is.na(p_value_omr))

        #Get models to remove, if necessary
        id_to_remove <- cand_final_updated$ID[p_value_omr > significance |
                                                is.na(p_value_omr)]

      } else { #If all models have already been calculated
        any_bad <- FALSE
      }

      #if there are non-significant models...
      if (any_bad) {
        if (verbose) {
          message("\nSome of the selected models have non-significant pROC values. Re-selecting models...\n")
        }

        #Update candidate_models with proc values calculated
        candidate_models[proc_values$ID, colnames(proc_values[, -ncpr])] <-
          proc_values[, -ncpr]
        #if there are not non-significant models...
      } else {
        #cand_final <- cand_final_updated
        insig_proc <- NULL
      }
    } #End of anybad
  } #End of calculate proc

  if (!compute_proc) {
    #### Initiate filtering if it's NOT necessary to calculate proc ####
    # Remove models with errors
    na_models <- candidate_models[is.na(candidate_models$is_concave), "ID"]

    if (verbose) {
      message("\nRemoving ", length(na_models), " model(s) that failed to fit.")
    }

    candidate_models <- candidate_models[!is.na(candidate_models$is_concave), ]

    # Remove concave curves if test_concave is TRUE
    if (test_concave) {
      concave_models <- candidate_models[candidate_models$is_concave, "ID"]

      if (verbose) {
        message("Removing ", length(concave_models),
                " model(s) with concave responses.")
      }

      candidate_models <- candidate_models[!candidate_models$is_concave, ]
    } else {
      concave_models <- 0
    }

    # Subset models with significant pROC
    insig_proc <- candidate_models[candidate_models[[proc_pval]] >= significance |
                                is.na(candidate_models[[proc_pval]]), "ID"]

    if (verbose) {
      message("Removing ", length(insig_proc),
              " model(s) with non-significant pROC values.")
    }

    candidate_models <- candidate_models[candidate_models[[proc_pval]] < significance &
                                 !is.na(candidate_models[[proc_pval]]), ]

    # Subset models by omission rate
    high_omr <- candidate_models[candidate_models[, om_thr] > omrat_threshold / 100, "ID"]
    cand_om <- candidate_models[candidate_models[, om_thr] <= omrat_threshold / 100, ]

    if (verbose) {
      message(nrow(cand_om), " models passed the ", omrat_threshold,
              "% omission criterion.")
    }

    # Stop if no models meet the omission rate threshold and allow_tolerance is FALSE
    if (nrow(cand_om) == 0 & !allow_tolerance) {
      stop("There were no models with omissions below ", omrat_threshold, "%.",
           "Try using the arguments 'allow_tolerance' and 'tolerance'.")
    }

    # Apply tolerance if no models meet the omission rate threshold and allow_tolerance is TRUE
    if (nrow(cand_om) == 0 & allow_tolerance) {
      min_thr <- min(candidate_models[, om_thr])
      cand_om <- subset(candidate_models, candidate_models[, om_thr] <= min_thr + tolerance)
      high_omr <- candidate_models[candidate_models[, om_thr] > min_thr + tolerance, "ID"]

      if (verbose) {
        message("Minimum value of omission in models (", round(min_thr * 100, 1),
                "%) > omission criterion (", omrat_threshold, "%).\n",
                "Applying tolerance: ", nrow(cand_om), " models with omission <",
                round((min_thr + tolerance) * 100, 1), "% were found.")
      }
    }

    # Calculate delta AIC and select models based on delta AIC
    cand_om$dAIC <- cand_om[, AIC_option] - min(cand_om[, AIC_option])
    high_aic <- cand_om[cand_om$dAIC > delta_aic, "ID"]
    cand_final_updated <- cand_om[cand_om$dAIC <= delta_aic, ]

    if (verbose) {
      message("Selecting ", nrow(cand_final_updated), " final model(s) with delta AIC <",
              delta_aic, ".")
    }
  }

  #Final message
  if (nrow(cand_final_updated) > 0 & verbose) {
    message("\nAll selected models have significant pROC values.")  #Do we need this?
  }

  if (nrow(cand_final_updated) == 0) {
    stop("It is impossible to select models with significant pROC values.")
  }

  # Final results
  sel_res <- list(selected_models = cand_final_updated,
                  summary = list(delta_AIC = delta_aic,
                                 omission_rate_thr = omrat_threshold,
                                 Errors = na_models,
                                 Concave = concave_models,
                                 Non_sig_pROC = insig_proc,
                                 High_omission_rate = high_omr,
                                 High_AIC = high_aic,
                                 Selected = cand_final_updated$ID))

  #Update calibration_results if necessary
  if (!is.null(calibration_results)) {
    calibration_results$selected_models <- sel_res$selected_models
    calibration_results$summary <- sel_res$summary
    calibration_results$omission_rate <- omrat_threshold
    return(calibration_results)
  } else {
    return(sel_res)
  }
}
