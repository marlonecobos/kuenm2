#' Select models that perform the best among all candidates
#'
#' @description
#' This function selects the best models according to user-defined criteria,
#' evaluating statistical significance (partial ROC), predictive ability
#' (omission rates), and model complexity (AIC).
#'
#' @usage
#' sel_best_models(calibration_results = NULL, cand_models = NULL,
#'                 algorithm = c("glmnet", "glm"), test_concave = TRUE,
#'                 omrat_threshold = 10, allow_tolerance = TRUE,
#'                 tolerance = 0.01, AIC_option = "ws", significance = 0.05,
#'                 delta_aic = 2, verbose = TRUE)
#'
#' @param calibration_results an object of class `calibration_results` returned
#' by the \code{\link{calibration}}() function. Default is NULL.
#' @param cand_models (data.frame) a summary of the evaluation metrics for each
#' candidate model. In the output of the \code{\link{calibration}}(), this
#' data.frame is located in `$calibration_results$Summary`. Default is NULL.
#' @param algorithm (character) model type, either "glm" or "glmnet".
#' @param test_concave (logical) whether to remove candidate models presenting
#' concave curves. Default is TRUE.
#' @param omrat_threshold (numeric) the maximum omission rate a candidate model
#' can have to be considered a best model. Default is 10. This value must match
#' one of the values specified in `omrat` in \code{\link{calibration}}().
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
#' and Kawano (2016). This is only applicable if algorithm = "glmnet".
#' Default is "ws". See References for details.
#' @param significance (numeric) the significance level to select models
#' based on the partial ROC (pROC). Default is 0.05. See Details.
#' @param delta_aic (numeric) the value of delta AIC used as a threshold to
#' select models. Default is 2.
#' @param verbose (logical) whether to display messages during processing.
#' Default is TRUE.
#'
#' @details
#' Partial ROC is calculated following Peterson et al.
#' (2008; http://dx.doi.org/10.1016/j.ecolmodel.2007.11.008).
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
#' @references
#' Ninomiya, Yoshiyuki, and Shuichi Kawano. "AIC for the Lasso in generalized
#' linear models." (2016): 2537-2560.
#'
#' Warren, D. L., & Seifert, S. N. (2011). Ecological niche modeling in Maxent:
#' the importance of model complexity and the performance of model selection
#' criteria. Ecological applications, 21(2), 335-342.
#'
#' @examples
#' # Import example of calibration results (output of calibration function)
#' ## GLM
#' data("calib_results_glm", package = "kuenm2")
#'
#' #Select new best models based on another value of omission rate
#' new_best_model <- sel_best_models(cand_models = calib_results_glm$calibration_results$Summary,
#'                                   algorithm = "glm",
#'                                   omrat_threshold = 5,
#'                                   delta_aic = 10)  # Higher value of delta AIC
#'
#' # Compare with best models selected previously, with omission rate of 10 and delta AIC of 2
#' calib_results_glm$summary$Selected  # Models 1, 2 and 5 selected
#' new_best_model$summary$Selected  # Models 1 and 5 selected
#'
#' # Replace selected models in calib_results
#' calib_results_glm$selected_models <- new_best_model$cand_final
#' calib_results_glm$summary <- new_best_model$summary

sel_best_models <- function(calibration_results = NULL,
                            cand_models = NULL,
                            algorithm = c("glmnet", "glm"),
                            test_concave = TRUE,
                            omrat_threshold = 10,
                            allow_tolerance = TRUE,
                            tolerance = 0.01,
                            AIC_option = "ws",
                            significance = 0.05,
                            delta_aic = 2,
                            verbose = TRUE) {
  #Check data
  if (is.null(calibration_results) & is.null(cand_models)) {
    stop("You must specified calibration_results or cand_models")
  }

  if (!is.null(calibration_results)) {
    cand_models <- calibration_results$calibration_results$Summary
  }

  # Adjust AIC column based on model type
  if (algorithm == "glmnet") {
    # Remove the unused AIC column in glmnet
    if (AIC_option == "nk") {
      AIC_option <- "AIC_nk"
      cand_models$AIC_ws <- NULL
    } else if (AIC_option == "ws") {
      AIC_option <- "AIC_ws"
      cand_models$AIC_nk <- NULL
    } else {
      stop("Unsupported AIC option. Please use 'nk' or 'ws'.")
    }
  } else if (algorithm == "glm") {
    AIC_option <- "AIC" # For glm models, we only use a single AIC column
  } else {
    stop("Unsupported algorithm. Please use 'glmnet' or 'glm'.")
  }

  # Omission rate column name
  om_thr <- paste0("Omission_rate_at_", omrat_threshold, "_mean")

  #proc-pval columns
  proc_pval <- paste0("pval_pROC_at_", omrat_threshold, "_mean")

  # Log the number of models being filtered
  if (verbose) {
    message("Selecting best among ", nrow(cand_models), " models")
  }

  # Remove models with errors
  na_models <- cand_models[is.na(cand_models$is_concave), "ID"]
  if (verbose) {
    message("Removing ", length(na_models), " model(s) that failed to fit")
  }
  cand_models <- cand_models[!is.na(cand_models$is_concave), ]

  # Remove concave curves if test_concave is TRUE
  if (test_concave) {
    concave_models <- cand_models[cand_models$is_concave, "ID"]
    if (verbose) {
      message("Removing ", length(concave_models), " model(s) with concave responses")
    }
    cand_models <- cand_models[!cand_models$is_concave, ]
  } else {
    concave_models <- 0
  }

  # Subset models with significant pROC
  insig_proc <- cand_models[cand_models[[proc_pval]] >= significance |
                              is.na(cand_models[[proc_pval]]), "ID"]
  if (verbose) {
    message("Removing ", length(insig_proc), " model(s) with non-significant pROC values")
  }
  cand_models <- cand_models[cand_models[[proc_pval]] < significance &
                               !is.na(cand_models[[proc_pval]]), ]

  # Subset models by omission rate
  high_omr <- cand_models[cand_models[, om_thr] > omrat_threshold / 100, "ID"]
  cand_om <- cand_models[cand_models[, om_thr] <= omrat_threshold / 100, ]
  if (verbose) {
    message(nrow(cand_om), " models passed the ", omrat_threshold, "% omission criterion")
  }

  # Stop if no models meet the omission rate threshold and allow_tolerance is FALSE
  if (nrow(cand_om) == 0 & !allow_tolerance) {
    stop("There were no models with omissions below ", omrat_threshold, "%.",
         "Try using the arguments 'allow_tolerance' and 'tolerance'.")
  }

  # Apply tolerance if no models meet the omission rate threshold and allow_tolerance is TRUE
  if (nrow(cand_om) == 0 & allow_tolerance) {
    min_thr <- min(cand_models[, om_thr])
    cand_om <- subset(cand_models, cand_models[, om_thr] <= min_thr + tolerance)
    high_omr <- cand_models[cand_models[, om_thr] > min_thr + tolerance, "ID"]
    if (verbose) {
      message("Minimum value of omission in models (", round(min_thr * 100, 1),
              "%) > omission criterion (", omrat_threshold, "%).\n",
              "Applying tolerance: ", nrow(cand_om), " models with omission <",
              round((min_thr + tolerance) * 100, 1), "% were found")
    }
  }

  # Calculate delta AIC and select models based on delta AIC
  cand_om$dAIC <- cand_om[, AIC_option] - min(cand_om[, AIC_option])
  high_aic <- cand_om[cand_om$dAIC > delta_aic, "ID"]
  cand_final <- cand_om[cand_om$dAIC <= delta_aic, ]

  if (verbose) {
    message("Selecting ", nrow(cand_final), " final model(s) with delta AIC <", delta_aic)
  }

  # Final results
  sel_res <- list(cand_final = cand_final,
                  summary = list("delta_AIC" = delta_aic,
                                 "omission_rate_thr" = omrat_threshold,
                                 "Errors" = na_models,
                                 "Concave" = concave_models,
                                 "Non_sig_pROC" = insig_proc,
                                 "High_omission_rate" = high_omr,
                                 "High_AIC" = high_aic,
                                 "Selected" = cand_final$ID))

  #Update calibration_results if necessary
  if(!is.null(calibration_results)) {
    calibration_results$selected_models <- sel_res$cand_final
    calibration_results$summary <- sel_res$summary
    calibration_results$omission_rate <- omrat_threshold
    return(calibration_results)
  } else {
  return(sel_res)}
}
