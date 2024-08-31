#' Calibration Function for Model Selection and Evaluation
#'
#' This function performs model calibration using a specified grid of formulas
#' and data of class `prepare_data`. It supports both Generalized Linear Models (GLM)
#' and GLMNET models, depending on the class of the `formula_grid` argument.
#'
#' @param data A dataset of class `prepare_data` containing the data to be used in model fitting.
#' @param formula_grid A grid of formulas of class `calibration_grid_glm` or `calibration_grid_glmnet`.
#' @param test_concave Logical. Should concave curves in quadratic models be tested? Default = TRUE.
#' @param addsamplestobackground Logical. Should samples be added to the background? Default = TRUE.
#' @param use_weights Logical. Should weights be used in the models? Default = FALSE.
#' @param parallel Logical. Should parallel processing be used? Default = TRUE.
#' @param ncores Integer. Number of cores to use for parallel processing. Default = 1.
#' @param progress_bar Logical. Should a progress bar be shown during parallel processing? Only works if `parallel_type` is "doSNOW". Default = TRUE.
#' @param write_summary Logical. Should a summary of each candidate evaluation be written to a file? Default = FALSE.
#' @param out_dir Character. The directory where summaries will be written if `write_summary` is TRUE. Default = NULL.
#' @param parallel_type Character. The type of parallel processing to use. Default = "doSNOW".
#' @param return_replicate Logical. Should replicate models be returned? Default = TRUE.
#' @param omrat_thr Numeric vector. Thresholds for omission rate. Default = c(5, 10).
#' @param omrat_threshold Numeric. Specific threshold for omission rate. Default = 10.
#' @param AIC Character. AIC criteria to use. Default = "ws".
#' @param delta_aic Numeric. The delta AIC threshold. Default = 2.
#' @param allow_tolerance Logical. If omission rate is higher than the set threshold, should the model with the minimum omission rate be selected? Default = TRUE.
#' @param tolerance Numeric. Tolerance level for model selection. Default = 0.01.
#' @param skip_existing_models Logical. Should existing models be skipped? Only works if `write_summary` is TRUE. Default = FALSE.
#' @param verbose Logical. Should detailed output be printed? Default = TRUE.
#'
#' @return A list or model object based on the specified model and parameters.
#'
#'
#' @export
model_calibration <- function(data,
                              formula_grid, # Grid with formulas
                              test_concave = TRUE, # Test concave curves in quadratic models?
                              addsamplestobackground = TRUE,
                              use_weights = FALSE,
                              parallel = TRUE,
                              ncores = 1,
                              progress_bar = TRUE, # Show progress bar? Only works if parallel_type = "doSNOW"
                              write_summary = FALSE, # Write summary of each candidate evaluation?
                              out_dir = NULL, # Name of the folder to write candidate evaluations
                              parallel_type = "doSNOW",
                              return_replicate = TRUE,
                              omrat_thr = c(5, 10),
                              omrat_threshold = 10,
                              AIC = "ws",
                              delta_aic = 2,
                              allow_tolerance = TRUE, # If omission rate is higher than set, select the model with minimum omission rate
                              tolerance = 0.01,
                              skip_existing_models = FALSE, # Only works if write_summary = TRUE
                              verbose = TRUE) {

  # Check if the 'data' argument is of class 'prepare_data'
  if (!inherits(data, "prepare_data")) {
    stop("The 'data' argument must be of class 'prepare_data'.")
  }

  # Check the class of 'formula_grid' and call the appropriate function
  if (inherits(formula_grid, "calibration_grid_glmnet")) {

    # Call the function to run GLMNET models
    result <- calibration_glmnetmx(data = data,
                                   formula_grid= formula_grid[[1]],
                                   test_concave = test_concave,
                                   addsamplestobackground = addsamplestobackground,
                                   use_weights = use_weights,
                                   parallel = parallel,
                                   ncores = ncores,
                                   progress_bar = progress_bar,
                                   write_summary = write_summary,
                                   out_dir = out_dir,
                                   parallel_type = parallel_type,
                                   return_replicate = return_replicate,
                                   omrat_thr = omrat_thr,
                                   omrat_threshold = omrat_threshold,
                                   AIC = AIC,
                                   delta_aic = delta_aic,
                                   allow_tolerance = allow_tolerance,
                                   tolerance =tolerance,
                                   skip_existing_models = skip_existing_models,
                                   verbose = verbose)

  } else if (inherits(formula_grid, "calibration_grid_glm")) {

    # Call the function to run GLM models
    stop("\n No implemented yet for GLM ... working on it :-)")

  } else {
    stop("The 'formula_grid' argument must be of class 'calibration_grid_glm' or 'calibration_grid_glmnet'.")
  }
  return(result)
}
