#' Fitting and evaluation of models, and selection of the best ones
#'
#' @description
#' This function fits and evaluates candidate models using the data and grid of
#' formulas prepared with [prepare_data()]. It supports both
#' algorithms `glm` and `maxnet`. The function then selects the best models
#' based on unimodality (optional), partial ROC, omission rate, and AIC values.
#'
#' @usage
#' calibration(data, error_considered, remove_concave = FALSE,
#'             proc_for_all = FALSE, omission_rate = NULL, delta_aic = 2,
#'             allow_tolerance = TRUE, tolerance = 0.01,
#'             addsamplestobackground = TRUE, use_weights = NULL,
#'             write_summary = FALSE, output_directory = NULL,
#'             skip_existing_models = FALSE, return_all_results = TRUE,
#'             parallel = FALSE, ncores = NULL, progress_bar = TRUE,
#'             verbose = TRUE)
#'
#' @param data an object of class `prepared_data` returned by the
#' [prepare_data()] function. It contains the calibration data,
#' formulas grid, kfolds, and model type.
#' @param proc_for_all (logical) whether to apply partial ROC tests to all
#' candidate models or only to the selected models. Default is FALSE, meaning
#' that tests are applied only to the selected models.
#' @param remove_concave (logical) whether to remove candidate models presenting
#' concave curves. Default is FALSE.
#' @param addsamplestobackground (logical) whether to add to the background any
#' presence sample that is not already there. Default is TRUE.
#' @param use_weights (logical) whether to apply the weights present in the
#' data. The default, NULL, uses weights provided in `data`. If they are not
#' present in `data`, NULL weights are 1 for presences and 100 for background.
#' If turned to FALSE, it uses NULL weights even if present in `data`.
#' @param parallel (logical) whether to fit the candidate models in parallel.
#' Default is FALSE.
#' @param ncores (numeric) number of cores to use for parallel processing.
#' Default is NULL and uses available cores - 1. This is only applicable if
#' `parallel = TRUE`.
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is TRUE.
#' @param write_summary (logical) whether to save the evaluation results for
#' each candidate model to disk. Default is FALSE.
#' @param output_directory (character) the file name, with or without a path, for saving
#' the evaluation results for each candidate model. This is only applicable if
#' `write_summary = TRUE`.
#' @param skip_existing_models (logical) whether to check for and skip candidate
#' models that have already been fitted and saved in `output_directory`. This is only
#' applicable if `write_summary = TRUE`. Default is FALSE.
#' @param return_all_results (logical) whether to return the evaluation results
#' for each replicate. Default is TRUE, meaning evaluation results for each
#' replicate will be returned.
#' @param error_considered (numeric) values from 0 to 100 representing the
#' percentage of potential error due to any source of uncertainty in your data.
#' This value is used to calculate omission rates and partial ROC. See details.
#' @param omission_rate (numeric) values from 0 - 100, the maximum omission rate
#' a candidate model can have to be considered as a potentially selected model.
#' The default, NULL, uses the value in `error_considered`. If more that one
#' value is used in `error_considered`, `omission_rate` must be defined.
#' @param delta_aic (numeric) the value of delta AIC used as a threshold to
#' select models. Default is 2.
#' @param allow_tolerance (logical) whether to allow selection of models with
#' minimum values of omission rates even if their omission rate surpasses the
#' `omission_rate`. This is only applicable if all candidate models have
#' omission rates higher than the `omission_rate`. Default is TRUE.
#' @param tolerance (numeric) The value added to the minimum omission rate if it
#' exceeds the `omission_rate`. If `allow_tolerance = TRUE`, selected models
#' will have an omission rate equal to or less than the minimum rate plus this
#' tolerance. Default is 0.01.
#' @param verbose (logical) whether to display messages during processing.
#' Default is TRUE.
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach `%dopar%`
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv head write.csv
#' @importFrom stats aggregate glm as.formula
#' @importFrom glmnet glmnet.control glmnet
#' @importFrom enmpa proc_enm
#'
#' @export
#'
#' @return
#' An object of class 'calibration_results' containing the following elements:
#' - species: a character string with the name of the species.
#' - calibration data: a data.frame containing a column (`pr_bg`) that
#' identifies occurrence points (1) and background points (0), along with the
#' corresponding values of predictor variables for each point.
#' - formula_grid: data frame containing the calibration grid with possible
#' formulas and parameters.
#' - kfolds: a list of vectors with row indices corresponding to each fold.
#' - data_xy: a data.frame with occurrence and background coordinates.
#' - continuous_variables: a character indicating the continuous variables.
#' - categorical_variables: a character, categorical variable names (if used).
#' - weights: a numeric vector specifying weights for data_xy (if used).
#' - pca: if a principal component analysis was performed with variables, a list
#' of class "prcomp". See \code{\link[stats]{prcomp}}() for details.
#' - algorithm: the model type (glm or maxnet)
#' - calibration_results: a list containing a data frame with all evaluation
#' metrics for all replicates (if `return_all_results = TRUE`) and a summary of
#' the evaluation metrics for each candidate model.
#' - omission_rate: The omission rate used to select models.
#' - addsamplestobackground: a logical value indicating whether any presence
#' sample not already in the background was added.
#' - selected_models: data frame with the ID and the summary of evaluation
#' metrics for the selected models.
#' - summary: A list containing the delta AIC values for model selection, and
#' the ID values of models that failed to fit, had concave curves,
#' non-significant pROC values, omission rates above the threshold, delta AIC
#' values above the threshold, and the selected models.
#'
#' @references
#' Ninomiya, Yoshiyuki, and Shuichi Kawano. "AIC for the Lasso in generalized
#' linear models." (2016): 2537-2560.
#'
#' Warren, D. L., & Seifert, S. N. (2011). Ecological niche modeling in Maxent:
#' the importance of model complexity and the performance of model selection
#' criteria. Ecological applications, 21(2), 335-342.
#'
#' @details
#' Partial ROC is calculated using the values defined in `error_considered`
#' following Peterson et al. (2008).
#'
#' Omission rates are calculated using separate testing data
#' subsets. Users can specify multiple values of `error_considered` to calculate
#' this metric (e.g., c(5, 10)), but only one can be used as the omission
#' rate for model selection.
#'
#' Model fitting and complexity (AICc) is assessed using models generated with
#' the complete set of occurrences. AICc values are computed as proposed by
#' Warren and Seifert (2011).
#'
#' @examples
#' # Import prepared data for maxnet models
#' data(sp_swd, package = "kuenm2")
#'
#' # Model calibration (maxnet)
#' m <- calibration(data = sp_swd, error_considered = 10)
#'
#' m
#'
#' # Import prepared data for GLM models
#' data(sp_swd_glm, package = "kuenm2")
#'
#' ## Model calibration (GLM)
#' m_glm <- calibration(data = sp_swd_glm, error_considered = 10)
#'
#' m_glm


calibration <- function(data,
                        error_considered,
                        remove_concave = FALSE,
                        proc_for_all = FALSE,
                        omission_rate = NULL,
                        delta_aic = 2,
                        allow_tolerance = TRUE,
                        tolerance = 0.01,
                        addsamplestobackground = TRUE,
                        use_weights = NULL,
                        write_summary = FALSE,
                        output_directory = NULL,
                        skip_existing_models = FALSE,
                        return_all_results = TRUE,
                        parallel = FALSE,
                        ncores = NULL,
                        progress_bar = TRUE,
                        verbose = TRUE) {

  #Check data
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (!inherits(data, "prepared_data")) {
    stop("'data' must be a 'prepared_data' object.")
  }
  if (missing(error_considered)) {
    stop("Argument 'error_considered' must be defined.")
  }
  if (!is.numeric(error_considered)) {
    stop("'error_considered' must be 'numeric' 0-100.")
  }
  if (!is.null(omission_rate)) {
    if (length(omission_rate) > 1) {
      stop("Argument 'omission_rate' must be defined of leght 1.")
    }
    if (!omission_rate %in% error_considered) {
      stop("None of the values in 'error_considered' matches 'omission_rate'.")
    }
  } else {
    if (length(error_considered) > 1) {
      stop("'omission_rate' must be defined if more than one 'error_considered' is used.")
    } else {
      omission_rate <- error_considered
    }
  }


   # Convert calibration data to dataframe if necessary
  if (is.matrix(data$calibration_data) || is.array(data$calibration_data)) {
    data$calibration_data <- as.data.frame(data$calibration_data)
  }

  # If write_summary = TRUE, create directory
  if (write_summary && !file.exists(output_directory)) {
    dir.create(output_directory)
  }

  # Define global vars
  algorithm <- data$algorithm
  formula_grid <- data$formula_grid
  weights <- data$weights

  # Skip existing models if requested
  if (skip_existing_models && write_summary) {
    ready_models <- list.files(path = output_directory, pattern = "summary",
                               full.names = TRUE)
    ready_models <- do.call("rbind",
                            lapply(seq_along(ready_models), function(i) {
                              utils::read.csv(ready_models[i])
    }))
    run_models <- setdiff(formula_grid$ID, ready_models$ID)

    if (length(run_models) == 0) {
      stop(paste("All models completed. Check the folder:", output_directory))
    } else {
      formula_grid <- formula_grid[formula_grid$ID %in% run_models, ]
    }
  }

  # Warning about samples added to background when weights are null
  if (!is.null(use_weights)) {
    if (verbose & addsamplestobackground & use_weights & is.null(weights)) {
      message("Weights for samples added to background are the same as in samples.")
    }

    # Validate weights
    if (use_weights && is.null(weights)) {
      message("'use_weights' = TRUE, but weights are not present in 'data'.\n",
              "Setting 'use_weights' = FALSE.")
      use_weights <- FALSE
    }

    if (use_weights && length(weights) != nrow(data$calibration_data)) {
      stop("In 'data', 'weights' length does not match rows in 'calibration_data'.")
    }
  }

  # Parallelization setup
  if (parallel) {
    if (is.null(ncores)) {
      ncores <- max(1, parallel::detectCores() - 1)
    }
    cl <- parallel::makeCluster(ncores)
  }

  # Task 1: Checking concave curves in quadratic models_________________________
  # ____________________________________________________________________________

  if (remove_concave) {
    if (verbose) {
      message("Task 1/2: checking for concave responses in models:")
    }

    q_grids <- formula_grid[grepl("q", formula_grid$Features), ]
    n_tot <- nrow(q_grids)

    if (n_tot == 0) {
      message("None of the models include quadratic terms.")
    } else {
      if (progress_bar) {
        pb <- txtProgressBar(min = 0, max = n_tot, style = 3)
        progress <- function(n) {
          setTxtProgressBar(pb, n)
        }
        opts <- list(progress = progress)
      } else {
        opts <- NULL
      }


      # Execute fit_eval_concave in parallel or sequentially
      if (parallel) {
        doSNOW::registerDoSNOW(cl)

        results_concave <- foreach::foreach(
          x = 1:n_tot,
          .packages = c("glmnet", "enmpa"),
          .options.snow = opts) %dopar% {
            fit_eval_concave(x = x, q_grids = q_grids, data = data,
                             formula_grid = formula_grid,
                             error_considered = error_considered,
                             omission_rate = omission_rate,
                             write_summary = write_summary,
                             addsamplestobackground = addsamplestobackground,
                             weights = weights,
                             return_all_results = return_all_results,
                             algorithm = algorithm,
                             proc_for_all = proc_for_all,
                             out_dir = output_directory)
          }
      } else {
        results_concave <- vector("list", length = n_tot)
        for (x in 1:n_tot) {
          results_concave[[x]] <-
            fit_eval_concave(x = x, q_grids = q_grids, data = data,
                             formula_grid = formula_grid,
                             error_considered = error_considered,
                             omission_rate = omission_rate,
                             write_summary = write_summary,
                             addsamplestobackground = addsamplestobackground,
                             weights = weights,
                             return_all_results = return_all_results,
                             algorithm = algorithm,
                             proc_for_all = proc_for_all,
                             out_dir = output_directory)
          # Sets the progress bar to the current state
          if (progress_bar) setTxtProgressBar(pb, x)
        }
      }
    } # End of if (n > 0)
  } # End of If remove_concave = TRUE

  # Update formula grid after concave test
  if (!remove_concave) {
    n_tot <- 0
  }

  if (remove_concave & n_tot > 0) {

    # Convert results to dataframe
    d_concave_rep <- do.call("rbind", lapply(results_concave,
                                             function(x) x$Replicates))
    row.names(d_concave_rep) <- NULL
    d_concave_sum <- do.call("rbind", lapply(results_concave,
                                             function(x) x$Summary))
    formula_grid <- formula_grid[!(formula_grid$ID %in% d_concave_sum$ID), ]
  }

  # Task 2: Fitting remaining models ___________________________________________
  # ____________________________________________________________________________

  n_tot <- nrow(formula_grid)

  if (n_tot == 0) {
    message("All candidate models have been tested in task 1.")
  } else {

    if (verbose) {
      if (remove_concave) {
        message("\n\nTask 2/2: fitting and evaluating models with no concave responses:")
      } else {
        message("Task 1/1: fitting and evaluating models:")
      }
    }

    if (progress_bar) {
      pb <- txtProgressBar(0, n_tot, style = 3)
      progress <- function(n) {
        setTxtProgressBar(pb, n)
      }
      opts <- list(progress = progress)
    } else {
      opts <- NULL
    }

    if (parallel) {
      doSNOW::registerDoSNOW(cl)

      results <- foreach(
        x = 1:n_tot,
        .packages = c("glmnet", "enmpa"),
        .options.snow = opts
      ) %dopar% {
          fit_eval_models(x = x, formula_grid = formula_grid, data = data,
                          error_considered = error_considered,
                          omission_rate = omission_rate,
                          write_summary = write_summary,
                          addsamplestobackground = addsamplestobackground,
                          weights = weights,
                          return_all_results = return_all_results,
                          algorithm = algorithm,
                          proc_for_all = proc_for_all,
                          out_dir = output_directory)
        }
    } else {
      results <- vector("list", length = n_tot)
      for (x in 1:n_tot) {
        results[[x]] <-
          fit_eval_models(x = x, formula_grid = formula_grid, data = data,
                          error_considered = error_considered,
                          omission_rate = omission_rate,
                          write_summary = write_summary,
                          addsamplestobackground =  addsamplestobackground,
                          weights = weights,
                          return_all_results = return_all_results,
                          algorithm = algorithm,
                          proc_for_all = proc_for_all,
                          out_dir = output_directory)

        if (progress_bar) {
          setTxtProgressBar(pb, x)
        }
      }
    }

    # Stop cluster
    if (parallel) {
      parallel::stopCluster(cl)
    }

    d_res_rep <- do.call("rbind", lapply(results, function(x) x$All_results))
    row.names(d_res_rep) <- NULL
    d_res_sum <- do.call("rbind", lapply(results, function(x) x$Summary))

    # Join results with results concave, if it exists
    if (remove_concave) {
      replicates_final <- rbind(d_concave_rep, d_res_rep)
      summary_final <- rbind(d_concave_sum, d_res_sum)
      res_final <- list(All_results = replicates_final,
                        Summary = summary_final)
    } else {
      res_final <- list(All_results =  d_res_rep, Summary = d_res_sum)
    }
  } # End of if (n == 0)

  if (n_tot == 0) {
    res_final <- list(All_results = d_concave_rep,
                      Summary = d_concave_sum)
  }

  # Select the best models
  if (verbose) {
    message("\n\nModel selection step:")
  }

  bm <- select_models(candidate_models = res_final$Summary,
                      remove_concave = remove_concave,
                      compute_proc = !proc_for_all,
                      data = data,
                      omission_rate = omission_rate,
                      allow_tolerance = allow_tolerance,
                      tolerance = tolerance,
                      significance = 0.05, verbose = verbose,
                      delta_aic = delta_aic,
                      parallel = parallel, ncores = ncores,
                      progress_bar = progress_bar)

  # Concatenate final results
  fm <- new_calibration_results(
    prepared_data = data, calibration_results = list(res_final),
    omission_rate = omission_rate,
    addsamplestobackground = addsamplestobackground,
    weights = weights, selected_models = list(bm$selected_models),
    summary = list(bm$summary)
  )

  class(fm) <- "calibration_results"
  return(fm)
}

