#' Fitting and evaluation of models, and selection of the best ones
#'
#' @description
#' This function fits and evaluates candidate models using the data and grid of
#' formulas prepared with \code{\link{prepare_data}}. It supports both
#' algorithms `glm` and `maxnet`. The function then selects the best models
#' based on unimodality (optional), partial ROC, omission rate, and AIC values.
#'
#' @usage
#' calibration(data, addsamplestobackground = TRUE,use_weights = FALSE,
#'             parallel = FALSE, ncores = 4, parallel_option = "doSNOW",
#'             progress_bar = TRUE, write_summary = FALSE,
#'             output_directory = NULL, skip_existing_models = FALSE,
#'             return_replicate = TRUE, test_concave = FALSE,
#'             omission_rate = 10, omrat_threshold = 10,
#'             AIC_option = "ws", delta_aic = 2, allow_tolerance = TRUE,
#'             tolerance = 0.01, verbose = TRUE)
#'
#' @param data an object of class `prepared_data` returned by the
#' \code{\link{prepare_data()}} function. It contains the calibration data,
#' formulas grid, kfolds, and model type.
#' @param test_concave (logical) whether to test for and remove candidate models
#' presenting concave curves. Default is FALSE.
#' @param addsamplestobackground (logical) whether to add to the background any
#' presence sample that is not already there. Default is TRUE.
#' @param use_weights (logical) whether to apply the weights present in the
#' data. Default is FALSE.
#' @param parallel (logical) whether to fit the candidate models in parallel.
#' Default is FALSE.
#' @param ncores (numeric) number of cores to use for parallel processing.
#' Default is 1. This is only applicable if `parallel = TRUE`.
#' @param parallel_option (character) the package to use for parallel processing:
#' "doParallel" or "doSNOW". Default is "doSNOW". This is only applicable if
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
#' @param return_replicate (logical) whether to return the evaluation results
#' for each replicate. Default is TRUE, meaning evaluation results for each
#' replicate will be returned.
#' @param omission_rate (numeric) values from 0 to 100 representing the
#' percentage of potential error due to any source of uncertainty. This value is
#' used to calculate the omission rate. Default is 10. See details.
#' @param omrat_threshold (numeric) the maximum omission rate a candidate model
#' can have to be considered a best model. This value must match one of the
#' values specified in omrat. Defaut is 10.
#' @param AIC_option (character) the type of AIC to be calculated: "ws" for AIC
#' proposed by Warren and Seifert (2011), or "nk" for AIC proposed by Ninomiya
#' and Kawano (2016). Default is "ws". See References for details.
#' @param delta_aic (numeric) the value of delta AIC used as a threshold to
#' select models. Default is 2.
#' @param allow_tolerance (logical) whether to allow selection of models with
#' minimum values of omission rates even if their omission rate surpasses the
#' `omrat_threshold`. This is only applicable if all candidate models have
#' omission rates higher than the `omrat_threshold`. Default is TRUE.
#' @param tolerance (numeric) The value added to the minimum omission rate if it
#' exceeds the `omrat_threshold`. If `allow_tolerance = TRUE`, selected models
#' will have an omission rate equal to or less than the minimum rate plus this
#' tolerance. Default is 0.01.
#' @param verbose (logical) whether to display messages during processing.
#' Default is TRUE.
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach `%dopar%`
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats aggregate
#' @importFrom glmnet glmnet.control glmnet
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
#' of class "prcomp". See ?stats::prcomp() for details.
#' - algorithm: the model type (glm or maxnet)
#' - calibration_results: a list containing a data frame with all evaluation
#' metrics for all replicates (if `return_replicate = TRUE`) and a summary of
#' the evaluation metrics for each candidate model.
#' - omission_rate: The omission rate determined by `omrat_threshold` for
#' selecting best models.
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
#' Partial ROC is calculated following Peterson et al.
#' (2008; http://dx.doi.org/10.1016/j.ecolmodel.2007.11.008).
#'
#' Omission rates are calculated using models trained with separate testing data
#' subsets. Users can specify multiple omission rates to be calculated
#' (e.g., c(5, 10)), though only one can be used as the threshold for selecting
#' the best models.
#'
#' Model complexity (AIC) is assessed using models generated with the complete
#' set of occurrences.
#'
#' @examples
#' # Import occurrences
#' data(occ_data, package = "kuenm2")
#'
#' # Import variables
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' # Use only variables 1, 2 and 3
#' var <- var[[1:3]]
#'
#' #### maxnet ####
#' # Prepare data for maxnet model
#' sp_swd <- prepare_data(algorithm = "maxnet", occ = occ_data,
#'                        species = occ_data[1, 1], x = "x", y = "y",
#'                        raster_variables = var,
#'                        n_background = 100,
#'                        features = c("l", "lq"),
#'                        reg_mult = 1)
#'
#' # Calibrate maxnet models
#' m <- calibration(data = sp_swd, omission_rate = c(5, 10))
#'
#' m
#'
#' #### GLM ####
#' # Prepare data for glm model
#' sp_swd_glm <- prepare_data(algorithm = "glm", occ = occ_data,
#'                            species = occ_data[1, 1], x = "x", y = "y",
#'                            raster_variables = var,
#'                            n_background = 100,
#'                            features = c("l", "lq", "q", "lqp"))
#'
#' # Calibrate glm models
#' m_glm <- calibration(data = sp_swd_glm, omission_rate = c(5, 10))
#'
#' m_glm


calibration <- function(data,
                        addsamplestobackground = TRUE,
                        use_weights = FALSE,
                        parallel = FALSE,
                        ncores = 4,
                        parallel_option = "doSNOW",
                        progress_bar = TRUE,
                        write_summary = FALSE,
                        output_directory = NULL,
                        skip_existing_models = FALSE,
                        return_replicate = TRUE,
                        test_concave = FALSE,
                        omission_rate = 10,
                        omrat_threshold = 10,
                        AIC_option = "ws",
                        delta_aic = 2,
                        allow_tolerance = TRUE,
                        tolerance = 0.01,
                        verbose = TRUE) {

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
    ready_models <- do.call("rbind", lapply(seq_along(ready_models), function(i) {
      read.csv(ready_models[i])
    }))
    run_models <- setdiff(formula_grid$ID, ready_models$ID)

    if (length(run_models) == 0) {
      stop(paste("All models completed. Check the folder:", output_directory))
    } else {
      formula_grid <- formula_grid[formula_grid$ID %in% run_models, ]
    }
  }

  # Warning about samples added to background when weights are null
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
    stop("Length of weights does not match number of rows in calibration_data.")
  }

  # Parallelization setup
  if (parallel) {
    if (!(parallel_option %in% c("doSNOW", "doParallel"))) {
      stop("Invalid parallel_option. Use 'doSNOW' or 'doParallel'.")
    }
    cl <- parallel::makeCluster(ncores)
  }

  # Task 1: Checking concave curves in quadratic models_________________________
  # ____________________________________________________________________________

  if (test_concave) {
    if (verbose) {
      message("Task 1/2: checking for concave responses in models:")
    }

    q_grids <- formula_grid[grepl("q", formula_grid$Features), ]
    n_tot <- nrow(q_grids)

    if(n_tot == 0) {
      message("None of the models include quadratic terms")
    } else {
      if (progress_bar) {
        pb <- txtProgressBar(min = 0, max = n_tot, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)} else {opts <- NULL}

      if (parallel & parallel_option == "doParallel") {
        doParallel::registerDoParallel(cl)
        opts <- NULL # Progress bar does not work with doParallel
      }

      if (parallel & parallel_option == "doSNOW") {
        doSNOW::registerDoSNOW(cl)
        if (isTRUE(progress_bar))
          opts <- list(progress = progress)
        else opts <- opts
      }

      # Execute fit_eval_concave in parallel or sequentially
      if(parallel){
        results_concave <- foreach::foreach(
          x = 1:n_tot,
          .packages = c("glmnet", "enmpa"),
          .options.snow = opts) %dopar% {
            fit_eval_concave(x = x, q_grids, data, formula_grid,
                             omission_rate = omission_rate,
                             omrat_thr = omrat_threshold,
                             write_summary = write_summary,
                             addsamplestobackground = addsamplestobackground,
                             weights = weights,
                             return_replicate = return_replicate,
                             algorithm = algorithm, AIC_option = AIC_option)
          }
      } else {
        results_concave <- vector("list", length = n_tot)
        for (x in 1:n_tot) {
          results_concave[[x]] <-
            fit_eval_concave(x = x, q_grids, data, formula_grid,
                             omission_rate = omission_rate,
                             omrat_thr = omrat_threshold,
                             write_summary = write_summary,
                             addsamplestobackground = addsamplestobackground,
                             weights = weights,
                             return_replicate = return_replicate,
                             algorithm = algorithm, AIC_option = AIC_option
            )
          # Sets the progress bar to the current state
          if(progress_bar) setTxtProgressBar(pb, x)
        }
      }
    } # End of if(n > 0)
  } # End of If test_concave = TRUE

  # Update formula grid after concave test
  if(!test_concave) {n_tot = 0}
  if(test_concave & n_tot > 0) {

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

  if(n_tot == 0) {
    message("All candidate models have been tested in task 1")
  } else {

    if(verbose) {
      if(test_concave) {
        message("\n\nTask 2/2: fitting and evaluating models with no concave responses:")
      } else {
        message("Task 1/1: fitting and evaluating models:")
      }
    }

    if (progress_bar) {
      pb <- txtProgressBar(0, n_tot, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n) }

    if (parallel_option == "doParallel") {
      doParallel::registerDoParallel(cl)
      opts <- NULL
    }

    if (parallel_option == "doSNOW") {
      doSNOW::registerDoSNOW(cl)
      if (isTRUE(progress_bar)) {
        opts <- list(progress = progress)
      } else {
        opts <- NULL
      }
    }

    if (parallel) {
      results <- foreach(
        x = 1:n_tot,
        .packages = c("glmnet", "enmpa"),
        .options.snow = opts ) %dopar% {
          fit_eval_models(x, formula_grid, data,
                          omission_rate = omission_rate,
                          omrat_thr = omrat_threshold,
                          write_summary = write_summary,
                          addsamplestobackground = addsamplestobackground,
                          weights = weights,
                          return_replicate = return_replicate, AIC_option = AIC_option,
                          algorithm = algorithm)
        }
    } else {
      results <- vector("list", length = n_tot)
      for (x in 1:n_tot) {
        results[[x]] <-
          fit_eval_models(x, formula_grid = formula_grid, data = data,
                          omission_rate = omission_rate,
                          omrat_thr = omrat_threshold,
                          addsamplestobackground =  addsamplestobackground,
                          weights = weights, write_summary,
                          return_replicate, AIC_option = AIC_option,
                          algorithm = algorithm)

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
    if (test_concave) {
      replicates_final <- rbind(d_concave_rep, d_res_rep)
      summary_final <- rbind(d_concave_sum, d_res_sum)
      res_final <- list(All_results = replicates_final,
                        Summary = summary_final)
    } else {
      res_final <- list(All_results =  d_res_rep, Summary = d_res_sum)
    }
  }# End of if(n == 0)

  if (n_tot == 0) {
    res_final <- list(All_results = d_concave_rep,
                      Summary = d_concave_sum)
  }

  # Select the best models
  if (verbose) {
    message("\n\nModel selection step:")
  }

  bm <- sel_best_models(cand_models = res_final$Summary,
                        test_concave = test_concave,
                        omrat_threshold = omrat_threshold,
                        allow_tolerance = allow_tolerance,
                        tolerance = tolerance, AIC_option = AIC_option,
                        significance = 0.05, verbose = verbose,
                        delta_aic = delta_aic,
                        algorithm = algorithm)

  # Concatenate final results
  fm <- new_calibration_results(
    prepared_data = data, calibration_results = list(res_final),
    omission_rate = omrat_threshold,
    addsamplestobackground = addsamplestobackground,
    weights = weights, selected_models = list(bm$cand_final),
    summary = list(bm$summary)
  )

  class(fm) <- "calibration_results"
  return(fm)
}

