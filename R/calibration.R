#' Fitting and evaluation of models, and selection of the best ones
#'
#' @description
#' This function fits and evaluates candidate models using the data and grid of
#' formulas prepared with [prepare_data()]. It supports both
#' algorithms `glm` and `maxnet`. The function then selects the best models
#' based on unimodality (optional), partial ROC, omission rate, and AIC values.
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
#' The default, NULL, uses the value in `error_considered`. If more than one
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
#' @param parallel (logical) whether to fit the candidate models in parallel.
#' Default is FALSE.
#' @param ncores (numeric) number of cores to use for parallel processing.
#' Default is NULL and uses available cores - 1. This is only applicable if
#' `parallel = TRUE`.
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is TRUE.
#' @param progress_type (character) type of the progress bar. See
#' `pbapply::pbtypes()` for available progress bar types depending on operating
#' system. Default is "timer".
#' @param use_lb (logical) whether to use load balancing when running in
#' parallel. Default is FALSE.
#' @param verbose (logical) whether to display messages during processing.
#' Default is TRUE.
#'
#' @importFrom parallel makeCluster stopCluster detectCores clusterCall
#' @importFrom utils read.csv head write.csv
#' @importFrom stats aggregate glm as.formula model.matrix sd
#' @importFrom glmnet glmnet.control glmnet
#' @importFrom enmpa predict_glm
#' @importFrom fpROC auc_metrics
#' @importFrom pbapply pboptions pblapply
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
#' metrics for all partitions (if `return_all_results = TRUE`) and a summary of
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
#' # Import occurrences
#' data(occ_data, package = "kuenm2")
#'
#' # Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#' # Use only continuous variables
#' var <- var[[c("bio_1", "bio_7", "bio_12", "bio_15")]]
#'
#' # Prepare data for maxnet model
#' sp_swd <- prepare_data(algorithm = "maxnet", occ = occ_data,
#'                        x = "x", y = "y",
#'                        raster_variables = var,
#'                        species = occ_data[1, 1],
#'                        n_background = 100,
#'                        features = c("l", "lq"),
#'                        r_multiplier = 1,
#'                        partition_method = "kfolds")
#' # Model calibration (maxnet)
#' m <- calibration(data = sp_swd, error_considered = 10)
#' m
#'
#' # Prepare data for glm model
#' sp_swd_glm <- prepare_data(algorithm = "glm", occ = occ_data,
#'                            x = "x", y = "y",
#'                            raster_variables = var,
#'                            species = occ_data[1, 1],
#'                            n_background = 100,
#'                            features = c("l", "lq"),
#'                            partition_method = "kfolds")
#' m_glm <- calibration(data = sp_swd_glm, error_considered = 10)
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
                        progress_type = "timer",
                        use_lb = FALSE,
                        verbose = TRUE) {

  # Check data
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
    if (length(omission_rate) > 1L) {
      stop("Argument 'omission_rate' must be defined of length 1.")
    }
    if (!omission_rate %in% error_considered) {
      stop("None of the values in 'error_considered' matches 'omission_rate'.")
    }
  } else {
    if (length(error_considered) > 1L) {
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
  if (write_summary) {
    if (is.null(output_directory)) {
      stop("If 'write_summary = TRUE', 'output_directory' must be defined.")
    }
    if (!dir.exists(output_directory)) {
      dir.create(output_directory, recursive = TRUE)
    }
  }

  # Define global vars
  algorithm <- data$algorithm
  formula_grid <- data$formula_grid
  weights <- data$weights

  # Skip existing models if requested
  if (skip_existing_models && write_summary) {
    ready_files <- list.files(path = output_directory, pattern = "summary",
                              full.names = TRUE)
    if (length(ready_files) > 0L) {
      ready_models <- do.call("rbind", lapply(ready_files, utils::read.csv))
      run_models <- setdiff(formula_grid$ID, ready_models$ID)

      if (length(run_models) == 0L) {
        stop(paste("All models completed. Check the folder:", output_directory))
      } else {
        formula_grid <- formula_grid[formula_grid$ID %in% run_models, ]
      }
    }
  }

  # Validate and apply use_weights properly
  if (!is.null(use_weights)) {
    if (!is.logical(use_weights) || length(use_weights) != 1L) {
      stop("Argument 'use_weights' must be NULL or a single 'logical'.")
    }
    if (!use_weights) {
      weights <- NULL
    } else {
      if (is.null(weights)) {
        if (verbose) {
          message("'use_weights' = TRUE, but weights are not present in 'data'.\n",
                  "Setting 'use_weights' = FALSE.")
        }
        use_weights <- FALSE
      } else {
        if (verbose && addsamplestobackground) {
          message("Weights for samples added to background are the same as in samples.")
        }
        if (length(weights) != nrow(data$calibration_data)) {
          stop("In 'data', 'weights' length does not match rows in 'calibration_data'.")
        }
      }
    }
  } else if (!is.null(weights)) {
    if (length(weights) != nrow(data$calibration_data)) {
      stop("In 'data', 'weights' length does not match rows in 'calibration_data'.")
    }
  }

  # Configure glmnet control once in main session
  if (algorithm == "maxnet") {
    glmnet::glmnet.control(pmin = 1.0e-8, fdev = 0)
  }

  # Store and set progress bar options safely for CRAN
  old_pbopt <- pbapply::pboptions(
    type = if (progress_bar) progress_type else "none",
    use_lb = use_lb
  )
  on.exit(pbapply::pboptions(old_pbopt), add = TRUE)

  # Setup parallel cluster once for Tasks 1 and 2
  if (parallel) {
    if (is.null(ncores)) {
      ncores <- max(1L, parallel::detectCores() - 1L, na.rm = TRUE)
      chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
      if (nzchar(chk) && chk == "TRUE") {
        ncores <- min(ncores, 2L)
      }
    }
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    if (algorithm == "maxnet") {
      parallel::clusterCall(cl, glmnet::glmnet.control, pmin = 1.0e-8, fdev = 0)
    }
  } else {
    cl <- NULL
  }

  # Initialize objects for concave results
  d_concave_rep <- NULL
  d_concave_sum <- NULL

  # Task 1: Checking concave curves in quadratic models ________________________
  if (remove_concave) {
    if (verbose) {
      message("Task 1/2: checking for concave responses in models:")
    }

    if (algorithm == "maxnet") {
      max_reg <- max(formula_grid$R_multiplier)
      q_grids <- formula_grid[grepl("q", formula_grid$Features) &
                                formula_grid$R_multiplier == max_reg, ]
    } else {
      q_grids <- formula_grid[grepl("q", formula_grid$Features), ]
    }

    n_tot_concave <- nrow(q_grids)

    if (n_tot_concave == 0L) {
      if (verbose) message("None of the models include quadratic terms.")
    } else {
      results_concave <- pbapply::pblapply(
        X = seq_len(n_tot_concave),
        FUN = fit_eval_concave,
        q_grids = q_grids,
        data = data,
        formula_grid = formula_grid,
        error_considered = error_considered,
        omission_rate = omission_rate,
        write_summary = write_summary,
        addsamplestobackground = addsamplestobackground,
        weights = weights,
        return_all_results = return_all_results,
        algorithm = algorithm,
        proc_for_all = proc_for_all,
        out_dir = output_directory,
        cl = cl
      )

      d_concave_rep <- do.call("rbind", lapply(results_concave,
                                               function(x) x$All_results))
      if (!is.null(d_concave_rep)) {
        row.names(d_concave_rep) <- NULL
      }
      d_concave_sum <- do.call("rbind", lapply(results_concave,
                                               function(x) x$Summary))
      formula_grid <- formula_grid[!(formula_grid$ID %in% d_concave_sum$ID), ]
    }
  }

  # Task 2: Fitting remaining models __________________________________________
  n_tot <- nrow(formula_grid)

  if (n_tot == 0L) {
    if (verbose) message("All candidate models have been tested in task 1.")
    res_final <- list(All_results = d_concave_rep,
                      Summary = d_concave_sum)
  } else {
    if (verbose) {
      if (remove_concave) {
        message("\n\nTask 2/2: fitting and evaluating models with no concave responses:")
      } else {
        message("Task 1/1: fitting and evaluating models:")
      }
    }

    if (algorithm == "maxnet") {
      unique_forms <- unique(formula_grid$Formulas)
      form_groups <- lapply(unique_forms, function(f_str) {
        formula_grid[formula_grid$Formulas == f_str, , drop = FALSE]
      })

      results <- pbapply::pblapply(
        X = form_groups,
        FUN = fit_eval_formula_group,
        data = data,
        error_considered = error_considered,
        omission_rate = omission_rate,
        write_summary = write_summary,
        addsamplestobackground = addsamplestobackground,
        weights = weights,
        return_all_results = return_all_results,
        proc_for_all = proc_for_all,
        out_dir = output_directory,
        cl = cl
      )
    } else {
      results <- pbapply::pblapply(
        X = seq_len(n_tot),
        FUN = fit_eval_models,
        formula_grid = formula_grid,
        data = data,
        error_considered = error_considered,
        omission_rate = omission_rate,
        write_summary = write_summary,
        addsamplestobackground = addsamplestobackground,
        weights = weights,
        return_all_results = return_all_results,
        algorithm = algorithm,
        proc_for_all = proc_for_all,
        out_dir = output_directory,
        cl = cl
      )
    }

    d_res_rep <- do.call("rbind", lapply(results, function(x) x$All_results))
    if (!is.null(d_res_rep)) {
      row.names(d_res_rep) <- NULL
    }
    d_res_sum <- do.call("rbind", lapply(results, function(x) x$Summary))

    replicates_final <- rbind(d_concave_rep, d_res_rep)
    summary_final <- rbind(d_concave_sum, d_res_sum)

    if (!is.null(summary_final) && "ID" %in% colnames(summary_final)) {
      summary_final <- summary_final[order(summary_final$ID), ]
      row.names(summary_final) <- NULL
    }
    if (!is.null(replicates_final) && "ID" %in% colnames(replicates_final)) {
      replicates_final <- replicates_final[order(replicates_final$ID,
                                                 replicates_final$Partition), ]
      row.names(replicates_final) <- NULL
    }

    res_final <- list(All_results = replicates_final,
                      Summary = summary_final)
  }

  # Select the best models
  if (verbose) {
    message("\n\nModel selection step:")
  }

  bm <- select_models(candidate_models = res_final$Summary,
                      data = data,
                      algorithm = algorithm,
                      compute_proc = !proc_for_all,
                      addsamplestobackground = addsamplestobackground,
                      weights = weights,
                      remove_concave = remove_concave,
                      omission_rate = omission_rate,
                      allow_tolerance = allow_tolerance,
                      tolerance = tolerance,
                      significance = 0.05,
                      delta_aic = delta_aic,
                      parallel = parallel,
                      ncores = ncores,
                      progress_bar = progress_bar,
                      progress_type = progress_type,
                      use_lb = use_lb,
                      verbose = verbose)

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
