#' Fit models selected after calibration
#'
#' @description
#' This function fits models selected after candidate model training and testing
#' using the function [calibration()].
#'
#' @usage
#' fit_selected(calibration_results, partition_method = "kfolds",
#'              n_partitions = 1, train_proportion = 0.7, type = "cloglog",
#'              write_models = FALSE,
#'              file_name = NULL, parallel = FALSE, ncores = NULL,
#'              progress_bar = TRUE, verbose = TRUE, seed = 1)
#'
#' @param calibration_results an object of class `calibration_results` returned
#' by the [calibration()] function.
#' @param partition_method (character) method used for data partitioning.
#' Available options are `"kfolds"`, `"subsample"`, and `"bootstrap"`.
#' See **Details** for more information.
#' @param n_partitions (numeric) number of partitions or folds to generate. If
#' `partition_method` is `"subsample"` or `"bootstrap"`, this defines the number
#' of replicates. If `"kfolds"`, it specifies the number of folds. Default is 4.
#' @param train_proportion (numeric) proportion of occurrence and background
#' points to be used for model training in each replicate. Only applicable when
#'  `partition_method` is `"subsample"` or `"bootstrap"`. Default is 0.7 (i.e.,
#'  70% for training and 30% for testing).
#' @param type (character) the format of prediction values for computing
#' thresholds. For maxnet models, valid options are "raw", "cumulative",
#' "logistic", and "cloglog". For glm models, valid options are "cloglog",
#' "response" and "raw". Default is "cloglog".
#' @param write_models (logical) whether to save the final fitted models to disk.
#' Default is FALSE.
#' @param file_name (character) the file name, with or without a path, for saving
#' the final models. This is only applicable if `write_models = TRUE`.
#' @param parallel (logical) whether to fit the final models in parallel. Default
#' is FALSE.
#' @param ncores (numeric) number of cores to use for parallel processing.
#' Default is NULL and uses available cores - 1. This is only applicable if
#' `parallel = TRUE`.
#' @param progress_bar (logical) whether to display a progress bar during processing.
#' Default is TRUE.
#' @param verbose (logical) whether to display detailed messages during processing.
#' Default is TRUE.
#' @param seed (numeric) integer value used to specify an initial seed to split
#' the data. Default is 1.
#'
#' @details
#' This function also computes model consensus (mean and median), the thresholds
#' to binarize model predictions based on the omission rate set during model
#' calibration to select models.
#'
#'
#' @return
#' An object of class 'fitted_models' containing the following elements:
#' \item{species}{a character string with the name of the species.}
#' \item{Models}{a list of fitted models, including partitions (trained with
#' the parts of the data) and full models (trained with all available records).}
#' \item{calibration_data}{a data.frame containing a column (`pr_bg`) that
#' identifies occurrence points (1) and background points (0), along with the
#' corresponding values of predictor variables for each point.}
#' \item{selected_models}{a data frame with the ID and summary of evaluation
#' metrics for the selected models.}
#' \item{weights}{a numeric vector specifying weights for the predictor variables
#' (if used).}
#' \item{pca}{a list of class \code{\link[terra]{prcomp}} representing the
#' result of principal component analysis (if performed).}
#' \item{addsamplestobackground}{a logical value indicating whether any presence
#' sample not already in the background was added.}
#' \item{omission_rate}{the omission rate determined during the calibration step.}
#' \item{thresholds}{the thresholds to binarize each partition and the consensus
#' (mean and median), calculated based on the omission rate set in
#' [calibration()].}
#'
#' @export
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom foreach foreach `%dopar%`
#' @importFrom glmnet glmnet.control glmnet
#' @importFrom enmpa predict_glm
#'
#' @examples
#' # An example with maxnet models
#' data(calib_results_maxnet, package = "kuenm2")
#'
#' # Fit models using calibration results
#' fm <- fit_selected(calibration_results = calib_results_maxnet,
#'                    n_partitions = 4)
#'
#' # Output the fitted models
#' fm
#'
#' # An example with GLMs
#' data(calib_results_glm, package = "kuenm2")
#'
#' # Fit models using calibration results
#' fm_glm <- fit_selected(calibration_results = calib_results_glm,
#'                        partition_method = "subsample",
#'                        n_partitions = 5)
#'
#' # Output the fitted models
#' fm_glm

fit_selected <- function(calibration_results,
                         partition_method = "kfolds",
                         n_partitions = 1,
                         train_proportion = 0.7,
                         type = "cloglog",
                         write_models = FALSE,
                         file_name = NULL,
                         parallel = FALSE,
                         ncores = NULL,
                         progress_bar = TRUE,
                         verbose = TRUE,
                         seed = 1) {

  if (missing(calibration_results)) {
    stop("Argument 'calibration_results' must be defined.")
  }
  if (!inherits(calibration_results, "calibration_results")) {
    stop("Argument 'calibration_results' must be a 'calibration_results' object.")
  }

  # Extract model IDs from selected models
  m_ids <- calibration_results$selected_models$ID
  algorithm <- calibration_results$algorithm

  # Fitting models over multiple partitions
  if (n_partitions > 1) {
    if (verbose) {
      message("Fitting partitions...")
    }

    # Create a grid of model IDs and partitions
    dfgrid <- expand.grid(models = m_ids, partitions = 1:n_partitions)
    n_tot <- nrow(dfgrid) # Total models * partitions

    #Prepare data (index) to partitions
    if (n_partitions > 1) {
      #Partitioning data
      rep_data <- part_data(data = calibration_results$calibration_data,
                            pr_bg = "pr_bg",
                            train_proportion = train_proportion,
                            n_partitions = n_partitions,
                            partition_method = partition_method, seed = seed)
    } else {
      rep_data <- NULL
    }


    # Adjust parallelization based on number of tasks and cores
    if (n_tot == 1 & parallel) {
      parallel <- FALSE
    }

    # Progress bar
    if (isTRUE(progress_bar)) {
      pb <- utils::txtProgressBar(0, n_tot, style = 3)
      progress <- function(n) utils::setTxtProgressBar(pb, n)
    }


    # Setup parallel cluster
    if (parallel) {
      if (is.null(ncores)) {
        ncores <- max(1, parallel::detectCores() - 1)
      }
      if (n_tot < ncores) {
        ncores <- n_tot
      }

      cl <- parallel::makeCluster(ncores)
      doSNOW::registerDoSNOW(cl)

      opts <- if (progress_bar) {
        list(progress = progress)
      } else {
        NULL
      }
    }

    # Fit the best models (either in parallel or sequentially)
    if (parallel) {
      best_models <- foreach::foreach(x = 1:n_tot,
                                      .packages = c("glmnet", "enmpa"),
                                      .options.snow = opts
      ) %dopar% {
        fit_best_model(x = x, dfgrid = dfgrid, cal_res = calibration_results,
                       n_partitions = n_partitions, rep_data = rep_data,
                       algorithm = algorithm)
      }
    } else {
      best_models <- vector("list", length = n_tot)
      for (x in 1:n_tot) {
        best_models[[x]] <- fit_best_model(
          x = x, dfgrid = dfgrid, cal_res = calibration_results,
          n_partitions = n_partitions, rep_data = rep_data,
          algorithm = algorithm
        )
        if (progress_bar) {
          utils::setTxtProgressBar(pb, x)
        }
      }
    }

    # Stop the cluster
    if (parallel) parallel::stopCluster(cl)

    # Split models by their respective IDs
    best_models <- split(best_models, dfgrid$models)
    best_models <- lapply(best_models, function(sublist) {
      names(sublist) <- paste0("Partition_", seq_along(sublist))
      return(sublist)
    })

    # Assign names to the models based on their IDs
    names(best_models) <- paste0("Model_", names(best_models))

  } else {
    best_models <- list()
  }

  # Fit full models ____________________________________________________________
  if (verbose) {
    message("\nFitting full models...")
  }
  # Full models grid setup
  n_models <- length(m_ids)
  dfgrid <- expand.grid(models = m_ids, partitions = 1)

  # Adjust parallelization for full models
  if (n_models == 1 & parallel) {
    parallel <- FALSE
  }

  # Progress bar setup for full models
  if (progress_bar) {
    pb <- utils::txtProgressBar(0, n_models, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
  }

  # Parallel cluster setup for full models
  if (parallel) {
    if (is.null(ncores)) {
      ncores <- max(1, parallel::detectCores() - 1)
    }
    if (n_models < ncores & parallel) {
      ncores <- n_models
    }
    cl <- parallel::makeCluster(ncores)
    doSNOW::registerDoSNOW(cl)

    opts <- if (progress_bar) {
      list(progress = progress)
    } else {
      NULL
    }
  }

  # Fit the full models (either in parallel or sequentially)
  if (parallel) {
    full_models <- foreach::foreach(x = 1:n_models,
                                    .packages = c("glmnet", "enmpa"),
                                    .options.snow = opts
    ) %dopar% {
      fit_best_model(x = x, dfgrid = dfgrid, cal_res = calibration_results,
                     n_partitions = 1, rep_data = rep_data,
                     algorithm = algorithm)
    }
  } else {
    full_models <- vector("list", length = n_models)
    for (x in 1:n_models) {
      full_models[[x]] <- fit_best_model(
        x = x, dfgrid = dfgrid, cal_res = calibration_results,
        n_partitions = 1, rep_data = rep_data, algorithm = algorithm
      )
      if (progress_bar) utils::setTxtProgressBar(pb, x)
    }
  }

  # Assign names to full models
  names(full_models) <- paste0("Model_", m_ids)

  # Combine partition models with full models
  for (i in names(full_models)) {
    best_models[[i]]$Full_model <- full_models[[i]]
  }

  # Stop the cluster for full models
  if (parallel) {
    parallel::stopCluster(cl)
  }

  # Compute thresholds for predictions
  occ <- calibration_results$calibration_data[
    calibration_results$calibration_data$pr_bg == 1, -1, drop = FALSE]

  # Predictions and consensus for occurrences
  p_occ <- lapply(names(best_models), function(x) {
    m_x <- best_models[[x]]
    if (any(grepl("Rep", names(m_x)))) {
      m_x$Full_model <- NULL
    }

    if (algorithm == "maxnet") {
      p_r <- sapply(m_x, function(i) predict.glmnet_mx(object = i,
                                                       newdata = occ,
                                                       type = type))
    } else if (algorithm == "glm") {
      p_r <- sapply(m_x, function(i) suppressWarnings(
        as.numeric(predict_glm_mx(model = i,
                       newdata = occ, type = type))
        ))
    }

    p_mean <- apply(p_r, 1, mean, na.rm = TRUE)
    p_median <- apply(p_r, 1, median, na.rm = TRUE)

    list(mean = p_mean, median = p_median)
  })

  names(p_occ) <- names(best_models)

  # Calculate consensus across models
  mean_consensus <- apply(sapply(p_occ, function(x) x$mean), 1,
                          mean, na.rm = TRUE)
  median_consensus <- apply(sapply(p_occ, function(x) x$median), 1,
                            median, na.rm = TRUE)
  consensus <- list(mean = mean_consensus, median = median_consensus)
  p_occ <- c(p_occ, list(consensus = consensus))

  # Calculate thresholds
  p_thr <- lapply(p_occ, function(model) {
    lapply(model, calc_thr,
           thr = calibration_results$summary$omission_rate_thr / 100)
  })
  #Append type of predictions
  p_thr$type <- type

  #Prepare final data
  if(partition_method %in% c("kfolds", "leave-one-out")){
    train_proportion <- NULL
  }

  # Prepare final results
  res <- new_fitted_models(
    species = calibration_results$species,
    Models = best_models,
    calibration_data = calibration_results$calibration_data,
    continuous_variables = calibration_results$continuous_variables,
    categorical_variables = calibration_results$categorical_variables,
    selected_models = calibration_results$selected_models,
    weights = calibration_results$weights,
    pca = calibration_results$pca,
    addsamplestobackground = calibration_results$addsamplestobackground,
    omission_rate = calibration_results$summary$omission_rate_thr,
    thresholds = p_thr,
    algorithm = algorithm,
    partition_method = partition_method,
    n_partitions = n_partitions,
    train_proportion = train_proportion
  )

  # Optionally save the fitted models
  if (write_models) {
    saveRDS(res, file = paste0(file_name, ".RDS"))
  }

  return(res)
}
