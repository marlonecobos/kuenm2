#' Fit models selected after calibration
#'
#' @description
#' This function fits models selected during model [calibration()].
#'
#' @usage
#' fit_selected(calibration_results, replicate_method = "kfolds",
#'              n_replicates = 1, sample_proportion = 0.7, type = "cloglog",
#'              write_models = FALSE,
#'              file_name = NULL, parallel = FALSE, ncores = NULL,
#'              progress_bar = TRUE, progress_type = "timer",
#'              use_lb = FALSE, verbose = TRUE, seed = 1)
#'
#' @param calibration_results an object of class `calibration_results` returned
#' by the [calibration()] function.
#' @param replicate_method (character) method used for producing replicates.
#' Available options are `"kfolds"`, `"subsample"`, and `"bootstrap"`.
#' See **Details** for more information.
#' @param n_replicates (numeric) number of replicates or folds to generate. If
#' `replicate_method` is `"subsample"` or `"bootstrap"`, this defines the number
#' of replicates. If `"kfolds"`, it specifies the number of folds. Default is 4.
#' @param sample_proportion (numeric) proportion of occurrence and background
#' points to be used to fit model replicates. Only applicable when
#'  `replicate_method` is `"subsample"` or `"bootstrap"`. Default is 0.7 (i.e.,
#'  70% of the data).
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
#' @param progress_type (character) type of the progress bar. See
#' `pbapply::pbtypes()` for available progress bar types depending on operating
#' system. Default is "timer".
#' @param use_lb (logical) whether to use load balancing when running in
#' parallel. Default is FALSE.
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
#' \item{Models}{a list of fitted models, including replicates (fitted with
#' part of the data) and full models (fitted with all data).}
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
#' \item{thresholds}{the thresholds to binarize each replicate and the consensus
#' (mean and median), calculated based on the omission rate set in
#' [calibration()].}
#'
#' @export
#'
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom glmnet glmnet.control glmnet
#' @importFrom enmpa predict_glm
#' @importFrom pbapply pboptions pblapply
#'
#' @examples
#' # An example with maxnet models
#' data(calib_results_maxnet, package = "kuenm2")
#'
#' # Fit models using calibration results
#' fm <- fit_selected(calibration_results = calib_results_maxnet,
#'                    n_replicates = 4)
#'
#' # Output the fitted models
#' fm
#'
#' # An example with GLMs
#' data(calib_results_glm, package = "kuenm2")
#'
#' # Fit models using calibration results
#' fm_glm <- fit_selected(calibration_results = calib_results_glm,
#'                        replicate_method = "subsample",
#'                        n_replicates = 5)
#'
#' # Output the fitted models
#' fm_glm

fit_selected <- function(calibration_results,
                         replicate_method = "kfolds",
                         n_replicates = 1,
                         sample_proportion = 0.7,
                         type = "cloglog",
                         write_models = FALSE,
                         file_name = NULL,
                         parallel = FALSE,
                         ncores = NULL,
                         progress_bar = TRUE,
                         progress_type = "timer",
                         use_lb = FALSE,
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

  # Fitting models if they have multiple replicates
  if (n_replicates > 1) {
    if (verbose) {
      message("Fitting replicates...")
    }

    # Create a grid of model IDs and replicates
    dfgrid <- expand.grid(models = m_ids, replicates = 1:n_replicates)
    n_tot <- nrow(dfgrid) # Total models * replicates

    #Prepare data (index) to replicates
    if (n_replicates > 1) {
      #Partitioning data
      rep_data <- part_data(data = calibration_results$calibration_data,
                            pr_bg = "pr_bg",
                            train_proportion = sample_proportion,
                            n_partitions = n_replicates,
                            partition_method = replicate_method, seed = seed)
      names(rep_data) <- gsub("Partition", "Replicate", rep_data)
    } else {
      rep_data <- NULL
    }

    if (n_tot == 1 & parallel) {
      parallel <- FALSE
    }

    # Store and set progress bar options safely for CRAN
    old_pbopt <- pbapply::pboptions(
      type = if (progress_bar) progress_type else "none",
      use_lb = use_lb
    )
    on.exit(pbapply::pboptions(old_pbopt), add = TRUE)

    if (parallel) {
      # Parallelization setup
      if (is.null(ncores)) {
        ncores <- max(1L, parallel::detectCores() - 1L)
        # Check if function is running in CRAN...
        chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
        if (nzchar(chk) && chk == "TRUE") {
          ncores <- min(ncores, 2L)
        }
      }
      cl <- parallel::makeCluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
    } else {
      cl <- NULL
    }


    best_models <- pbapply::pblapply(
      X = seq_len(n_tot),
      FUN = fit_best_model,
      dfgrid = dfgrid,
      cal_res = calibration_results,
      n_replicates = n_replicates,
      rep_data = rep_data,
      algorithm = algorithm,
      cl = cl)


    # Split models by their respective IDs
    best_models <- split(best_models, dfgrid$models)
    best_models <- lapply(best_models, function(sublist) {
      names(sublist) <- paste0("Replicate_", seq_along(sublist))
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
  dfgrid <- expand.grid(models = m_ids, replicates = 1)

  if (n_tot == 1 & parallel) {
    parallel <- FALSE
  }

  # Store and set progress bar options safely for CRAN
  old_pbopt <- pbapply::pboptions(
    type = if (progress_bar) progress_type else "none",
    use_lb = use_lb
  )
  on.exit(pbapply::pboptions(old_pbopt), add = TRUE)

  if (parallel) {
    # Parallelization setup
    if (is.null(ncores)) {
      ncores <- max(1L, parallel::detectCores() - 1L)
      # Check if function is running in CRAN...
      chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
      if (nzchar(chk) && chk == "TRUE") {
        ncores <- min(ncores, 2L)
      }
    }
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else {
    cl <- NULL
  }

  # Fit the full models
  full_models <- pbapply::pblapply(
    X = seq_len(n_models),
    FUN = fit_best_model,
    dfgrid = dfgrid, cal_res = calibration_results,
    n_replicates = 1, rep_data = rep_data,
    algorithm = algorithm,
    cl = cl)

  # Assign names to full models
  names(full_models) <- paste0("Model_", m_ids)

  # Combine replicate models with full models
  for (i in names(full_models)) {
    best_models[[i]]$Full_model <- full_models[[i]]
  }


  # Compute thresholds for predictions
  occ <- calibration_results$calibration_data[
    calibration_results$calibration_data$pr_bg == 1, -1, drop = FALSE]

  # Predictions and consensus for occurrences
  p_occ <- lapply(names(best_models), function(x) {
    m_x <- best_models[[x]]
    if (n_replicates > 1) {
      m_x$Full_model <- NULL
    }

    if (algorithm == "maxnet") {
      p_r <- sapply(m_x, function(i) predict.glmnet_mx(object = i,
                                                       newdata = occ,
                                                       type = type))
    } else if (algorithm == "glm") {
      p_r <- sapply(m_x, function(i) suppressWarnings(
        as.numeric(predict_glm_mx(model = i, newdata = occ, type = type))
        ))
    }

    if (n_replicates > 1) {
      p_mean <- apply(p_r, 1, mean, na.rm = TRUE)
      p_median <- apply(p_r, 1, median, na.rm = TRUE)

      list(mean = p_mean, median = p_median, rep = p_r)
    } else {
      list(Full_model = p_r[, 1])
    }
  })

  names(p_occ) <- names(best_models)

  # Calculate consensus across models
  if (length(p_occ) == 1) {
    if (n_replicates > 1) {
      mean_consensus <- p_occ[[1]]$mean
      median_consensus <- p_occ[[1]]$median
    } else {
      mean_consensus <- p_occ[[1]]$Full_model
      median_consensus <- p_occ[[1]]$Full_model
    }

  } else {
    if (n_replicates > 1) {
      mean_consensus <- apply(sapply(p_occ, function(x) x$mean), 1,
                              mean, na.rm = TRUE)
      median_consensus <- apply(do.call(cbind, lapply(p_occ, `[[`, "rep")), 1,
                                median, na.rm = TRUE)
    } else {
      mean_consensus <- apply(sapply(p_occ, function(x) x$Full_model), 1,
                              mean, na.rm = TRUE)
      median_consensus <- apply(sapply(p_occ, function(x) x$Full_model), 1,
                                median, na.rm = TRUE)
    }
  }

  p_occ <- lapply(p_occ, function(x) x[names(x) != "rep"])

  consensus <- list(mean = mean_consensus, median = median_consensus)

  p_occ <- c(p_occ, list(consensus = consensus))

  # Calculate thresholds
  p_thr <- lapply(p_occ, function(model) {
    lapply(model, calc_thr,
           thr = calibration_results$summary$omission_rate_thr / 100)
  })

  # Append type of predictions
  p_thr$type <- type

  #Prepare final data
  if(replicate_method %in% c("kfolds", "leave-one-out")){
    sample_proportion <- NULL
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
    replicate_method = replicate_method,
    n_replicates = n_replicates,
    sample_proportion = sample_proportion
  )

  # Optionally save the fitted models
  if (write_models) {
    saveRDS(res, file = paste0(file_name, ".RDS"))
  }

  return(res)
}
