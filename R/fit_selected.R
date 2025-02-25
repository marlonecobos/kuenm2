#' Fit selected models
#'
#' @description
#' This function fits and calculates the consensus (mean and median) of models
#' that were evaluated and selected using the \code{\link{calibration}}()
#' function. It also calculates the thresholds based on the omission rate set in
#'  \code{\link{calibration}}().The function supports parallelization for faster
#'  model fitting.
#'
#' @param calibration_results an object of class `calibration_results` returned
#' by the \code{\link{calibration}}() function.
#' @param n_replicates (numeric) the number of model replicates. Default is 5.
#' @param rep_type (character) the replicate type. It can be: "kfold", "bootstrap",
#' or "subsample". Default is "kfold".
#' @param train_portion (numeric) the proportion of occurrence records used to
#' train the model in each replicate. This parameter is applicable only when
#' `rep_type` is set to "bootstrap" or "subsample". Default is 0.7.
#' @param write_models (logical) whether to save the final fitted models to disk.
#' Default is FALSE.
#' @param file_name (character) the file name, with or without a path, for saving
#' the final models. This is only applicable if `write_models = TRUE`.
#' @param parallel (logical) whether to fit the final models in parallel. Default
#' is FALSE.
#' @param ncores (numeric) the number of cores to use for parallel processing.
#' Default is 2. This is only applicable if `parallel = TRUE`.
#' @param parallel_option (character) the parallelization package to use:
#' either "doParallel" or "doSNOW". Default is "doSNOW". This is only applicable
#' if `parallel = TRUE`.
#' @param progress_bar (logical) whether to display a progress bar during processing.
#' Default is TRUE.
#' @param verbose (logical) whether to display detailed messages during processing.
#' Default is TRUE.
#' @param seed (numeric) integer value used to specify an initial seed to split
#' the data. Default is 42.
#'
#' @return
#' An object of class 'fitted_models' containing the following elements:
#' \item{species}{a character string with the name of the species.}
#' \item{Models}{a list of fitted models, including replicates (trained with
#' the training data) and full models (trained with all available records).}
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
#' \code{\link{calibration}}()}
#'
#' @export
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom foreach foreach `%dopar%`
#' @importFrom glmnet glmnet.control glmnet
#' @importFrom enmpa predict_glm
#'
#' @usage fit_selected(calibration_results, n_replicates = 5, rep_type = "kfold",
#'                     train_portion = 0.7, write_models = FALSE, file_name = NULL,
#'                     parallel = FALSE, ncores = 2, parallel_option = "doSNOW",
#'                     progress_bar = TRUE, verbose = TRUE, seed = 42)
#'
#' @examples
#' # Import example of calibration results (output of calibration function)
#' ## maxnet
#' data("calib_results_glmnet", package = "kuenm2")
#'
#' # Fit models using calibration results
#' fm <- fit_selected(calibration_results = calib_results_glmnet,
#'                    n_replicates = 2,
#'                    rep_type = "kfold",
#'                    train_portion = 0.7,
#'                    write_models = FALSE,
#'                    file_name = NULL,
#'                    parallel = FALSE,
#'                    ncores = 1,
#'                    parallel_option = "doSNOW",
#'                    progress_bar = TRUE,
#'                    verbose = TRUE,
#'                    seed = 42)
#'
#' # Output the fitted models
#' print(fm)
#'
#' ##GLM
#' data("calib_results_glm", package = "kuenm2")
#'
#' # Fit models using calibration results
#' fm_glm <- fit_selected(calibration_results = calib_results_glm,
#'                        n_replicates = 2,
#'                        rep_type = "kfold",
#'                        train_portion = 0.7,
#'                        write_models = FALSE,
#'                        file_name = NULL,
#'                        parallel = FALSE,
#'                        ncores = 1,
#'                        parallel_option = "doSNOW",
#'                        progress_bar = TRUE,
#'                        verbose = TRUE,
#'                        seed = 42)
#'
#' # Output the fitted models
#' print(fm_glm)

fit_selected <- function(calibration_results, n_replicates = 5,
                         rep_type = "kfold", train_portion = 0.7,
                         write_models = FALSE, file_name = NULL,
                         parallel = FALSE, ncores = 2, parallel_option = "doSNOW",
                         progress_bar = TRUE, verbose = TRUE, seed = 42) {

  # Extract model IDs from selected models
  m_ids <- calibration_results$selected_models$ID
  algorithm <- calibration_results$algorithm

  # Fitting models over multiple replicates_____________________________________
  if(n_replicates > 1){
    if(verbose){
      message("Fitting replicates...")
    }

    # Create a grid of model IDs and replicates
    dfgrid <- expand.grid(models = m_ids, replicates = 1:n_replicates)
    n_tot <- nrow(dfgrid) # Total models * replicates

    #Prepare data (index) to replicates
    if(n_replicates > 1) {
      #Partitioning data
      rep_data <- part_data(data = calibration_results$calibration_data,
                            pr_bg = "pr_bg",
                            train_portion = train_portion,
                            n_replicates = n_replicates,
                            method = rep_type, seed = seed)
    } else {
      rep_data <- NULL
    }

    # Progress bar setup
    if (isTRUE(progress_bar)) {
      pb <- utils::txtProgressBar(0, n_tot, style = 3)
      progress <- function(n) utils::setTxtProgressBar(pb, n)
    }

    # Adjust parallelization based on number of tasks and cores
    if (n_tot == 1 & parallel) {
      parallel <- FALSE
    }
    if (n_tot < ncores & parallel) {
      ncores <- n_tot
    }

    # Setup parallel cluster
    if (parallel) {
      cl <- parallel::makeCluster(ncores)
      if (parallel_option == "doParallel") {
        doParallel::registerDoParallel(cl)
        opts <- NULL
      } else if (parallel_option == "doSNOW") {
        doSNOW::registerDoSNOW(cl)
        opts <- if (progress_bar) list(progress = progress) else NULL
      }
    } else {
      opts <- NULL
    }

    # Fit the best models (either in parallel or sequentially)
    if (parallel) {
      best_models <- foreach::foreach(x = 1:n_tot,
                                      .packages = c("glmnet", "enmpa"),
                                      .options.snow = opts
      ) %dopar% {
        fit_best_model(x, dfgrid, calibration_results, n_replicates,
                       rep_data, algorithm)
      }
    } else {
      best_models <- vector("list", length = n_tot)
      for (x in 1:n_tot) {
        best_models[[x]] <- fit_best_model(x, dfgrid, calibration_results,
                                           n_replicates, rep_data, algorithm)
        if (progress_bar) utils::setTxtProgressBar(pb, x)
      }
    }

    # Stop the cluster
    if (parallel) parallel::stopCluster(cl)

    # Split models by their respective IDs
    best_models <- split(best_models, dfgrid$models)
    best_models <- lapply(best_models, function(sublist) {
      names(sublist) <- paste0("Rep_", seq_along(sublist))
      return(sublist)
    })

    # Assign names to the models based on their IDs
    names(best_models) <- paste0("Model_", names(best_models))

  } else {
    best_models <- list()
  }

  # Fit full models ____________________________________________________________
  if(verbose){
    message("\nFitting full models...")
  }
  # Full models grid setup
  n_models <- length(m_ids)
  dfgrid <- expand.grid(models = m_ids, replicates = 1)

  # Adjust parallelization for full models
  if (n_models == 1 & parallel) {
    parallel <- FALSE
  }
  if (n_models < ncores & parallel) {
    ncores <- n_models
  }

  # Progress bar setup for full models
  if (progress_bar) {
    pb <- utils::txtProgressBar(0, n_models, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
  }

  # Parallel cluster setup for full models
  if (parallel) {
    cl <- parallel::makeCluster(ncores)
    if (parallel_option == "doParallel") {
      doParallel::registerDoParallel(cl)
      opts <- NULL
    } else if (parallel_option == "doSNOW") {
      doSNOW::registerDoSNOW(cl)
      opts <- if (progress_bar) list(progress = progress) else NULL
    }
  } else {
    opts <- NULL
  }

  # Fit the full models (either in parallel or sequentially)
  if (parallel) {
    full_models <- foreach::foreach(x = 1:n_models,
                                    .packages = c("glmnet", "enmpa"),
                                    .options.snow = opts
    ) %dopar% {
      fit_best_model(x, dfgrid, calibration_results, n_replicates = 1, rep_data,
                     algorithm)
    }
  } else {
    full_models <- vector("list", length = n_models)
    for (x in 1:n_models) {
      full_models[[x]] <- fit_best_model(x, dfgrid, calibration_results,
                                         n_replicates = 1, rep_data, algorithm)
      if (progress_bar) utils::setTxtProgressBar(pb, x)
    }
  }

  # Assign names to full models
  names(full_models) <- paste0("Model_", m_ids)

  # Combine replicate models with full models
  for (i in names(full_models)) {
    best_models[[i]]$Full_model <- full_models[[i]]
  }

  # Stop the cluster for full models
  if (parallel) parallel::stopCluster(cl)

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
                                                       type = "cloglog"))
    } else if (algorithm == "glm") {
      p_r <- sapply(m_x, function(i) suppressWarnings(
        enmpa::predict_glm(model = i,
                           data = calibration_results$calibration_data,
                           newdata = occ, type = "response"))
        )
    }

    p_mean <- apply(p_r, 1, mean, na.rm = TRUE)
    p_median <- apply(p_r, 1, median, na.rm = TRUE)

    list(mean = p_mean, median = p_median)
  })

  names(p_occ) <- names(best_models)

  # Calculate consensus across models
  mean_consensus <- apply(sapply(p_occ, function(x) x$mean), 1, mean, na.rm = TRUE)
  median_consensus <- apply(sapply(p_occ, function(x) x$median), 1, median, na.rm = TRUE)
  consensus <- list(mean = mean_consensus, median = median_consensus)
  p_occ <- c(p_occ, list(consensus = consensus))

  # Calculate thresholds
  p_thr <- lapply(p_occ, function(model) {
    lapply(model, calc_thr,
           thr = calibration_results$summary$omission_rate_thr / 100)
  })

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
    algorithm = algorithm
  )

  # Optionally save the fitted models
  if (write_models) {
    saveRDS(res, file = paste0(file_name, ".RDS"))
  }

  return(res)
}
