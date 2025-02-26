#' Partial ROC calculation for candidate models
#'
#' @description
#' Applies partial ROC tests to model predictions
#'
#' @usage partial_roc(formula_grid, data, omission_rate = 10,
#'                    addsamplestobackground = TRUE, weights = NULL,
#'                    algorithm = maxnet, parallel = FALSE, ncores = 1,
#'                    parallel_type = "doSNOW", progress_bar = TRUE)
#' @param formula_grid (data.frame) the formula grid of the candidate models
#' to test.
#' @param data an object of class `prepared_data` returned by the
#' \code{\link{prepare_data()}} function or an object of class
#' calibration_results returned by the \code{\link{calibration()}} function.
#' It contains the calibration data and kfolds.
#' @param omission_rate (numeric) values from 0 to 100 representing the
#' percentage of potential error due to any source of uncertainty. This value is
#' used to calculate the omission rate. Default is 10. See details.
#' @param addsamplestobackground (logical) whether to add to the background any
#' presence sample that is not already there. Default is TRUE.
#' @param weights (numeric) a numeric vector specifying weights for the
#' occurrence records. Default is NULL.
#' @param algorithm (character) type algorithm, either "glm" or "maxnet".
#' Default is "maxnet".
#' @param parallel (logical) whether to fit the candidate models in parallel.
#' Default is FALSE.
#' @param ncores (numeric) number of cores to use for parallel processing.
#' Default is 1. This is only applicable if `parallel = TRUE`.
#' @param parallel_type (character) the package to use for parallel processing:
#' "doParallel" or "doSNOW". Default is "doSNOW". This is only applicable if
#' `parallel = TRUE`.
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is TRUE.
#'
#' @return
#' A data frame with summary statistics of the and AUC ratios and significance
#' calculated from the replicates of each candidate model. Specifically, it
#' includes the mean and standard deviation of these metrics for each model.
#'
#' @importFrom parallel makeCluster
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach `%dopar%`
#' @importFrom enmpa predict_glm proc_enm
#'
#' @details
#' Partial ROC is calculated following Peterson et al. (2008) doi:10.1016/j.ecolmodel.2007.11.008.
#'
#' @examples
#' # Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' # Import occurrences
#' data(occ_data, package = "kuenm2")
#'
#'
#' # Prepare data for maxnet model
#' sp_swd <- prepare_data(algorithm = "maxnet", occ = occ_data,
#'                        species = occ_data[1, 1], x = "x", y = "y",
#'                        raster_variables = var,
#'                        categorical_variables = "SoilType",
#'                        n_background = 500,
#'                        features = c("l", "q", "p", "lq", "lqp"),
#'                        reg_mult = c(0.1, 1, 2, 3, 5))
#'
#' # Calculate proc for the first 5 candidate models
#' res_proc <- partial_roc(formula_grid = sp_swd$formula_grid[1:5,],
#'                         data = sp_swd,
#'                         omission_rate = 10)
#'
partial_roc <- function(formula_grid, data, omission_rate = 10,
                        addsamplestobackground = TRUE, weights = NULL,
                        algorithm = "maxnet", parallel = FALSE, ncores = 1,
                        parallel_type = "doSNOW", progress_bar = TRUE) {

  # Parallelization setup
  n_tot <- nrow(formula_grid)

  if (parallel) {
    if (!(parallel_type %in% c("doSNOW", "doParallel"))) {
      stop("Invalid parallel_type. Use 'doSNOW' or 'doParallel'.")
    }

    if(ncores > n_tot){
      ncores <- n_tot
    }

    cl <- parallel::makeCluster(ncores)
  }

  #Progress bar setup

  if (progress_bar) {
    pb <- utils::txtProgressBar(min = 0, max = n_tot, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)} else {opts <- NULL}

  if (parallel & parallel_type == "doParallel") {
    doParallel::registerDoParallel(cl)
    opts <- NULL # Progress bar does not work with doParallel
  }

  if (parallel & parallel_type == "doSNOW") {
    doSNOW::registerDoSNOW(cl)
    if (isTRUE(progress_bar))
      opts <- list(progress = progress)
    else opts <- opts
  }

  # Execute in parallel or sequentially
  if(parallel){
    res_proc <- foreach::foreach(
      x = 1:n_tot,
      .packages = c("glmnet", "enmpa"),
      .options.snow = opts) %dopar% {
        proc(x, formula_grid, data, omission_rate,
             addsamplestobackground, weights,
             algorithm)
      }
  } else {
    res_proc <- vector("list", length = n_tot)
    for (x in 1:n_tot) {
      res_proc[[x]] <-
        proc(x, formula_grid, data, omission_rate,
             addsamplestobackground, weights,
             algorithm)
      # Sets the progress bar to the current state
      if(progress_bar) setTxtProgressBar(pb, x)
    }
  }
  #Bind results
  return(do.call("rbind", res_proc))
}
