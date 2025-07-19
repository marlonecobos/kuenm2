#' Partial ROC calculation for multiple candidate models
#'
#' @description
#' Computes partial ROC tests for multiple candidate models.
#'
#' @usage
#' partial_roc(formula_grid, data, omission_rate = 10,
#'             addsamplestobackground = TRUE, weights = NULL,
#'             algorithm = "maxnet", parallel = FALSE, ncores = NULL,
#'             progress_bar = TRUE)
#'
#' @param formula_grid a data.frame with the grid of formulas defining the
#' candidate models to test.
#' @param data an object of class `prepared_data` returned by the
#' [prepare_data()] function or an object of class
#' calibration_results returned by the [calibration()] function.
#' It contains the calibration data and k-folds.
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
#' Default is NULL and uses available cores - 1. This is only applicable if
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
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach `%dopar%`
#' @importFrom enmpa predict_glm
#' @importFrom fpROC auc_metrics
#'
#' @export
#'
#' @details
#' Partial ROC is calculated following Peterson et al. (2008)
#' <doi:10.1016/j.ecolmodel.2007.11.008>.
#'
#' @examples
#' # Import prepared data to get model formulas
#' data(sp_swd, package = "kuenm2")
#'
#' # Calculate proc for the first 5 candidate models
#' res_proc <- partial_roc(formula_grid = sp_swd$formula_grid[1:2,],
#'                         data = sp_swd, omission_rate = 10,
#'                         algorithm = "maxnet")

partial_roc <- function(formula_grid, data, omission_rate = 10,
                        addsamplestobackground = TRUE, weights = NULL,
                        algorithm = "maxnet", parallel = FALSE, ncores = NULL,
                        progress_bar = TRUE) {

  # error check
  if (missing(formula_grid)) {
    stop("Argument 'formula_grid' must be defined.")
  }
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (!inherits(formula_grid, "data.frame")) {
    stop("Argument 'formula_grid' must be a 'data.frame' object.")
  }
  if (!inherits(data, "prepared_data")) {
    stop("Argument 'data' must be a 'prepared_data' object.")
  }

  n_tot <- nrow(formula_grid)

  if (n_tot == 1 & parallel) {
    parallel <- FALSE
  }

  if (parallel) {
    if (is.null(ncores)) {
      ncores <- max(1, parallel::detectCores() - 1)
    }
    if (ncores > n_tot) {
      ncores <- n_tot
    }

    cl <- parallel::makeCluster(ncores)
  }

  #Progress bar setup

  if (progress_bar) {
    pb <- utils::txtProgressBar(min = 0, max = n_tot, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else {
    opts <- NULL
  }

  # Execute in parallel or sequentially
  if (parallel) {
    doSNOW::registerDoSNOW(cl)

    res_proc <- foreach::foreach(
      x = 1:n_tot, .packages = c("glmnet", "enmpa", "fpROC"),
      .options.snow = opts) %dopar% {
        proc(x, formula_grid, data, omission_rate, addsamplestobackground,
             weights, algorithm)
      }
  } else {
    res_proc <- vector("list", length = n_tot)
    for (x in 1:n_tot) {
      res_proc[[x]] <- proc(x, formula_grid, data, omission_rate,
                            addsamplestobackground, weights,
                            algorithm)

      # Sets the progress bar to the current state
      if (progress_bar) {
        setTxtProgressBar(pb, x)
      }
    }
  }
  #Bind results
  return(do.call("rbind", res_proc))
}
