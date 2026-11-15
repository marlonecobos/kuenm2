#' Partial ROC calculation for multiple candidate models
#'
#' @description
#' Computes partial ROC tests for multiple candidate models.
#'
#' @usage
#' partial_roc(formula_grid, data, omission_rate = 10,
#'             addsamplestobackground = TRUE, weights = NULL,
#'             algorithm = "maxnet", parallel = FALSE, ncores = NULL,
#'             progress_bar = TRUE, progress_type = "timer",
#'             use_lb = FALSE)
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
#' @param progress_type (character) type of the progress bar. See
#' `pbapply::pbtypes()` for available progress bar types depending on operating
#' system. Default is "timer".
#' @param use_lb (logical) whether to use load balancing when running in
#' parallel. Default is FALSE.
#'
#' @return
#' A data frame with summary statistics of the and AUC ratios and significance
#' calculated from the replicates of each candidate model. Specifically, it
#' includes the mean and standard deviation of these metrics for each model.
#'
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom enmpa predict_glm
#' @importFrom fpROC auc_metrics
#' @importFrom pbapply pboptions pblapply
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
                        progress_bar = TRUE, progress_type = "timer",
                        use_lb = FALSE) {

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

  # Run
  res_proc <- pbapply::pblapply(
    X = seq_len(n_tot),
    FUN = proc,
    formula_grid = formula_grid,
    data = data,
    error_considered = omission_rate,
    addsamplestobackground = addsamplestobackground,
    weights = weights,
    algorithm = algorithm,
    cl = cl)

  #Bind results
  return(do.call("rbind", res_proc))
}
