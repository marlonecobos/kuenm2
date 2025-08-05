#' Response curves for selected models according to training/testing partitions
#'
#' @description
#' Variable responses in models selected after model calibration. Responses
#' are based on training partitions and points are testing presence records.
#'
#' @usage
#' partition_response_curves(calibration_results, modelID, n = 100,
#'                           averages_from = "pr_bg", col = "darkblue",
#'                           ylim = NULL, las = 1, parallel = FALSE,
#'                           ncores = NULL, ...)
#'
#' @param calibration_results an object of class `calibration_results` returned
#' by the [calibration()] function.
#' @param modelID (character or numeric) number of the Model (its ID) to be
#' considered for plotting.
#' @param n (numeric) an integer guiding the number of breaks. Default = 100
#' @param averages_from (character) specifies how the averages or modes of the
#' variables are calculated. Available options are "pr" (to calculate averages
#' from the presence localities) or "pr_bg" (to use the combined set of presence
#' and background localities). Default is "pr_bg". See details.
#' @param col (character) color for lines. Default = "darkblue".
#' @param ylim (numeric) vector of length two indicating minimum and maximum
#' limits for the y axis. The default, NULL, uses `c(0, 1)`.
#' @param las (numeric) the stile of axis tick labels; options are: 0, 1, 2, 3.
#' Default = 1.
#' @param parallel (logical) whether to fit the models in parallel. Default
#' is FALSE.
#' @param ncores (numeric) number of cores to use for parallel processing.
#' Default is NULL and uses available cores - 1. This is only applicable if
#' `parallel = TRUE`.
#' @param ... additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @details
#' Response curves are generated using training portions of the data and points
#' showed are the ones left out for testing. The partition labeled in plot
#' panels indicates the portion left out for testing.
#'
#' The response curves are generated with all other variables set to their mean
#' values (or mode for categorical variables), calculated either from the
#' presence localities (if averages_from = "pr") or from the combined set of
#' presence and background localities (if averages_from = "pr_bg").
#'
#' For categorical variables, a bar plot is generated with error bars showing
#' variability across models (if multiple models are included).
#'
#' @return
#' A plot with response curves for all `variables` used in the selected model
#' corresponding to `modelID`. Each row in the plot shows response curves
#' produced with training data that leaves out the partition laveled. The points
#' represent the records left out for testing.
#'
#' @export
#'
#' @seealso
#' [response_curve()]
#'
#' @examples
#' # Example with maxnet
#' # Import example of calibration results
#' data(calib_results_maxnet, package = "kuenm2")
#'
#' # Options of models that can be tested
#' calib_results_maxnet$selected_models$ID
#'
#' # Response curves
#' partition_response_curves(calibration_results = calib_results_maxnet,
#'                           modelID = 192)

partition_response_curves <- function(calibration_results,
                                      modelID,
                                      n = 100,
                                      averages_from = "pr_bg",
                                      col = "darkblue",
                                      ylim = NULL,
                                      las = 1,
                                      parallel = FALSE,
                                      ncores = NULL,
                                      ...) {

  if (missing(calibration_results)) {
    stop("Argument 'calibration_results' must be defined.")
  }
  if (missing(modelID)) {
    stop("Argument 'modelID' must be defined.")
  }

  if (!inherits(calibration_results, "calibration_results")) {
    stop("Argument 'calibration_results' must be a 'calibration_results' object.")
  }


  # base arguments
  m_id <- as.numeric(gsub("Model_", "", modelID))
  algorithm <- calibration_results$algorithm

  # Create a grid of model IDs and partitions
  tparts <- length(calibration_results$part_data)

  dfgrid <- expand.grid(models = m_id, replicates = 1:tparts)
  n_tot <- nrow(dfgrid) # Total models * partitions

  if (n_tot == 1 & parallel) {
    parallel <- FALSE
  }

  # set up parallel processing
  if (parallel) {
    if (is.null(ncores)) {
      ncores <- max(1, parallel::detectCores() - 1)
    }
    if (n_tot < ncores) {
      ncores <- n_tot
    }

    cl <- parallel::makeCluster(ncores)
    doSNOW::registerDoSNOW(cl)
  }

  # fit selected models
  if (parallel) {
    best_models <- foreach::foreach(x = 1:n_tot,
                                    .packages = c("glmnet", "enmpa")
    ) %dopar% {
      fit_best_model(x = x, dfgrid = dfgrid, cal_res = calibration_results,
                     n_replicates = tparts,
                     rep_data = calibration_results$part_data,
                     algorithm = algorithm)
    }
  } else {
    best_models <- vector("list", length = n_tot)
    for (x in 1:n_tot) {
      best_models[[x]] <- fit_best_model(
        x = x, dfgrid = dfgrid, cal_res = calibration_results,
        n_replicates = tparts, rep_data = calibration_results$part_data,
        algorithm = algorithm
      )
    }
  }

  # Stop the cluster
  if (parallel) parallel::stopCluster(cl)

  # get all variables to plot
  ## patterns to remove
  pattern <- "I\\(|\\^2\\)|categorical\\(|\\)|thresholds\\(|hinge\\("

  ## get all variables
  coefss <- lapply(best_models, function(x) {
    sort(unique(unlist(strsplit(
      gsub(pattern, "",
           unname(
             coefs <- if (inherits(x, "glmnet")) {
               names(x$betas)
             } else if (inherits(x, "glm")) {
               names(coef(x)[-1])
             }
           )),
      ":"
    ))))
  })

  # plot response curves
  ## par settings
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))

  ## arrangement
  allvarun <- sort(unique(unlist(coefss)))
  mfrow <- c(n_tot, length(allvarun))
  graphics::par(mfrow = mfrow, mar = c(2, 4, 1.5, 0.5))

  ## modify coefss to include NA when a variable is not present
  acoefss <- lapply(coefss, function(x) {
    unname(sapply(allvarun, function(v) {if (v %in% x) v else NA}))
  })

  ## other plotting parameters
  if (is.null(ylim)) {
    ylim <- c(0, 1)
  }
  p_col <- adjustcolor("black", alpha.f = 0.5)

  ## plot curves
  for (h in 1:length(coefss)) {
    datat <- calibration_results$calibration_data[calibration_results$part_data[[h]], ]

    ## new limits for plot
    newdata <- apply(calibration_results$calibration_data[, coefss[[h]], drop = FALSE],
                     2, range)
    extension <- 0.1 * apply(newdata, 2, diff)
    newdata[1, ] <- newdata[1, ] - extension
    newdata[2, ] <- newdata[2, ] + extension

    if (h == 1) {
      main <- paste("Suitability:" , allvarun)
    } else {
      main <- rep(NULL, length(allvarun))
    }

    for (i in 1:length(acoefss[[h]])) {
      if (i == 1) {
        ylab <- names(calibration_results$part_data)[h]
      } else {
        ylab <- ""
      }

      if (!is.na(acoefss[[h]][i])) {
        response_curve_consmx(model_list = list(best_models[[h]]),
                              variable = acoefss[[h]][i],
                              data = calibration_results$calibration_data[-calibration_results$part_data[[h]], ],
                              show_lines = FALSE, n = n, new_data = newdata,
                              xlab = "", ylab = ylab, col = col,
                              categorical_variables = calibration_results$categorical_variables,
                              averages_from = averages_from, ylim = ylim,
                              main = main[i], cex.main = 1, font.main = 1,
                              las = las, ...)

        ## add testing points
        points(datat[datat$pr_bg == 1, c(acoefss[[h]][i], "pr_bg")],
               bg = p_col, pch = 21, cex = 0.6)
      } else {
        plot.new()
        title(main = main[i], ylab = ylab, cex.main = 1, font.main = 1)
      }
    }
  }
}
