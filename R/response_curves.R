#' Variable response curves
#'
#' @description
#' A view of variable responses in models. Responses based on single or multiple
#' models can be provided.
#'
#' @usage
#' response_curve(fitted, variable, modelID = NULL, n = 100,
#'                by_replicates = FALSE, data = NULL, new_data = NULL,
#'                extrapolate = TRUE, xlab = NULL, ylab = "Suitability",
#'                col = "darkblue" ...)
#
#' @param fitted object.
#' @param variable (character) name of the variables to be plotted.
#' @param data data.frame or matrix of data used in the model calibration step.
#' Default = NULL.
#' @param modelID (character) vector of ModelID(s) to be considered in the
#' fitted object. By default all models are included.Default = NULL.
#' @param n (numeric) an integer guiding the number of breaks. Default = 100
#' @param by_replicates (logical) whether use replicates or full_model to
#' estimate the model's response curve.
#' @param new_data a `SpatRaster`, data.frame, or  matrix of variables
#' representing the range of variable values in an area of interest.
#' Default = NULL. It must be defined in case the model entered does not
#' explicitly include a data component.
#' @param extrapolate (logical) whether to allow extrapolation to study the
#' behavior of the response outside the calibration limits. Ignored if
#' `new_data` is defined. Default = TRUE.
#' @param xlab (character) a label for the x axis. The default, NULL, uses the
#' name defined in `variable`.
#' @param ylab (character) a label for the y axis. Default = "Suitability".
#' @param col (character) color for lines. Default = "darkblue".
#' @param ... additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return
#' A plot with the response curve for a `variable`.
#'
#' @export
#'
#' @importFrom stats predict coef
#' @importFrom graphics abline polygon
#' @importFrom terra minmax
#' @importFrom mgcv gam

response_curve <- function(fitted, variable, modelID = NULL, n = 100,
                           by_replicates = FALSE, data = NULL,
                           new_data = NULL, extrapolate = TRUE,
                           xlab = NULL, ylab = "Suitability",
                           col = "darkblue", ...) {

  # initial tests
  if (missing(fitted) | missing(variable) ) {
    stop("Arguments 'fitted', and 'variable' must be defined.")
  }

  if (!is.null(new_data)) {
    if (!class(new_data)[1] %in% c("matrix", "data.frame", "SpatRaster")) {
      stop("'new_data' must be of class 'matrix', 'data.frame', 'SpatRaster'")
    }
  }

  # if data is not defined it is extratec from the fitted kuenm2 object
  if (is.null(data)){
    data <- fitted$calibration_data
  }

  ## add a warming message indicating that the are not replicates.


  if (!is.null(modelID)){

    if (!modelID %in% names(fitted[["Models"]])){
      stop(paste0(
        "The 'ModelID' is not correct, check the following: [",
        paste(names(fitted[["Models"]]), collapse = ", ")),
        "]"
           )
    }

    if (by_replicates){
      glmnet_list <- fitted[["Models"]][[modelID]]
      glmnet_list$Full_model <- NULL


      coefs <- names(glmnet_list[[1]]$betas)
      c1 <- any(c(variable, paste0("I(", variable, "^2)")) %in% coefs)
      c2 <- any(grepl(paste0(variable, ":"), coefs))
      c3 <- any(grepl(paste0(":", variable,"$"), coefs))

      if (any(c1, c2, c3) == FALSE){
        stop("Defined 'variable' is not present in the fitted model.")
      }

    } else {
      glmnet_list <- fitted[["Models"]][[modelID]]["Full_model"]
    }

  } else {

    # extract the slot Full_model for each model
    glmnet_list <- lapply(fitted[["Models"]], function(x){x$Full_model})
  }

  # Response curve for all selected models
  response_curve_consmx(glmnet_list, data = data, variable = variable, n = n,
                        new_data = new_data, extrapolate = extrapolate,
                        xlab = xlab, ylab = ylab,
                        col = col, ...)
  }

