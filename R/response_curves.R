#' Variable response curves
#'
#' @description
#' A view of variable responses in fitted models. Responses based on single or multiple
#' models can be provided.
#'
#' @usage response_curve(models, variable, modelID = NULL, n = 100,
#'                      by_replicates = FALSE, data = NULL, new_data = NULL,
#'                      extrapolate = TRUE, xlab = NULL, ylab = "Suitability",
#'                      col = "darkblue", ...)
#
#' @param models an object of class `fitted_models` returned by the
#' \code{\link{fit_selected}}() function.
#' @param variable (character) name of the variables to be plotted.
#' @param data data.frame or matrix of data used in the model calibration step.
#' Default = NULL.
#' @param modelID (character) vector of ModelID(s) to be considered in the
#' models object. By default all models are included.Default = NULL.
#' @param n (numeric) an integer guiding the number of breaks. Default = 100
#' @param by_replicates (logical) whether use replicates or full_model to
#' estimate the model's response curve. Default = FALSE.
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
#'
#' @examples
#' ##Example with glmnet
#' # Import example of fitted_models (output of fit_selected())
#' data("fitted_model_glmnet", package = "kuenm2")
#'
#' #Response curves
#' response_curve(models = fitted_model_glmnet,
#'                variable = "bio_1", by_replicates = T)
#' response_curve(models = fitted_model_glmnet,
#'                variable = "bio_7", by_replicates = T)
#' response_curve(models = fitted_model_glmnet, variable = "bio_1",
#'                modelID = "Model_13", by_replicates = T)
#' response_curve(models = fitted_model_glmnet, variable = "bio_1",
#'                modelID = "Model_1", by_replicates = T)
#'
#' ##Example with glm
#' # Import example of fitted_models (output of fit_selected())
#' data("fitted_model_glm", package = "kuenm2")
#'
#' #Response curves
#' response_curve(models = fitted_model_glm,
#'                variable = "bio_1", by_replicates = T)
#' response_curve(models = fitted_model_glm,
#'                variable = "bio_7", by_replicates = T)
#' response_curve(models = fitted_model_glm, variable = "bio_1",
#'                modelID = "Model_4", by_replicates = T)
#' response_curve(models = fitted_model_glm, variable = "bio_12",
#'                modelID = "Model_4", by_replicates = T)
#'
response_curve <- function(models, variable, modelID = NULL, n = 100,
                           by_replicates = FALSE, data = NULL,
                           new_data = NULL, extrapolate = TRUE,
                           xlab = NULL, ylab = "Suitability",
                           col = "darkblue", ...) {

  # initial tests
  if (missing(models) | missing(variable) ) {
    stop("Arguments 'models' and 'variable' must be defined.")
  }

  if (!inherits(models, "fitted_models")) {
    stop(paste0("Argument models must be a fitted_models object, not ",
                class(models)))
  }

  if (!inherits(variable, "character")) {
    stop(paste0("Argument variable must be a character, not ",
                class(variable)))
  }

  if (!is.null(modelID) & !inherits(modelID, "character")) {
    stop(paste0("Argument modelID must be a character, not ",
                class(modelID)))
  }

  if (!inherits(n, "numeric")) {
    stop(paste0("Argument n must be numeric, not ",
                class(n)))
  }

  if (!inherits(by_replicates, "logical")) {
    stop(paste0("Argument by_replicates must be logical, not ",
                class(by_replicates)))
  }

  if (!is.null(new_data)) {
    if (!class(new_data)[1] %in% c("matrix", "data.frame", "SpatRaster")) {
      stop("'new_data' must be of class 'matrix', 'data.frame', 'SpatRaster'")
    }
  }

  if (!inherits(extrapolate, "logical")) {
    stop(paste0("Argument extrapolate must be logical, not ",
                class(extrapolate)))
  }

  if (!is.null(xlab) & !inherits(xlab, "character")) {
    stop(paste0("Argument xlab must be NULL or a character, not ",
                class(xlab)))
  }

  if (!is.null(ylab) & !inherits(ylab, "character")) {
    stop(paste0("Argument ylab must be NULL or a character, not ",
                class(ylab)))
  }

  if (!inherits(col, "character")) {
    stop(paste0("Argument col must be a character, not ",
                class(col)))
  }

  if (!is.null(models$categorical_variables)) {
    if (variable %in% models$categorical_variables) {
      stop(paste0("Response curves for categorical variables are not yet",
                  " implemented. We're actively working on it."))
    }
  }

  # if data is not defined it is extratec from the models kuenm2 object
  if (is.null(data)){
    data <- models$calibration_data
  }

  ## add a warming message indicating that the are not replicates.
  if (!is.null(modelID)){

    if (!modelID %in% names(models[["Models"]])){
      stop(paste0(
        "The 'ModelID' is not correct, check the following: [",
        paste(names(models[["Models"]]), collapse = ", ")),
        "]"
      )
    }

    # Handling replicates or the full model
    if (by_replicates){
      model_list <- models[["Models"]][[modelID]]
      model_list$Full_model <- NULL

      # Check if the variable is present in any of the replicates
      coefs <- if (inherits(model_list[[1]], "glmnet")) {
        names(model_list[[1]]$betas)
      } else if (inherits(model_list[[1]], "glm")) {
        names(coef(model_list[[1]])[-1])
      }

      c1 <- any(c(variable, paste0("I(", variable, "^2)")) %in% coefs)
      c2 <- any(grepl(paste0(variable, ":"), coefs))
      c3 <- any(grepl(paste0(":", variable,"$"), coefs))

      if (any(c1, c2, c3) == FALSE){
        stop("Defined 'variable' is not present in the models model.")
      }
    } else {
      model_list <- models[["Models"]][[modelID]]["Full_model"]
    }
  } else {
    # extract the slot Full_model for each model
    model_list <- lapply(models[["Models"]], function(x){x$Full_model})
  }

  # Response curve for all selected models
  response_curve_consmx(model_list, data = data, variable = variable, n = n,
                        new_data = new_data, extrapolate = extrapolate,
                        xlab = xlab, ylab = ylab,
                        col = col,
                        categorical_variables = models$categorical_variables,
                        ...)
}

### Helpers functions for response curves

# Consensus response curve
response_curve_consmx <- function(model_list, variable, data, n = 100,
                                  extrapolate = FALSE, new_data = NULL,
                                  xlab = NULL, ylab = NULL, col = "darkblue",
                                  categorical_variables = NULL,
                                  ...) {

  # initial tests
  if (missing(model_list) | missing(variable) | missing(data)) {
    stop("Arguments 'model_list', 'data', and 'variable' must be defined.")
  }

  if (!is.null(new_data)) {
    if (!class(new_data)[1] %in% c("matrix", "data.frame", "SpatRaster")) {
      stop("'new_data' must be of class 'matrix', 'data.frame', 'SpatRaster'")
    }
  }

  if (!variable %in% colnames(data)) {
    stop("The name of the 'variable' was not defined correctly.")
  }


  if (length(model_list) == 1 ){
    model <- model_list[[1]]

    # Handling glmnet and glm models differently
    coefs <- if (inherits(model, "glmnet")) {
      names(model$betas)
    } else if (inherits(model, "glm")) {
      names(coef(model)[-1])
    }

    c1 <- any(c(variable, paste0("I(", variable, "^2)")) %in% coefs) # check linear or quadratic term
    c2 <- any(grepl(paste0(variable, ":"), coefs))                   # check product terms 1st position
    c3 <- any(grepl(paste0(":", variable,"$"), coefs))               # check product terms 2nd position

    if (any(c1, c2, c3) == FALSE){
      stop("Defined 'variable' is not present in the models model.")
    }

    response_out <- response(model = model_list[[1]], data = data,
                             variable = variable, n = n,
                             extrapolate = extrapolate, new_data = new_data,
                             categorical_variables = categorical_variables)

    limits <- range(data[,variable])

    ## Plotting curve
    if (is.null(xlab)) {xlab <- variable}
    if (is.null(ylab)) {ylab <- "Suitability"}

    plotcurve_args <- list(x = response_out[, variable],
                           y = response_out$predicted,
                           type = "l",
                           xlab= xlab,
                           ylab = ylab,
                           col = col,
                           ...)

    # plot using do.call()
    do.call(plot, plotcurve_args)

    abline(v = limits,
           col = c("black", "black"),
           lty = c(2, 2),
           lwd = c(1, 1)
    )

  } else {
    # for multiple models in  model_list

    # extract the response of the variable for each models
    response_out <- lapply(model_list, function(x) {

      coefs <- if (inherits(x, "glmnet")) {
        names(x$betas)
      } else if (inherits(x, "glm")) {
        names(coef(x)[-1])
      }
      c1 <- any(c(variable, paste0("I(", variable, "^2)")) %in% coefs)
      c2 <- any(grepl(paste0("^", variable, ":"), coefs))
      c3 <- any(grepl(paste0(":", variable, "$"), coefs))

      if (any(c1, c2, c3)){

        out <- response(x, data, variable, new_data = new_data,
                        extrapolate = extrapolate,
                        categorical_variables = categorical_variables)
        return(out)

      } else {
        return(NULL)
      }
    })

    response_out <- do.call(rbind, response_out)
    limits <- range(data[, variable])


    x <- response_out[, variable]
    y <- response_out$predicted

    # Fit GAM model
    fitting <- mgcv::gam(y ~ s(x, bs = "cs"))

    # Generate predicted values and standard error.
    x_seq <- seq(min(x), max(x), length.out = 100)
    pred <- predict(fitting, newdata = data.frame(x = x_seq), se = T)

    # Extract predicted values, confidence intervals (95%), and standard errors
    y_pred <- pred$fit
    lower_ci <- y_pred - 1.96 * pred$se.fit
    upper_ci <- y_pred + 1.96 * pred$se.fit


    ## Plotting curve
    if (is.null(xlab)) {xlab <- variable}
    if (is.null(ylab)) {ylab <- "Suitability"}

    # Create a list of arguments to pass to the plot function
    plotcurve_args <- list(x = x, y = y, type = "n", xlab= xlab, ylab = ylab,
                           ...)

    # plot using do.call()
    do.call(plot, plotcurve_args)

    # Create shading interval using polygon
    x_polygon <- c(x_seq, rev(x_seq))
    y_polygon <- c(lower_ci, rev(upper_ci))
    polygon(x_polygon, y_polygon, col = "lightgrey", border = NA)

    # Add the regression curve
    lines(x_seq, y_pred, col = col)

    # It adds the calibration limits
    abline(v = limits,
           col = c("black", "black"),
           lty = c(2, 2),
           lwd = c(1, 1)
    )
  }
}

# It gets the response from an individual model
response <- function(model, data, variable, type = "cloglog", n = 100,
                     new_data = NULL, extrapolate = FALSE,
                     categorical_variables = NULL) {

  # initial tests
  if (missing(model) | missing(variable) | missing(data)) {
    stop("Arguments 'model', 'data', and 'variable' must be defined.")
  }

  if (!is.null(new_data)) {
    if (!class(new_data)[1] %in% c("matrix", "data.frame", "SpatRaster")) {
      stop("'new_data' must be of class 'matrix', 'data.frame', 'SpatRaster'")
    }
  }

  # Extract variable names used in the model
  vnames <- if (inherits(model, "glmnet")) {
    colSums(sapply(colnames(data), grepl, names(model$betas))) > 0
  } else if (inherits(model, "glm")) {
    colSums(sapply(colnames(data), grepl, names(coef(model)))) > 0
  }

  if (any(!variable %in% names(vnames))) {
    stop("The name of the 'variable' was not defined correctly.")
  }

  # Extract calibration data from the model object
  cal_data <- data[, vnames, drop = FALSE]

  # Extract the limits of the calibration data
  cal_maxs <-  apply(cal_data, 2, FUN = max)
  cal_mins <-  apply(cal_data, 2, FUN = min)

  # Get the average of all variables

  ####Check - For deal with categorical variables####
  means <- colMeans(cal_data[sapply(cal_data, is.numeric)])
  if(!is.null(categorical_variables) & vnames[categorical_variables]){
    mode_cat <- sapply(categorical_variables, function(x){
      as.numeric(names(which.max(table(cal_data[, x]))))
    })
    means <- c(means, mode_cat)
  }
  #########################################

  if (is.null(new_data)) {
    if (extrapolate) {

      rr <- range(cal_data[, variable]) # range of the calibration data
      extension <- 0.11 * diff(rr)

      l_limit <- rr[1] - extension
      u_limit <- rr[2] + extension

      rangev <- c(l_limit, u_limit)

    } else {
      rangev <- c(cal_mins[variable], cal_maxs[variable])
    }

  } else {

    if (class(new_data)[1] == "SpatRaster") {
      rangev <- terra::minmax(new_data[[variable]])
    } else {
      rangev <- range(new_data[, variable])
    }

  }

  newvar <- seq(rangev[1], rangev[2], length = n)

  m <- data.frame(matrix(means, n, length(means), byrow = T))
  colnames(m) <- names(means)

  m[, variable] <- newvar

  # Predict using glmnet or glm
  m$predicted <- if (inherits(model, "glmnet")) {
    predict.glmnet_mx(model, m, type = type)
  } else if (inherits(model, "glm")) {
    predict.glm(model, newdata = m, type = "response")
  }

  return(m[,c(variable, "predicted")])

}
