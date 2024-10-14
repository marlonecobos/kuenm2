### Helpers functions for response curves

# Consensus response curve
response_curve_consmx <- function(model_list, variable, data, n = 100,
                                  extrapolate = FALSE, new_data = NULL,
                                  xlab = NULL, ylab = NULL, col = "darkblue",
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
      names(coef(model))
    }

    c1 <- any(c(variable, paste0("I(", variable, "^2)")) %in% coefs) # check linear or quadratic term
    c2 <- any(grepl(paste0(variable, ":"), coefs))                   # check product terms 1st position
    c3 <- any(grepl(paste0(":", variable,"$"), coefs))               # check product terms 2nd position

    if (any(c1, c2, c3) == FALSE){
      stop("Defined 'variable' is not present in the fitted model.")
    }

    response_out <- response(model = model_list[[1]], data = data,
                             variable = variable, n = n,
                             extrapolate = extrapolate, new_data = new_data)

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

      coefs <- if (inherits(model, "glmnet")) {
        names(model$betas)
      } else if (inherits(model, "glm")) {
        names(coef(model))
      }
      c1 <- any(c(variable, paste0("I(", variable, "^2)")) %in% coefs)
      c2 <- any(grepl(paste0("^", variable, ":"), coefs))
      c3 <- any(grepl(paste0(":", variable, "$"), coefs))

      if (any(c1, c2, c3)){

        out <- response(x, data, variable, new_data = new_data,
                        extrapolate = extrapolate)
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
                     new_data = NULL, extrapolate = FALSE) {

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
  means <- colMeans(cal_data)

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



