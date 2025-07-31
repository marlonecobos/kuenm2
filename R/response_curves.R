#' Variable response curves for fitted models
#'
#' @description
#' A view of variable responses in fitted models. Responses based on single or
#' multiple models can be plotted.
#'
#' @usage
#' response_curve(models, variable, modelID = NULL, n = 100,
#'                show_variability = FALSE, show_lines = FALSE, data = NULL,
#'                new_data = NULL, averages_from = "pr_bg", extrapolate = TRUE,
#'                extrapolation_factor = 0.1, add_points = FALSE, p_col = NULL,
#'                l_limit = NULL, u_limit = NULL,
#'                xlab = NULL, ylab = "Suitability",
#'                col = "darkblue", ...)
#
#' @param models an object of class `fitted_models` returned by the
#' [fit_selected()] function.
#' @param variable (character) name of the variable to be plotted.
#' @param data data.frame or matrix of data used in the model calibration step.
#' Default = NULL.
#' @param modelID (character) vector of ModelID(s) to be considered in the
#' models object. By default all models are included. Default = NULL.
#' @param n (numeric) an integer guiding the number of breaks. Default = 100
#' @param show_variability (logical) if `modelID` is defined, shows variability
#' in response curves considering replicates. If `modelID` is not defined, the
#' default, FALSE, always shows variability from multiple models if present in
#' `models`.
#' @param show_lines (logical) whether to show variability by plotting lines for
#' all models or replicates. The default = FALSE, uses a GAM to characterize a
#' median trend and variation among modes or replicates. Ignored if
#' `show_variability` = FALSE.
#' @param new_data a `SpatRaster`, data.frame, or  matrix of variables
#' representing the range of variable values in an area of interest.
#' Default = NULL. It must be defined in case the model entered does not
#' explicitly include a data component.
#' @param averages_from (character) specifies how the averages or modes of the
#' variables are calculated. Available options are "pr" (to calculate averages
#' from the presence localities) or "pr_bg" (to use the combined set of presence
#' and background localities). Default is "pr_bg". See details.
#' @param extrapolate (logical) whether to allow extrapolation to study the
#' behavior of the response outside the calibration limits. Ignored if
#' `new_data` is defined. Default = TRUE.
#' @param extrapolation_factor (numeric) a multiplier used to calculate the
#' extrapolation range. Larger values allow broader extrapolation beyond the
#' observed data range. Default is 0.1.
#' @param l_limit (numeric) specifies the lower limit for the variable. Default
#' is \code{NULL}, meaning the lower limit will be calculated based on the
#' data's minimum value and the \code{extrapolation_factor}
#' (if \code{extrapolation = TRUE}).
#' @param u_limit (numeric) specifies the upper limit for the variable. Default
#' is \code{NULL}, meaning the upper limit will be calculated based on the
#' data's minimum value and the \code{extrapolation_factor}
#' (if \code{extrapolation = TRUE}).
#' @param xlab (character) a label for the x axis. The default, NULL, uses the
#' name defined in `variable`.
#' @param ylab (character) a label for the y axis. Default = "Suitability".
#' @param col (character) color for lines. Default = "darkblue".
#' @param col (character) color for lines. Default = "darkblue".
#' @param add_points (logical) if \code{TRUE}, adds the original observed
#' points (0/1) to the plot. This also sets `ylim = c(0, 1)`, unless these
#' limits are defined as part of `...`. Default = \code{FALSE}.
#' @param p_col (character) color for the observed points when
#' \code{add_points = TRUE}. Any valid R color name or hexadecimal code.
#' Default = "black".
#' @param ... additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @details
#' The response curves are generated with all other variables set to their mean
#' values (or mode for categorical variables), calculated either from the
#' presence localities (if averages_from = "pr") or from the combined set of
#' presence and background localities (if averages_from = "pr_bg").
#'
#' For categorical variables, a bar plot is generated with error bars showing
#' variability across models (if multiple models are included).
#'
#' @return
#' A plot with the response curve for a `variable`.
#'
#' @export
#'
#' @seealso
#' [bivariate_response()]
#'
#' @importFrom stats predict coef approx binomial complete.cases formula
#' @importFrom stats model.matrix predict.glm sd
#' @importFrom graphics abline polygon arrows image layout lines points
#' @importFrom grDevices adjustcolor
#' @importFrom terra minmax
#' @importFrom mgcv gam
#'
#' @examples
#' # Example with maxnet
#' # Import example of fitted_models (output of fit_selected())
#' data(fitted_model_maxnet, package = "kuenm2")
#'
#' #Response curves
#' response_curve(models = fitted_model_maxnet, variable = "bio_1",
#'                show_variability = TRUE)
#' response_curve(models = fitted_model_maxnet, variable = "bio_1",
#'                show_variability = TRUE, add_points = TRUE)
#' response_curve(models = fitted_model_maxnet, variable = "bio_1",
#'                show_variability = TRUE, show_lines = TRUE)
#' response_curve(models = fitted_model_maxnet, variable = "bio_1",
#'                modelID = "Model_192", show_variability = TRUE)
#' response_curve(models = fitted_model_maxnet, variable = "bio_1",
#'                modelID = "Model_192", show_variability = TRUE,
#'                show_lines = TRUE)
#'
#' # Example with GLM
#' # Import example of fitted_models (output of fit_selected())
#' data(fitted_model_glm, package = "kuenm2")
#'
#' #Response curves
#' response_curve(models = fitted_model_glm,
#'                variable = "bio_1", show_variability = TRUE)
#' response_curve(models = fitted_model_glm, variable = "bio_1",
#'                modelID = "Model_85", show_variability = TRUE)

response_curve <- function(models, variable, modelID = NULL, n = 100,
                           show_variability = FALSE, show_lines = FALSE,
                           data = NULL, new_data = NULL, averages_from = "pr_bg",
                           extrapolate = TRUE, extrapolation_factor = 0.1,
                           add_points = FALSE, p_col = NULL,
                           l_limit = NULL, u_limit = NULL,
                           xlab = NULL, ylab = "Suitability",
                           col = "darkblue", ...) {

  # initial tests
  if (missing(models) | missing(variable) ) {
    stop("Arguments 'models' and 'variable' must be defined.")
  }

  if (!inherits(models, "fitted_models")) {
    stop("Argument 'models' must be a fitted_models object.")
  }

  if (!inherits(variable, "character")) {
    stop("Argument 'variable' must be a 'character'.")
  }

  if (!is.null(modelID) & !inherits(modelID, "character")) {
    stop("Argument 'modelID' must be a 'character'.")
  }

  if (!inherits(n, "numeric")) {
    stop("Argument 'n' must be 'numeric'.")
  }

  if (!is.null(new_data)) {
    if (!class(new_data)[1] %in% c("matrix", "data.frame", "SpatRaster")) {
      stop("'new_data' must be of class 'matrix', 'data.frame', 'SpatRaster'")
    }
  }

  if (!inherits(extrapolation_factor, "numeric")) {
    stop("Argument 'extrapolation_factor' must be 'numeric'.")
  }

  if (!averages_from %in% c("pr_bg", "pr")) {
    stop("Argument 'averages_from' must be 'pr_bg' or 'pr'.")
  }

  if (!is.null(xlab) & !inherits(xlab, "character")) {
    stop("Argument 'xlab' must be NULL or a 'character'.")
  }

  if (!is.null(ylab) & !inherits(ylab, "character")) {
    stop("Argument 'ylab' must be NULL or a 'character'.")
  }



  # if data is not defined it is extracted from the models kuenm2 object
  if (is.null(data)) {
    data <- models$calibration_data
  }

  ## add a warming message indicating that the are not replicates.
  if (!is.null(modelID)) {

    if (!modelID %in% names(models[["Models"]])) {
      stop("'ModelID' is not correct, check the following: ",
           paste(names(models[["Models"]]), collapse = ", "))
    }

    # Handling replicates or the full model
    lmods <- length(models[["Models"]][[modelID]])
    if (show_variability & lmods > 1) {
      model_list <- models[["Models"]][[modelID]]
      model_list$Full_model <- NULL
    } else {
      model_list <- models[["Models"]][[modelID]]["Full_model"]
    }
  } else {
    # extract the slot Full_model for each model
    model_list <- lapply(models[["Models"]], function(x) {x$Full_model})
  }

  # Response curve for all selected models
  response_curve_consmx(model_list, variable = variable, data = data,
                        show_lines = show_lines, n = n,
                        new_data = new_data, extrapolate = extrapolate,
                        extrapolation_factor = extrapolation_factor,
                        xlab = xlab, ylab = ylab, col = col,
                        categorical_variables = models$categorical_variables,
                        averages_from = averages_from,
                        l_limit = l_limit, u_limit = u_limit,
                        add_points = add_points, p_col = p_col,
                        ...)
}

### Helpers functions for response curves

# Consensus response curve
response_curve_consmx <- function(model_list, variable, data,
                                  show_lines = FALSE,
                                  n = 100,
                                  extrapolate = FALSE,
                                  extrapolation_factor = 0.11,
                                  new_data = NULL,
                                  xlab = NULL, ylab = NULL, col = "darkblue",
                                  categorical_variables = NULL,
                                  averages_from = "pr_bg",
                                  l_limit = NULL,
                                  u_limit = NULL,
                                  add_points = FALSE, p_col = NULL,
                                  ...) {

  # initial tests
  if (missing(model_list) | missing(variable) | missing(data)) {
    stop("Arguments 'model_list', 'data', and 'variable' must be defined.")
  }

  if (!is.null(new_data)) {
    if (!class(new_data)[1] %in% c("matrix", "data.frame", "SpatRaster")) {
      stop("'new_data' must be of class 'matrix', 'data.frame', 'SpatRaster'.")
    }
  }

  if (!variable %in% colnames(data)) {
    stop("The name of the 'variable' was not defined correctly.")
  }


  if (length(model_list) == 1) {
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
    c4 <- any(grepl(paste0("categorical(", variable, ")"), coefs, fixed = TRUE)) # check categorical variables


    if (any(c1, c2, c3, c4) == FALSE) {
      stop("Defined 'variable' is not present in the model(s).")
    }

    response_out <- response(model = model_list[[1]], data = data,
                             variable = variable, n = n,
                             extrapolate = extrapolate,
                             extrapolation_factor = extrapolation_factor,
                             new_data = new_data,
                             categorical_variables = categorical_variables,
                             averages_from = averages_from,
                             l_limit = l_limit, u_limit = u_limit)

    # Plot response for categorical variable for a single model
    if (!is.null(categorical_variables) && variable %in% categorical_variables) {
      x <- response_out[, variable]
      y <- c(response_out$predicted)

      if (is.null(xlab)) {
        xlab <- variable
      }
      if (is.null(ylab)) {
        ylab <- "Suitability"
      }

      # Create a list of arguments to pass to the plot function
      plotcurve_args <- list(height = y, names.arg = x, col = "lightblue",
                             xlab= xlab, ylab = ylab, ...)

      # plot using do.call()
      do.call(barplot, plotcurve_args)

    } else {
      # Plot response for continuous variable for a single model
      limits <- range(data[, variable])

      ## Plotting curve
      if (is.null(xlab)) {
        xlab <- variable
      }
      if (is.null(ylab)) {
        ylab <- "Suitability"
      }

      plotcurve_args <- list(x = response_out[, variable],
                             y = response_out$predicted,
                             type = "l", xlab= xlab, ylab = ylab,
                             col = col, ...)

      if (add_points){
        if (is.null(p_col)){
          p_col <- adjustcolor("black", alpha.f = 0.5)
        }

        if (!"ylim" %in% names(plotcurve_args)) {
          plotcurve_args <- c(plotcurve_args, list(ylim = c(0, 1)))
        }
      }
      #print(plotcurve_args)

      # plot using do.call()
      do.call(plot, plotcurve_args)

      abline(v = limits, # It adds the calibration limits
             col = c("black", "black"), lty = c(2, 2), lwd = c(1, 1))

      if (add_points) { # Add points to the plot
        points(data[, c(variable, "pr_bg")], bg = p_col, pch = 21, cex = 0.6)
      }
    }

  } else {
    # for multiple models in model_list

    # extract the response of the variable for each models
    response_out <- lapply(model_list, function(x) {

      coefs <- if (inherits(x, "glmnet")) {
        names(x$betas)
      } else if (inherits(x, "glm")) {
        names(coef(x)[-1])
      }

      c1 <- any(c(variable, paste0("I(", variable, "^2)")) %in% coefs) # check linear or quadratic term
      c2 <- any(grepl(paste0(variable, ":"), coefs))                   # check product terms 1st position
      c3 <- any(grepl(paste0(":", variable,"$"), coefs))               # check product terms 2nd position
      c4 <- any(grepl(paste0("categorical(", variable, ")"), coefs, fixed = TRUE)) # check categorical variables

      if (any(c1, c2, c3, c4)) {

        out <- response(x, data = data, variable = variable,
                        new_data = new_data,
                        extrapolate = extrapolate,
                        extrapolation_factor = extrapolation_factor,
                        categorical_variables = categorical_variables,
                        averages_from = averages_from,
                        l_limit = l_limit, u_limit = u_limit)
        return(out)

      } else {
        return(NULL)
      }
    })

    if (show_lines) {
      response_out0 <- response_out
    }

    ## summary of responses
    response_out <- do.call(rbind, response_out)

    if (all(sapply(response_out, is.null))) {
      stop("Defined 'variable' is not present in the model(s).")
    }

    if (!is.null(categorical_variables) && variable %in% categorical_variables) {
      # Plot response for categorical variable for multiple models
      # Function to calculate mean and standard error
      calc_stats <- function(x) {
        n <- length(x)  # Number of observations
        mean_x <- mean(x)  # Mean
        sd_x <- sd(x)  # Standard deviation
        se_x <- sd(x) / sqrt(n)  # Standard error
        return(c(mean = mean_x, sd = sd_x, se = se_x))
      }

      # Compute statistics: mean and standard error
      stats <- aggregate(as.formula(paste("predicted ~", variable)),
                         data = response_out, FUN = calc_stats)
      stats <- do.call(data.frame, stats)
      stats <- stats[order(as.numeric(levels(stats[, variable]))), ]

      x <- stats[, variable] # Variable
      y <- stats[, 2] # Mean
      z <- stats[, 3] # Standard deviation
      #z <- stats[, 4] # Standard error

      if (is.null(xlab)) {xlab <- variable}
      if (is.null(ylab)) {ylab <- "Suitability"}

      # Create a list of arguments to pass to the plot function
      plotcurve_args <- list(height = y, names.arg = x, col = "lightblue",
                             xlab= xlab, ylab = ylab,
                             ...)

      # plot using do.call()
      bar_positions <- do.call(barplot, plotcurve_args)

      # Add error bars
      error_args <- list(x0 = bar_positions, y0 = y - z,
                         x1 = bar_positions, y1 = y + z,
                         angle = 90, code = 3, length = 0.05, col = "black",
                         lwd = 1.5)

      do.call(arrows, error_args)
    } else {
      # Plot response for continuous variable for multiple models
      limits <- range(data[, variable])

      x <- response_out[, variable]
      y <- response_out$predicted

      ## Plotting curve
      if (is.null(xlab)) {
        xlab <- variable
      }
      if (is.null(ylab)) {
        ylab <- "Suitability"
      }

      # Create a list of arguments to pass to the plot function
      plotcurve_args <- list(x = x, y = y, type = "n", xlab= xlab, ylab = ylab,
                             ...)

      if (add_points){
        if (is.null(p_col)){
          p_col <- adjustcolor("black", alpha.f = 0.5)
        }

        if (!"ylim" %in% names(plotcurve_args)) {
          plotcurve_args <- c(plotcurve_args, list(ylim = c(0, 1)))
        }
      }
      #print(plotcurve_args)

      # getting basic line plottinf arguments
      ltys <- plotcurve_args$lty
      lwds <- plotcurve_args$lwd

      if (show_lines) {
        # plot using do.call()
        do.call(plot, plotcurve_args)

        # adding lines
        col <- adjustcolor(col, alpha.f = 0.5)
        for (i in response_out0) {
          lines(i[, variable], i[, "predicted"], col = col, lty = ltys,
                lwd = lwds)
        }

      } else {
        # Fit GAM model
        fitting <- mgcv::gam(y ~ s(x, bs = "cs"))

        # Generate predicted values and standard error.
        x_seq <- seq(min(x), max(x), length.out = 100)
        pred <- predict(fitting, newdata = data.frame(x = x_seq), se = T)

        # Extract predicted values, confidence intervals (95%), and standard errors
        y_pred <- pred$fit
        lower_ci <- y_pred - 1.96 * pred$se.fit
        upper_ci <- y_pred + 1.96 * pred$se.fit

        # plot using do.call()
        do.call(plot, plotcurve_args)

        # Create shading interval using polygon
        x_polygon <- c(x_seq, rev(x_seq))
        y_polygon <- c(lower_ci, rev(upper_ci))
        polygon(x_polygon, y_polygon, col = "lightgrey", border = NA)

        # adding GAM line
        lines(x_seq, y_pred, col = col, lty = ltys, lwd = lwds)
      }

      # It adds the calibration limits
      abline(v = limits, col = c("black", "black"), lty = c(2, 2),
             lwd = c(1, 1))

      if (add_points) { # Add points to the plot
        points(data[, c(variable, "pr_bg")], bg = p_col, pch = 21, cex = 0.6)
      }
    }
  }
}


# It gets the response from an individual model
response <- function(model, data, variable, type = "cloglog", n = 100,
                     new_data = NULL, extrapolate = FALSE,
                     extrapolation_factor = 0.11,
                     categorical_variables = NULL,
                     averages_from = "pr_bg",
                     l_limit = NULL, u_limit = NULL) {

  # initial tests
  if (missing(model) | missing(variable) | missing(data)) {
    stop("Arguments 'model', 'data', and 'variable' must be defined.")
  }

  if (!is.null(new_data)) {
    if (!class(new_data)[1] %in% c("matrix", "data.frame", "SpatRaster")) {
      stop("'new_data' must be of class 'matrix', 'data.frame', 'SpatRaster'.")
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
  if (averages_from == "pr") {
    cal_data <- data[data$pr_bg == 1, vnames, drop = FALSE]
  } else if (averages_from == "pr_bg") {
    cal_data <- data[, vnames, drop = FALSE]
  }

  # Extract the limits of the calibration data
  cal_maxs <-  apply(cal_data, 2, FUN = max)
  cal_mins <-  apply(cal_data, 2, FUN = min)

  # Get the average of all variables

  # Dealing with categorical variables
  means <- colMeans(cal_data[sapply(cal_data, is.numeric)])
  if (!is.null(categorical_variables) && vnames[categorical_variables]) {
    mode_cat <- sapply(categorical_variables, function(x) {
      as.numeric(names(which.max(table(cal_data[, x]))))
    })
    means <- c(means, mode_cat)
  }

  if (is.null(new_data)) {
    if (extrapolate && !variable %in% categorical_variables) {

      rr <- range(cal_data[, variable]) # range of the calibration data
      extension <- extrapolation_factor * diff(rr)

      if (is.null(l_limit)) {
        l_limit <- rr[1] - extension
      }
      if (is.null(u_limit)) {
        u_limit <- rr[2] + extension
      }

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

  if (variable %in% categorical_variables) {
    newvar <- factor(levels(data[[variable]]))
  } else {
    newvar <- seq(rangev[1], rangev[2], length = n)
  }


  m <- data.frame(matrix(means, length(newvar), length(means), byrow = T))
  colnames(m) <- names(means)

  m[, variable] <- newvar

  # Predict using glmnet or glm
  m$predicted <- if (inherits(model, "glmnet")) {
    predict.glmnet_mx(model, m, type = type)
  } else if (inherits(model, "glm")) {
    predict.glm(model, newdata = m, type = "response")
  }

  return(m[, c(variable, "predicted")])

}
