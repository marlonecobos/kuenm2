#' Bivariate response plot for fitted models
#'
#' @description
#' A plot of suitability prediction in a two-dimensional environmental space.
#'
#' @usage
#' bivariate_response(models, variable1 , variable2, modelID = NULL, n = 500,
#'                    new_data = NULL, extrapolate = FALSE, add_bar = TRUE ,
#'                    add_limits = TRUE, color_palette	= NULL,
#'                    xlab = NULL, ylab = NULL, ...)
#
#' @param models an object of class `fitted_models` returned by the
#' \code{\link{fit_selected}}() function.
#' @param variable1 (character) name of the variable to be plotted in x axis.
#' @param variable2 (character) name of the variable to be plotted in y axis.
#' @param modelID (character) name of the ModelID presents in the fitted
#' object. Default = NULL.
#' @param n (numeric) the number of breaks for plotting grid. Default = 500
#' @param new_data a `SpatRaster`, data.frame, or  matrix of variables
#' representing an area of interest. Default = NULL.
#' @param extrapolate (logical) whether to allow extrapolation to study the
#' behavior of the response outside the calibration limits. Ignored if
#' `new_data` is defined. Default = TRUE.
#' @param add_bar (logical) whether to add bar legend. Default = TRUE.
#' @param add_limits (logical) whether to add calibration limits if
#' `extrapolate = TRUE`. Default = TRUE.
#' @param color_palette (function) a color palette function to be used to assign
#' colors in the plot. The default, NULL uses `rev(hcl.colors(n, "terrain"))`.
#' @param xlab (character) a label for the x axis. The default, NULL, uses the
#' name defined in `variable1`.
#' @param ylab (character) a label for the y axis. The default, NULL, uses the
#' name defined in `variable2`.
#' @param ... additional arguments passed to \code{\link[graphics]{image}}.
#'
#' @return
#' A bivariate plot considering `variable1` and `variable2`.
#'
#' @export
#'
#' @importFrom stats predict coef
#' @importFrom graphics filled.contour par title rect axis mtext abline plot.new
#' @importFrom graphics plot.window
#' @importFrom terra minmax
#' @importFrom grDevices hcl.colors
#'
#' @seealso
#' [response_curve()]
#'
#' @examples
#' # Example with glmnet
#' # Import example of fitted_models (output of fit_selected())
#' data(fitted_model_maxnet, package = "kuenm2")
#'
#' # Response curve (notice response affected by covariance)
#' bivariate_response(models = fitted_model_maxnet, modelID = "Model_13",
#'                    variable1 = "bio_1", variable2 = "bio_12")
#'
#' # Example with glm
#' # Import example of fitted_models (output of fit_selected())
#' data(fitted_model_glm, package = "kuenm2")
#'
#' # Response curve
#' bivariate_response(models = fitted_model_glm, modelID = "Model_1",
#'                    variable1 = "bio_1", variable2 = "bio_7")

bivariate_response <- function(models,
                               variable1 ,
                               variable2,
                               modelID = NULL,
                               n = 500,
                               new_data = NULL,
                               extrapolate = FALSE,
                               add_bar = TRUE,
                               add_limits = TRUE,
                               color_palette = NULL,
                               xlab = NULL,
                               ylab = NULL,
                               ...) {

  # initial tests
  if (missing(models) | missing(variable1) | missing(variable2)) {
    stop("Argument 'models', 'variables' must be defined.")
  }

  if (!inherits(models, "fitted_models")) {
    stop("Argument 'models' must be a 'fitted_models' object.")
  }
  if (!inherits(variable1, "character")) {
    stop("Argument 'variable1' must be a 'character'.")
  }
  if (!inherits(variable2, "character")) {
    stop("Argument 'variable2' must be a 'character'.")
  }
  if (!inherits(n, "numeric")) {
    stop("Argument 'n' must be 'numeric'.")
  }

  if (!is.null(new_data)) {
    if (!class(new_data)[1] %in% c("matrix", "data.frame", "SpatRaster")) {
      stop("'new_data' must be of class 'matrix', 'data.frame', 'SpatRaster'.")
    }
  }

  if (!is.null(xlab) & !inherits(xlab, "character")) {
    stop("Argument 'xlab' must be NULL or a 'character'.")
  }
  if (!is.null(ylab) & !inherits(ylab, "character")) {
    stop("Argument 'ylab' must be NULL or a 'character'.")
  }

  if (!is.null(models$categorical_variables)) {
    if (any(c(variable1, variable2) %in% models$categorical_variables)) {
      stop("This function is designed exclusively for continuous variables.")
    }
  }

  if (is.null(xlab)) {
    xlab <- variable1
  }
  if (is.null(ylab)) {
    ylab <- variable2
  }

  if (is.null(color_palette)) {
    color_palette <- function(n) {
      rev(hcl.colors(n, "terrain"))
    }
  }

  model <- models[["Models"]][[modelID]][["Full_model"]]
  data <- models[["calibration_data"]]
  variables <- c(variable1, variable2)


  # Check if variables are defined in the model's formula as:
  # lineal, quadratic or product.
  coefs <- if (inherits(model, "glmnet")) {
    names(model$betas)
  } else if (inherits(model, "glm")) {
    names(coef(model)[-1])
  }

  c11 <- any(c(variable1, paste0("I(", variable1, "^2)")) %in% coefs)
  c12 <- any(grepl(paste0(variable1, ":"), coefs))
  c13 <- any(grepl(paste0(":", variable1,"$"), coefs))

  c21 <- any(c(variable2, paste0("I(", variable2, "^2)")) %in% coefs)
  c22 <- any(grepl(paste0(variable2, ":"), coefs))
  c23 <- any(grepl(paste0(":", variable2,"$"), coefs))

  if (any(c11, c12, c13) == FALSE) {
    stop("Defined 'variable1' is not present in the models model.")
  } else if (any(c21, c22, c23) == FALSE) {
    stop("Defined 'variable2' is not present in the models model.")
  }

  # Extract calibration data from the model object
  # It gets only the variable names used in the models model
  vnames <- if (inherits(model, "glmnet")) {
    colSums(sapply(colnames(data), grepl, names(model$betas))) > 0
  } else {
    colSums(sapply(colnames(data), grepl, names(coef(model)))) > 0
  }
  cal_data <- data[, vnames]

  # Extract the limits of the calibration data
  cal_maxs <-  apply(cal_data, 2, FUN = max)
  cal_mins <-  apply(cal_data, 2, FUN = min)

  # # Get the average of all variables
  # means <- apply(cal_data, 2, FUN = mean)

  ####Check - For deal with categorical variables####
  means <- colMeans(cal_data[sapply(cal_data, is.numeric)])

  if (!is.null(models$categorical_variables)) {
    if (vnames[models$categorical_variables]) {
      mode_cat <- sapply(models$categorical_variables, function(x) {
        as.numeric(names(which.max(table(cal_data[, x]))))
      })
      means <- c(means, mode_cat)
    }
  }

  #########################################
  if (is.null(new_data)) {
    if (extrapolate) {
      rr1 <- range(cal_data[, variables[1]]) # range of the calibration data
      rr2 <- range(cal_data[, variables[2]]) # range of the calibration data

      l_limit1 <- rr1[1] - 0.11 * diff(rr1)
      u_limit1 <- rr1[2] + 0.11 * diff(rr1)

      l_limit2 <- rr2[1] - 0.11 * diff(rr2)
      u_limit2 <- rr2[2] + 0.11 * diff(rr2)

      rangev1 <- c(l_limit1, u_limit1)
      rangev2 <- c(l_limit2, u_limit2)
    } else {

      rangev1 <- range(cal_data[, variables[1]])
      rangev2 <- range(cal_data[, variables[2]])
    }
  } else {
    if (class(new_data)[1] == "SpatRaster") {
      rangev1 <- terra::minmax(new_data[[variables[1]]])
      rangev2 <- terra::minmax(new_data[[variables[2]]])

    } else {
      rangev1 <- range(new_data[, variables[1]])
      rangev2 <- range(new_data[, variables[2]])
    }
  }

  newvar <- expand.grid(x = seq(rangev1[1], rangev1[2], length = n),
                        y = seq(rangev2[1], rangev2[2], length = n))

  m <- data.frame(matrix(means, nrow(newvar), length(means), byrow = TRUE))
  colnames(m) <- names(means)

  m[, variables[1]] <- newvar[, 1]
  m[, variables[2]] <- newvar[, 2]

  # Response of the variables
  predicted <- if (inherits(model, "glmnet")) {
    predict.glmnet_mx(model, m, type = "cloglog")
  } else {
    predict(model, newdata = m, type = "response")
  }

  # Arguments for filled.contour
  # "x,y" locations of grid lines at which the values in z are measured
  # "z" a numeric matrix containing the values to be plotted
  x <- unique(m[, variables[1]])
  y <- unique(m[, variables[2]])
  z <- matrix(data = predicted, nrow = n, ncol = n)

  colors <- color_palette(10) # colors

  # Save the original par settings
  original_par <- par(no.readonly = TRUE)

  if (add_bar == FALSE) {
    # Plot the main image
    image(x, y, z, zlim = c(0, 1), col = colors,
          xlab = xlab, ylab = ylab, useRaster = FALSE, ...)

    if (extrapolate & add_limits ) {
      abline(v = c(cal_mins[variables[1]], cal_maxs[variables[1]]),
             h = c(cal_mins[variables[2]], cal_maxs[variables[2]]),
             lty = 2)
    }
  } else {

    layout(matrix(1:2, ncol = 2), widths = c(4, 1))  # Adjust widths to allocate space for the legend
    par(mar = c(5, 4, 4, 2) + 0.1)  # Set margins for the main plot

    # Plot the main image
    image(x, y, z, zlim = c(0,1), col = colors,
          xlab = xlab, ylab = ylab, useRaster = TRUE, ...)

    if (extrapolate & add_limits ) {
      abline(v = c(cal_mins[variables[1]], cal_maxs[variables[1]]),
             h = c(cal_mins[variables[2]], cal_maxs[variables[2]]),
             lty = 2)
    }

    # Add the color bar legend
    par(mar = c(5, 0, 4, 2) + 0.1)  # Set margins for the legend
    color_levels <- seq(0, 1, length.out = length(colors) + 1)
    legend_y <- seq(min(y), max(y), length.out = length(colors) + 1)

    # Plot the legend
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(legend_y))

    # Draw rectangles for the legend
    for (i in seq_along(colors)) {
      rect(0.25, legend_y[i], 0.5, legend_y[i + 1], col = colors[i],
           border = "transparent")
    }
    rect(0.25, min(legend_y), 0.5, max(legend_y), col = NA, border = "black")

    # Add text labels to the legend
    axis(4, at = legend_y, labels = round(color_levels, 2), las = 1,
         tick = FALSE, pos = 0.5  )
    mtext("Suitability", side = 3, line = 0, cex = 1.2)

    # Restore the original par settings
    par(original_par)
  }
}
