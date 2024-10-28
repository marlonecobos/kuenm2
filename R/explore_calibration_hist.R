#' Explore variable distribution for occurrence and background points
#'
#' @description
#' This function generates overlaid histograms to visualize the distribution of
#' predictor variables for occurrence (presence) and background points.
#'
#' @param data an object of class `prepare_data` returned by the
#'            \code{\link{prepare_data}}() function.
#' @param color_background (character) color used to fill the histogram bars for
#'        background data. Default is "#0000FF80".
#' @param color_presence (character) color used to fill the histogram bars for
#'        presence data. Default is "#FF000080".
#' @param mfrow (numeric) a vector specifying the number of rows and columns in
#'        the plot layout, e.g., c(rows, columns). Default is c(1, 1), meaning
#'        one plot per variable.
#' @param plot_median (logical) whether to plot a vertical dashed line
#'        representing the median of the presence points on the histogram.
#'        Default is TRUE.
#' @param ... further arguments and graphical parameters passed to
#'        \code{\link{hist}}(), allowing customization of titles, axes, and
#'        other elements.
#'
#' @importFrom stats median
#' @importFrom graphics hist abline par
#' @export
#'
#' @usage explore_calibration_hist(data, color_background = "#0000FF80",
#'                              color_presence = "#FF000080",
#'                              mfrow = c(1, 1), plot_median = TRUE, ...)
#' @examples
#' #Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#' #Import occurrences
#' data(occ_data, package = "kuenm2")
#' #Prepare data for glmnet model
#' sp_swd <- prepare_data(model_type = "glmnet", occ = occ_data,
#'                        species = occ_data[1, 1], x = "x", y = "y",
#'                        spat_variables = var, mask = NULL,
#'                        categorical_variables = "SoilType",
#'                        do_pca = FALSE, deviance_explained = 95,
#'                        min_explained = 5, center = TRUE, scale = TRUE,
#'                        write_pca = FALSE, output_pca = NULL, nbg = 2000,
#'                        kfolds = 4, weights = NULL, min_number = 2,
#'                        min_continuous = NULL,
#'                        features = c("l", "q", "p", "lq", "lqp"),
#'                        regm = c(0.1, 1, 2, 3, 5),
#'                        include_xy = TRUE,
#'                        write_file = FALSE, file_name = NULL,
#'                        seed = 1)
#' #Explore calibration data
#' explore_calibration_hist(data = sp_swd,
#'                          color_background = "#0000FF80",
#'                          color_presence = "#FF000080",
#'                          mfrow = c(2, 3), plot_median = TRUE,
#'                          breaks = "Scott")
#'
explore_calibration_hist <- function(data, color_background = "#0000FF80",
                                     color_presence = "#FF000080",
                                     mfrow = c(1, 1), plot_median = TRUE, ...){

  #Extract calibration data
  d <- data$calibration_data

  #Convert categorical variables to numeric
  if(!is.null(data$categorical_variables)){
    d[, data$categorical_variables] <- as.numeric(d[, data$categorical_variables])
  }

  # Split background and presence data
  bg_data <- d[d$pr_bg == 0, ]
  pr_data <- d[d$pr_bg == 1, ]

  #Get variables
  v <- setdiff(colnames(d), "pr_bg")

  # Define layout of histograms
  graphics::par(mfrow = mfrow)

  for(i in v){
    graphics::hist(bg_data[[i]], col = color_background, main = i, xlab = i,
                   xlim = range(d[[i]]), ...)
    graphics::hist(pr_data[[i]], col = color_presence, add = TRUE, ...) # Overlap

    if(plot_median & !(i %in% data$categorical_variables)){
    #Plot median
    m <- stats::median(pr_data[[i]])
    graphics::abline(v = m, col = "red", lwd = 2, lty = 2)
    }
    }

  # Reset layout
  graphics::par(mfrow = c(1, 1))
  return(invisible(NULL))
}

