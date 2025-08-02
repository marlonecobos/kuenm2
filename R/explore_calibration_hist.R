#' Explore variable distribution for occurrence and background points
#'
#' @description
#' This function prepares data to generate overlaid histograms to visualize the
#' distribution of predictor variables for occurrence (presence) and background
#' points.
#'
#' @usage
#' explore_calibration_hist(data, include_m = FALSE, raster_variables = NULL,
#'                          magnify_occurrences = 2, breaks = 15)
#'
#' @param data an object of class `prepared_data` returned by the
#' \code{\link{prepare_data}}() function.
#' @param include_m (logical) whether to include data for plotting the histogram
#' of the entire area from which background points were sampled. Default is
#' FALSE, meaning only background and presence information will be plotted.
#' @param raster_variables (SpatRaster) predictor variables used to prepared the
#' data with `prepared_data`. Only applicable if `include_m` is TRUE.
#' @param magnify_occurrences (numeric) factor by which the frequency of
#' occurrences is magnified for better visualization. Default is 2, meaning
#' occurrence frequencies in the plot will be doubled.
#' @param breaks (numeric) a single number giving the desired number of
#' intervals in the histogram.
#'
#' @importFrom stats na.omit quantile
#' @importFrom graphics hist
#' @importFrom terra as.data.frame
#'
#' @export
#'
#' @seealso
#' [plot_explore_calibration()]
#'
#' @returns
#' A list of with information to plot informative histograms to explore data
#' to be used in the modeling process. Histogram plots can be plotted with
#' the function [plot_explore_calibration()].
#'
#' @examples
#' # Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' # Import occurrences
#' data(sp_swd, package = "kuenm2")
#'
#' # Explore calibration data
#' calib_hist <- explore_calibration_hist(data = sp_swd,
#'                                        raster_variables = var,
#'                                        include_m = TRUE)
#'
#' # To visualize results use the function plot_explore_calibration()

explore_calibration_hist <- function(data, include_m = FALSE,
                                     raster_variables = NULL,
                                     magnify_occurrences = 2,
                                     breaks = 15) {
  #Check data
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (!inherits(data, "prepared_data")) {
    stop("'data' must be a 'prepared_data' object.")
  }
  if (!inherits(include_m, "logical")) {
    stop("'include_m' must be 'logical'.")
  }
  if (!inherits(magnify_occurrences, "numeric")) {
    stop("'magnify_occurrences' must be 'numeric'.")
  }
  if (!inherits(breaks, "numeric")) {
    stop("'breaks' must be 'numeric'.")
  }


  #Subset raster_variables
  if (include_m) {
    if (is.null(raster_variables)) {
      stop("If 'include_m' is TRUE, 'raster_variables' can not be NULL.")
    }
    if (!inherits(raster_variables, "SpatRaster")) {
      stop("'raster_variables' must be a 'SpatRaster'.")
    }

    #Check pca
    if(!is.null(data$pca)){
      if (!("vars_out" %in% names(data$pca))) {
        raster_variables <- terra::predict(raster_variables[[data$pca$vars_in]],
                                           data$pca)
      } else {
        raster_variables <- c(terra::predict(raster_variables[[data$pca$vars_in]],
                                             data$pca),
                              raster_variables[[data$pca$vars_out]])
      }
      raster_variables <- raster_variables[[colnames(data$calibration_data)[-1]]]
    }

    r_out <- setdiff(terra::names(raster_variables),
                     colnames(data$calibration_data)[-1])
    if (length(r_out) >= 1) {
      stop("Variables in 'data' don't match exactly those in 'raster_variables'.")
    }
  }

  #Extract calibration data
  d <- data$calibration_data

  #Convert categorical variables to numeric
  if (!is.null(data$categorical_variables)) {
    d[, data$categorical_variables] <- as.numeric(levels(
      d[, data$categorical_variables]))[d[, data$categorical_variables]]
  }

  # Split background and presence data
  bg_data <- d[d$pr_bg == 0, ]
  pr_data <- d[d$pr_bg == 1, ]

  #Get variables
  v <- setdiff(colnames(d), "pr_bg") #All variables
  cont_v <- setdiff(v, data$categorical_variables) #Continuous variables
  cat_v <- data$categorical_variables #Categorical variables

  if (include_m) {
  #Get values of entire M, if necessary
  m_data <- terra::as.data.frame(raster_variables)
  #Change orders
  m_data <- m_data[, v]}

  # preparing histogram information
  var_hist <- lapply(v, function(i) {
    #print(i)

    #Get means, cl and range
    if (!(i %in% data$categorical_variables)) {

      #Pretty breaks
      pr_bg <- stats::na.omit(c(pr_data[[i]], bg_data[[i]]))
      if (include_m) {
        pr_bg <- stats::na.omit(c(pr_bg, m_data[[i]]))
      }
      ax <- pretty((min(pr_bg) - 0.1):max(pr_bg), n = breaks)
      df <- diff(ax[1:2])
      ax <- c(ax[1] - df, ax, ax[length(ax)] + df)
      #Check if ax cover all values in the variable
      if(min(ax) > min(pr_bg)){
        ax <- c(min(ax) - df, ax)
      }
      if(max(ax) < max(pr_bg)){
        ax <- c(ax, max(ax) + df)
      }

      #M
      if (include_m) {
        m_m <- mean(m_data[[i]], na.rm = TRUE)
        l_m <- range(m_data[[i]], na.rm = TRUE)
        cl_m <- stats::quantile(m_data[[i]], c(0.05, 0.95), na.rm = TRUE)
        hist_m <- graphics::hist(m_data[[i]], breaks = ax, plot = FALSE)
      } else {
        m_m <- NA
        l_m <- NA
        cl_m <- NA
        hist_m <- NA
      }
      #Background
      m_bg <- mean(bg_data[[i]], na.rm = TRUE)
      l_bg <- range(bg_data[[i]], na.rm = TRUE)
      cl_bg <- stats::quantile(bg_data[[i]], c(0.05, 0.95), na.rm = TRUE)
      hist_bg <- graphics::hist(bg_data[[i]], breaks = ax, plot = FALSE)

      #Presence
      m_pr <- mean(pr_data[[i]], na.rm = TRUE)
      l_pr <- range(pr_data[[i]], na.rm = TRUE)
      cl_pr <- stats::quantile(pr_data[[i]], c(0.05, 0.95), na.rm = TRUE)
      hist_pr <- graphics::hist(rep(pr_data[[i]], magnify_occurrences),
                     breaks = ax, plot = FALSE)
    } else {
      #If categorical, get mode instead of mean
      #M
      if (include_m) {
        table_m <- table(m_data[[i]])
        m_m <- as.numeric(names(table_m[table_m == max(table_m)])[[1]])
        l_m <- range(m_data[[i]], na.rm = TRUE)
        cl_m <- NA
        hist_m <- table_m
      } else {
        m_m <- NA
        l_m <- NA
        cl_m <- NA
        hist_m <- NA
      } #End of plot M
      #Background
      table_bg <- table(bg_data[[i]])
      m_bg <- as.numeric(names(table_bg[table_bg == max(table_bg)])[[1]])
      l_bg <- range(bg_data[[i]], na.rm = TRUE)
      cl_bg <- NA
      hist_bg <- table_bg
      #Presence
      table_pr <- table(pr_data[[i]])
      m_pr <- as.numeric(names(table_pr[table_pr == max(table_pr)])[[1]])
      l_pr <- range(pr_data[[i]], na.rm = TRUE)
      cl_pr <- NA
      hist_pr <- table_pr

    } #End of else (categorical)

    new_explore_back(hist_m = hist_m, hist_bg = hist_bg, hist_pr = hist_pr,
                     mean_m = m_m, mean_bg = m_bg, mean_pr = m_pr,
                     cl_m = cl_m, cl_bg = cl_bg, cl_pr = cl_pr,
                     range_m = l_m, range_bg = l_bg, range_pr = l_pr)
  })

  names(var_hist) <- v

  #Summarize results
  if (include_m) {
    n_all <- nrow(stats::na.omit(m_data))
  } else {n_all <- NA}


  summa <- data.frame(n_all = n_all,
                      n_background = nrow(bg_data),
                      n_occurrences = nrow(pr_data),
                      sample_percentage = nrow(bg_data) / n_all * 100)

  # results
  return(new_explore_list(summary = summa, exploration_stats = var_hist,
         continuous_variables = cont_v,
         categorical_variables = cat_v))
} #end of function
