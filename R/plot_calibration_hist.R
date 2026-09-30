#' Histograms to visualize data from explore_calibration objects
#'
#' @description
#' Plots histograms to visualize data from an `explore_calibration` object
#' generated with the `explore_calibration_hist` function.
#'
#' @usage
#' plot_calibration_hist(explore_calibration, variables = NULL,
#'                       color_m = "grey",
#'                       color_background = "#56B4E9",
#'                       color_presence = "#009E73", alpha = 0.4,
#'                       lines = FALSE, which_lines = c("cl", "mean"),
#'                       lty_range = 1, lty_cl = 2, lty_mean = 3,
#'                       lwd_range = 3, lwd_cl = 2, lwd_mean = 2,
#'                       xlab = NULL, ylab = NULL, mfrow = NULL,
#'                       show_legend = FALSE, legend_position = c("bottom"),
#'                       legend_offset = 0.08)
#'
#' @param explore_calibration an object of class `explore_calibration` generated
#'        by the `explore_calibration_hist` function.
#' @param variables (character) vector specifying the variables for which to
#'        plot histograms. Defaults to NULL, which plots histograms for all
#'        variables.
#' @param color_m (character) color used to fill the histogram bars for the
#'        entire area (M). Default is "grey".
#' @param color_background (character) color used to fill the histogram bars for
#'        background data. Default is "#56B4E9".
#' @param color_presence (character) color used to fill the histogram bars for
#'        presence data. Default is "#009E73".
#' @param alpha (numeric) opacity factor to fill the bars, typically in the
#'        range 0-1. Default is 0.4.
#' @param lines (logical) whether to add vertical lines to the plot representing
#'        the range, confidence interval, and mean of variables. Default = FALSE.
#' @param which_lines (character) a vector indicating which lines to plot.
#'        Available options are "range", "cl" (confidence interval), and "mean".
#'        Default is c("range", "cl", "mean").
#' @param lty_range (numeric) line type for plotting the ranges of variables.
#'        Default is 1, meaning a solid line.
#' @param lty_cl (numeric) line type for plotting the confidence interval of
#'        variables. Default is 2, meaning a dashed line.
#' @param lty_mean (numeric) line type for plotting the mean of variables.
#'        Default is 3, meaning a dotted line.
#' @param lwd_range (numeric) line width for the line representing the range.
#'        Default is 3.
#' @param lwd_cl (numeric) line width for the line representing the confidence
#'        interval. Default is 2.
#' @param lwd_mean (numeric) line width for the line representing the mean.
#'        Default is 2.
#' @param xlab (character) a vector of names for labeling the x-axis. It must
#'        have the same length as the number of variables. Default is NULL,
#'        meaning the labels will be extracted from the `explore_calibration`
#'        object.
#' @param ylab (character) the label for the y-axis. Default is NULL, meaning
#'        the y-axis will be labeled as "Frequency".
#' @param mfrow (numeric) a vector specifying the number of rows and columns in
#'        the plot layout, e.g., c(rows, columns). Default is NULL, meaning
#'        the grid will be arranged automatically based on the number of plots.
#' @param show_legend (logical) whether to display a shared legend for the
#'        histogram bar colors. Default is FALSE.
#' @param legend_position (character) position of the shared legend when
#'        `show_legend = TRUE`. Available options are "bottom", "right", and
#'        "top". Default is "bottom".
#' @param legend_offset (numeric) distance between the legend and the outer edge
#'        of the plotting area. Larger values move the legend closer to the
#'        plotting area. Default is 0.08.
#'
#' @importFrom grDevices adjustcolor
#' @importFrom graphics par abline box barplot plot
#' @importFrom stats na.omit setNames
#' @importFrom grid grid.rect gpar grid.text
#'
#' @return
#' No return value, called for side effects (plots histograms).
#'
#' @export
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
#' # Plot histograms
#' plot_calibration_hist(explore_calibration = calib_hist)

plot_calibration_hist <- function(explore_calibration,
                                  variables = NULL,
                                  color_m = "grey",
                                  color_background = "#56B4E9",
                                  color_presence = "#009E73",
                                  alpha = 0.4,
                                  lines = FALSE,
                                  which_lines = c("cl", "mean"),
                                  lty_range = 1,
                                  lty_cl = 2,
                                  lty_mean = 3,
                                  lwd_range = 3,
                                  lwd_cl = 2,
                                  lwd_mean = 2,
                                  xlab = NULL,
                                  ylab = NULL,
                                  mfrow = NULL,
                                  show_legend = FALSE,
                                  legend_position = c("bottom"),
                                  legend_offset = 0.08) {
  #Check errors####
  if (missing(explore_calibration)) {
    stop("Argument 'explore_calibration' must be defined.")
  }
  if (!inherits(explore_calibration, "explore_calibration")) {
    stop("'explore_calibration' must be a 'explore_calibration' object.")
  }
  if (!inherits(color_m, "character")) {
    stop("'color_m' must be a 'character'.")
  }
  if (!inherits(color_background, "character")) {
    stop("'color_background' must be a 'character'.")
  }
  if (!inherits(color_presence, "character")) {
    stop("'color_presence' must be a 'character'.")
  }
  if (!inherits(alpha, "numeric")) {
    stop("'alpha' must be 'numeric'.")
  }

  out_lines <- setdiff(which_lines, c("range", "cl", "mean"))
  if (length(out_lines) > 0) {
    stop("'which_lines' specified are not valid.")
  }

  if (!inherits(lty_range, "numeric")) {
    stop("'lty_range' must be 'numeric'.")
  }
  if (!inherits(lty_cl, "numeric")) {
    stop("'lty_cl' must be 'numeric'.")
  }
  if (!inherits(lty_mean, "numeric")) {
    stop("'lty_mean' must be 'numeric'.")
  }
  if (!inherits(lwd_range, "numeric")) {
    stop("'lwd_range' must be 'numeric'.")
  }
  if (!inherits(lwd_cl, "numeric")) {
    stop("'lwd_cl' must be 'numeric'.")
  }
  if (!inherits(lwd_mean, "numeric")) {
    stop("'lwd_mean' must be 'numeric'.")
  }
  if (!is.logical(show_legend) || length(show_legend) != 1L ||
      is.na(show_legend)) {
    stop("'show_legend' must be TRUE or FALSE.")
  }
  if(show_legend){
    legend_unmatch <- setdiff(legend_position, c("bottom", "right", "top"))
    if(length(legend_unmatch) > 0 | length(legend_position) > 1){
      stop("If show_legend is TRUE, legend_position must be 'bottom', 'right', or 'top'")
    }
  }

  #### End of checking errors ####

  #Adjust bar colors
  color_m_b <- grDevices::adjustcolor(color_m, alpha)
  color_background_b <- grDevices::adjustcolor(color_background, alpha)
  color_presence_b <- grDevices::adjustcolor(color_presence, alpha)

  #Get variables
  v <- names(explore_calibration$exploration_stats)

  if(!is.null(variables)){
    # Check variables
    v_out <- setdiff(variables, v)
    if(length(v_out) > 0){
      stop("The following variables are absent from 'explore_partition':\n",
           paste(v_out, collapse = ", "))
    }
    v <- intersect(v, variables)
  }

  #y Labels
  ylab <- ifelse(is.null(ylab), "Frequency", ylab)
  #Check x labels
  if (!is.null(xlab)) {
    if (length(xlab) != length(v)) {
      stop("If 'xlab' is not NULL, 'xlab' lenght must match the number of variables in 'explore_calibration'")
    }
  }

  #Par settings
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))

  # If show legend...
  if (show_legend) {
    oma <- graphics::par("oma")

    if (legend_position == "bottom") oma[1] <- oma[1] + 2.5
    if (legend_position == "top")    oma[3] <- oma[3] + 2.5
    if (legend_position == "right")  oma[4] <- oma[4] + 9

    graphics::par(oma = oma)
  }

  #Set mfrow
  if(is.null(mfrow)){ #If NULL, arrange automatically
    nl <- length(v)
    mfrow <- c(ceiling(nl / ceiling(sqrt(nl))), ceiling(sqrt(nl)))
    graphics::par(mfrow = mfrow)
  } else {
    graphics::par(mfrow = mfrow)
  }

  #Loop variables
  for(i in v) {

    #Get xlab
    xlab_i <- ifelse(is.null(xlab), i, xlab[[i]])

    #Get variable histogram
    var_res <- explore_calibration$exploration_stats[[i]]

    # plot


    #Continuous variables
    if (i %in% explore_calibration$continuous_variables) {
      if (all(!is.na(var_res$hist_m))) {
        graphics::plot(var_res$hist_m, col = color_m_b, main = "",
                       xlab = xlab_i, border = color_m, freq = TRUE,
                       ylab = ylab)
        add_next <- TRUE
      } else {
        add_next <- FALSE
      }
      graphics::plot(var_res$hist_bg, col = color_background_b, main = "",
                     xlab = xlab_i, add = add_next, border = color_background,
                     ylab = ylab)
      graphics::plot(var_res$hist_pr, col = color_presence_b, add = TRUE,
                     border = color_presence)
    } #End of is continuous

    if (i %in% explore_calibration$categorical_variables) {
      #Create comum x-axis for all
      all_categories <- stats::na.omit(sort(unique(c(names(var_res$hist_m),
                                                     names(var_res$hist_bg),
                                                     names(var_res$hist_pr)))))
      #Reorder
      all_categories <- sort(as.numeric(all_categories))


      #Add absent categories in bg and pr
      freq_bg <- stats::setNames(rep(0, length(all_categories)), all_categories)
      freq_bg[names(var_res$hist_bg)] <- var_res$hist_bg
      freq_pr <- stats::setNames(rep(0, length(all_categories)), all_categories)
      freq_pr[names(var_res$hist_pr)] <- var_res$hist_pr
      freq_m <- stats::setNames(rep(0, length(all_categories)), all_categories)
      freq_m[names(var_res$hist_m)] <- var_res$hist_m

      if (all(!is.na(var_res$hist_m))) {
        graphics::barplot(freq_m, col = color_m_b, main = "", xlab = xlab_i,
                          border = color_m, ylab = ylab)
        add_next <- TRUE
      } else {
        add_next <- FALSE
      }
      graphics::barplot(freq_bg, col = color_background_b, main = "",
                        xlab = xlab_i, ylab = ylab, add = add_next,
                        border = color_background)
      graphics::barplot(freq_pr, col = color_presence_b, add = TRUE,
                        border = color_presence)
    } #End of is categorical

    #Add lines?
    if (lines) {
      #For continuous variables
      if (i %in% explore_calibration$continuous_variables) {
        if ("range" %in% which_lines) {
          graphics::abline(v = var_res$range_m, col = color_m, lwd = lwd_range,
                           lty = lty_range)
          graphics::abline(v = var_res$range_bg, col = color_background,
                           lwd = lwd_range, lty = lty_range)
          graphics::abline(v = var_res$range_pr, col = color_presence,
                           lwd = lwd_range, lty = lty_range)
        }

        if ("cl" %in% which_lines) {
          graphics::abline(v = var_res$cl_m, col = color_m, lwd = lwd_cl,
                           lty = lty_cl)
          graphics::abline(v = var_res$cl_bg, col = color_background,
                           lwd = lwd_cl, lty = lty_cl)
          graphics::abline(v = var_res$cl_pr, col = color_presence,
                           lwd = lwd_cl, lty = lty_cl)
        }

        if ("mean" %in% which_lines) {
          graphics::abline(v = var_res$mean_m, col = color_m, lwd = lwd_mean,
                           lty = lty_mean)
          graphics::abline(v = var_res$mean_bg, col = color_background,
                           lwd = lwd_mean, lty = lty_mean)
          graphics::abline(v = var_res$mean_pr, col = color_presence,
                           lwd = lwd_mean, lty = lty_mean)
        }
      } #End of continuous

    } #End of lines

    #Add box
    graphics::box(bty = "l")
  } #End of for in

  if (show_legend) {
    cores <- c(color_m, color_background, color_presence)
    preenchimentos <- grDevices::adjustcolor(cores, alpha.f = alpha)
    rotulos <- c("Calibration area", "Background", "Presences")

    # Dimensões e margens externas, em polegadas
    omi <- graphics::par("omi")
    din <- graphics::par("din")

    if (legend_position == "right") {
      centro_x <- 1 - omi[4] / (2 * din[1])
      x_caixas <- rep(centro_x - 0.04, 3)
      x_textos <- rep(centro_x - 0.025, 3)
      y <- c(0.60, 0.50, 0.40)
    } else {
      y_centro <- if (legend_position == "top") {
        1 - omi[3] / (2 * din[2])
      } else {
        omi[1] / (2 * din[2])
      }

      x_caixas <- c(0.10, 0.42, 0.74)
      x_textos <- c(0.12, 0.44, 0.76)
      y <- rep(y_centro, 3)
    }

    # Control legend position
    if (legend_position == "top") {
      y <- y - legend_offset / din[2]
    } else if (legend_position == "bottom") {
      y <- y + legend_offset / din[2]
    } else {
      x_caixas <- x_caixas - legend_offset / din[1]
      x_textos <- x_textos - legend_offset / din[1]
    }

    # Add legend
    grid::grid.rect(
      x = x_caixas, y = y,
      width = 0.025, height = 0.025,
      gp = grid::gpar(fill = preenchimentos, col = cores)
    )
    grid::grid.text(
      rotulos, x = x_textos, y = y,
      just = "left", gp = grid::gpar(fontsize = 9)
    )
  }

  return(invisible(NULL))
}
