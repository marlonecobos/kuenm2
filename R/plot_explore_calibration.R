#' Plot histograms from an explore_calibration object
#'
#' @description
#' Customizes and plots histograms stored in an `explore_calibration` object
#' generated with the `explore_calibration_hist` function.
#'
#' @usage plot_explore_calibration(explore_calibration, color_m = "grey",
#'        color_background = "#56B4E9", color_presence = "#009E73", alpha = 0.4,
#'        lines = TRUE, which_lines = c("range", "cl", "mean"), lty_range = 1,
#'        lty_cl = 2, lty_mean = 3, lwd_range = 3, lwd_cl = 2, lwd_mean = 2,
#'        xlab = NULL, ylab = NULL, mfrow = c(1, 1))
#' @param explore_calibration an object of class `explore_calibration` generated
#'        by the `explore_calibration_hist` function.
#' @param color_m (character) color used to fill the histogram bars for the
#'        entire area (M). Default is "grey".
#' @param color_background (character) color used to fill the histogram bars for
#'        background data. Default is "#56B4E9".
#' @param color_presence (character) color used to fill the histogram bars for
#'        presence data. Default is "#009E73".
#' @param alpha (numeric) opacity factor to fill the bars, typically in the
#'        range [0,1]. Default is 0.4.
#' @param lines (logical) whether to add vertical lines to the plot representing
#'        the range, confidence interval, and mean of variables.
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
#'        the plot layout, e.g., c(rows, columns). Default is c(1, 1), meaning
#'        one plot per variable.
#'
#' @importFrom grDevices adjustcolor
#' @importFrom graphics par abline box barplot
#' @importFrom stats na.omit setNames
#' @export
#'
#' @examples
#' #Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#' #Import occurrences
#' data(occ_data, package = "kuenm2")
#' #Prepare data for glmnet model
#' sp_swd <- prepare_data(algorithm = "maxnet", occ = occ_data,
#'                        species = occ_data[1, 1], x = "x", y = "y",
#'                        raster_variables = var,
#'                        categorical_variables = "SoilType",
#'                        n_background = 500)
#' #Explore calibration data
#' calib_hist <- explore_calibration_hist(data = sp_swd, raster_variables = var,
#'                                        include_m = TRUE,
#'                                        magnify_occurrences = 2)
#' #Plot histograms
#' plot_explore_calibration(explore_calibration = calib_hist, mfrow = c(2,3))
#'
plot_explore_calibration <- function(explore_calibration, color_m = "grey",
                                    color_background = "#56B4E9",
                                    color_presence = "#009E73",
                                    alpha = 0.4, lines = TRUE,
                                    which_lines = c("range", "cl", "mean"),
                                    lty_range = 1, lty_cl = 2, lty_mean = 3,
                                    lwd_range = 3, lwd_cl = 2, lwd_mean = 2,
                                    xlab = NULL, ylab = NULL,
                                    mfrow = c(1, 1)){
  #Check errors####
  if (!inherits(explore_calibration, "explore_calibration")) {
    stop("'explore_calibration' must be of class 'explore_calibration', not ",
         class(explore_calibration))
  }
  if(!inherits(color_m, "character"))
    stop("'color_m' must be a 'character', not ", class(color_m))
  if(!inherits(color_background, "character"))
    stop("'color_background' must be a 'character', not ", class(color_background))
  if(!inherits(color_presence, "character"))
    stop("'color_presence' must be a 'character', not ", class(color_presence))
  if(!inherits(alpha, "numeric"))
    stop("'alpha' must be 'numeric', not ", class(alpha))
  if(!inherits(lines, "logical"))
    stop("'lines' must be 'logical', not ", class(lines))

  out_lines <- setdiff(which_lines, c("range", "cl", "mean"))
  if(length(out_lines) > 0)
    stop("'which_lines' specified are not valid")

  if(!inherits(lty_range, "numeric"))
    stop("'lty_range' must be 'numeric', not ", class(lty_range))
  if(!inherits(lty_cl, "numeric"))
    stop("'lty_cl' must be 'numeric', not ", class(lty_cl))
  if(!inherits(lty_mean, "numeric"))
    stop("'lty_mean' must be 'numeric', not ", class(lty_mean))
  if(!inherits(lwd_range, "numeric"))
    stop("'lwd_range' must be 'numeric', not ", class(lwd_range))
  if(!inherits(lwd_cl, "numeric"))
    stop("'lwd_cl' must be 'numeric', not ", class(lwd_cl))
  if(!inherits(lwd_mean, "numeric"))
    stop("'lwd_mean' must be 'numeric', not ", class(lwd_mean))
  #### End of checking errors ####

  #Adjust bar colors
  color_m_b <- grDevices::adjustcolor(color_m, alpha)
  color_background_b <- grDevices::adjustcolor(color_background, alpha)
  color_presence_b <- grDevices::adjustcolor(color_presence, alpha)


  #Get variables
  v <- names(explore_calibration$exploration_stats)

  #y Labels
  ylab <- ifelse(is.null(ylab), "Frequency", ylab)
  #Check x labels
  if(!is.null(xlab)){
    if(length(xlab) != length(v)){
      stop("If xlab is not NULL, the lenght of xlab must be the same as the number of variables in explore_calibration")
    }
  }

  #Set mfrow
  graphics::par(mfrow = mfrow)

  #Loop variables
  for(i in v){

    #Get xlab
    xlab_i <- ifelse(is.null(xlab), i, xlab[[i]])

    #Get variable histogram
    var_res <- explore_calibration$exploration_stats[[i]]

    # plot
    #Continuous variables
    if(i %in% explore_calibration$continuous_variables){
    if(all(!is.na(var_res$hist_m))){
    plot(var_res$hist_m, col = color_m_b, main = "", xlab = xlab_i,
         border = color_m,
         freq = TRUE, ylab = ylab)
      add_next <- TRUE} else {add_next <- FALSE}
    plot(var_res$hist_bg, col = color_background_b, xlab = xlab_i,
         add = add_next, border = color_background)
    plot(var_res$hist_pr, col = color_presence_b, add = TRUE,
         border = color_presence)
    } #End of is continuous

    if(i %in% explore_calibration$categorical_variables){
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

      if(all(!is.na(var_res$hist_m))){
        graphics::barplot(freq_m, col = color_m_b, main = "", xlab = xlab_i,
             border = color_m, ylab = ylab)
        add_next <- TRUE} else {add_next <- FALSE}
      graphics::barplot(freq_bg, col = color_background_b, xlab = xlab_i,
           add = add_next, border = color_background)
      graphics::barplot(freq_pr, col = color_presence_b, add = TRUE,
           border = color_presence)
    } #End of is categorical

    #Add lines?
    if(lines){
      #For continuous variables
      if(i %in% explore_calibration$continuous_variables){
        if("range" %in% which_lines){
          graphics::abline(v = var_res$range_m, col = color_m, lwd = lwd_range,
                           lty = lty_range)
          graphics::abline(v = var_res$range_bg, col = color_background,
                           lwd = lwd_range, lty = lty_range)
          graphics::abline(v = var_res$range_pr, col = color_presence,
                           lwd = lwd_range, lty = lty_range)
        }

        if("cl" %in% which_lines){
          graphics::abline(v = var_res$cl_m, col = color_m, lwd = lwd_cl,
                           lty = lty_cl)
          graphics::abline(v = var_res$cl_bg, col = color_background,
                           lwd = lwd_cl, lty = lty_cl)
          graphics::abline(v = var_res$cl_pr, col = color_presence,
                           lwd = lwd_cl, lty = lty_cl)
        }

        if("mean" %in% which_lines){
          graphics::abline(v = var_res$mean_m, col = color_m, lwd = lwd_mean,
                           lty = lty_mean)
          graphics::abline(v = var_res$mean_bg, col = color_background,
                           lwd = lwd_mean, lty = lty_mean)
          graphics::abline(v = var_res$mean_pr, col = color_presence,
                           lwd = lwd_mean, lty = lty_mean)
        }
        } #End of continuous

        #For categorical variables
        if(i %in% explore_calibration$categorical_variables){
          #Find position of values in barplot
          bp <- as.numeric(barplot(freq_bg, plot = FALSE))
          names(bp) <- all_categories


          if("range" %in% which_lines){
            #Find positions
            range_m <- bp[names(bp) %in% var_res$range_m] - 0.075
            range_bg <- bp[names(bp) %in% var_res$range_bg]
            range_pr <- bp[names(bp) %in% var_res$range_pr] + 0.075
            #Plot
            graphics::abline(v = range_m, col = color_m, lwd = lwd_range,
                             lty = lty_range)
            graphics::abline(v = range_bg, col = color_background,
                             lwd = lwd_range, lty = lty_range)
            graphics::abline(v = range_pr, col = color_presence, lwd = lwd_range,
                             lty = lty_range)
          }

          if("mean" %in% which_lines){
            #Find positions
            mean_m <- bp[names(bp) %in% var_res$mean_m] - 0.075
            mean_bg <- bp[names(bp) %in% var_res$mean_bg]
            mean_pr <- bp[names(bp) %in% var_res$mean_pr] + 0.075
            #Plot
            graphics::abline(v = mean_m, col = color_m, lwd = lwd_mean,
                             lty = lty_mean)
            graphics::abline(v = mean_bg, col = color_background, lwd = lwd_mean,
                             lty = lty_mean)
            graphics::abline(v = mean_pr, col = color_presence, lwd = lwd_mean,
                             lty = lty_mean)
          }
        }#End of categorical


      } #End of lines

    #Add box
    graphics::box(bty = "l")
  } #End of for in

  #Reset layout
  graphics::par(mfrow = c(1, 1))
  return(invisible(NULL))
}
