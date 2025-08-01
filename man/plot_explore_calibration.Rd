% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_explore_calibration.R
\name{plot_explore_calibration}
\alias{plot_explore_calibration}
\title{Histograms to visualize data from explore_calibration objects}
\usage{
plot_explore_calibration(explore_calibration, color_m = "grey",
                         color_background = "#56B4E9",
                         color_presence = "#009E73", alpha = 0.4,
                         lines = FALSE, which_lines = c("cl", "mean"),
                         lty_range = 1, lty_cl = 2, lty_mean = 3,
                         lwd_range = 3, lwd_cl = 2, lwd_mean = 2,
                         xlab = NULL, ylab = NULL, mfrow = NULL)
}
\arguments{
\item{explore_calibration}{an object of class \code{explore_calibration} generated
by the \code{explore_calibration_hist} function.}

\item{color_m}{(character) color used to fill the histogram bars for the
entire area (M). Default is "grey".}

\item{color_background}{(character) color used to fill the histogram bars for
background data. Default is "#56B4E9".}

\item{color_presence}{(character) color used to fill the histogram bars for
presence data. Default is "#009E73".}

\item{alpha}{(numeric) opacity factor to fill the bars, typically in the
range 0-1. Default is 0.4.}

\item{lines}{(logical) whether to add vertical lines to the plot representing
the range, confidence interval, and mean of variables. Default = FALSE.}

\item{which_lines}{(character) a vector indicating which lines to plot.
Available options are "range", "cl" (confidence interval), and "mean".
Default is c("range", "cl", "mean").}

\item{lty_range}{(numeric) line type for plotting the ranges of variables.
Default is 1, meaning a solid line.}

\item{lty_cl}{(numeric) line type for plotting the confidence interval of
variables. Default is 2, meaning a dashed line.}

\item{lty_mean}{(numeric) line type for plotting the mean of variables.
Default is 3, meaning a dotted line.}

\item{lwd_range}{(numeric) line width for the line representing the range.
Default is 3.}

\item{lwd_cl}{(numeric) line width for the line representing the confidence
interval. Default is 2.}

\item{lwd_mean}{(numeric) line width for the line representing the mean.
Default is 2.}

\item{xlab}{(character) a vector of names for labeling the x-axis. It must
have the same length as the number of variables. Default is NULL,
meaning the labels will be extracted from the \code{explore_calibration}
object.}

\item{ylab}{(character) the label for the y-axis. Default is NULL, meaning
the y-axis will be labeled as "Frequency".}

\item{mfrow}{(numeric) a vector specifying the number of rows and columns in
the plot layout, e.g., c(rows, columns). Default is NULL, meaning
the grid will be arranged automatically based on the number of plots.}
}
\description{
Plots histograms to visualize data from an \code{explore_calibration} object
generated with the \code{explore_calibration_hist} function.
}
\examples{
# Import raster layers
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2"))

# Import occurrences
data(sp_swd, package = "kuenm2")

# Explore calibration data
calib_hist <- explore_calibration_hist(data = sp_swd,
                                       raster_variables = var,
                                       include_m = TRUE)

# Plot histograms
plot_explore_calibration(explore_calibration = calib_hist)
}
