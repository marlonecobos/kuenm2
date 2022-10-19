
explore_background <- function(calibration_variables, occurrence_data,
                               longitude_column, latitude_column,
                               magnify_occurrences = 3, breaks = 15,
                               sample_size = NULL, replace = FALSE,
                               sampling_type = c("random", "biased"),
                               bias = NULL, bias_effect = c("direct", "inverse"),
                               set_seed = 1) {

  # error checking


  # getting values for the three sets of data
  mvalues <- as.data.frame(calibration_variables)

  bvalues <- sample_background(mvalues, sample_size, replace, sampling_type,
                               bias, bias_effect, set_seed)

  xy <- c(longitude_column, latitude_column)
  ovalues <- terra::extract(calibration_variables,
                            as.matrix(occurrence_data[, xy]))

  # preparing histogram information
  var_hist <- lapply(1:ncol(mvalues), function(y) {
    m <- mvalues[, y]
    b <- bvalues$background[, y]
    o <- ovalues[, y]

    mbo <- c(m, b, o)
    ax <- pretty((min(mbo) - 0.1):max(mbo), n = breaks)
    df <- diff(ax[1:2])
    ax <- c(ax[1] - df, ax, ax[length(ax)] + df)

    new_explore_back(hist_m = hist(m, breaks = ax, plot = FALSE),
                     hist_b = hist(b, breaks = ax, plot = FALSE),
                     hist_o = hist(rep(o, magnify_occurrences),
                                   breaks = ax, plot = FALSE),
                     mean_m = mean(m), mean_b = mean(b), mean_o = mean(o),
                     cl_m = quantile(m, c(0.05, 0.95)),
                     cl_b = quantile(b, c(0.05, 0.95)),
                     cl_o = quantile(o, c(0.05, 0.95)), range_m = range(m),
                     range_b = range(b), range_o = range(o))
  })

  names(var_hist) <- colnames(mvalues)

  summa <- data.frame(n_all = nrow(mvalues),
                      n_background = nrow(bvalues$background),
                      n_occurrences = nrow(ovalues),
                      bvalues$summary)

  # results
  return(new_explore_list(summary = summa, exploration_stats = var_hist))
}







plot_explore <- function(exploration_list, variable = 1, col_all = "gray70",
                         col_background = "gray40", col_occurrences = "gray5",
                         alpha = 0.4, xlab = NULL, ylab = NULL, lines = TRUE,
                         which_lines = c("range", "cl", "mean"),
                         lty_all = 1, lty_background = 2, lty_occurrences = 3,
                         lwd_all = 3, lwd_background = 2, lwd_occurrences = 1) {

  # error checking
  if (missing(exploration_list)) {
    stop("Argument 'exploration_list' must be defined.")
  }
  if (class(exploration_list)[1] != "explore_list") {
    stop("'exploration_list' must be of class 'explore_list'.")
  }

  var_res <- exploration_list$exploration_stats[[variable]]

  if (class(var_res)[1] != "explore_back") {
    stop("Each element of 'exploration_stats' in 'exploration_list' must be\nof class 'explore_back'.")
  }

  # color
  mcol <- scales::alpha(col_all, alpha)
  bcol <- scales::alpha(col_background, alpha)
  ocol <- scales::alpha(col_occurrences, alpha)

  # labels
  vname <- ifelse(
    is.numeric(variable),
    gsub("_", " ", names(exploration_list$exploration_stats)[variable]),
    gsub("_", " ", varaible)
  )

  ylab <- ifelse(is.null(ylab), "Frequency", ylab)
  xlab <- ifelse(is.null(xlab), vname, xlab)

  # preparing data for plotting
  lins <- paste0(which_lines[1], "_", c("m", "b", "o"))

  # plot
  plot(var_res$hist_m, col = mcol, main = "", xlab = xlab, border = mcol,
       freq = TRUE, ylab = ylab)
  plot(var_res$hist_b, col = bcol, add = TRUE, border = ocol)
  plot(var_res$hist_o, col = ocol, add = TRUE, border = ocol)

  abline(v = var_res[[lins[1]]], col = col_all, lwd = lwd_all, lty = lty_all)
  abline(v = var_res[[lins[2]]], col = col_background, lwd = lwd_background,
         lty = lty_background)
  abline(v = var_res[[lins[3]]], col = col_occurrences, lwd = lwd_occurrences,
         lty = lty_occurrences)

  box(bty = "l")
}
