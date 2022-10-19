new_explore_back <- function(hist_m, hist_b, hist_o, mean_m, mean_b, mean_o,
                             cl_m, cl_b, cl_o, range_m, range_b, range_o) {
  # error checking

  val <- list(hist_m = hist_m, hist_b = hist_b, hist_o = hist_o,
              mean_m = mean_m, mean_b = mean_b, mean_o = mean_o,
              cl_m = cl_m, cl_b = cl_b, cl_o = cl_o, range_m = range_m,
              range_b = range_b, range_o = range_o)

  class(val) <- "explore_back"
  return(val)
}


new_explore_list <- function(summary, exploration_stats) {
  # error checking

  val <- list(summary = summary, exploration_stats = exploration_stats)

  class(val) <- "explore_list"
  return(val)
}


new_occ_cal <- function(summary, occurrences, block_limits = NULL,
                        calibration_regime = NULL) {
  # error checking

  val <- list(summary = summary, occurrences = occurrences,
              block_limits = block_limits,
              calibration_regime = calibration_regime)
  class(val) <- "occ_cal"
  return(val)
}



new_back_cal <- function(summary, background, background_bias = NULL,
                         training_regime = NULL) {
  # error checking

  val <- list(summary = summary, background = background,
              background_bias = background_bias,
              training_regime = training_regime)
  class(val) <- "back_cal"
  return(val)
}
