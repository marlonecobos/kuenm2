# background preparation helper functions

# helper to sample background points
sample_background <- function(initial_points, sample_size = NULL, replace = FALSE,
                              sampling_type = c("random", "biased"),
                              bias = NULL, bias_effect = c("direct", "inverse"),
                              set_seed = 1) {
  # error checking
  sampling_type <- sampling_type[1]
  bias_effect <- bias_effect[1]



  # preparations
  ipoints <- nrow(initial_points)

  # sampling
  if (is.null(sample_size)) {
    ## if NULL size
    message("'sample_size' is NULL, using all points")

    summa <- data.frame(initial_n = ipoints, sample_size = NA, replace = NA,
                        sampling_type = "none", bias_effect = NA,
                        rare_percent_insample = NA, rare_threshold = NA,
                        seed = NA)

    return(list(summary = summa, background = initial_points, bias = NULL))

  } else {
    if (sample_size >= ipoints) {
      ## if size equal to all points
      message("'sample_size' >= initial number of points, using all points")

      summa <- data.frame(initial_n = ipoints, sample_size = NA, replace = NA,
                          sampling_type = "none", bias_effect = NA,
                          rare_percent_insample = NA, rare_threshold = NA,
                          seed = NA)

      return(list(summary = summa, background = initial_points, bias = NULL))

    } else {
      if (sampling_type == "random") {
        ## random sample
        set.seed(set_seed)
        initial_points <- initial_points[sample(ipoints, size = sample_size,
                                                replace = replace), ]

        bias_effect <- NA
      }

      if (sampling_type == "biased") {
        ## biased sample
        set.seed(set_seed)

        if (bias_effect == "inverse") {
          initial_points <- initial_points[sample(ipoints, size = sample_size,
                                                  replace = replace,
                                                  prob = max(bias) - bias), ]
        } else {
          initial_points <- initial_points[sample(ipoints, size = sample_size,
                                                  replace = replace,
                                                  prob = bias), ]
        }
      }

      summa <- data.frame(initial_n = ipoints, sample_size = sample_size,
                          replace = replace, sampling_type = sampling_type,
                          bias_effect = bias_effect, seed = set_seed)

      return(list(summary = summa, background = initial_points, bias = bias))
    }
  }
}




