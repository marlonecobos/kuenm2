#' Analysis of extrapolation risks using the MOP metric
#'
#' @description mop calculates a the mobility-oriented parity metric and other
#' sub-products to represent conditions in a projection area that are
#' non-analogous to those in a calibration area.
#'
#' @param m SpatRaster or matrix of variables representing the original area
#' of interest (where a model was calibrated). If a matrix is used, each column
#' represents a variable.
#' @param g a SpatRaster or matrix of variables representing the area or
#' scenario where dissimilarity and non-analogous conditions will be detected
#' (where a model is projected). Variable names must match between \code{m} and
#' \code{g}.
#' @param mop_type (character) type of MOP analyses to be performed. Options
#' are: "basic", "simple", and "detailed".
#' @param distance (character) one of the options: "euclidean" or "mahalanobis".
#' Default = "euclidean".
#' @param scale scaling options (logical or numeric-alike) as in
#' \code{\link{scale}}. Default = FALSE.
#' @param center (logic or numeric-alike) center options as in
#' \code{\link{scale}}. Default = FALSE.
#' @param fix_NA (logical) whether to fix the layers so all cells with NA values
#' coincide among all layers.
#' @param percent (numeric) percentage of points of m (the closest ones) used to
#' derive mean environmental distances to each g point.
#' @param comp_each (numeric) number of points of the g matrix to be used for
#' distance calculations at a time (default = 2000). Increasing this number
#' requires more RAM.
#' @param rescale_mop (logical) whether or not to re-scale distances 0-1.
#' Default = TRUE. Re-scaling complicates comparisons of dissimilarity values
#' obtained from exercises in which the percentage is changed.
#' @param parallel (logical) if TRUE, calculations will be performed in parallel
#' using \code{n_cores} of the computer. Using this option will speed up the
#' analysis  but will demand more RAM. Default = FALSE.
#' @param n_cores (numeric) number of cores to be used in parallel processing.
#' Default = NULL, in which case all CPU cores on current host - 1 will be used.
#'
#' @return **change to include matrix results**
#' A list of results containing:
#' - summary.- a list with details on the data used in the analysis.
#' - mop_basic.- a SpatRaster for the projection area (\code{g}) with
#' dissimilarity values (0-1), where 1 means that values of at least one of the
#' variables in the projection area are non-analogous to those in the
#' calibration area.
#' - mop_simple.- a SpatRaster for the projection area (\code{g}) with values
#' representing how many variables in the projection area are non-analogous to
#' those in the calibration area.
#' - mop_detailed.- a list containing:
#'     - interpretation_combined.- a data.frame to help identify the variables
#'     that are non-analogous to \code{m} in the layers towards_low_combined and
#'     towards_high_combined.
#'     - towards_low_end.- a SpatRaster or matrix with layers for all variables
#'     representing where non-analogous conditions were found towards low values
#'     of the variable.
#'     - towards_high_end.- a SpatRaster or matrix with layers for all variables
#'     representing where non-analogous conditions were found towards high values
#'     of the variable.
#'     - towards_low_combined.- a SpatRaster or vector with values representing
#'     the identity of the variables found to have non-analogous conditions
#'     towards low values. Interpretation requires the use of the table in
#'     interpretation_combined.
#'     - towards_high_combined.- a SpatRaster or vector with values representing
#'     the identity of the variables found to have non-analogous conditions
#'     towards high values. Interpretation requires the use of the table in
#'     interpretation_combined.
#'
#' @details
#' The options for the argument \code{mop_type} return results that differ in
#' the detail of how non-analogous conditions are identified. The option "basic"
#' makes calculation as proposed by Owens et al. (2013;
#' \url{https://doi.org/10.1016/j.ecolmodel.2013.04.011}). Other options perform
#' further analyses that help to identify non-analogous conditions in more
#' detail (see description of results returned).
#'
#' When the variables used to represent environmental conditions in the areas of
#' interest differ considerably in the way they were measured, using the
#' arguments \code{scale} and \code{center} is recommended.
#'
#' @usage
#' mop(m, g, mop_type = "basic", distance = "euclidean", scale = FALSE,
#'     center = FALSE, fix_NA = FALSE, percent = 5, comp_each = 2000,
#'     parallel = FALSE, n_cores = NULL)

#if (!require(raster)) {
#  install.packages("raster")
#}
#
#if (!require(foreach)) {
#  install.packages("foreach")
#}
#
#if (!require(Kendall)) {
#  install.packages("Kendall")
#}
#
#if (!require(snow)) {
#  install.packages("snow")
#}
#
#if (!require(doSNOW)) {
#  install.packages("doSNOW")
#}

mop <- function(m, g, mop_type = "basic", distance = "euclidean",
                scale = FALSE, center = FALSE, fix_NA = FALSE, percent = 5,
                comp_each = 2000, rescale_mop = TRUE, parallel = FALSE,
                n_cores = NULL) {

  # initial tests
  if (missing(m) | missing(g)) {
    stop("Arguments 'm' and 'g' must be defined.")
  }

  clasm <- class(m)[1]
  clasg <- class(g)[1]

  if (!clasm %in% c("SpatRaster", "matrix", "data.frame")) {
    stop("'m' must be a 'SpatRaster', 'matrix', or 'data.frame'.")
  } else {
    if (clasm == "SpatRaster") {
      var_names <- terra::varnames(m)
    } else {
      var_names <- colnames(m)
    }
  }
  if (!clasg %in% c("SpatRaster", "matrix")) {
    stop("'g' must be a 'SpatRaster', 'matrix', or 'data.frame'.")
  } else {
    if (clasg == "SpatRaster") {
      gnames <- terra::varnames(g)
    } else {
      gnames <- colnames(g)
    }
  }

  if (!identical(var_names, gnames)) {
    stop("Variables in 'm' and 'g' must be named identically.")
  }

  mop_type <- mop_type[1]


  # preprocessing
  ## layer for mop results
  if (clasg == "SpatRaster") {
    ### fix NA mismatches across layers
    if (fix_NA == TRUE) {
      g <- fix_na_rast(raster_layers = g)
    }

    ### layer
    mop <- g[[1]]
    names(mop) <- "mop"
  }

  ## getting values
  m <- m[]
  g <- g[]

  ## variables and test
  nvar <- ncol(m)

  ## exclude NAs and number of cells
  m <- na.omit(m)
  g <- na.omit(g)

  nm <- nrow(m)
  ng <- nrow(g)

  if (clasg == "matrix") {
    mop <- rep(NA, 1:nrow(g))
  }

  ## scaling options
  if (distance == "euclidean") {
    if (scale == TRUE | center == TRUE) {
      m <- scale(rbind(m, g), center = center, scale = scale)
      g <- m[(nm + 1):nrow(m), ]
      m <- m[1:nm, ]
    }
  }


  # identifying dimensions where g is out of m box
  if (!mop_type %in% c("basic", "simple", "detailed")) {
    stop("Option for 'mop_type' not valid")
  } else {
    out_ranges <- out_m(m_matrix = m, g_matrix = g, type = mop_type)
  }


  # distance calculation
  ## only relevant points
  reduced <- out_ranges$general == 0
  g_red <- g[reduced, ]

  ## analysis
  mop1 <- mop_distance(m, g_red, distance, percent, comp_each,
                       parallel, n_cores)


  # re-assigning values to rasters
  ## na values in mop layer
  if (clasg == "SpatRaster") {
    nona <- which(!is.na(mop[]))
  } else {
    nona <- rep(TRUE, length(mop))
  }

  ## simple results
  if (mop_type != "basic") {
    mop2 <- mop
    out_ranges$general[reduced] <- NA
    mop2[nona] <- out_ranges$general
  } else {
    mop2 <- NULL
  }

  ## detailed results
  if (mop_type == "detailed") {
    mop3 <- mop
    mop4 <- mop

    ### results for combined variables
    mop3[nona] <- out_ranges$low_combined
    mop4[nona] <- out_ranges$high_combined

    ### results for independent variables
    if (clasg == "SpatRaster") {
      mop5 <- terra::rast(lapply(1:nvar, function(x) {
        mop0 <- mop
        mop0[nona] <- out_ranges$low_all[, x]
        mop0
      }))
      terra::varnames(mop5) <- var_names
    } else {
      mop5 <- out_ranges$low_all
    }

    if (clasg == "SpatRaster") {
      mop6 <- terra::rast(lapply(1:nvar, function(x) {
        mop0 <- mop
        mop0[nona] <- out_ranges$high_all[, x]
        mop0
      }))
      terra::varnames(mop6) <- var_names
    } else {
      mop6 <- out_ranges$high_all
    }

  } else {
    mop3 <- NULL
    mop4 <- NULL
    mop5 <- NULL
    mop6 <- NULL
  }

  ## MOP results
  maxmop <- max(mop1)

  ### re-scaling if needed
  if (rescale_mop == TRUE) {
    out_ranges$general[reduced] <- mop1 / (maxmop * 1.1)
    out_ranges$general[!reduced] <- 1
  } else {
    out_ranges$general[reduced] <- mop1
    out_ranges$general[!reduced] <- maxmop * 1.1
  }

  mop[nona] <- out_ranges$general


  # returning results
  results <- list(
    summary = list(
      variables = var_names, percent = percent, mop_type = mop_type,
      fix_NA = fix_NA, cells_m = nm, cells_g = ng
    ),
    mop_basic = mop, mop_simple = mop2,
    mop_detailed = list(
      interpretation_combined = out_ranges$interpretation,
      towards_low_end = mop5,
      towards_high_end = mop6,
      towards_low_combined = mop3,
      towards_high_combined = mop4
    )
  )

  return(results)
}




# helper to detect values outside range of m conditions
out_m <- function(m_matrix, g_matrix, type = c("basic", "simple", "detailed")) {

  type <- type[1]
  var_names <- colnames(m_matrix)
  nvar <- length(var_names)

  if (!type %in% c("basic", "simple", "detailed")) {
    stop("Option for 'type' not valid")
  } else {
    ### defining range of what is inside M realms
    m_range <- apply(m_matrix, 2, range)

    ### what is out
    out <- sapply(1:nvar, function(x) {
      gx <- g_matrix[, x]
      gx < m_range[1, x] | gx > m_range[2, x]
    })

    out <- rowSums(out)

    ### what is out in more detail
    if (type == "detailed") {
      mul <- 10^(1:nvar)

      ### which lower end
      outl <- sapply(1:nvar, function(x) {
        (g_matrix[, x] < m_range[1, x]) * mul[x]
      })
      colnames(outl) <- var_names

      outl1 <- rowSums(outl)
      outl[outl > 0] <- 1

      ## which high end
      outh <- sapply(1:nvar, function(x) {
        (g_matrix[, x] > m_range[2, x]) * mul[x]
      })
      colnames(outh) <- var_names

      outh1 <- rowSums(outh)
      outh[outh > 0] <- 1

      ## fixing NAs
      outh1[outh1 == 0] <- NA
      outl1[outl1 == 0] <- NA
      outh[outh == 0] <- NA
      outl[outl == 0] <- NA

      ## interpretation table
      inter_table <- ext_interpret(var_names, var_codes = mul)
    } else {
      outh <- NULL
      outl <- NULL
      outh1 <- NULL
      outl1 <- NULL
      inter_table <- NULL
    }
  }

  return(list(general = out, high_all = outh, low_all = outl,
              high_combined = outh1, low_combined = outl1,
              interpretation = inter_table))
}



# helper function to create table for mop detailed interpretation
ext_interpret <- function (var_names, var_codes) {
  var_comb <- lapply(1:length(var_names), function(x) {
    apply(combn(var_names, m = x), 2, paste, collapse = ", ")
  })
  var_comb <- unlist(var_comb)

  var_cod <- lapply(1:length(var_codes), function(x) {
    apply(combn(var_codes, m = x), 2, sum)
  })
  var_cod <- unlist(var_cod)

  return(data.frame(values = var_cod, extrapolation_variables = var_comb))
}



# helper to calculate mop distances
mop_distance <- function(m_matrix, g_matrix, distance = "euclidean", percent = 5,
                         comp_each = 2000, parallel = FALSE, n_cores = NULL) {
  ### preparing groups of points
  nred <- nrow(g_matrix)
  per_out <- round(nrow(m_matrix) * (percent / 100))

  if (nred <= comp_each) {
    comp_each <- ceiling(nred / 3)
  }

  groups <- c(seq(1, nred, comp_each),  nred + 1)
  nprocess <- (length(groups) - 1)

  ### running analysis
  if (parallel == FALSE) {
    pb <- txtProgressBar(min = 1, max = nprocess, style = 3)

    mop1 <- lapply(1:nprocess, function(x) {
      Sys.sleep(0.1)
      setTxtProgressBar(pb, x)

      #### defining sets and all distances
      seq_rdist <- groups[x]:(groups[x + 1] - 1)

      if (distance == "euclidean") {
        cdist <- fields::rdist(g_matrix[seq_rdist, ], m_matrix)
      } else {
        cv <- cov(g_matrix)
        cdist <- lapply(seq_rdist, function(y) {
          stats::mahalanobis(x = m_matrix, center = g_matrix[y, ], cov = cv)
        })
        cdist <- do.call(rbind, cdist)
      }

      #### getting mean of closer distances
      apply(cdist, 1, function(y) {mean(sort(y)[1:per_out])})
    })

    close(pb)
    mop1 <- unlist(mop1)

  } else {
    # packages used internally
    suppressPackageStartupMessages(library(Kendall))
    suppressPackageStartupMessages(library(foreach))

    ### parallel processing preparation
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
    }

    cl <- snow::makeSOCKcluster(n_cores)
    doSNOW::registerDoSNOW(cl)

    ### progress bar preparation
    pb <- txtProgressBar(min = 1, max = nprocess, style = 3)
    progress <- function(n) {
      setTxtProgressBar(pb, n)
    }
    opts <- list(progress = progress)

    ### parallel running
    mop1 <- foreach(i = 1:nprocess, .packages = "Kendall", .inorder = TRUE,
                    .options.snow = opts, .combine = "c") %dopar% {

                      #### defining sets and all distances
                      seq_rdist <- groups[i]:(groups[i + 1] - 1)

                      if (distance == "euclidean") {
                        cdist <- fields::rdist(g_matrix[seq_rdist, ], m_matrix)
                      } else {
                        cv <- cov(g_matrix)
                        cdist <- lapply(seq_rdist, function(y) {
                          stats::mahalanobis(x = m_matrix, center = g_matrix[y, ],
                                             cov = cv)
                        })
                        cdist <- do.call(rbind, cdist)
                      }

                      #### getting mean of closer distances
                      return(
                        apply(cdist, 1, function(y) {mean(sort(y)[1:per_out])})
                      )
                    }

    close(pb)
    snow::stopCluster(cl)
  }

  return(mop1)
}




# helper to get information to plot raster files and their legends properly
plot_raster_inf <- function(mop_result, result = "basic",
                            col_nonanalogous = "#000000",
                            color_palette = heat.colors, reverse = FALSE) {
  if (result %in% c("basic", "simple")) {
    result1 <- paste0("mop_", result)
    raster_layer <- mop_result[[result1]]
  } else {
    raster_layer <- mop_result$mop_detailed[[result]]
    mop_detail <- mop_result$mop_detailed$interpretation_combined
  }

  vals <- na.omit(raster_layer[])
  fvals <- as.factor(vals)
  fac <- levels(fvals)
  colsa <- color_palette(length(fac))

  if (result == "basic") {
    colsa[1] <- col_nonanalogous
    colsa <- rev(colsa)
  }

  if (reverse == TRUE) {
    colsa <- rev(colsa)
  }

  facn <- as.numeric(fac)
  breaksa <- c(0, facn)

  if (result %in% c("basic", "simple")) {
    leg_col <- colsa
    if (result == "basic") {
      leg_text <- c("Low", "High")
    } else {
      leg_text <- fac
    }
  } else {
    leg_col <- colsa
    leg_text <- mop_detail[mop_detail$values %in% facn, ]
    ovals <- leg_text$values
    colsa <- colsa[order(ovals)]
    leg_text <- leg_text[, 2]
  }

  return(list(col = colsa, breaks = breaksa, legend_col = leg_col,
              legend_text = leg_text, val_factor = fvals))
}
