#' Maxent-like glmnet models
#'
#' @description
#' This function fits Maxent-like models using the \code{glmnet} package, designed
#' for presence-background data.
#'
#' @usage
#' glmnet_mx(p, data, f, regmult = 1.0, regfun = maxnet.default.regularization,
#'           addsamplestobackground = TRUE, weights = NULL, only_fit = FALSE,
#'           precomp = NULL, ...)
#'
#' @param p A vector of binary presence-background labels, where 1 indicates
#' presence and 0 indicates background.
#' @param data A \code{data.frame} containing the predictor variables for the
#' model. This must include the same number of rows as the length of \code{p}.
#' @param f A formula specifying the model to be fitted, in the format used by
#' \code{\link[stats]{model.matrix}}.
#' @param regmult (numeric) Regularization multiplier, default is 1.0.
#' @param regfun A function that calculates regularization penalties. Default is
#' \code{maxnet.default.regularization}.
#' @param addsamplestobackground (logical) Whether to add presence points not in
#' the background to the background data. Default is \code{TRUE}.
#' @param weights (numeric) A numeric vector of weights for each observation.
#' Default is \code{NULL}, which sets weights to 1 for presence points
#' and 100 for background points.
#' @param only_fit (logical) Whether to skip computing auxiliary metadata
#' (`entropy`, `featuremins`, `featuremaxs`, `varmin`, `varmax`, `samplemeans`,
#' `levels`) that are not required during cross-validation fold evaluation.
#' Default is \code{FALSE}.
#' @param precomp (list) Optional precomputed design list returned by
#' `prep_design()`, avoiding redundant `model.matrix()` and regularization
#' calculations across multiple regularization multipliers. Default is \code{NULL}.
#' @param ... Additional arguments to pass to \code{\link[glmnet]{glmnet}}.
#'
#' @return
#' A fitted Maxent-like model object of class \code{glmnet_mx}, which
#' includes model coefficients, AIC (if requested), and other elements
#' such as feature mins and maxes, sample means, and entropy.
#'
#' @details
#' This function is modified from the package maxnet and fits a Maxent-like
#' model using regularization to avoid over-fitting. Regularization weights
#' are computed using a provided function (which can be changed) and can be
#' multiplied by a regularization multiplier (\code{regmult}).
#'
#' @importFrom glmnet glmnet.control glmnet
#'
#' @export

glmnet_mx <- function(p,
                      data,
                      f,
                      regmult = 1.0,
                      regfun = maxnet.default.regularization,
                      addsamplestobackground = TRUE,
                      weights = NULL,
                      only_fit = FALSE,
                      precomp = NULL,
                      ...) {

  if (is.null(precomp)) {
    if (missing(p)) {
      stop("Argument 'p' must be defined.")
    }
    if (missing(data)) {
      stop("Argument 'data' must be defined.")
    }
    if (missing(f)) {
      stop("Argument 'f' must be defined.")
    }
    if (!inherits(p, "numeric")) {
      stop("Argument 'p' must be a 'numeric' vector.")
    }
    if (anyNA(data)) {
      stop("NA values in data table. Please remove them and rerun.")
    }

    iniweigth <- is.null(weights)
    if (iniweigth) {
      weights <- ifelse(p == 1, 1, 100)
    }

    if (addsamplestobackground) {
      is_p <- (p == 1)
      pdata <- data[is_p, , drop = FALSE]
      ndata <- data[!is_p, , drop = FALSE]

      wadd <- !(do.call(paste, pdata) %in% do.call(paste, ndata))
      if (any(wadd)) {
        n_add <- sum(wadd)
        p <- c(p, rep(0, n_add))
        data <- rbind(data, pdata[wadd, , drop = FALSE])

        if (iniweigth) {
          weights <- c(weights, rep(100, n_add))
        } else {
          pweight <- weights[is_p]
          weights <- c(weights, pweight[wadd])
        }
      }
    }

    mm <- stats::model.matrix(f, data)
    if (ncol(mm) < 2L) {
      mm <- cbind(0, mm)
    }

    reg <- regfun(p, mm) * regmult
    y_fac <- structure(as.integer(p) + 1L, levels = c("0", "1"), class = "factor")
  } else {
    p <- precomp$p
    y_fac <- precomp$y_fac
    data <- precomp$data
    weights <- precomp$weights
    mm <- precomp$mm
    reg <- precomp$reg * regmult
  }

  lambdas <- 10^(seq(4, 0, length.out = 200)) * sum(reg) / length(reg) *
    sum(p) / sum(weights)

  glmnet::glmnet.control(pmin = 1.0e-8, fdev = 0)
  model <- glmnet::glmnet(x = mm, y = y_fac, family = "binomial",
                          standardize = FALSE, penalty.factor = reg,
                          lambda = lambdas, weights = weights, ...)

  class(model) <- c("glmnet_mx", class(model))
  if (length(model$lambda) < 200L) {
    msg <- "glmnet failed to complete regularization path. Model may be infeasible."
    if (!addsamplestobackground) {
      msg <- paste(msg, "\n\tTry re-running with 'addsamplestobackground' = TRUE.")
    }
    stop(msg)
  }

  b_last <- model$beta[, 200L]
  filter <- (b_last != 0)
  betas <- b_last[filter]
  model$betas <- betas

  bg_rows <- (p == 0)
  if (length(betas) > 0L) {
    rr <- exp(as.numeric(mm[bg_rows, names(betas), drop = FALSE] %*% betas))
  } else {
    rr <- rep(1, sum(bg_rows))
  }

  rr <- rr + .Machine$double.eps
  model$alpha <- -log(sum(rr))

  if (only_fit) {
    return(model)
  }

  model$AIC <- NA
  raw <- rr / sum(rr)
  model$entropy <- -sum(raw * log(raw))
  model$penalty.factor <- reg

  model$featuremins <- apply(mm, 2L, min)
  model$featuremaxs <- apply(mm, 2L, max)

  vv <- (sapply(data, class) != "factor")
  model$varmin <- apply(data[, vv, drop = FALSE], 2L, min)
  model$varmax <- apply(data[, vv, drop = FALSE], 2L, max)

  means <- colMeans(data[p == 1, vv, drop = FALSE])
  majorities <- sapply(names(data)[!vv], function(n) {
    which.max(table(data[p == 1, n, drop = FALSE]))
  })
  names(majorities) <- names(data)[!vv]
  model$samplemeans <- unlist(c(means, majorities))

  model$levels <- lapply(data, levels)

  return(model)
}
