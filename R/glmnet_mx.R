#' Maxent-like glmnet models
#'
#' @description
#' This function fits Maxent-like models using the \code{glmnet} package, designed
#' for presence-background data. It includes options for regularization and
#' calculating the Akaike Information Criterion (AIC).
#' The model can automatically add presence points to the background if needed.
#'
#' @usage
#' glmnet_mx(p, data, f, regmult = 1.0, regfun = maxnet.default.regularization,
#'           addsamplestobackground = TRUE, weights = NULL,
#'           calculate_AIC = FALSE, AIC_option = "ws", ...)
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
#' @param calculate_AIC (logical) Whether to calculate AIC. Default is \code{FALSE}.
#' @param AIC_option (character) Method for calculating AIC, either "nk" or "ws".
#' Default is "ws".
#' @param ... Additional arguments to pass to \code{\link[glmnet]{glmnet}}.
#'
#' @return
#' A fitted Maxent-like model object of class \code{glmnet_mx}, which
#' includes model coefficients, AIC (if requested), and other elements
#' such as feature mins and maxes, sample means, and entropy.
#'
#' @details
#' This function is modified from the package maxnet and fits a Maxent-like
#' model using regularization to avoid overfitting. Regularization weights
#' are computed using a user-provided function and can be multiplied by
#' a regularization multiplier (\code{regmult}). The function also
#' includes an option to calculate AIC.
#'
#' @importFrom glmnet glmnet.control glmnet
#' @export


glmnet_mx <- function(p,
                      data,
                      f,
                      regmult = 1.0,
                      regfun = maxnet.default.regularization,
                      addsamplestobackground = TRUE,
                      weights = NULL,
                      calculate_AIC = FALSE,
                      AIC_option = "ws",
                      ...) {
  if (anyNA(data)) {
    stop("NA values in data table. Please remove them and rerun.")
  }
  if (!is.vector(p)) {
    stop("p must be a vector.")
  }
  iniweigth <- is.null(weights)

  if (iniweigth) {
    weights <- ifelse(p == 1, 1, 100)
  }

  if (addsamplestobackground) {
    pdata <- data[p == 1, ]
    ndata <- data[p == 0, ]

    # add to background any presence data that isn't there already
    wadd <- !do.call(paste, pdata) %in% do.call(paste, ndata)
    if (sum(wadd) > 0) {
      p <- c(p, rep(0, sum(wadd)))
      data <- rbind(data, pdata[wadd, ])

      ## adding extra weights if required
      if (iniweigth) {
        weights <- c(weights, rep(100, sum(wadd)))
      } else {
        pweight <- weights[p == 1]
        weights <- c(weights, pweight[wadd])
      }
    }
  }

  mm <- model.matrix(f, data)

  reg <- regfun(p, mm) * regmult
  lambdas <- 10^(seq(4, 0, length.out = 200)) * sum(reg) / length(reg) *
    sum(p) / sum(weights)

  glmnet::glmnet.control(pmin = 1.0e-8, fdev = 0)
  model <- glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",
                          standardize = FALSE, penalty.factor = reg,
                          lambda = lambdas, weights = weights, ...)

  bb <- coef(model)[, 200]

  class(model) <- c("glmnet_mx", class(model))
  if (length(model$lambda) < 200) {
    msg <- "glmnet failed to complete regularization path. Model may be infeasible."
    if (!addsamplestobackground) {
      msg <- paste(msg, "\n\tTry re-running with 'addsamplestobackground' = TRUE.")
    }
    stop(msg)
  }

  # AIC calculation
  filter <- bb[-1] != 0
  bb <- c(bb[1], bb[-1][filter])

  if (calculate_AIC & sum(filter) != 0 & AIC_option == "nk") {
    model$AIC <- aic_nk(x = as.matrix(mm[, filter]), y = p, beta = bb)
  } else {
    model$AIC <- NA
  }

  # returning other elements
  model$betas <- bb[-1]
  model$alpha <- 0
  rr <- predict.glmnet_mx(model, data[p == 0, , drop = FALSE],
                          type = "exponent", clamp = FALSE)
  #.Machine$double.eps added to deal with rr = 0
  rr <- rr + .Machine$double.eps
  raw <- rr / sum(rr)
  model$entropy <- -sum(raw * log(raw))
  model$alpha <- -log(sum(rr))

  model$penalty.factor <- reg
  model$featuremins <- apply(mm, 2, min)
  model$featuremaxs <- apply(mm, 2, max)

  vv <- (sapply(data, class) != "factor")
  model$varmin <- apply(data[, vv, drop = FALSE], 2, min)
  model$varmax <- apply(data[, vv, drop = FALSE], 2, max)

  means <- apply(data[p == 1, vv, drop = FALSE], 2, mean)
  majorities <- sapply(names(data)[!vv], function(n) {
    which.max(table(data[p == 1, n, drop = FALSE]))
  })
  names(majorities) <- names(data)[!vv]
  model$samplemeans <- unlist(c(means, majorities))

  model$levels <- lapply(data, levels)

  return(model)
}

