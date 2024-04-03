#### New functions for maxnet ####
#Modified from maxnet ‘0.1.4’
glmnet_mx <- function(p, data, f, regmult = 1.0,
                      regfun = maxnet.default.regularization,
                      addsamplestobackground = TRUE,
                      weights = NULL,
                      calculate_AIC = FALSE, ...) {
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

  if (calculate_AIC & sum(filter) != 0) {
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

