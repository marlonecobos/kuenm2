####Omission Rate####
omrat <- function(threshold, pred_train, pred_test) {
  if (missing(threshold)) {
    stop("Argument 'threshold' must be defined.")
  }
  if (missing(pred_train)) {
    stop("Argument 'pred_train' must be defined.")
  }
  if (missing(pred_test)) {
    stop("Argument 'pred_test' must be defined.")
  }
  om_rate <- vector("numeric", length = length(threshold))
  for (i in 1:length(threshold)) {
    val <- ceiling(length(pred_train) * threshold[i] / 100) + 1
    omi_val_suit <- sort(pred_train)[val]
    om_rate[i] <- as.numeric(
      length(pred_test[pred_test < omi_val_suit]) / length(pred_test)
    )
  }
  names(om_rate) <- paste0("Omission_rate_at_", threshold)
  return(om_rate)
}



####AIC - Warren and Seifert ####
aic_ws <- function(pred_occs, ncoefs) {
  if (missing(pred_occs)) {
    stop("Argument 'pred_occs' must be defined.")
  }
  if (missing(ncoefs)) {
    stop("Argument 'ncoefs' must be defined.")
  }
  LL <- sum(log(pred_occs + .Machine$double.eps))
  AICc <- ((2 * ncoefs) - (2 * LL)) + (2 * ncoefs * (ncoefs + 1) /
                                         (length(pred_occs) - ncoefs - 1))

  return(AICc)
}


#### AIC - Ninomiya and Kawano ####
aic_nk <- function(x, y, beta) {
  if (missing(x)) {
    stop("Argument 'x' must be defined.")
  }
  if (missing(y)) {
    stop("Argument 'y' must be defined.")
  }
  if (missing(beta)) {
    stop("Argument 'beta' must be defined.")
  }
  if (!is.matrix(x)) {
    stop("'x' must be a 'matrix'.")
  }
  if (mode(x) != "numeric") {
    stop("'x' must be 'numeric'.")
  }

  x_glm <- glm(y ~ x, family = "binomial")
  predict_glm <- predict(x_glm)

  pr0 <- 1 - 1 / (1 + exp(predict_glm))
  pr <- 1 - 1 / (1 + exp(cbind(1, x) %*% beta))
  #.Machine$double.eps added to deal with pr = 0
  pr <- pr + .Machine$double.eps
  cc <- as.vector(beta)
  jc <- abs(cc) > 1e-05
  xj <- cbind(1, x)[, jc]
  Pi <- diag(as.vector(pr * (1 - pr)))
  Pi0 <- diag(as.vector(pr0 * (1 - pr0)))
  j22 <- t(xj) %*% Pi %*% xj
  i22 <- t(xj) %*% Pi0 %*% xj

  aic <- -2 * sum(y * log(pr) + (1 - y) * log(1 - pr)) +
    2 * sum(diag(solve(j22) %*% i22))

  return(aic)
}
