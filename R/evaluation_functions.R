####Omission Rate####
omrat <- function(threshold = 5, pred_train, pred_test) {
  om_rate <- vector("numeric", length = length(threshold))
  for (i in 1:length(threshold)) {
    val <- ceiling(length(pred_train) * threshold[i]/100) + 1
    omi_val_suit <- sort(pred_train)[val]
    om_rate[i] <- as.numeric(length(pred_test[pred_test <
                                                omi_val_suit])/length(pred_test))
  }
  names(om_rate) <- paste("Omission_rate_at_", threshold, sep = "")
  return(om_rate)
}
# t2 <- omrat_maxnet(threshold = 5,
#                    pred_train =  suit_val_cal,
#                    pred_test = suit_val_eval)


####AIC - Warren and Seifert ####
aic_ws <- function(pred_occs, ncoefs) {
  LL <- sum(log(pred_occs + .Machine$double.eps))
  AICc <- ((2 * ncoefs) - (2 * LL)) + (2 * ncoefs * (ncoefs +
                                                   1)/(length(pred_occs) - ncoefs - 1))

  return(AICc)
}

#### AIC - Ninomiya and Kawano ####
aic_nk <- function(x, y, beta) {
  if (!is.matrix(x)) {
    stop("x must be a matrix.")
  }
  if (mode(x) != "numeric") {
    stop("x must be numeric.")
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
