####Helpers glm.net####
#From maxnet ‘0.1.4’

####Functions to fit hinge features####
hinge <- function(x, nknots=50) {
  min <- min(x)
  max <- max(x)
  k <- seq(min, max, length=nknots)
  lh <- outer(x, utils::head(k,-1), function(w,h) hingeval(w, h, max))
  rh <- outer(x, k[-1], function(w,h) hingeval(w, min, h))
  colnames(rh) <- paste("", min, k[-1], sep=":")
  cbind(lh, rh)
}

hingeval <- function (x, min, max) {
  pmin(1, pmax(0, (x - min)/(max - min)))
}

####Functions to fit thresholds features####
thresholds <- function (x, nknots = 50) {
  min <- min(x)
  max <- max(x)
  k <- seq(min, max, length = nknots + 2)[2:nknots + 1]
  f <- outer(x, k, function(w, t) ifelse(w >= t, 1, 0))
  colnames(f) <- paste("", k, sep = ":")
  f
}

thresholdval <- function (x, knot) {
  ifelse(x >= knot, 1, 0)
}

####Functions to fit categorical features####
categorical <- function (x) {
  f <- outer(x, levels(x), function(w, f) ifelse(w == f, 1,
                                                 0))
  colnames(f) <- paste("", levels(x), sep = ":")
  f
}
categoricalval <- function (x, category) {
  ifelse(x == category, 1, 0)
}

#### maxnet default regularization ####
maxnet.default.regularization <- function (p, m) {
  isproduct <- function(x) {grepl(":", x) & !grepl("\\(", x)}
  isquadratic <- function(x) {grepl("^I\\(.*\\^2\\)", x)}
  ishinge <- function(x) {grepl("^hinge\\(", x)}
  isthreshold <- function(x) {grepl("^thresholds\\(", x)}
  iscategorical <- function(x) {grepl("^categorical\\(", x)}

  regtable <- function(name, default) {
    if (ishinge(name)) {
      return(list(c(0, 1), c(0.5, 0.5)))
    }
    if (iscategorical(name)) {
      return(list(c(0, 10, 17), c(0.65, 0.5, 0.25)))
    }
    if (isthreshold(name)) {
      return(list(c(0, 100), c(2, 1)))
    }
    default
  }

  lregtable <- list(c(0, 10, 30, 100), c(1, 1, 0.2, 0.05))
  qregtable <- list(c(0, 10, 17, 30, 100), c(1.3, 0.8, 0.5, 0.25, 0.05))
  pregtable <- list(c(0, 10, 17, 30, 100), c(2.6, 1.6, 0.9, 0.55, 0.05))

  mm <- m[p == 1, , drop = FALSE]
  np <- nrow(mm)
  lqpreg <- lregtable

  if (sum(isquadratic(colnames(mm)))) {
    lqpreg <- qregtable
  }
  if (sum(isproduct(colnames(mm)))) {
    lqpreg <- pregtable
  }

  classregularization <- sapply(colnames(mm), function(n) {
    t <- regtable(n, lqpreg)
    approx(t[[1]], t[[2]], np, rule = 2)$y
  })
  classregularization <- classregularization / sqrt(np)


  ishinge <- grepl("^hinge\\(", colnames(mm))

  hmindev <- sapply(1:ncol(mm), function(i) {
    if (!ishinge[i]) {
      return(0)
    }
    avg <- mean(mm[, i])
    std <- max(sd(mm[, i]), 1 / sqrt(np))
    std * 0.5 / sqrt(np)
  })

  tmindev <- sapply(1:ncol(mm), function(i) {
    ifelse(isthreshold(colnames(mm)[i]) && (sum(mm[, i]) ==
                                              0 || sum(mm[, i]) == nrow(mm)), 1, 0)
  })

  pmax(0.001 * (apply(m, 2, max) - apply(m, 2, min)), hmindev,
       tmindev, apply(as.matrix(mm), 2, sd) * classregularization)
}


