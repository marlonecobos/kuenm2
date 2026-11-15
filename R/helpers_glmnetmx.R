####Helpers glm.net####
#From maxnet ‘0.1.4’

####Functions to fit hinge features####
hinge <- function(x, nknots=50) {
  min <- min(x)
  max <- max(x)
  k <- seq(min, max, length=nknots)
  lh <- outer(x, utils::head(k,-1), function(w,h) hingeval(w, h, max))
  rh <- outer(x, k[-1], function(w,h) hingeval(w, min, h))
  colnames(lh) <- paste("", utils::head(k, -1), max, sep = ":") #Fix issue #118
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

# Helper to predict glms
predict_glm_mx <- function(model, newdata, type = "response"){

  if(type == "response"){
    return(stats::predict.glm(model, newdata, type))
  }

  if(type == "cloglog"){
    exp_p <- stats::predict.glm(model, newdata, type = "link")
    return(1 - exp(-exp(exp_p)))
  }

  if(type == "raw"){
    return(exp(stats::predict.glm(model, newdata, type = "link")))
  }

  if (type == "cumulative") {
    return(cumulative_predictions(
      exp(stats::predict.glm(model, newdata, type = "link"))))
  }

  if(type == "link"){
    return(stats::predict.glm(model, newdata, type = "link"))
  }

}

# Extract variable names from coefs
extract_var_names <- function(terms) {
  # Remove hinge/threshold numeric suffixes (e.g., :13.14:22.53)
  clean_terms <- gsub("(:[0-9.-]+)+$", "", terms)

  # Extract content inside function parentheses: func(content)
  # This matches "hinge(bio_1)", "thresholds(bio_2)", "categorical(soil)", "I(bio_1^2)"
  clean_terms <- gsub("^(hinge|thresholds|categorical|I)\\(([^\\)]+)\\)(.*)$", "\\2", clean_terms)

  # Aggressively remove the ^2 suffix if it wasn't caught by the step above
  clean_terms <- gsub("\\^2$", "", clean_terms)

  # Split by ':' (for product interactions like bio_1:bio_7)
  vars <- unlist(strsplit(clean_terms, ":"))

  # Return unique clean variable names
  return(unique(trimws(vars)))
}

# Helper for fast row-wise consensus on numeric matrices
calc_mat_consensus <- function(mat, cons_types) {
  out <- list()
  if ("median" %in% cons_types) {
    out$median <- if (ncol(mat) == 1L) mat[, 1L] else apply(mat, 1L, stats::median)
  }
  if ("mean" %in% cons_types) {
    out$mean <- rowMeans(mat)
  }
  if ("stdev" %in% cons_types) {
    out$stdev <- if (ncol(mat) == 1L) {
      rep(0, nrow(mat))
    } else {
      m <- rowMeans(mat)
      sqrt(rowSums((mat - m)^2) / (ncol(mat) - 1L))
    }
  }
  if ("range" %in% cons_types) {
    out$range <- if (ncol(mat) == 1L) {
      rep(0, nrow(mat))
    } else {
      do.call(pmax, as.data.frame(mat)) - do.call(pmin, as.data.frame(mat))
    }
  }
  return(out)
}

# Helper to map numeric vector back to SpatRaster or full vector
vec_to_output <- function(v, is_raster, valid_cells, template_rast = NULL,
                          n_cells = NULL, n_rows = NULL) {
  if (is_raster) {
    r_out <- terra::rast(template_rast)
    full_v <- rep(NA_real_, n_cells)
    full_v[valid_cells] <- v
    terra::values(r_out) <- full_v
    return(r_out)
  } else {
    full_v <- rep(NA_real_, n_rows)
    full_v[valid_cells] <- v
    return(full_v)
  }
}

# Helper to stack a named list of vectors into SpatRaster or data.frame
list_to_container <- function(v_list, is_raster, valid_cells,
                              template_rast = NULL, n_cells = NULL,
                              n_rows = NULL) {
  mapped <- lapply(
    X = v_list,
    FUN = vec_to_output,
    is_raster = is_raster,
    valid_cells = valid_cells,
    template_rast = template_rast,
    n_cells = n_cells,
    n_rows = n_rows
  )

  if (is_raster) {
    r_stack <- terra::rast(mapped)
    names(r_stack) <- names(v_list)
    return(r_stack)
  } else {
    df_out <- as.data.frame(mapped)
    colnames(df_out) <- names(v_list)
    return(df_out)
  }
}
