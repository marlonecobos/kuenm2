#### Helpers glm.net ####
# From maxnet '0.1.4'

#### Functions to fit hinge features ####
hinge <- function(x, nknots = 50) {
  min <- min(x)
  max <- max(x)
  k <- seq(min, max, length = nknots)
  lh <- outer(x, utils::head(k, -1), function(w, h) hingeval(w, h, max))
  rh <- outer(x, k[-1], function(w, h) hingeval(w, min, h))
  colnames(lh) <- paste("", utils::head(k, -1), max, sep = ":")
  colnames(rh) <- paste("", min, k[-1], sep = ":")
  cbind(lh, rh)
}

hingeval <- function(x, min, max) {
  pmin(1, pmax(0, (x - min) / (max - min)))
}

#### Functions to fit thresholds features ####
thresholds <- function(x, nknots = 50) {
  min <- min(x)
  max <- max(x)
  k <- seq(min, max, length = nknots + 2)[2:nknots + 1]
  f <- +(outer(x, k, ">="))
  colnames(f) <- paste("", k, sep = ":")
  f
}

thresholdval <- function(x, knot) {
  +(x >= knot)
}

#### Functions to fit categorical features ####
categorical <- function(x) {
  levs <- levels(x)
  f <- +(outer(x, levs, "=="))
  colnames(f) <- paste("", levs, sep = ":")
  f
}

categoricalval <- function(x, category) {
  +(x == category)
}

#### maxnet default regularization ####
maxnet.default.regularization <- function(p, m) {
  cnames <- colnames(m)
  isproduct <- grepl(":", cnames, fixed = TRUE) & !grepl("(", cnames, fixed = TRUE)
  isquadratic <- grepl("^I\\(.*\\^2\\)", cnames)
  ishinge <- grepl("^hinge\\(", cnames)
  isthreshold <- grepl("^thresholds\\(", cnames)
  iscategorical <- grepl("^categorical\\(", cnames)

  mm <- m[p == 1, , drop = FALSE]
  np <- nrow(mm)

  if (any(isproduct)) {
    lqpreg <- list(c(0, 10, 17, 30, 100), c(2.6, 1.6, 0.9, 0.55, 0.05))
  } else if (any(isquadratic)) {
    lqpreg <- list(c(0, 10, 17, 30, 100), c(1.3, 0.8, 0.5, 0.25, 0.05))
  } else {
    lqpreg <- list(c(0, 10, 30, 100), c(1, 1, 0.2, 0.05))
  }

  default_reg <- stats::approx(lqpreg[[1]], lqpreg[[2]], np, rule = 2)$y
  hinge_reg <- stats::approx(c(0, 1), c(0.5, 0.5), np, rule = 2)$y
  cat_reg <- stats::approx(c(0, 10, 17), c(0.65, 0.5, 0.25), np, rule = 2)$y
  thr_reg <- stats::approx(c(0, 100), c(2, 1), np, rule = 2)$y

  classregularization <- rep(default_reg, length(cnames))
  classregularization[ishinge] <- hinge_reg
  classregularization[iscategorical] <- cat_reg
  classregularization[isthreshold] <- thr_reg
  classregularization <- classregularization / sqrt(np)

  if (np > 1L) {
    mm_means <- colMeans(mm)
    mm_vars <- (colSums(mm^2) - np * (mm_means^2)) / (np - 1L)
    mm_sds <- sqrt(pmax(0, mm_vars))
  } else {
    mm_sds <- rep(0, ncol(mm))
  }

  hmindev <- ifelse(ishinge, pmax(mm_sds, 1 / sqrt(np)) * 0.5 / sqrt(np), 0)
  col_sums_mm <- colSums(mm)
  tmindev <- ifelse(isthreshold & (col_sums_mm == 0 | col_sums_mm == np), 1, 0)

  m_min <- apply(m, 2L, min)
  m_max <- apply(m, 2L, max)

  pmax(0.001 * (m_max - m_min), hmindev, tmindev, mm_sds * classregularization)
}

# Helper to predict glms
predict_glm_mx <- function(model, newdata, type = "response") {
  if (type == "response") {
    return(stats::predict.glm(model, newdata, type))
  }

  if (type == "cloglog") {
    exp_p <- stats::predict.glm(model, newdata, type = "link")
    return(1 - exp(-exp(exp_p)))
  }

  if (type == "raw") {
    return(exp(stats::predict.glm(model, newdata, type = "link")))
  }

  if (type == "cumulative") {
    return(cumulative_predictions(
      exp(stats::predict.glm(model, newdata, type = "link"))))
  }

  if (type == "link") {
    return(stats::predict.glm(model, newdata, type = "link"))
  }
}

# Extract variable names from coefs
extract_var_names <- function(terms) {
  clean_terms <- gsub("(:[0-9.-]+)+$", "", terms)
  clean_terms <- gsub("^(hinge|thresholds|categorical|I)\\(([^\\)]+)\\)(.*)$", "\\2", clean_terms)
  clean_terms <- gsub("\\^2$", "", clean_terms)
  vars <- unlist(strsplit(clean_terms, ":"))
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
