#' Maxent-like Generalized Linear Models (GLM)
#'
#' @description
#' This function fits a Generalized Linear Model (GLM) to binary
#' presence-background data. It allows for the specification of custom weights,
#' with a default in which presences have a weight of 1 and background 100.
#'
#' @usage
#' glm_mx(formula, family = binomial(link = "cloglog"), data,
#'        weights = NULL, ...)
#'
#' @param formula A formula specifying the model to be fitted, in the format
#' used by \code{\link[stats]{glm}}.
#' @param family A description of the error distribution and link function to be
#' used in the model. Defaults to \code{binomial(link = "cloglog")},
#' which is commonly used for presence-background data.
#' @param data A \code{data.frame} containing the variables in the model. Must
#' include a column named \code{pr_bg} that indicates whether a record is
#' a presence (1) or background (0), and at least another column with an
#' independent variable (predictor).
#' @param weights Optional. A numeric vector of weights for each observation. If
#' not provided, default weights of 1 for presences and 100 for background
#' are used.
#' @param ... Additional arguments to be passed to \code{\link[stats]{glm}}.
#'
#' @details
#' For more details about glms using presence and background emulating what
#' Maxent does, see Fithian and Hastie (2013) DOI: 10.1214/13-AOAS667.
#'
#'
#' @return
#' A fitted \code{\link[stats]{glm}} object. The model object includes
#' the minimum and maximum values of the non-factor variables in the
#' dataset, stored as \code{model$varmin} and \code{model$varmax}.
#'
#' @importFrom stats glm
#' @export


glm_mx <- function(formula,
                   family = binomial(link = "cloglog"),
                   data,
                   weights = NULL,
                   ...) {

  # Check for missing values in the data
  if (anyNA(data)) stop("NA values in data. Please remove them and rerun.")

  if (!is.data.frame(data)){
    stop("data must be a data.frame")
  }

  # Initialize weights: if weights is NULL, assign default weights based on pr_bg
  if (is.null(weights)) {
    weights <- ifelse(data$pr_bg == 1, 1, 100)
  }
  # Ensure that the weights are a numeric vector
  if (!is.numeric(weights)) {
    stop("'weights' must be a numeric vector.")
  }

  # Fit the model with the provided or default weights
  model <- suppressWarnings(
    glm(formula = formula, family = family, data = data,
        weights = NULL, ...))

  vv <- (sapply(data, class) != "factor")
  model$varmin <- apply(data[, vv, drop = FALSE], 2, min)
  model$varmax <- apply(data[, vv, drop = FALSE], 2, max)

  return(model)
}
