#' Maxent-like GLM models
#'
#' for more details see Fithian and Hastie (2013) DOI: 10.1214/13-AOAS667

#'
#' @importFrom stats glm
#' @export

glm_mx <- function(formula, family = binomial(link = "cloglog"), data = data,
                   weights = NULL, ...) {

  if(anyNA(data)) stop("NA values in data. Please remove them and rerun.")

  if (is.null(weights))
    weights <- ifelse(data$calibration_data$pr_bg == 1, 1, 10000)

  model <- suppressWarnings(
    glm(formula = formula, family = family, data = data, weights = weights))

  return(model)
}

