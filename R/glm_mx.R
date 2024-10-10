#' Maxent-like GLM models
#'
#' for more details see Fithian and Hastie (2013) DOI: 10.1214/13-AOAS667

#'
#' @importFrom stats glm
#' @export

glm_mx <- function(formula, family = binomial(link = "cloglog"), data,
                   weights = NULL, ...) {

  # Check for missing values in the data
  if (anyNA(data)) stop("NA values in data. Please remove them and rerun.")

  if (!is.data.frame(data)){
    stop("data must be a data.frame")
  }

  # Initialize weights: if weights is NULL, assign default weights based on pr_bg
  if (is.null(weights)) {
    weights <- ifelse(data$pr_bg == 1, 1, 10000)
  }
  # Ensure that the weights are a numeric vector
  if (!is.numeric(weights)) {
    stop("'weights' must be a numeric vector.")
  }

  # print(length(weights))
  # print(data)

  # Fit the model with the provided or default weights
  model <- suppressWarnings(
    glm(formula = formula, family = family, data = data,
        weights = NULL, ...))

  return(model)
}
