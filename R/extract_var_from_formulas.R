#' Extract predictor names from formulas
#'
#' @usage
#' extract_var_from_formulas(formulas, ...)
#'
#' @param formulas (character or formula) model formulas.
#' @param ... Arguments to pass to [all.vars()]
#'
#' @return
#' A character vector or a list of the same length as `formulas`,
#' containing the names of the predictors each formula.
#'
#' @export
#'
#' @importFrom stats terms as.formula
#'
#' @examples
#' # Import an example of calibration results
#' data(calib_results_maxnet, package = "kuenm2")
#'
#' # Extract predictor names
#' vars <- extract_var_from_formulas(calib_results_maxnet$formula_grid$Formulas)

extract_var_from_formulas <- function(formulas, ...) {
  if (missing(formulas)) {
    stop("Argument 'formulas' must be defined.")
  }

  res <- sapply(formulas, function(x) {
    if (inherits(x, "formula")) {
      x <- stats::terms(x)
    } else {
      x <- stats::terms(stats::as.formula(x))
    }
    y_index <- attr(x, "response")

    ## If there is something on the lhs of the formula,
    ## remove it and get vars
    if (y_index != 0) {
      x[[2]] <- NULL
      x <- stats::terms(x)
    }
    all.vars(x, ...)
  })
  names(res) <- NULL

  if(length(res) == 1) {
    return(unlist(res))
  } else {
    return(res)
  }
}
