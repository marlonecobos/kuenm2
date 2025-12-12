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

  if(inherits(res, "matrix")) {
    return(res[,1])
  } else {
    return(res)
  }
}

extract_categorical <- function(formulas){
  res <- sapply(formulas, function(formula_str) {

    # Try to find categorical
    match_location <- regexpr("categorical\\(([^)]+)\\)", formula_str, perl = TRUE)

    if (match_location > 0) {
      full_match <- regmatches(formula_str, match_location)

      variable_name <- gsub(
        pattern = "categorical\\(([^)]+)\\)",
        replacement = "\\1",
        x = full_match
      )
      return(unique(variable_name))

    } else {
      return(NULL)
    }
  }, USE.NAMES = FALSE)
  return(unlist(res))
}
