#' Extract Predictor Names from Formulas
#'
#' @usage extract_var_from_formulas(formulas, ...)
#' @param formulas (character or formula) model formulas.
#' @param ... Arguments to pass to [all.vars()]
#'
#' @return A character vector or a list of the same length as the formulas,
#' containing the variable names.
#'
#' @export
#'
#' @importFrom stats terms as.formula
#'
#' @examples
#' # Extract variables from the formula grid
#' # Import an example of calibration results (output from the calibration function)
#' data("calib_results_glmnet", package = "kuenm2")
#' extract_var_from_formulas(formulas = calib_results_glmnet$formula_grid$Formulas)
#'
extract_var_from_formulas <- function(formulas, ...) {
  res <- sapply(formulas, function(object) {
    if (inherits(object, "formula")) {
      object <- stats::terms(object)
    } else {
      object <- stats::terms(stats::as.formula(object))
    }
    y_index <- attr(object, "response")

    ## If there is something on the lhs of the formula,
    ## remove it and get vars
    if (y_index != 0) {
      object[[2]] <- NULL
      object <- stats::terms(object)
    }
    all.vars(object, ...)
  })
  names(res) <- NULL
  if(length(res) == 1) {
    return(unlist(res))} else {
  return(res)}
}
