#' Extracts Environmental Variables for Occurrences
#'
#' @description
#' This function extracts values from environmental or predictor variables
#' (`SpatRaster`) for georeferenced occurrence points. It also adds a column
#' indicating that these are presence points(pr_bg = 1).
#'
#' @param occ A data.frame containing occurrence data. It must include columns
#' with longitude (x) and latitude (y) coordinates.
#' @param x (character) a string specifying the name of the column in occ that
#' contains the longitude values.
#' @param y (character) a string specifying the name of the column in occ that
#' contains the latitude values.
#' @param raster_variables (SpatRaster) predictor variables used to calibrate
#' the models.
#'
#' @return
#' A data.frame containing the original x and y coordinates of the occurrence
#' points (`x` and `y`), the values of the variables extracted
#' from `raster_variables`, and a new column `pr_bg` with a value of 1 for all
#' occurrences.
#' @importFrom terra extract
#' @examples
#' # Import occurrences
#' data(occ_data, package = "kuenm2")
#'
#' # Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' # Extracts environmental variables for occurrences
#' occ_var <- extract_occurrence_variables(occ = occ_data, x = "x", y = "y",
#'                                         raster_variables = var)
#' @export
extract_occurrence_variables <- function(occ, x, y, raster_variables) {
  xy <- as.matrix(occ[, c(x, y)])
  colnames(xy) <- c("x", "y")
  occ_var <- cbind("pr_bg" = 1, xy, terra::extract(x = raster_variables,
                                                   y = xy))
  return(occ_var)
}
