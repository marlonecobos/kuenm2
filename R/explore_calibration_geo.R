#' Explore the spatial distribution of occurrence and background points
#'
#' @usage
#' explore_calibration_geo(data, raster_variables, plot = TRUE)
#'
#' @param data an object of class `prepared_data` returned by the
#' \code{\link{prepare_data}}() function.
#' @param raster_variables (SpatRaster) predictor variables used for model
#' calibration.
#' @param plot (logical) wheter to plot the SpatRaster. Default is TRUE.
#'
#' @return
#' A categorical `SpatRaster` with four factor values representing:
#' \describe{
#'   \item{1 - Background cells}{}
#'   \item{2 - Presence cells}{}
#'   \item{3 - Cells with both presence and background}{}
#'   \item{4 - Non-used cells}{}
#' }
#'
#' @export
#'
#' @importFrom terra extract NAflag levels trim crop plot
#'
#' @examples
#' # Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' # Import occurrences
#' data(sp_swd, package = "kuenm2")
#'
#' # Explore calibration data
#' pbg <- explore_calibration_geo(data = sp_swd, raster_variables = var[[1]])

explore_calibration_geo <- function(data,
                                    raster_variables,
                                    plot = TRUE) {

  #Check data
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (missing(raster_variables)) {
    stop("Argument 'raster_variables' must be defined.")
  }
  if(!inherits(data, "prepared_data")){
    stop("'data' must be a 'prepared_data' object.")
  }

  if(!inherits(raster_variables, "SpatRaster"))
    stop("'raster_variables' must be a 'SpatRaster' object.")

  if(!inherits(plot, "logical"))
    stop("'plot' must be logical.")

  #Extract cell values
  bg <- data$data_xy[data$calibration_data$pr_bg == 0,]
  p <- data$data_xy[data$calibration_data$pr_bg == 1,]

  #Get cells
  bg_cells <- terra::extract(raster_variables[[1]], bg,
                             cells = TRUE, ID = FALSE)$cell
  p_cells <- terra::extract(raster_variables[[1]], p,
                            cells = TRUE, ID = FALSE)$cell
  p_bg_cells <- intersect(bg_cells, p_cells)

  #Fill rasters
  r <- (!is.na(raster_variables[[1]])) * 1
  terra::NAflag(r) <- 0  # Set NA
  r[bg_cells] <- 2
  r[p_cells] <- 3
  r[p_bg_cells] <- 4

  #Set levels
  levels(r) <- data.frame(id = 1:4,
                          category = c("Unused data", "Background", "Presence",
                                       "Presence & Background"))

  #Crop and trim final raster
  r <- terra::trim(terra::crop(r, raster_variables[[1]], mask = TRUE))

  if (plot) {
    terra::plot(r)
  }

  return(r)
}
