#' Explore the spatial distribution of occurrence and background points
#'
#' @param data an object of class `prepare_data` returned by the
#'            \code{\link{prepare_data}}() function.
#' @param raster_variables (SpatRaster) predictor variables used for model calibration.
#' @param plot (logical) wheter to plot the SpatRaster. Default is TRUE.
#'
#' @return
#' A categorical `SpatRaster` with four factors representing:
#' \describe{
#'   \item{1 - Background cells}{}
#'   \item{2 - Presence cells}{}
#'   \item{3 - Cells with both presence and background}{}
#'   \item{4 - Non-used cells}{}
#' }
#'
#' @export
#' @importFrom terra extract NAflag levels trim crop plot
#' @usage explore_calibration_geo(data, raster_variables, plot = TRUE)
#' @examples
#' #Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#' #Import occurrences
#' data(occ_data, package = "kuenm2")
#' #Prepare data
#' sp_swd <- prepare_data(model_type = "glmnet", occ = occ_data,
#'                        species = occ_data[1, 1], x = "x", y = "y",
#'                        raster_variables = var, mask = NULL,
#'                        categorical_variables = "SoilType",
#'                        do_pca = FALSE, deviance_explained = 95,
#'                        min_explained = 5, center = TRUE, scale = TRUE,
#'                        write_pca = FALSE, output_pca = NULL, nbg = 500,
#'                        kfolds = 4, weights = NULL, min_number = 2,
#'                        min_continuous = NULL,
#'                        features = c("l", "q", "p", "lq", "lqp"),
#'                        regm = c(0.1, 1, 2, 3, 5),
#'                        include_xy = TRUE,
#'                        write_file = FALSE, file_name = NULL,
#'                        seed = 1)
#' #Explore calibration data
#' pbg <- explore_calibration_geo(data = sp_swd, raster_variables = var[[1]],
#'                                plot = FALSE)
explore_calibration_geo <- function(data,
                                    raster_variables,
                                    plot = TRUE){

  #Check data
  if(!inherits(data, "prepared_data")){
    stop("'data' must be a 'prepared_data' object, not ", class(data))
  }

  if(!inherits(raster_variables, "SpatRaster"))
    stop("'raster_variables' must be a 'SpatRaster' object, not ", class(raster_variables))

  if(!inherits(plot, "logical"))
    stop("'plot' must be logical, not ", class(plot))

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
  terra::NAflag(r) <- 0 #Set NA
  r[bg_cells] <- 2
  r[p_cells] <- 3
  r[p_bg_cells] <- 4
  #Set levels
  levels(r) <- data.frame(id=1:4,
                          category=c("Unused data", "Background", "Presence",
                                     "Presence & Background"))
  #Crop and trim final raster
  r <- terra::trim(terra::crop(r, raster_variables[[1]], mask = TRUE))

  if(plot){
  terra::plot(r)}
  return(r)
}
