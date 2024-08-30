#' Explore the spatial distribution of occurrence and background points
#'
#' @param data an object of class `???` returned by the prepare_data() function
#' @param spat_variables (SpatRaster) predictor variables used for model calibration.
#' @param plot (logical) wheter to plot the SpatRaster. Default is TRUE.
#'
#' @return
#' A categorical `SpatRaster` with four factors representing:
#' \describe{
#'   \item{1}{Background cells}
#'   \item{2}{Presence cells}
#'   \item{3}{Cells with both presence and background}
#'   \item{4}{Non-used cells}
#' }
#'
#' @export
#' @importFrom terra extract NAflag levels trim crop plot
#' @usage explore_calibration_geo(data, spat_variables, plot = TRUE)
#' @examples
#' #Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#' #Import occurrences
#' data(occ_data, package = "kuenm2")
#' #Prepare data
#' sp_swd <- prepare_data(occ = occ_data,
#'                        species = occ_data[1, 1],
#'                        x = "x",
#'                        y = "y",
#'                        spat_variables = var,
#'                        mask = NULL,
#'                        categorical_variables = "SoilType",
#'                        do_pca = TRUE,
#'                        exclude_from_pca = NULL,
#'                        nbg = 500,
#'                        kfolds = 4,
#'                        weights = NULL,
#'                        include_xy = TRUE,
#'                        write_file = FALSE,
#'                        file_name = NULL,
#'                        seed = 1)
#' #Explore calibration data
#' pbg <- explore_calibration_geo(data = sp_swd, spat_variables = var[[1]],
#'                                plot = FALSE)
#'
explore_calibration_geo <- function(data,
                                    spat_variables,
                                    plot = TRUE){

  #Extract cell values
  bg <- data$data_xy[data$calibration_data$pr_bg == 0,]
  p <- data$data_xy[data$calibration_data$pr_bg == 1,]

  #Get cells
  bg_cells <- terra::extract(spat_variables[[1]], bg, cells = TRUE, ID = FALSE)$cell
  p_cells <- terra::extract(spat_variables[[1]], p, cells = TRUE, ID = FALSE)$cell
  p_bg_cells <- intersect(bg_cells, p_cells)

  #Fill rasters
  r <- (!is.na(spat_variables[[1]])) * 1
  terra::NAflag(r) <- 0 #Set NA
  r[bg_cells] <- 2
  r[p_cells] <- 3
  r[p_bg_cells] <- 4
  #Set levels
  levels(r) <- data.frame(id=1:4, category=c("Unused data", "Background", "Presence",
                                             "Presence & Background"))
  #Crop and trim final raster
  r <- terra::trim(terra::crop(r, spat_variables[[1]], mask = TRUE))

  if(plot){
  terra::plot(r)}
  return(r)
}
