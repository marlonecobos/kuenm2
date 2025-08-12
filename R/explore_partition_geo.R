#' Explore the spatial distribution of partitions for occurrence and background points
#'
#' @usage
#' explore_partition_geo(data, raster_variables, mask = NULL,
#'                       show_partitions = TRUE, partition_palette = "cols25",
#'                       custom_partition_palette = NULL, pr_col = "#D55E00",
#'                       bg_col = "#0072B2", pr_bg_col = "#CC79A7",
#'                       calibration_area_col = "gray80", ...)
#'
#' @param data an object of class `prepared_data` returned by the
#' [prepare_data()] function.
#' @param raster_variables (SpatRaster) predictor variables used for model
#' calibration.
#' @param mask (SpatRaster, SpatVector, or SpatExtent) spatial object used to
#' mask `raster_variables` to the area where the model will be calibrated.
#' Preferably the same object used in `prepare_data` (if applicable). Default is
#' NULL.
#' @param show_partitions (logical) whether to return `SpatRaster` showing the
#' spatial distribution of each partition for presence and background points.
#' Default is TRUE.
#' @param partition_palette (character) the color palette used to color the
#' different partitions. See `?kuenm2_discrete_palettes` to check available
#' options. Default is `"cols25"`. Only applicable if `show_partitions = TRUE`.
#' @param custom_partition_palette (character) a character vector defining
#' custom colors for the different partitions. The number of values must match
#' the number of partitions in `data`. Default is NULL, meaning the palette
#' defined in `partition_palette` will be used.
#' @param pr_col (character) the color used for cells with presence records.
#' Default is "#D55E00".
#' @param bg_col (character) the color used for cells with background points.
#' Default is "#0072B2".
#' @param pr_bg_col (character) the color used for cells with presences and
#' background points. Default is "#CC79A7".
#' @param calibration_area_col (character) the color used for cells without
#' presences or background points. Default is "gray80".
#' @param ... additional arguments passed to `terra::plot()`.
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
#' If show_partitions = TRUE, it also returns `SpatRaster` showing the spatial
#' distribution of each partition for presence and background points.
#'
#' @export
#'
#' @importFrom terra extract NAflag levels trim crop mask coltab
#'
#' @examples
#' # Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' # Import prepared_data
#' data(sp_swd, package = "kuenm2")
#'
#' # Explore partitions in the geographic space
#' pbg <- explore_partition_geo(data = sp_swd, raster_variables = var[[1]])
#' terra::plot(pbg)

explore_partition_geo <- function(data,
                                  raster_variables,
                                  mask = NULL,
                                  show_partitions = TRUE,
                                  partition_palette = "cols25",
                                  custom_partition_palette = NULL,
                                  pr_col = "#D55E00",
                                  bg_col = "#0072B2",
                                  pr_bg_col = "#CC79A7",
                                  calibration_area_col = "gray80",
                                  ...) {

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

  if (!is.null(mask) &
      !inherits(mask, c("SpatRaster", "SpatVector", "SpatExtent"))) {
    stop("Argument 'mask' must be a 'SpatVector', 'SpatExtent' or 'SpatRaster'.")
  }

  if(is.null(data$data_xy)){
    stop("'data' must have the 'xy' element.
Make sure you run 'prepare_data()' with 'include_xy' set to TRUE")
  }

  if(is.null(custom_partition_palette) &&
     !(partition_palette %in% names(kuenm2::kuenm2_discrete_palletes)))
    stop("Invalid 'partition_palette'. Check the available options in '?kuenm2_palletes'")

  if(!is.null(custom_partition_palette) &&
     length(custom_partition_palette) < length(data$part_data))
    stop("Insufficient number of colors in 'custom_partition_palette'.
Provide at least ", length(data$part_data), " colors.")

  if(!inherits(pr_col, "character")){
    stop("'pr_col' must be a character")
  }
  if(!inherits(bg_col, "character")){
    stop("'bg_col' must be a character")
  }
  if(!inherits(pr_bg_col, "character")){
    stop("'pr_bg_col' must be a character")
  }
  if(!inherits(calibration_area_col, "character")){
    stop("'calibration_area_col' must be a character")
  }

  #Extract cell values
  bg <- data$data_xy[data$calibration_data$pr_bg == 0,]
  p <- data$data_xy[data$calibration_data$pr_bg == 1,]

  #Get cells
  bg_cells <- terra::extract(raster_variables[[1]], bg,
                             cells = TRUE, ID = FALSE)$cell
  p_cells <- terra::extract(raster_variables[[1]], p,
                            cells = TRUE, ID = FALSE)$cell
  p_bg_cells <- intersect(bg_cells, p_cells)

  if(!is.null(mask)){
    raster_variables <- terra::mask(raster_variables[[1]], mask)
  }

  #Fill rasters
  r <- (!is.na(raster_variables[[1]])) * 1
  terra::NAflag(r) <- 0  # Set NA
  r[bg_cells] <- 2
  r[p_cells] <- 3
  r[p_bg_cells] <- 4

  #Set levels
  levels(r) <- data.frame(id = 1:4,
                          category = c("Calibration area", "Background", "Presence",
                                       "Presence & Background"))

  #Crop and trim final raster
  r <- terra::trim(terra::mask(r, raster_variables[[1]]))
  names(r) <- "Presence and background"

  #Set colors
  r_col <- data.frame(values = 1:4,
                      col = c(calibration_area_col, bg_col, pr_col, pr_bg_col))
  terra::coltab(r) <- r_col

  #If get results showing partition
  if(show_partitions){
    # Get partitions
    np <- names(data$part_data)

    #Create empty rasters to fill
    res_presence <- (!is.na(raster_variables[[1]])) * 0
    res_background <- res_presence

    #For each partition, fill the R with the values
    for(i in 1:length(np)) {
      pi <- data$part_data[[i]]
      pi_pr <- intersect(which(data$calibration_data$pr_bg == 1), pi)
      presence_i <- data$data_xy[pi_pr,]
      pi_bg <- intersect(which(data$calibration_data$pr_bg == 0), pi)
      background_i <- data$data_xy[pi_bg,]

      #Get cells
      pi_pr_cells <- terra::extract(res_presence, presence_i, cells = TRUE)[["cell"]]
      pi_bg_cells <- terra::extract(res_presence, background_i, cells = TRUE)[["cell"]]

      #Fill
      res_presence[pi_pr_cells] <- i
      res_background[pi_bg_cells] <- i
    }

    #Mask rasters
    res_presence <- terra::trim(terra::mask(res_presence, raster_variables[[1]]))
    res_background <- terra::trim(terra::mask(res_background, raster_variables[[1]]))

    #Set levels
    l <- data.frame(id = 0:length(np),
                    Partition = c("Calibration area", np))
    levels(res_presence) <- l
    levels(res_background) <- l

    #Set colors
    if(is.null(custom_partition_palette)){
      cores <- kuenm2::kuenm2_discrete_palletes[[partition_palette]]
    } else {
      cores <- custom_partition_palette
    }
    #Create dataframe with colors
    df_colors <- data.frame(values = 0:length(np),
                            col = c(calibration_area_col, cores[1:length(np)]))

    terra::coltab(res_presence) <- df_colors
    terra::coltab(res_background) <- df_colors

    #Join results
    res_partition <- c(res_presence, res_background)
    names(res_partition) <- c("Presence", "Background")

  } else {
    res_partition <- NULL
  }

  return(c(r, res_partition))


}
