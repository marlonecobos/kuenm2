# makes NA values in all layers match
fix_na_rast <- function(raster_layers, return_mask = FALSE) {
  if (missing(raster_layers)) {
    stop("Argument 'raster_layers' must be defined.")
  }
  if (class(raster_layers)[1] != "SpatRaster") {
    stop("Argument 'raster_layers' must of class 'SpatRaster'.")
  }

  # processing
  if (terra::nlyr(raster_layers) <= 1) {
    warning("'raster_layers' has only one layer, processing is not needed.")
    msk <- NULL
  } else {
    msk <- rowSums(raster_layers[])
    msk[!is.na(msk)] <- 1
    raster_layers[] <- raster_layers[] * msk
  }

  # results
  if (return_mask == TRUE) {
    return(list(raster_layers = raster_layers, mask = msk))
  } else {
    return(raster_layers)
  }
}
