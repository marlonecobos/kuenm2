advanced_cleaning <- function(data, longitude_column, latitude_column,
                              raster_layer, cell_duplicates = TRUE,
                              move_points_inside = TRUE,
                              move_limit_distance = NULL, verbose = TRUE) {
  # error checking


  # remove cell duplicates
  if (cell_duplicates == TRUE) {
    data <- remove_cell_duplicates(data, longitude_column, latitude_column,
                                   raster_layer)
  }

  # move to closest pixel
  if (move_points_inside == TRUE) {
    data <- move_2closest_cell(data, longitude_column, latitude_column,
                               raster_layer, move_limit_distance,
                               verbose)
  }

  # return results (metadata)
  return(data)
}




remove_cell_duplicates <- function(data, longitude_column, latitude_column,
                                   raster_layer) {
  # detecting potential errors
  if (missing(data)) {
    stop("Argument 'data' is necessary to perform the analysis.")
  }
  if (missing(longitude_column)) {
    stop("Argument 'longitude_column' is not defined.")
  }
  if (missing(latitude_column)) {
    stop("Argument 'latitude_column' is not defined.")
  }

  # preparing data
  xy <- c(longitude_column, latitude_column)
  cells <- terra::extract(raster_layer, as.matrix(data[, xy]), cells = TRUE)

  # returning data (metadata?)
  return(data[!duplicated(cells[, "cell"]), ])
}



move_2closest_cell <- function(data, longitude_column, latitude_column,
                               raster_layer, move_limit_distance,
                               verbose = TRUE) {

  # detecting potential errors
  if (missing(data)) {
    stop("Argument 'data' is necessary to perform the analysis.")
  }
  if (missing(longitude_column)) {
    stop("Argument 'longitude_column' is not defined.")
  }
  if (missing(latitude_column)) {
    stop("Argument 'latitude_column' is not defined.")
  }
  if (missing(raster_layer)) {
    stop("Argument 'raster_layer' is not defined.")
  }
  if (missing(move_limit_distance)) {
    stop("Argument 'move_limit_distance' is not defined.")
  }

  # preparing data
  xy <- data[, c(longitude_column, latitude_column)]
  vals <- terra::extract(raster_layer, xy)

  tomove <- which(is.na(vals[, 1]))

  # finding pixels to move in
  if (length(tomove) > 0) {
    xyout <- data[tomove, ]

    ## buffer from out
    limdist <- move_limit_distance * 1000

    xyvec <- terra::vect(xyout, crs = "+proj=longlat")
    xyvec <- terra::buffer(xyvec, width = limdist)

    ## relevant pixels to move
    raster_layer <- terra::crop(raster_layer, xyvec, mask = TRUE)

    xyras <-  as.data.frame(raster_layer, xy = TRUE)[, 1:2]
    dists <- terra::distance(xyout, xyras, lonlat = TRUE)

    condition <- rep("Correct", nrow(data))
    distss <- rep(0, nrow(data))

    # running process
    if (verbose == TRUE) {
      message("Moving occurrences to closest pixels...")
    }

    no <- nrow(xyout)

    for (i in 1:no) {
      mindis <- min(dists[i, ])

      if (mindis <= limdist) {
        xyin <- xyras[dists[i, ] == mindis, ]

        data[tomove[i], longitude_column] <- xyin[1, 1]
        data[tomove[i], latitude_column] <- xyin[1, 2]
        condition[tomove[i]] <- "Moved"
        distss[tomove[i]] <- mindis / 1000
      } else {
        condition[tomove[i]] <- "Not_moved"
        distss[tomove[i]] <- mindis / 1000
      }
    }
    data <- data.frame(data, condition = condition, distance_km = distss,
                       initial_lon = xy[, 1], initial_lat = xy[, 2])
  } else {
    if (verbose == TRUE) {
      message("All points are inside valid raster boundaries.")
    }
  }

  return(data)
}



