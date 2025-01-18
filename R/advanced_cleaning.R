#' Advanced occurrence data cleaning
#'
#' @description Advanced processes of data cleaning involving duplicate removal
#' and movement of records.
#'
#' @param data data.frame with occurrence records. Rows with NA values will be
#' omitted.
#' @param longitude_column (character) name of the column in \code{data}
#' containing longitude values.
#' @param latitude_column (character) name of the column in \code{data}
#' containing latitude values.
#' @param raster_layer a raster layer (object of class
#' \code{\link[terra]{SpatRaster}}).
#' @param cell_duplicates (logical) whether to remove duplicate coordinates
#' considering raster cells. Default = TRUE.
#' @param move_points_inside (logical) whether to move records outside of raster
#' cells with valid values to the closest cell with values. Default = TRUE.
#' @param move_limit_distance maximum distance to move records outside cells
#' with valid values. Default = NULL. Must be defined if \code{move_points_inside}
#' = TRUE.
#' @param verbose (logical) whether to print messages of progress. Default =
#' TRUE.
#'
#' @details
#' Data used in this functions should have gone through initial processes of
#' cleaning and filtering.
#'
#' @return
#' A data.frame with occurrence records resulting from advanced cleaning
#' procedures. Other columns will be added to describe changes made in the
#' original data.
#'
#' @seealso
#' \code{\link{initial_cleaning}}
#'
#' @export
#'
#' @rdname advanced_cleaning
#'
#' @usage
#' advanced_cleaning(data, longitude_column, latitude_column, raster_layer,
#'                   cell_duplicates = TRUE, move_points_inside = TRUE,
#'                   move_limit_distance = NULL, verbose = TRUE)

advanced_cleaning <- function(data, longitude_column, latitude_column,
                              raster_layer, cell_duplicates = TRUE,
                              move_points_inside = TRUE,
                              move_limit_distance = NULL, verbose = TRUE) {
  # error checking
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (missing(longitude_column)) {
    stop("Argument 'longitude_column' must be defined.")
  }
  if (missing(latitude_column)) {
    stop("Argument 'latitude_column' must be defined.")
  }
  if (missing(raster_layer)) {
    stop("Argument 'raster_layer' must be defined.")
  }
  if (!class(data)[1] %in% c("data.frame", "matrix")) {
    stop("'data' must be of class 'data.frame' or 'matrix'.")
  }
  if (class(raster_layer)[1] != "SpatRaster") {
    stop("'raster_layer' must be of class 'SpatRaster'")
  }

  data <- na.omit(data)

  # remove cell duplicates
  if (cell_duplicates == TRUE) {
    data <- remove_cell_duplicates(data, longitude_column, latitude_column,
                                   raster_layer)
  }

  # move to closest pixel
  if (move_points_inside == TRUE) {
    if (is.null(move_limit_distance)) {
      stop("'move_limit_distance' must be defined if 'move_points_inside' = TRUE.")
    }
    data <- move_2closest_cell(data, longitude_column, latitude_column,
                               raster_layer, move_limit_distance,
                               verbose)
  }

  # return results (metadata)
  return(data)
}



#' @export
#' @importFrom terra extract
#' @rdname advanced_cleaning
#' @usage
#' remove_cell_duplicates(data, longitude_column, latitude_column,
#'                        raster_layer)

remove_cell_duplicates <- function(data, longitude_column, latitude_column,
                                   raster_layer) {
  # detecting potential errors
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (missing(longitude_column)) {
    stop("Argument 'longitude_column' must be defined.")
  }
  if (missing(latitude_column)) {
    stop("Argument 'latitude_column' must be defined.")
  }
  if (missing(raster_layer)) {
    stop("Argument 'raster_layer' must be defined.")
  }
  if (class(raster_layer)[1] != "SpatRaster") {
    stop("'raster_layer' must be of class 'SpatRaster'")
  }

  # preparing data
  xy <- c(longitude_column, latitude_column)
  cells <- terra::extract(raster_layer, as.matrix(data[, xy]), cells = TRUE)

  # returning data (metadata?)
  return(data[!duplicated(cells[, "cell"]), ])
}


#' @export
#' @importFrom terra extract vect buffer crop distance
#' @rdname advanced_cleaning
#' @usage
#' move_2closest_cell(data, longitude_column, latitude_column, raster_layer,
#'                    move_limit_distance, verbose = TRUE)

move_2closest_cell <- function(data, longitude_column, latitude_column,
                               raster_layer, move_limit_distance,
                               verbose = TRUE) {

  # detecting potential errors
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (missing(longitude_column)) {
    stop("Argument 'longitude_column' must be defined.")
  }
  if (missing(latitude_column)) {
    stop("Argument 'latitude_column' must be defined.")
  }
  if (missing(raster_layer)) {
    stop("Argument 'raster_layer' must be defined.")
  }
  if (missing(move_limit_distance)) {
    stop("Argument 'move_limit_distance' must be defined.")
  }
  if (class(raster_layer)[1] != "SpatRaster") {
    stop("'raster_layer' must be of class 'SpatRaster'")
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



