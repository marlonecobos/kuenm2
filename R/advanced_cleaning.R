#' Advanced occurrence data cleaning
#'
#' @description Advanced processes of data cleaning involving duplicate removal
#' and movement of records.
#'
#' @usage
#' advanced_cleaning(data, x, y, raster_layer, cell_duplicates = TRUE,
#'                   move_points_inside = FALSE, move_limit_distance = NULL,
#'                   verbose = TRUE)
#'
#' @param data data.frame with occurrence records. Rows with NA values will be
#' omitted.
#' @param x (character) name of the column in \code{data}
#' containing longitude values.
#' @param y (character) name of the column in \code{data}
#' containing latitude values.
#' @param raster_layer a raster layer (object of class
#' \code{\link[terra]{SpatRaster}}).
#' @param cell_duplicates (logical) whether to remove duplicate coordinates
#' considering raster cells. Default = TRUE.
#' @param move_points_inside (logical) whether to move records outside of raster
#' cells with valid values to the closest cell with values. Default = FALSE.
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
#' [initial_cleaning()]
#'
#' @export
#'
#' @rdname advanced_cleaning
#'
#' @examples
#' # Import occurrences
#' data(occ_data_noclean, package = "kuenm2")
#'
#' # Import raster layers
#' var <- rast(system.file("extdata", "Current_variables.tif", package = "kuenm2"))
#'
#' # Keep only one layer
#' var <- var$bio_1
#'
#' # all basic cleaning steps
#' clean_init <- initial_cleaning(data = occ_data_noclean, species = "species",
#'                                x = "x", y = "y", remove_na = TRUE,
#'                                remove_empty = TRUE, remove_duplicates = TRUE,
#'                                by_decimal_precision = TRUE,
#'                                decimal_precision = 2)
#'
#' # Advanced cleaning steps
#' # exclude duplicates based on raster cell (pixel)
#' celldup <- remove_cell_duplicates(data = clean_init, x = "x", y = "y",
#'                                   raster_layer = var)
#'
#' # move records to valid pixels
#' moved <- move_2closest_cell(data = celldup, x = "x", y = "y",
#'                             raster_layer = var, move_limit_distance = 10)
#'
#' # the steps at a time
#' clean_data <- advanced_cleaning(data = clean_init, x = "x", y = "y",
#'                                 raster_layer = var, cell_duplicates = TRUE,
#'                                 move_points_inside = TRUE,
#'                                 move_limit_distance = 10)

advanced_cleaning <- function(data, x, y,
                              raster_layer, cell_duplicates = TRUE,
                              move_points_inside = FALSE,
                              move_limit_distance = NULL, verbose = TRUE) {
  # error checking
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (missing(x)) {
    stop("Argument 'x' must be defined.")
  }
  if (missing(y)) {
    stop("Argument 'y' must be defined.")
  }
  if (missing(raster_layer)) {
    stop("Argument 'raster_layer' must be defined.")
  }
  if (!class(data)[1] %in% c("data.frame", "matrix")) {
    stop("'data' must be of class 'data.frame' or 'matrix'.")
  }
  if (class(raster_layer)[1] != "SpatRaster") {
    stop("'raster_layer' must be of class 'SpatRaster'.")
  }

  data <- na.omit(data)

  # remove cell duplicates
  if (cell_duplicates == TRUE) {
    data <- remove_cell_duplicates(data, x, y, raster_layer)
  }

  # move to closest pixel
  if (move_points_inside == TRUE) {
    if (is.null(move_limit_distance)) {
      stop("'move_limit_distance' must be defined if 'move_points_inside' = TRUE.")
    }
    data <- move_2closest_cell(data, x, y, raster_layer, move_limit_distance,
                               verbose)

    # remove cell duplicates again after moving points
    if (cell_duplicates == TRUE) {
      data <- remove_cell_duplicates(data, x, y, raster_layer)
    }
  }

  # return results (metadata)
  return(data)
}



#' @export
#' @importFrom terra extract
#' @rdname advanced_cleaning
#' @usage
#' remove_cell_duplicates(data, x, y,
#'                        raster_layer)

remove_cell_duplicates <- function(data, x, y, raster_layer) {
  # detecting potential errors
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (missing(x)) {
    stop("Argument 'x' must be defined.")
  }
  if (missing(y)) {
    stop("Argument 'y' must be defined.")
  }
  if (missing(raster_layer)) {
    stop("Argument 'raster_layer' must be defined.")
  }
  if (class(raster_layer)[1] != "SpatRaster") {
    stop("'raster_layer' must be of class 'SpatRaster'.")
  }

  # preparing data
  xy <- c(x, y)
  cells <- terra::extract(raster_layer, as.matrix(data[, xy]), cells = TRUE)

  # returning data (metadata?)
  return(data[!duplicated(cells[, "cell"]), ])
}


#' @export
#' @importFrom terra extract vect buffer crop distance
#' @rdname advanced_cleaning
#' @usage
#' move_2closest_cell(data, x, y, raster_layer,
#'                    move_limit_distance, verbose = TRUE)

move_2closest_cell <- function(data, x, y, raster_layer, move_limit_distance,
                               verbose = TRUE) {

  # detecting potential errors
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (missing(x)) {
    stop("Argument 'x' must be defined.")
  }
  if (missing(y)) {
    stop("Argument 'y' must be defined.")
  }
  if (missing(raster_layer)) {
    stop("Argument 'raster_layer' must be defined.")
  }
  if (missing(move_limit_distance)) {
    stop("Argument 'move_limit_distance' must be defined.")
  }
  if (class(raster_layer)[1] != "SpatRaster") {
    stop("'raster_layer' must be of class 'SpatRaster'.")
  }

  # preparing data
  xy <- data[, c(x, y)]
  vals <- terra::extract(raster_layer, xy, ID = FALSE)
  tomove <- which(is.na(vals[, 1]))

  # finding pixels to move in
  if (length(tomove) > 0) {
    xyout <- data[tomove, ]
    no <- nrow(xyout)

    ## basic conditions to fill in table
    condition <- rep("Correct", nrow(data))
    distss <- rep(0, nrow(data))

    ## buffer from out
    limdist <- move_limit_distance * 1000
    xyvec <- terra::vect(xyout, crs = "+proj=longlat",
                         geom = c(x, y))
    xyvec <- terra::buffer(xyvec, width = limdist)

    ## relevant pixels to move
    raster_layer <- terra::crop(raster_layer, xyvec, mask = TRUE)

    xyras <- as.data.frame(raster_layer, xy = TRUE)[, 1:2]

    if (nrow(xyras) >= 1) {
      dists <- terra::distance(xyout[, c(x, y)],
                               xyras, lonlat = TRUE)

      # running process
      if (verbose == TRUE) {
        message("Moving occurrences to closest pixels...")
      }

      for (i in 1:no) {
        mindis <- min(dists[i, ])
        if (mindis <= limdist) {
          xyin <- xyras[dists[i, ] == mindis, ]
          data[tomove[i], x] <- xyin[1, 1]
          data[tomove[i], y] <- xyin[1, 2]
          condition[tomove[i]] <- "Moved"
          distss[tomove[i]] <- mindis/1000
        } else {
          condition[tomove[i]] <- "Not_moved"
          distss[tomove[i]] <- mindis/1000
        }
      }

    } else {
      for (i in 1:no) {
        condition[tomove[i]] <- "Not_moved"
        distss[tomove[i]] <- NA
      }

      if (verbose == TRUE) {
        message("Occurrences could not be moved due to long distances.")
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



