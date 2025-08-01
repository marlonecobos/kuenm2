#' Initial occurrence data cleaning steps
#'
#' @description Simple occurrence data cleaning procedures.
#'
#' @usage
#' initial_cleaning(data, species, x, y,
#'                  other_columns = NULL, keep_all_columns = TRUE,
#'                  sort_columns = TRUE, remove_na = TRUE, remove_empty = TRUE,
#'                  remove_duplicates = TRUE, by_decimal_precision = FALSE,
#'                  decimal_precision = 0, longitude_precision = NULL,
#'                  latitude_precision = NULL)
#'
#' @param data data.frame with occurrence records.
#' @param species (character) name of the column in \code{data}
#' containing species name.
#' @param x (character) name of the column in \code{data}
#' containing longitude values.
#' @param y (character) name of the column in \code{data}
#' containing latitude values.
#' @param other_columns (character) vector of other column name(s) in
#' \code{data} to be considered while performing cleaning steps, default = NULL.
#' @param keep_all_columns (logical) whether to keep all columns in \code{data}.
#' Default = TRUE.
#' @param sort_columns (logical) whether to sort species, longitude, and
#' latitude columns in \code{data}. Default = TRUE.
#' @param remove_na (logical) whether to remove NA values in the columns
#' considered. Default = TRUE.
#' @param remove_empty (logical) whether to remove empty (missing) values in
#' the columns considered. Default = TRUE.
#' @param remove_duplicates (logical) whether to remove duplicates in
#' the columns considered. Default = TRUE.
#' @param by_decimal_precision (logical) whether to remove certain records with
#' coordinate precision lower than that of the following three parameters.
#' Default = FALSE
#' @param decimal_precision (numeric) decimal precision threshold for
#' coordinates. Default = 0. Ignored if the following two parameters are defined.
#' @param longitude_precision (numeric) decimal precision threshold for
#' longitude. Default = NULL.
#' @param latitude_precision (numeric) decimal precision threshold for
#' latitude. Default = NULL.
#'
#' @details
#' Function \code{initial_cleaning} helps to perform all simple steps of data
#' cleaning.
#'
#' @return
#' A data.frame with resulting occurrence records.
#'
#' @seealso
#' \code{\link{advanced_cleaning}}
#'
#' @export
#'
#' @rdname initial_cleaning


initial_cleaning <- function(data, species, x,
                             y, other_columns = NULL,
                             keep_all_columns = TRUE,
                             sort_columns = TRUE, remove_na = TRUE,
                             remove_empty = TRUE, remove_duplicates = TRUE,
                             by_decimal_precision = FALSE,
                             decimal_precision = 0, longitude_precision = NULL,
                             latitude_precision = NULL) {

  # error checking
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (missing(species)) {
    stop("Argument 'species' must be defined.")
  }
  if (missing(x)) {
    stop("Argument 'x' must be defined.")
  }
  if (missing(y)) {
    stop("Argument 'y' must be defined.")
  }
  if (class(data)[1] != "data.frame") {
    stop("'data' must be of class 'data.frame'.")
  }

  clnames <- colnames(data)
  if (!is.null(other_columns)) {
    if (all(!other_columns %in% clnames)) {
      stop("None of the 'other_columns' is in 'data'.")
    }
    if (any(!other_columns %in% clnames)) {
      other_columns <- other_columns[other_columns %in% clnames]
      warning("Some 'other_columns' were not in 'data'.")
    }
  }

  # preparing arguments
  if (!is.null(other_columns)) {
    columns <- c(species, x, y, other_columns)
  } else {
    columns <- c(species, x, y)
  }

  # cleaning steps
  if (sort_columns == TRUE) {
    data <- sort_columns(data, species, x, y, keep_all_columns)
  }

  if (remove_na == TRUE | remove_empty == TRUE) {
    data <- remove_missing(data, columns, remove_na, remove_empty,
                           keep_all_columns)
  }

  if (remove_duplicates == TRUE) {
    data <- remove_duplicates(data, columns, keep_all_columns)
  }

  if (by_decimal_precision == TRUE) {
    data <- filter_decimal_precision(data, x, y, decimal_precision,
                                     longitude_precision, latitude_precision)
  }

  # return results (metadata?)
  return(data)
}



#' @export
#' @rdname initial_cleaning
#' @usage
#' sort_columns(data, species, x, y, keep_all_columns = FALSE)

sort_columns <- function(data, species, x, y, keep_all_columns = FALSE) {
  # error checking
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (missing(species)) {
    stop("Argument 'species' must be defined.")
  }
  if (missing(x)) {
    stop("Argument 'x' must be defined.")
  }
  if (missing(y)) {
    stop("Argument 'y' must be defined.")
  }
  if (class(data)[1] != "data.frame") {
    stop("'data' must be of class 'data.frame'.")
  }

  # format data
  ## required columns
  required <- c(species, x, y)

  ## returning formatted columns
  if (keep_all_columns == TRUE) {
    allcol <- sort(colnames(data))
    retcol <- c(required, allcol[!allcol %in% required])

    return(data[, retcol])
  } else {
    return(data[, required])
  }
}



#' @export
#' @rdname initial_cleaning
#' @param columns (character) vector of additional column name(s) in
#' \code{data} to be considered while removing missing or duplicate records,
#' default = NULL.
#' @usage
#' remove_missing(data, columns = NULL, remove_na = TRUE,
#'                remove_empty = TRUE, keep_all_columns = TRUE)

remove_missing <- function(data, columns = NULL, remove_na = TRUE,
                           remove_empty = TRUE, keep_all_columns = TRUE) {
  # error checking
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (class(data)[1] != "data.frame") {
    stop("'data' must be of class 'data.frame'.")
  }

  clnames <- colnames(data)
  if (!is.null(columns)) {
    if (all(!columns %in% clnames)) {
      stop("None of the 'columns' is in 'data'.")
    }
    if (any(!columns %in% clnames)) {
      columns <- columns[columns %in% clnames]
      warning("Some 'columns' were not in 'data'.")
    }
  }

  # conditions
  if (is.null(columns)) {
    columns <- clnames
  }

  # rows to be removed
  if (remove_na == TRUE) {
    nas <- which(is.na(data[, columns]), arr.ind = TRUE)
  }

  if (remove_empty == TRUE) {
    miss <- which(data[, columns] == "", arr.ind = TRUE)
  }

  to_remove <- unique(c(nas[, 1], miss[, 1]))

  # returning results (metadata?)
  if (length(to_remove) >= 1) {
    if (keep_all_columns == TRUE) {
      return(data[-to_remove, ])
    } else {
      return(data[-to_remove, columns])
    }
  } else {
    if (keep_all_columns == TRUE) {
      return(data)
    } else {
      return(data[, columns])
    }
  }
}



#' @export
#' @rdname initial_cleaning
#' @usage
#' remove_duplicates(data, columns = NULL, keep_all_columns = TRUE)

remove_duplicates <- function(data, columns = NULL, keep_all_columns = TRUE) {
  # error checking
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (class(data)[1] != "data.frame") {
    stop("'data' must be of class 'data.frame'.")
  }

  clnames <- colnames(data)
  if (!is.null(columns)) {
    if (all(!columns %in% clnames)) {
      stop("None of the 'columns' is in 'data'.")
    }
    if (any(!columns %in% clnames)) {
      columns <- columns[columns %in% clnames]
      warning("Some 'columns' were not in 'data'.")
    }
  }

  # conditions
  if (is.null(columns)) {
    columns <- clnames
  }

  # exclude duplicates
  if (length(columns) < 2) {
    data <- data[!duplicated(data[, columns]), ]
  } else {
    data <- data[!duplicated(do.call(paste, data[, columns])), ]
  }

  # returning results (metadata?)
  if (keep_all_columns == TRUE) {
    return(data)
  } else {
    return(data[, columns])
  }
}



#' @export
#' @rdname initial_cleaning
#' @usage
#' remove_corrdinates_00(data, x, y)

remove_corrdinates_00 <- function(data, x, y) {
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
  if (class(data)[1] != "data.frame") {
    stop("'data' must be of class 'data.frame'.")
  }

  # filter
  to_remove <- data[, x] == 0 & data[, y] == 0

  # return result
  return(data[!to_remove, ])
}




#' @export
#' @rdname initial_cleaning
#' @usage
#' filter_decimal_precision(data, x,
#'                          y, decimal_precision = 0,
#'                          longitude_precision = NULL,
#'                          latitude_precision = NULL)

filter_decimal_precision <- function(data, x, y, decimal_precision = 0,
                                     longitude_precision = NULL,
                                     latitude_precision = NULL) {
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
  if (class(data)[1] != "data.frame") {
    stop("'data' must be of class 'data.frame'.")
  }

  # conditions
  longitude_precision <- ifelse(is.null(longitude_precision), decimal_precision,
                                longitude_precision)
  latitude_precision <- ifelse(is.null(latitude_precision), decimal_precision,
                               latitude_precision)

  # filter
  lon_dec <- vapply(data[, x], decimal_places, numeric(1))
  lat_dec <- vapply(data[, y], decimal_places, numeric(1))

  to_remove <- unique(c(which(lon_dec < longitude_precision),
                        which(lat_dec < latitude_precision)))

  # return results (metadata?)
  if (length(to_remove) >= 1) {
    return(data[-to_remove, ])
  } else {
    return(data)
  }
}



# helper to get decimal places
decimal_places <- function(x) {
  if (missing(x)) {stop("Argument 'x' must be defined.")}
  if (abs(x - round(x)) > (.Machine$double.eps^0.5)) {
    nchar(strsplit(sub("0+$", "", as.character(x)), ".",
                   fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
