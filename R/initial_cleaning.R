# data cleaning helper functions

initial_cleaning <- function(data, species_column, longitude_column,
                             latitude_column, other_columns = NULL,
                             keep_all_columns = TRUE,
                             sort_columns = TRUE, remove_na = TRUE,
                             remove_empty = TRUE, remove_duplicates = TRUE,
                             by_decimal_precission = TRUE,
                             decimal_precision = 0, longitude_precision = NULL,
                             latitude_precision = NULL) {

  # error checking

  # preparing arguments
  if (!is.null(other_columns)) {
    columns <- c(species_column, longitude_column, latitude_column,
                 other_columns)
  } else {
    columns <- c(species_column, longitude_column, latitude_column)
  }

  # cleaning steps
  if (sort_columns == TRUE) {
    data <- sort_columns(data, species_column, longitude_column,
                         latitude_column, keep_all_columns)
  }

  if (remove_na == TRUE | remove_empty == TRUE) {
    data <- remove_missing(data, columns, remove_na, remove_empty,
                           keep_all_columns)
  }

  if (remove_duplicates == TRUE) {
    data <- remove_duplicates(data, columns, keep_all_columns)
  }

  if (by_decimal_precission == TRUE) {
    data <- filter_decimal_precision(data, longitude_column,
                                     latitude_column, decimal_precision,
                                     longitude_precision, latitude_precision)
  }

  # return results (metadata?)
  return(data)
}

# ordering columns
sort_columns <- function(data, species_column, longitude_column,
                          latitude_column, keep_all_columns = FALSE) {
  # error checking

  # format data
  ## required columns
  required <- c(species_column, longitude_column, latitude_column)

  ## returning formatted columns
  if (keep_all_columns == TRUE) {
    allcol <- colnames(data)
    retcol <- c(required, allcol[!allcol %in% required])

    return(data[, retcol])
  } else {
    return(data[, required])
  }
}


# remove missing data
remove_missing <- function(data, columns = NULL, remove_na = TRUE,
                           remove_empty = TRUE, keep_all_columns = TRUE) {
  # error checking
  ## missing
  ## columns in data

  # conditions
  if (is.null(columns)) {
    columns <- colnames(data)
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
  if (keep_all_columns == TRUE) {
    return(data[to_remove, ])
  } else {
    return(data[!to_remove, columns])
  }
}



# remove duplicates
remove_duplicates <- function(data, columns = NULL, keep_all_columns = TRUE) {
  # error checking
  ## missing
  ## is data.frame
  ## columns in data

  # conditions
  if (is.null(columns)) {
    columns <- colnames(data)
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



# remove 0, 0
remove_corrdinates_00 <- function(data, longitude_column, latitude_column) {
  # error checking

  # filter
  to_remove <- data[, longitude_column] == 0 & data[, latitude_column] == 0

  # return result
  return(data[!to_remove, ])
}



# filter by coordinate precision
filter_decimal_precision <- function(data, longitude_column,
                                     latitude_column, decimal_precision = 0,
                                     longitude_precision = NULL,
                                     latitude_precision = NULL) {
  # error checking

  # conditions
  longitude_precision <- ifelse(is.null(longitude_precision), decimal_precision,
                                longitude_precision)
  latitude_precision <- ifelse(is.null(latitude_precision), decimal_precision,
                               latitude_precision)

  # filter
  lon_dec <- vapply(data[, longitude_column], decimal_places, numeric(1))
  lat_dec <- vapply(data[, latitude_column], decimal_places, numeric(1))

  to_remove <- unique(c(which(lon_dec < longitude_precision),
                        which(lat_dec < latitude_precision)))

  # return results (metadata?)
  return(data <- data[!to_remove, ])
}



# helper to get decimal places
decimal_places <- function(x) {
  if (abs(x - round(x)) > (.Machine$double.eps^0.5)) {
    nchar(strsplit(sub("0+$", "", as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
