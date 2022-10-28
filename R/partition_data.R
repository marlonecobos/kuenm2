#' Partition of data according to geographic blocks
#'
#' @description analysis to assign geographic blocks to occurrence data.
#'
#' @param data data.frame with occurrence records.
#' @param longitude_column (character) name of the column in \code{data}
#' containing longitude values.
#' @param latitude_column (character) name of the column in \code{data}
#' containing latitude values.
#' @param n_columns (numeric) number of columns to be considered for the matrix
#' that generates blocks.
#' @param n_rows (numeric) number of rows to be considered for the matrix
#' that generates blocks. The default, NULL, uses value in \code{n_columns}.
#'
#' @details
#' Blocks are assigned using a matrix produced based on the values defined in
#' \code{n_columns} and \code{n_rows}. The algorithm creates the matrix of
#' blocks so all blocks have similar number of occurrences.
#'
#' @return
#' A list of results including:
#' - blocked_data.- original data.frame with an extra column indicating block
#' per each record.
#' - block_limits.- data.frame indicating limits used to assign blocks in data.
#'
#' @export
#' @importFrom stats quantile
#' @rdname block
#'
#' @usage
#' block(data, longitude_column, latitude_column, n_columns, n_rows = NULL)


block <- function(data, longitude_column, latitude_column,
                  n_columns, n_rows = NULL) {
  # error checking
  if (missing(data)) {
    stop("Argument 'data' needs to be defined.")
  }
  if (missing(longitude_column) | missing(latitude_column)) {
    stop("'longitude_column' and 'latitude_column' must be defined.")
  }
  if (missing(n_columns)) {
    stop("Argument 'n_columns' must be defined.")
  }
  if (is.null(n_rows)) {
    n_rows <- n_columns
  }
  coln <- colnames(data)
  if (!longitude_column %in% coln) {
    stop(longitude_column, " is not a column in 'data'.")
  }
  if (!latitude_column %in% coln) {
    stop(latitude_column, " is not a column in 'data'.")
  }

  # detecting ranges and intervals of blocks
  ## column for block id
  data <- cbind(data, block = NA)

  ## longitude
  xlb <- seq(0, 1, (1 / n_columns))
  q1 <- quantile(data[, longitude_column], xlb)
  xlims <- data.frame(x_start = c(q1[1], q1[2:n_columns] + .Machine$double.eps),
                      x_end = q1[-1], row.names = NULL)

  blocking_limits <- lapply(1:n_columns, function(x) {
    ## data by longitude blocks
    xid <- which(data[, longitude_column] >= xlims[x, 1] &
                   data[, longitude_column] <= xlims[x, 2])

    pd <- data[xid, ]

    ## latitude block assignment
    ylb <-  seq(0, 1, (1 / n_rows))
    q2 <- quantile(pd[, latitude_column], ylb)
    yli <- data.frame(y_start = c(q2[1], q2[2:n_rows] + .Machine$double.eps),
                        y_end = q2[-1], row.names = NULL)

    for (y in 1:n_rows) {
      yid <- which(pd[, latitude_column] >= yli[y, 1] &
                     pd[, latitude_column] <= yli[y, 2])
      nb <- ifelse(x == 1, y, ((x - 1) * n_rows + y))

      pd[yid, ncol(pd)] <- rep(nb, length(yid))
    }

    return(list(latitude_limits = yli, blocked = pd))
  })

  # assigning block numbers
  data <- do.call(rbind, lapply(blocking_limits, function(x) {x$blocked}))


  # all limits
  ylim_list <- lapply(blocking_limits, function(x) {x$latitude_limits})
  ylims <- do.call(rbind, ylim_list)

  bl <- data.frame(xlims[expand.grid(1:n_rows, 1:n_columns)[, 2], ],
                   ylims, block = sort(unique(data$block)), row.names = NULL)


  # returning results
  return(list(blocked_data = data, block_limits = bl))
}



#' @param block_limits data.frame with details used to assign blocks in
#' previous processes.
#' @param new_data data.frame containing new data for blocks to be assigned.
#' @param return_all (logical) whether to return the complete resulting
#' data.frame. Default = TRUE. FALSE returns only block numbers corresponding
#' to each row of \code{new_data}.
#'
#' @rdname block
#' @export
#' @usage
#' assign_block(block_limits, new_data, longitude_column, latitude_column,
#'              return_all = TRUE)

assign_block <- function(block_limits, new_data, longitude_column,
                         latitude_column, return_all = TRUE) {
  # error checking
  if (missing(block_limits)) {
    stop("Argument 'block_limits' needs to be defined.")
  }
  if (missing(new_data)) {
    stop("Argument 'new_data' needs to be defined.")
  }
  if (missing(longitude_column)) {
    stop("Argument 'longitude_column' needs to be defined.")
  }
  if (missing(latitude_column)) {
    stop("Argument 'latitude_column' needs to be defined.")
  }
  coln <- colnames(data)
  if (!longitude_column %in% coln) {
    stop(longitude_column, " is not a column in 'data'.")
  }
  if (!latitude_column %in% coln) {
    stop(latitude_column, " is not a column in 'data'.")
  }

  # fixing x and y limits
  ## limits of new data and original data combined
  new_lims <- apply(new_data[, c(longitude_column, latitude_column)], 2, range)

  old_lims <- cbind(c(min(block_limits$x_start), max(block_limits$x_end)),
                    c(min(block_limits$y_start), max(block_limits$y_end)))

  colnames(old_lims) <- c(longitude_column, latitude_column)

  new_lims <- apply(rbind(new_lims, old_lims), 2, range)

  ## fixing
  block_limits[block_limits[, 1] == old_lims[1, 1], 1] <- new_lims[1, 1]
  block_limits[block_limits[, 2] == old_lims[2, 1], 2] <- new_lims[2, 1]

  ulongs <- unique(block_limits$x_start)

  for (i in ulongs) {
    ix <- block_limits[, 1] == i

    ylm <- range(c(block_limits[ix, 3:4]))

    block_limits[ix & block_limits[, 3] == ylm[ 1], 3] <- new_lims[1, 2]
    block_limits[ix & block_limits[, 4] == ylm[2], 4] <- new_lims[2, 2]
  }


  # assigning blocks
  ## new column for blocks
  new_data <- cbind(new_data, block = NA)

  ## per block
  for (i in 1:nrow(block_limits)) {
    ids <- which(new_data[, longitude_column] >= block_limits[i, 1] &
                   new_data[, longitude_column] <= block_limits[i, 2] &
                   new_data[, latitude_column] >= block_limits[i, 3] &
                   new_data[, latitude_column] <= block_limits[i, 4])

    new_data[ids, "block"] <- i
  }

  # return results
  if (return_all == TRUE) {
    return(new_data)
  } else {
    return(new_data$block)
  }
}
