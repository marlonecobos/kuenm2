% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partition_data.R
\name{block}
\alias{block}
\alias{assign_block}
\title{Partition of data according to geographic blocks}
\usage{
block(data, longitude_column, latitude_column, n_columns, n_rows = NULL)

assign_block(block_limits, new_data, longitude_column, latitude_column,
             return_all = TRUE)
}
\arguments{
\item{data}{data.frame with occurrence records.}

\item{longitude_column}{(character) name of the column in \code{data}
containing longitude values.}

\item{latitude_column}{(character) name of the column in \code{data}
containing latitude values.}

\item{n_columns}{(numeric) number of columns to be considered for the matrix
that generates blocks.}

\item{n_rows}{(numeric) number of rows to be considered for the matrix
that generates blocks. The default, NULL, uses value in \code{n_columns}.}

\item{block_limits}{data.frame with details used to assign blocks in
previous processes.}

\item{new_data}{data.frame containing new data for blocks to be assigned.}

\item{return_all}{(logical) whether to return the complete resulting
data.frame. Default = TRUE. FALSE returns only block numbers corresponding
to each row of \code{new_data}.}
}
\value{
A list of results including:
\itemize{
\item blocked_data.- original data.frame with an extra column indicating block
per each record.
\item block_limits.- data.frame indicating limits used to assign blocks in data.
}
}
\description{
analysis to assign geographic blocks to occurrence data.
}
\details{
Blocks are assigned using a matrix produced based on the values defined in
\code{n_columns} and \code{n_rows}. The algorithm creates the matrix of
blocks so all blocks have similar number of occurrences.
}
