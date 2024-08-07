% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_data.R
\name{prepare_data}
\alias{prepare_data}
\title{Prepare data for model calibration}
\usage{
prepare_data(
  occ,
  species = NULL,
  x,
  y,
  spat_variables,
  mask = NULL,
  categorical_variables = NULL,
  do_pca = FALSE,
  deviance_explained = 95,
  min_explainded = 5,
  exclude_from_pca = NULL,
  center = TRUE,
  scale = TRUE,
  write_pca = FALSE,
  output_pca = NULL,
  nbg = 10000,
  kfolds = 4,
  weights = NULL,
  include_xy = TRUE,
  write_files = FALSE,
  file_name = NULL,
  seed = 1
)
}
\description{
Prepare data for model calibration
}
