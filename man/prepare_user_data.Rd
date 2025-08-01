% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_user_data.R
\name{prepare_user_data}
\alias{prepare_user_data}
\title{Prepare data for model calibration with user-prepared calibration data}
\usage{
prepare_user_data(algorithm, user_data, pr_bg, species = NULL, x = NULL,
                  y = NULL, features = c("lq", "lqp"),
                  r_multiplier = c(0.1, 0.5, 1, 2, 3),
                  partition_method = "subsample", n_partitions = 4,
                  train_proportion = 0.7, user_part = NULL,
                  categorical_variables = NULL,
                  do_pca = FALSE, center = TRUE, scale = TRUE,
                  exclude_from_pca = NULL, variance_explained = 95,
                  min_explained = 5, min_number = 2, min_continuous = NULL,
                  weights = NULL, include_xy = TRUE, write_pca = FALSE,
                  pca_directory = NULL, write_file = FALSE, file_name = NULL,
                  seed = 1)
}
\arguments{
\item{algorithm}{(character) modeling algorithm, either "glm" or "maxnet".}

\item{user_data}{(data frame) A data.frame with a column with presence (1)
and background (0) records, together with variable values (one variable per
column). See an example with \code{data("user_data", package = "kuenm2")}.}

\item{pr_bg}{(character) the name of the column in \code{user_data} that contains
the presence/background records.}

\item{species}{(character) string specifying the species name (optional).
Default is NULL.}

\item{x}{(character) a string specifying the name of the column in \code{user_data}
that contains the longitude values. Default is NULL. Must be defined if
present in \code{user_data} otherwise it will be considered as another predictor
variable.}

\item{y}{(character) a string specifying the name of the column in \code{user_data}
that contains the latitude values. Default is NULL. Must be defined if
present in \code{user_data} otherwise it will be considered as another predictor
variable.}

\item{features}{(character) a vector of feature classes. Default is c("q",
"lq", "lp", "qp", "lqp").}

\item{r_multiplier}{(numeric) a vector of regularization parameters for maxnet.
Default is c(0.1, 1, 2, 3, 5).}

\item{partition_method}{(character) method used for data partitioning.
Available options are \code{"kfolds"}, \code{"subsample"}, and \code{"bootstrap"}.
See \strong{Details} for more information.}

\item{n_partitions}{(numeric) number of partitions to generate. If
\code{partition_method} is \code{"subsample"} or \code{"bootstrap"}, this defines the number
of partitions. If \code{"kfolds"}, it specifies the number of folds. Default is 4.}

\item{train_proportion}{(numeric) proportion of occurrence and background
points to be used for model training in each partition. Only applicable when
\code{partition_method} is \code{"subsample"} or \code{"bootstrap"}. Default is 0.7 (i.e.,
70\% for training and 30\% for testing).}

\item{user_part}{a user provided list with replicates or folds for
cross-validation to be used in model calibration. Each element of the list
should contain a vector of indices indicating the test points, which will be
used to split \code{use_data} into training and testing sets.}

\item{categorical_variables}{(character) names of the variables that are
categorical. Default is NULL.}

\item{do_pca}{(logical) whether to perform a principal component analysis
(PCA) with the set of variables. Default is FALSE.}

\item{center}{(logical) whether the variables should be zero-centered. Default
is TRUE.}

\item{scale}{(logical) whether the variables should be scaled to have unit
variance before the analysis takes place. Default is FALSE.}

\item{exclude_from_pca}{(character) variable names within raster_variables
that should not be included in the PCA transformation. Instead, these
variables will be added directly to the final set of output variables
without being modified. The default is NULL, meaning all variables
will be used unless specified otherwise.}

\item{variance_explained}{(numeric) the cumulative percentage of total
variance that must be explained by the selected principal components.
Default is 95.}

\item{min_explained}{(numeric) the minimum percentage of total variance that
a principal component must explain to be retained. Default is 5.}

\item{min_number}{(numeric) the minimum number of variables to be included in
the model formulas to be generated.}

\item{min_continuous}{(numeric) the minimum number of continuous variables
required in a combination. Default is NULL.}

\item{weights}{(numeric) a numeric vector specifying weights for the
occurrence records. Default is NULL.}

\item{include_xy}{(logical) whether to include the coordinates (longitude and
latitude) in the results from preparing data. Default is TRUE.}

\item{write_pca}{(logical) whether to save the PCA-derived raster layers
(principal components) to disk. Default is FALSE.}

\item{pca_directory}{(character) the path or name of the folder where the PC
raster layers will be saved. This is only applicable if \code{write_pca = TRUE}.
Default is NULL.}

\item{write_file}{(logical) whether to write the resulting prepared_data list
in a local directory. Default is FALSE.}

\item{file_name}{(character) the path or name of the folder where the
resulting list will be saved. This is only applicable if \code{write_file = TRUE}. Default is NULL.}

\item{seed}{(numeric) integer value to specify an initial seed to split the
data. Default is 1.}
}
\value{
An object of class \code{prepared_data} containing all elements to run a model
calibration routine. The elements include: species, calibration data,
a grid of model parameters, indices of k-folds for cross validation,
xy coordinates, names of continuous and categorical variables, weights,
results from PCA, and modeling algorithm.
}
\description{
This function prepares data for model calibration using user-prepared
calibration data. It includes optional PCA, k-fold partitioning, and the
creation of a grid parameter combinations, including distinct regularization
multiplier values, various feature classes, and different sets of
environmental variables.
}
\details{
The available data partitioning methods are:
\itemize{
\item \strong{"kfolds"}: Splits the dataset into \emph{K} subsets (folds) of approximately equal size. In each partition, one fold is used as the test set, while the remaining folds are combined to form the training set.
\item \strong{"bootstrap"}: Creates the training dataset by sampling observations from the original dataset \emph{with replacement} (i.e., the same observation can be selected multiple times). The test set consists of the observations that were not selected in that specific replicate.
\item \strong{"subsample"}: Similar to bootstrap, but the training set is created by sampling \emph{without replacement} (i.e., each observation is selected at most once). The test set includes the observations not selected for training.
\item \strong{"leave-one-out"}: A special case of kfolds where the number of folds equals the number of presence records. In each replicate, a single presence is left out to serve as the test set, while the remaining observations are used for training.
}
}
\examples{
# Import user-prepared data
data("user_data", package = "kuenm2")

# Prepare data for maxnet model
maxnet_swd_user <- prepare_user_data(algorithm = "maxnet",
                                     user_data = user_data, pr_bg = "pr_bg",
                                     species = "Myrcia hatschbachii",
                                     categorical_variables = "SoilType",
                                     features = c("l", "q", "p", "lq", "lqp"),
                                     r_multiplier = c(0.1, 1, 2, 3, 5))
maxnet_swd_user

# Prepare data for glm model
glm_swd_user <- prepare_user_data(algorithm = "glm",
                                  user_data = user_data, pr_bg = "pr_bg",
                                  species = "Myrcia hatschbachii",
                                  categorical_variables = "SoilType",
                                  features = c("l", "q", "p", "lq", "lqp"))
glm_swd_user
}
