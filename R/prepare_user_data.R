#' Prepare data for model calibration with user-prepared calibration data
#'
#' @description
#' This function prepares data for model calibration using user-prepared
#' calibration data. It includes optional PCA, k-fold partitioning, and the
#' creation of a grid parameter combinations, including distinct regularization
#' multiplier values, various feature classes, and different sets of
#' environmental variables.
#'
#' @usage
#' prepare_user_data(algorithm, user_data, pr_bg, species = NULL, x = NULL,
#'                   y = NULL, features = c("lq", "lqp"),
#'                   r_multiplier = c(0.1, 0.5, 1, 2, 3),
#'                   partition_method = "subsample", n_partitions = 4,
#'                   train_proportion = 0.7, user_part = NULL,
#'                   categorical_variables = NULL,
#'                   do_pca = FALSE, center = TRUE, scale = TRUE,
#'                   exclude_from_pca = NULL, variance_explained = 95,
#'                   min_explained = 5, min_number = 2, min_continuous = NULL,
#'                   weights = NULL, include_xy = TRUE, write_pca = FALSE,
#'                   pca_directory = NULL, write_file = FALSE, file_name = NULL,
#'                   seed = 1)
#'
#' @param algorithm (character) modeling algorithm, either "glm" or "maxnet".
#' @param user_data (data frame) A data.frame with a column with presence (1)
#' and background (0) records, together with variable values (one variable per
#' column). See an example with \code{data("user_data", package = "kuenm2")}.
#' @param pr_bg (character) the name of the column in `user_data` that contains
#' the presence/background records.
#' @param species (character) string specifying the species name (optional).
#' Default is NULL.
#' @param x (character) a string specifying the name of the column in `user_data`
#' that contains the longitude values. Default is NULL. Must be defined if
#' present in `user_data` otherwise it will be considered as another predictor
#' variable.
#' @param y (character) a string specifying the name of the column in `user_data`
#' that contains the latitude values. Default is NULL. Must be defined if
#' present in `user_data` otherwise it will be considered as another predictor
#' variable.
#' @param categorical_variables (character) names of the variables that are
#' categorical. Default is NULL.
#' @param do_pca (logical) whether to perform a principal component analysis
#' (PCA) with the set of variables. Default is FALSE.
#' @param exclude_from_pca (character) variable names within raster_variables
#' that should not be included in the PCA transformation. Instead, these
#' variables will be added directly to the final set of output variables
#' without being modified. The default is NULL, meaning all variables
#' will be used unless specified otherwise.
#' @param variance_explained (numeric) the cumulative percentage of total
#' variance that must be explained by the selected principal components.
#' Default is 95.
#' @param min_explained (numeric) the minimum percentage of total variance that
#' a principal component must explain to be retained. Default is 5.
#' @param center (logical) whether the variables should be zero-centered. Default
#' is TRUE.
#' @param scale (logical) whether the variables should be scaled to have unit
#' variance before the analysis takes place. Default is FALSE.
#' @param write_pca (logical) whether to save the PCA-derived raster layers
#' (principal components) to disk. Default is FALSE.
#' @param pca_directory (character) the path or name of the folder where the PC
#' raster layers will be saved. This is only applicable if `write_pca = TRUE`.
#' Default is NULL.
#' @param user_part a user provided list with replicates or folds for
#' cross-validation to be used in model calibration. Each element of the list
#' should contain a vector of indices indicating the test points, which will be
#' used to split `use_data` into training and testing sets.
#' @param weights (numeric) a numeric vector specifying weights for the
#' occurrence records. Default is NULL.
#' @param min_number (numeric) the minimum number of variables to be included in
#' the model formulas to be generated.
#' @param min_continuous (numeric) the minimum number of continuous variables
#' required in a combination. Default is NULL.
#' @param features (character) a vector of feature classes. Default is c("q",
#' "lq", "lp", "qp", "lqp").
#' @param r_multiplier (numeric) a vector of regularization parameters for maxnet.
#' Default is c(0.1, 1, 2, 3, 5).
#' @param partition_method (character) method used for data partitioning.
#' Available options are `"kfolds"`, `"subsample"`, and `"bootstrap"`.
#' See **Details** for more information.
#' @param n_partitions (numeric) number of partitions to generate. If
#' `partition_method` is `"subsample"` or `"bootstrap"`, this defines the number
#' of partitions. If `"kfolds"`, it specifies the number of folds. Default is 4.
#' @param train_proportion (numeric) proportion of occurrence and background
#' points to be used for model training in each partition. Only applicable when
#'  `partition_method` is `"subsample"` or `"bootstrap"`. Default is 0.7 (i.e.,
#'  70% for training and 30% for testing).
#' @param include_xy (logical) whether to include the coordinates (longitude and
#' latitude) in the results from preparing data. Default is TRUE.
#' @param write_file (logical) whether to write the resulting prepared_data list
#' in a local directory. Default is FALSE.
#' @param file_name (character) the path or name of the folder where the
#' resulting list will be saved. This is only applicable if `write_file =
#' TRUE`. Default is NULL.
#' @param seed (numeric) integer value to specify an initial seed to split the
#' data. Default is 1.
#'
#' @details
#' The available data partitioning methods are:
#'
#' - **"kfolds"**: Splits the dataset into *K* subsets (folds) of approximately equal size. In each partition, one fold is used as the test set, while the remaining folds are combined to form the training set.
#' - **"bootstrap"**: Creates the training dataset by sampling observations from the original dataset *with replacement* (i.e., the same observation can be selected multiple times). The test set consists of the observations that were not selected in that specific replicate.
#' - **"subsample"**: Similar to bootstrap, but the training set is created by sampling *without replacement* (i.e., each observation is selected at most once). The test set includes the observations not selected for training.
#' - **"leave-one-out"**: A special case of kfolds where the number of folds equals the number of presence records. In each replicate, a single presence is left out to serve as the test set, while the remaining observations are used for training.
#'
#' @return
#' An object of class `prepared_data` containing all elements to run a model
#' calibration routine. The elements include: species, calibration data,
#' a grid of model parameters, indices of k-folds for cross validation,
#' xy coordinates, names of continuous and categorical variables, weights,
#' results from PCA, and modeling algorithm.
#'
#' @importFrom enmpa aux_var_comb
#' @importFrom terra crop prcomp extract as.data.frame nlyr
#' @importFrom utils combn
#' @importFrom stats prcomp predict
#'
#' @export
#'
#' @examples
#' # Import user-prepared data
#' data("user_data", package = "kuenm2")
#'
#' # Prepare data for maxnet model
#' maxnet_swd_user <- prepare_user_data(algorithm = "maxnet",
#'                                      user_data = user_data, pr_bg = "pr_bg",
#'                                      species = "Myrcia hatschbachii",
#'                                      categorical_variables = "SoilType",
#'                                      features = c("l", "q", "p", "lq", "lqp"),
#'                                      r_multiplier = c(0.1, 1, 2, 3, 5))
#' maxnet_swd_user
#'
#' # Prepare data for glm model
#' glm_swd_user <- prepare_user_data(algorithm = "glm",
#'                                   user_data = user_data, pr_bg = "pr_bg",
#'                                   species = "Myrcia hatschbachii",
#'                                   categorical_variables = "SoilType",
#'                                   features = c("l", "q", "p", "lq", "lqp"))
#' glm_swd_user

prepare_user_data <- function(algorithm,
                              user_data,
                              pr_bg,
                              species = NULL,
                              x = NULL,
                              y = NULL,
                              features = c("lq", "lqp"),
                              r_multiplier = c(0.1, 0.5, 1, 2, 3),
                              partition_method = "subsample",
                              n_partitions = 4,
                              train_proportion = 0.7,
                              user_part = NULL,
                              categorical_variables = NULL,
                              do_pca = FALSE,
                              center = TRUE,
                              scale = TRUE,
                              exclude_from_pca = NULL,
                              variance_explained = 95,
                              min_explained = 5,
                              min_number = 2,
                              min_continuous = NULL,
                              weights = NULL,
                              include_xy = TRUE,
                              write_pca = FALSE,
                              pca_directory = NULL,
                              write_file = FALSE,
                              file_name = NULL,
                              seed = 1) {
  #Check data
  if (missing(algorithm)) {
    stop("Argument 'algorithm' must be defined.")
  }
  if (missing(user_data)) {
    stop("Argument 'user_data' must be defined.")
  }
  if (missing(pr_bg)) {
    stop("Argument 'pr_bg' must be defined.")
  }

  if (length(algorithm) != 1) {
    stop("'algorithm' must be a single value: either 'glm' or 'maxnet'.")
  } else {
    if (!algorithm %in% c("glm", "maxnet")) {
      stop("'algorithm' must be 'glm' or 'maxnet'.")
    }
  }

  if (!inherits(user_data, "data.frame")) {
    stop("Argument 'user_data' must be a 'data.frame'.")
  }

  if (!is.null(species) & !inherits(species, "character")) {
    stop("Argument 'species' must be a character.")
  }

  if (!is.null(x) | !is.null(y)) {
    if (!inherits(x, "character") | !inherits(y, "character")) {
      stop("Arguments 'x' and 'y' must be of class character.")

    if (!x %in% colnames(user_data) | !y %in% colnames(user_data)) {
        stop("'x' and/or 'y' are not columns in 'user_data'.")
      }
  }}


  if (!is.null(categorical_variables) &
      !inherits(categorical_variables, "character")) {
    stop("Argument 'categorical_variables' must be a 'character' vector.")
  }

  if (!is.null(exclude_from_pca) & !inherits(exclude_from_pca, "character")) {
    stop("Argument 'exclude_from_pca' must be a 'character'.")
  }

  if (write_pca & is.null(pca_directory)) {
    stop("If 'write_pca' = TRUE, 'pca_directory' must be specified.")
  }

  if (write_pca & !inherits(pca_directory, "character")) {
    stop("Argument 'pca_directory' must be a 'character'.")
  }

  if (!is.null(weights) && nrow(user_data) != length(weights)) {
    stop("Length of 'weights' must match the number of occurrences in 'user_data'.")
  }

  check_features <- unique(unlist(strsplit(features, "")))
  check_features <- setdiff(check_features, c("l", "q", "p", "h", "t"))

  if (length(check_features) > 0) {
    stop("'features' defined are not valid:\n",
         "'features' can be l, q, p, t, h, or combinations of these.")
  }

  if (write_file & is.null(file_name)) {
    stop("If 'write_file' = TRUE, 'file_name' must be specified.")
  }

  if (write_file & !is.character(file_name)) {
    stop("Argument 'file_name' must be a 'character'.")
  }

  #Remove NAs
  user_data <- handle_missing_data(user_data, weights)

  #Change name of pr_bg column if necessary
  colnames(user_data)[which(colnames(user_data) == pr_bg)] <- "pr_bg"

  #Get variables
  vars <- setdiff(colnames(user_data), c("pr_bg", x, y))
  user_vars <- user_data[, vars]

  if (!is.null(categorical_variables)) {
    if (any(!categorical_variables %in% vars)) {
      stop("Defined 'categorical_variables' must be a column in 'user_data'.")
    }
    #Make sure it is a factor
    user_data[, categorical_variables] <- as.factor(user_data[, categorical_variables])
    continuous_variable_names <- setdiff(vars, categorical_variables)
  } else {
    continuous_variable_names <- vars
  }

  sp_name <- species

  if (do_pca) {
    if (!is.null(categorical_variables)) {
      exclude_from_pca <- c(categorical_variables, exclude_from_pca)
    } else {
      exclude_from_pca <- exclude_from_pca
    }

  #Get vars to do pca
    vars_in <- setdiff(continuous_variable_names, exclude_from_pca)


    pca <- stats::prcomp(user_vars[,vars_in], center = center, scale = scale)
    pca$x <- NULL
    pca$vars_in <- vars_in
    pca$vars_out <- exclude_from_pca
    d_exp <- cumsum(pca$sdev^2/sum(pca$sdev^2)) * 100
    d_exp <- d_exp[(pca$sdev^2/sum(pca$sdev^2) * 100) > min_explained]
    ind_exp <- if (max(d_exp) > variance_explained) {
      min(which(d_exp >= variance_explained))
    } else {
      length(d_exp)
    }

    pca_variables <- stats::predict(pca, user_vars[,vars_in])[, 1:ind_exp]
    if (ind_exp == 1) {
      pca_variables <- data.frame("PC1" = pca_variables)
    }

    pca$projection_directory <- NULL #Remove projection directory
    #Update user data
    user_data <- cbind(user_data[, c("pr_bg", x, y, exclude_from_pca)],
                       pca_variables)

  } else {
    pca <- NULL
  }

  #Partition
  if (is.null(user_part)) {
    #Partitione data
    pd <- part_data(data = user_data, pr_bg = "pr_bg",
                    train_proportion = train_proportion,
                    n_partitions = n_partitions,
                    partition_method = partition_method,
                    seed = seed)
  } else {
    pd <- user_part
    partition_method <- "User defined"
  }

  #Check min_number and min_continuous
  if (min_number > ncol(user_data) - 1) {
    stop("'min_number' can't be greater than the number of variables in 'user_data'.")
  }

  if (!is.null(min_continuous)) {
    c_var <- setdiff(colnames(user_data), c("pr_bg", categorical_variables))
    if (!is.numeric(min_continuous)) {
      stop("Argument 'min_continuous' must be 'numeric'.")
    } else {
      if (min_continuous > length(c_var)) {
        stop("'min_continuous' can't be greater than the number of continuous variables in 'user_data'.")
      }
    }
  }

  #Formula grid
  formula_grid <- calibration_grid(user_data, min_number, min_continuous,
                                   categorical_var = categorical_variables,
                                   features, algorithm, r_multiplier)

  #Extract xy?
  if (!is.null(x) & !is.null(y)) {
    data_xy <- data.frame(x = x, t = y)
  } else {
    data_xy <- NULL
  }


  #Prepare final data
  if(partition_method %in% c("kfolds", "leave-one-out")){
    train_proportion <- NULL
  }


  data <- new_prepared_data(species = species, calibration_data = user_data,
                            formula_grid = formula_grid,
                            part_data = pd, partition_method = partition_method,
                            n_partitions = n_partitions,
                            train_proportion = train_proportion,
                            data_xy = data_xy,
                            continuous_variables = continuous_variable_names,
                            categorical_variables = categorical_variables,
                            weights = weights, pca = pca$pca,
                            algorithm = algorithm)


  if (write_file) {
    saveRDS(data, paste0(file_name, ".RDS"))
  }

  return(data)
}
