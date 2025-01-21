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
#' prepare_user_data(algorithm, user_data, pr_bg = "pr_bg", species = NULL,
#'                   x = NULL,y = NULL, categorical_variables = NULL,
#'                   do_pca = FALSE, exclude_from_pca = NULL,
#'                   variance_explained = 95, min_explained = 5, center = TRUE,
#'                   scale = FALSE, write_pca = FALSE, pca_directory = NULL,
#'                   kfolds = 4, weights = NULL, min_number = 2,
#'                   min_continuous = NULL,
#'                   features = c("q", "lq", "lp", "qp", "lqp"),
#'                   reg_mult = c(0.1, 0.5, 1, 2, 3), include_xy = TRUE,
#'                   write_file = FALSE, file_name = NULL, seed = 1)
#'
#' @param algorithm (character) type algorithm, either "glm" or "glmnet".
#' @param user_data (data frame) A data.frame with columns indicating presence (1)
#' and background (0) data, as well as the corresponding predictor values. For
#' an example, see \code{data("user_data", package = "kuenm2")}.
#' @param pr_bg (character) a string specifying the name of the column in
#' user_data that contains the presence/background.
#' @param species (character) string specifying the species name (optional).
#' Default is NULL.
#' @param x (character) a string specifying the name of the column in `occ` that
#' contains the longitude values.
#' @param y (character) a string specifying the name of the column in `occ` that
#' contains the latitude values.
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
#' @param kfolds (numeric) the number of groups (folds) the occurrence
#' data will be split into for cross-validation. Default is 4.
#' @param weights (numeric) a numeric vector specifying weights for the
#' occurrence records. Default is NULL.
#' @param min_number (numeric) the minimum number of variables to be included in
#' the model formulas to be generated.
#' @param min_continuous (numeric) the minimum number of continuous variables
#' required in a combination. Default is NULL.
#' @param features (character) a vector of feature classes. Default is c("q",
#' "lq", "lp", "qp", "lqp").
#' @param reg_mult (numeric) a vector of regularization parameters for glmnet.
#' Default is c(0.1, 1, 2, 3, 5).
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
#' @return
#' An object of class `prepared_data` containing all elements to run a model
#' calibration routine. The elements include: species, calibration data,
#' a grid of model parameters, indices of k-folds for cross validation,
#' xy coordinates, names of continuous and categorical variables, weights,
#' results from PCA, and modeling algorithm.
#'
#' @importFrom enmpa kfold_partition aux_var_comb
#' @importFrom terra crop prcomp extract as.data.frame nlyr
#' @importFrom utils combn
#' @importFrom stats prcomp predict
#'
#' @export
#'
#' @examples
#' # Import user-prepared calibration data
#' data("user_data", package = "kuenm2")
#' # Prepare data for glmnet model
#' glmnet_swd_user <- prepare_user_data(algorithm = "glmnet", user_data = user_data,
#'                                      species = "Myrcia hatschbachii",
#'                                      categorical_variables = "SoilType",
#'                                      features = c("l", "q", "p", "lq", "lqp"),
#'                                      reg_mult = c(0.1, 1, 2, 3, 5))
#' print(glmnet_swd_user)
#'
#' # Prepare data for glm model
#' glm_swd_user <- prepare_user_data(algorithm = "glm", user_data = user_data,
#'                                   species = "Myrcia hatschbachii",
#'                                   categorical_variables = "SoilType",
#'                                   features = c("l", "q", "p", "lq", "lqp"),
#'                                   reg_mult = c(0.1, 1, 2, 3, 5))
#' print(glm_swd_user)



prepare_user_data <- function(algorithm, user_data, pr_bg = "pr_bg",
                         species = NULL, x = NULL,y = NULL,
                         categorical_variables = NULL, do_pca = FALSE,
                         exclude_from_pca = NULL, variance_explained = 95,
                         min_explained = 5, center = TRUE, scale = FALSE,
                         write_pca = FALSE, pca_directory = NULL,
                         kfolds = 4, weights = NULL,
                         min_number = 2, min_continuous = NULL,
                         features = c("q", "lq", "lp", "qp", "lqp"),
                         reg_mult = c(0.1, 0.5, 1, 2, 3), include_xy = TRUE,
                         write_file = FALSE, file_name = NULL, seed = 1) {
  #Check data
  if(length(algorithm) != 1) {
    stop("'algorithm' must be a single value: either 'glm' or 'glmnet'")
  } else {
    if (!algorithm %in% c("glm", "glmnet")) {
      stop("'algorithm' must be 'glm' or 'glmnet'")
    }
  }

  if (!inherits(user_data, "data.frame")) {
    stop("Argument 'user_data' must be a 'data.frame'")
  }

  if (!is.null(species) & !inherits(species, "character")) {
    stop("Argument 'species' must be a character")
  }

  if(!is.null(x) | !is.null(y)){
    if (!inherits(x, "character") | !inherits(y, "character")) {
      stop("Arguments 'x' and 'y' must be of class character")

    if(!x %in% colnames(occ) | !y %in% colnames(occ)) {
        stop("'x' and/or 'y' are not column names in 'occ'")
      }
  }}


  if (!is.null(categorical_variables) &
      !inherits(categorical_variables, "character")) {
    stop("Argument 'categorical_variables' must be a 'character' vector")
  }

  if (!is.null(exclude_from_pca) & !inherits(exclude_from_pca, "character")) {
    stop("Argument 'exclude_from_pca' must be a 'character'")
  }

  if (write_pca & is.null(pca_directory)) {
    stop("If 'write_pca' = TRUE, 'pca_directory' must be specified")
  }

  if (write_pca & !inherits(pca_directory, "character")) {
    stop("Argument 'pca_directory' must be a 'character'")
  }

  if (!is.null(weights) && nrow(occ) != length(weights)) {
    stop("Length of 'weights' must match the number of occurrences in 'occ'")
  }

  check_features <- unique(unlist(strsplit(features, "")))
  check_features <- setdiff(check_features, c("l", "q", "p", "h", "t"))

  if(length(check_features) > 0){
    stop(paste(check_features, collapse = ", "), " is not a valid feature\n",
         "Features can be l, q, p, t, h, or combinations of these")
  }

  if (write_file & is.null(file_name)) {
    stop("If 'write_file' = TRUE, 'file_name' must be specified")
  }

  if (write_file & !is.character(file_name)) {
    stop("Argument 'file_name' must be a 'character'")
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
      stop("Defined 'categorical_variables' must be a column in 'user_data'")
    }
    #Make sure it is a factor
    user_data[,categorical_variables] <- as.factor(user_data[,categorical_variables])
    continuous_variable_names <- setdiff(vars,
                                         categorical_variables)
  } else {
    continuous_variable_names <- vars
  }

  sp_name <- species

    if (do_pca) {
    if (!is.null(categorical_variables)){
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
    ind_exp <- if (max(d_exp) > variance_explained) min(which(d_exp >= variance_explained)) else length(d_exp)

    pca_variables <- stats::predict(pca, user_vars[,vars_in])[,1:ind_exp]
    if(ind_exp == 1){
      pca_variables <- data.frame("PC1" = pca_variables)
    }

    pca$projection_directory <- NULL #Remove projection directory
    #Update user data
    user_data <- cbind(user_data[,c("pr_bg", x, y, exclude_from_pca)],
                       pca_variables)

  } else {
    pca <- NULL
  }

  #Partition
  k_f <- enmpa::kfold_partition(data = user_data, dependent = "pr_bg", k = kfolds,
                                seed = seed)

  #Check min_number and min_continuous
  if (min_number > ncol(user_data) - 1) {
    stop("'min_number' can't be greater than the number of variables in 'user_data'")
  }

  if (!is.null(min_continuous)) {
    c_var <- setdiff(colnames(user_data), c("pr_bg", categorical_variables))
    if (!is.numeric(min_continuous)) {
      stop("Argument 'min_continuous' must be 'numeric'")
    } else {
      if (min_continuous > length(c_var)) {
        stop("'min_continuous' can't be greater than the number of continuous variables in 'user_data'")
      }
    }
  }

  #Formula grid
  formula_grid <- calibration_grid(user_data, min_number, min_continuous,
                                   categorical_var = categorical_variables,
                                   features, algorithm, reg_mult)

  #Extract xy?
  if(!is.null(x) & !is.null(y)){
    data_xy <- data.frame(x = x, t = y)
  } else {
    data_xy <- NULL
  }

  data <- new_prepared_data(species, calibration_data = user_data, formula_grid,
                            kfolds = k_f, data_xy = data_xy,
                            continuous_variables = continuous_variable_names,
                            categorical_variables, weights, pca = pca,
                            algorithm)

  if (write_file) {
    saveRDS(data, paste0(file_name, ".RDS"))
  }

  return(data)
}
