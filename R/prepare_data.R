#' Prepare data for model calibration
#'
#' @description
#' This function prepares data for model calibration, including optional PCA,
#' background point generation, k-fold partitioning, and the creation of a grid
#' of parameter combinations, including regularization multiplier values,
#' feature classes, and sets of environmental variables.
#'
#' @usage
#' prepare_data(algorithm, occ, x, y, raster_variables, species = NULL,
#'              n_background = 1000, features = c("lq", "lqp"),
#'              r_multiplier = c(0.1, 0.5, 1, 2, 3), partition_method,
#'              n_partitions = 4, train_proportion = 0.7,
#'              categorical_variables = NULL,
#'              do_pca = FALSE, center = TRUE, scale = TRUE,
#'              exclude_from_pca = NULL, variance_explained = 95,
#'              min_explained = 5, min_number = 2, min_continuous = NULL,
#'              bias_file = NULL, bias_effect = NULL, weights = NULL,
#'              include_xy = TRUE, write_pca = FALSE, pca_directory = NULL,
#'              write_file = FALSE, file_name = NULL, seed = 1)
#'
#' @param algorithm (character) modeling algorithm, either "glm" or "maxnet".
#' @param occ (data frame) a data.frame containing the coordinates (longitude
#' and latitude) of the occurrence records.
#' @param x (character) a string specifying the name of the column in `occ` that
#' contains the longitude values.
#' @param y (character) a string specifying the name of the column in `occ` that
#' contains the latitude values.
#' @param raster_variables (SpatRaster) predictor variables from which
#' environmental values will be extracted using `occ` and a background will be
#' sampled. Must correspond geographically with the area where model is
#' calibrated.
#' @param species (character) string specifying the species name (optional).
#' Default is NULL.
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
#' @param n_background (numeric) number of points to represent the background
#' for the model. Default is 1000.
#' @param bias_file (SpatRaster) a raster containing bias values (probability
#' weights) that influence the selection of background points. It must have the
#' same extent, resolution, and number of cells as the raster variables.
#' Default is NULL.
#' @param bias_effect (character) a string specifying how the values in the
#' `bias_file` should be interpreted. Options are "direct" or "inverse". If
#' "direct", higher values in bias file increase the likelihood of selecting
#' background points. If "inverse", higher values decrease the likelihood.
#' Default = NULL. Must be defined if `bias_file` is provided.
#' @param partition_method (character) method used for data partitioning.
#' Available options are `"kfolds"`, `"subsample"`, and `"bootstrap"`.
#' See **Details** for more information.
#' @param n_partitions (numeric) number of partitions to generate. If
#' `partition_method` is `"subsample"` or `"bootstrap"`, this defines the number
#' of replicates. If `"kfolds"`, it specifies the number of folds. Default is 4.
#' @param train_proportion (numeric) proportion of occurrence and background
#' points to be used for model training in each partition. Only applicable when
#'  `partition_method` is `"subsample"` or `"bootstrap"`. Default is 0.7 (i.e.,
#'  70% for training and 30% for testing).
#' @param weights (numeric) a numeric vector specifying weights for the
#' occurrence records. Default is NULL.
#' @param min_number (numeric) the minimum number of variables to be included in
#' the model formulas to be generated. Default = 2.
#' @param min_continuous (numeric) the minimum number of continuous variables
#' required in a combination. Default is NULL.
#' @param features (character) a vector of feature classes. Default is c("q",
#' "lq", "lp", "qp", "lqp").
#' @param r_multiplier (numeric) a vector of regularization parameters for maxnet.
#' Default is c(0.1, 1, 2, 3, 5).
#' @param include_xy (logical) whether to include the coordinates (longitude and
#' latitude) in the results from preparing data. Columns containing coordinates
#' will be renamed as "x" and "y". Default is TRUE.
#' @param write_file (logical) whether to write the resulting prepared_data list
#' in a local directory. Default is FALSE.
#' @param file_name (character) name of file (no extension needed) to write
#' resulting object in a local directory. Only needed if `write_file = TRUE`.
#' Default is NULL.
#' @param seed (numeric) integer value to specify an initial seed to split the
#' data and extract background. Default is 1.
#'
#' @details
#' The available data partitioning methods are:
#'
#' - **"kfolds"**: Splits the dataset into *K* subsets (folds) of approximately
#'                 equal size. In each partition, one fold is used as the test
#'                 set, while the remaining folds are combined to form the
#'                 training set.
#' - **"bootstrap"**: Creates the training dataset by sampling observations
#'                    from the original dataset *with replacement* (i.e., the
#'                    same observation can be selected multiple times). The test
#'                    set consists of the observations that were not selected in
#'                    that specific replicate.
#' - **"subsample"**: Similar to bootstrap, but the training set is created by
#'                    sampling *without replacement* (i.e., each observation is
#'                    selected at most once). The test set includes the
#'                    observations not selected for training.
#'
#' @return
#' An object of class `prepared_data` containing all elements to run a model
#' calibration routine. The elements include: species, calibration data,
#' a grid of model parameters, indices of test data for cross validation,
#' xy coordinates, names of continuous and categorical variables, weights,
#' results from PCA, and modeling algorithm.
#'
#' @importFrom enmpa aux_var_comb
#' @importFrom terra crop prcomp extract as.data.frame nlyr ncell ext res
#' @importFrom utils combn
#'
#' @export
#'
#' @examples
#' # Import occurrences
#' data(occ_data, package = "kuenm2")
#'
#' # Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' # Import a bias file
#' bias <- terra::rast(system.file("extdata", "bias_file.tif",
#'                                 package = "kuenm2"))
#'
#' # Prepare data for maxnet model
#' sp_swd <- prepare_data(algorithm = "maxnet", occ = occ_data,
#'                        x = "x", y = "y",
#'                        raster_variables = var,
#'                        species = occ_data[1, 1],
#'                        categorical_variables = "SoilType",
#'                        n_background = 500, bias_file = bias,
#'                        bias_effect = "direct",
#'                        features = c("l", "q", "p", "lq", "lqp"),
#'                        r_multiplier = c(0.1, 1, 2, 3, 5),
#'                        partition_method = "kfolds")
#' print(sp_swd)
#'
#' # Prepare data for glm model
#' sp_swd_glm <- prepare_data(algorithm = "glm", occ = occ_data,
#'                            x = "x", y = "y",
#'                            raster_variables = var,
#'                            species = occ_data[1, 1],
#'                            categorical_variables = "SoilType",
#'                            n_background = 500, bias_file = bias,
#'                            bias_effect = "direct",
#'                            features = c("l", "q", "p", "lq", "lqp"),
#'                            partition_method = "kfolds")
#' print(sp_swd_glm)

prepare_data <- function(algorithm,
                         occ,
                         x,
                         y,
                         raster_variables,
                         species = NULL,
                         n_background = 1000,
                         features = c("lq", "lqp"),
                         r_multiplier = c(0.1, 0.5, 1, 2, 3),
                         partition_method,
                         n_partitions = 4,
                         train_proportion = 0.7,
                         categorical_variables = NULL,
                         do_pca = FALSE,
                         center = TRUE,
                         scale = TRUE,
                         exclude_from_pca = NULL,
                         variance_explained = 95,
                         min_explained = 5,
                         min_number = 2,
                         min_continuous = NULL,
                         bias_file = NULL,
                         bias_effect = NULL,
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
  if (missing(occ)) {
    stop("Argument 'occ' must be defined.")
  }
  if (missing(x)) {
    stop("Argument 'x' must be defined.")
  }
  if (missing(y)) {
    stop("Argument 'y' must be defined.")
  }
  if (missing(raster_variables)) {
    stop("Argument 'raster_variables' must be defined.")
  }
  if (length(algorithm) != 1) {
    stop("'algorithm' must be a single value: either 'glm' or 'maxnet'.")
  } else {
    if (!algorithm %in% c("glm", "maxnet")) {
      stop("'algorithm' must be 'glm' or 'maxnet'.")
    }
  }

  if (!inherits(occ, "data.frame")) {
    stop("Argument 'data' must be a 'data.frame'.")
  }

  if (!is.null(species) & !inherits(species, "character")) {
    stop("Argument 'species' must be a 'character'.")
  }

  if (!inherits(x, "character") | !inherits(y, "character")) {
    stop("Arguments 'x' and 'y' must be of class 'character'.")
  }

  if (!x %in% colnames(occ) | !y %in% colnames(occ)) {
    stop("'x' and/or 'y' are not column names in 'occ'.")
  }

  if (!inherits(raster_variables, "SpatRaster")) {
    stop("Argument 'raster_variables' must be a 'SpatRaster'.")
  }

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

  if (!is.null(weights) && nrow(occ) != length(weights)) {
    stop("Length of 'weights' must match the number of occurrences in 'occ'.")
  }

  if (min_number > terra::nlyr(raster_variables)) {
    stop("'min_number' can't be greater than the number of layers in 'raster_variables'.")
  }

  if (!is.null(min_continuous)) {
    if (!is.numeric(min_continuous)) {
      stop("Argument 'min_continuous' must be 'numeric'.")
    } else {
      if (min_continuous > terra::nlyr(raster_variables)) {
        stop("'min_continuous' can't be greater than the number of layers in 'raster_variables'.")
      }
    }
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

  if (!is.null(categorical_variables)) {
    if (any(!categorical_variables %in% names(raster_variables))) {
      stop("Defined 'categorical_variables' must be layers in 'raster_variables'.")
    }
    continuous_variable_names <- setdiff(names(raster_variables),
                                         categorical_variables)
  } else {
    continuous_variable_names <- names(raster_variables)
  }

  if (!is.null(bias_file)) {
    if (!inherits(bias_file, "SpatRaster")) {
      stop("Argument 'bias_file' must be a 'SpatRaster'.")
    }
    if (is.null(bias_effect)) {
      stop("Argument 'bias_effect' must be defined if 'bias_file' is provided.")
    }
    if (!inherits(bias_effect, "character")) {
      stop("Argument 'bias_effect' must be a 'character'.")
    }
    if (!(bias_effect %in% c("direct", "inverse")) | length(bias_effect) > 1) {
      stop("Argument 'bias_effect' must be 'direct' or 'inverse'.")
    }
  }

  if(!(partition_method %in% c("kfolds", "subsample", "bootstrap"))){
    stop("Invalid 'partition_method'. Available options include 'kfolds', 'subsample', and 'bootstrap'")
  }

  if(!(n_partitions %% 1 == 0) || n_partitions <= 0){
    stop("'n_partitions' must be a positive numeric integer (e.g., 1, 2, 3...)")
  }

  if(partition_method %in% c("bootstrap", "subsample")){
    if(train_proportion > 1 || train_proportion <= 0 ||
       is.na(train_proportion) || is.null(train_proportion)){
      stop("'train_proportion' must be a positive numeric between 0 and 1")
    }
  }

  sp_name <- species

  if (do_pca) {
    if (!is.null(categorical_variables)) {
      exclude_from_pca <- c(categorical_variables, exclude_from_pca)
    } else {
      exclude_from_pca <- exclude_from_pca
    }

    pca <- perform_pca(raster_variables = raster_variables,
                       exclude_from_pca = exclude_from_pca,
                       project = FALSE, projection_data = NULL,
                       out_dir = pca_directory,
                       overwrite = FALSE, progress_bar = FALSE,
                       center = center, scale = scale,
                       variance_explained = variance_explained,
                       min_explained = min_explained)

    pca$projection_directory <- NULL #Remove projection directory
    raster_variables <- pca$env
  } else {
    pca <- NULL
  }

  occ_var <- extract_occurrence_variables(occ, x, y, raster_variables)
  bg_var <- generate_background_variables(raster_variables, n_background,
                                          bias_file, bias_effect, seed = seed)
  bg_var <- bg_var[,colnames(occ_var)]

  # combine occurrence and background data
  occ_bg <- rbind(occ_var, bg_var)
  # occ_bg <- occ_bg[, c(1, 2, which(names(occ_bg) == "pr_bg"),
  #                      (3:(ncol(occ_bg)-1))[-which(names(occ_bg) == "pr_bg")])]
  occ_bg <- handle_missing_data(occ_bg, weights)

  occ_bg_xy <- if (include_xy) {occ_bg[, c("x", "y")]} else {NULL}
  occ_bg <- subset(occ_bg, select = -c(x, y))

  if (!is.null(categorical_variables)) {
    occ_bg[categorical_variables] <- lapply(occ_bg[categorical_variables],
                                            factor)
  }

  #Fix row.names
  row.names(occ_bg) <- NULL

  # k_f <- enmpa::kfold_partition(data = occ_bg, dependent = "pr_bg", k = kfolds,
  #                               seed = seed)

  #Partitione data
  pd <- part_data(data = occ_bg, pr_bg = "pr_bg",
                  train_proportion = train_proportion,
                  n_partitions = n_partitions,
                  partition_method = partition_method,
                  seed = seed)

  #Formula grid
  formula_grid <- calibration_grid(occ_bg = occ_bg, min_number = min_number,
                                   min_continuous = min_continuous,
                                   categorical_var = categorical_variables,
                                   features = features, algorithm = algorithm,
                                   r_multiplier = r_multiplier)

  #Prepare final data
  if(partition_method == "kfolds"){
    train_proportion <- NULL
  }


  data <- new_prepared_data(species = species, calibration_data = occ_bg,
                            formula_grid = formula_grid,
                            part_data = pd, partition_method = partition_method,
                            n_partitions = n_partitions,
                            train_proportion = train_proportion,
                            data_xy = occ_bg_xy,
                            continuous_variables = continuous_variable_names,
                            categorical_variables = categorical_variables,
                            weights = weights, pca = pca$pca,
                            algorithm = algorithm)

  if (write_file) {
    saveRDS(data, paste0(file_name, ".RDS"))
  }

  return(data)
}
