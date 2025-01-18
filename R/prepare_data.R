#' Prepare Data for Model Calibration
#'
#' @description
#' This function prepares data for model calibration, including optional PCA,
#' background point generation, k-fold partitioning, and the creation of a grid
#' parameter combinations, including distinct regularization multiplier values,
#' various feature classes, and different sets of environmental variables.
#'
#' @usage
#' prepare_data(algorithm, occ, species = NULL, x, y,
#'              raster_variables, mask = NULL, categorical_variables = NULL,
#'              do_pca = FALSE, exclude_from_pca = NULL, variance_explained = 95,
#'              min_explained = 5, center = TRUE, scale = FALSE,
#'              write_pca = FALSE, pca_directory = NULL,
#'              n_background = 10000, kfolds = 4, weights = NULL,
#'              min_number = 2, min_continuous = NULL,
#'              features = c("q", "lq", "lp", "qp", "lqp"),
#'              reg_mult = c(0.1, 0.5, 1, 2, 3), include_xy = TRUE,
#'              write_file = FALSE, file_name = NULL, seed = 1)
#'
#' @param algorithm (character) type algorithm, either "glm" or "glmnet".
#' @param occ (data frame) a data.frame containing the coordinates (longitude
#' and latitude) of the occurrence records.
#' @param species (character) string specifying the species name (optional).
#' Default is NULL.
#' @param x (character) a string specifying the name of the column in `occ` that
#' contains the longitude values.
#' @param y (character) a string specifying the name of the column in `occ` that
#' contains the latitude values.
#' @param raster_variables (SpatRaster) predictor variables used to calibrate the
#' models.
#' @param mask (SpatRaster, SpatVector, or SpatExtent) spatial object used to
#' mask the variables to the area where the model will be calibrated.
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
#' for the model. Default is 10000.
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
#'
#' @export
#'
#' @examples
#' #Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#' #Import occurrences
#' data(occ_data, package = "kuenm2")
#' #Prepare data for glmnet model
#' sp_swd <- prepare_data(algorithm = "glmnet", occ = occ_data,
#'                        species = occ_data[1, 1], x = "x", y = "y",
#'                        raster_variables = var, mask = NULL,
#'                        categorical_variables = "SoilType",
#'                        do_pca = FALSE, variance_explained = 95,
#'                        min_explained = 5, center = TRUE, scale = TRUE,
#'                        write_pca = FALSE, pca_directory = NULL,
#'                        exclude_from_pca = NULL, n_background = 500,
#'                        kfolds = 4, weights = NULL, min_number = 2,
#'                        min_continuous = NULL,
#'                        features = c("l", "q", "p", "lq", "lqp"),
#'                        reg_mult = c(0.1, 1, 2, 3, 5),
#'                        include_xy = TRUE,
#'                        write_file = FALSE, file_name = NULL,
#'                        seed = 1)
#' print(sp_swd)
#' #Prepare data for glm model
#' sp_swd_glm <- prepare_data(algorithm = "glm", occ = occ_data,
#'                            species = occ_data[1, 1], x = "x", y = "y",
#'                            raster_variables = var, mask = NULL,
#'                            categorical_variables = "SoilType",
#'                            do_pca = FALSE, variance_explained = 95,
#'                            min_explained = 5, center = TRUE, scale = TRUE,
#'                            write_pca = FALSE, pca_directory = NULL,
#'                            exclude_from_pca = NULL, n_background = 500,
#'                            kfolds = 4, weights = NULL, min_number = 2,
#'                            min_continuous = NULL, features = c("l", "q", "lq", "lpq"),
#'                            reg_mult = c(0.1, 1, 2), include_xy = TRUE,
#'                            write_file = F, file_name = FALSE, seed = 1)
#' print(sp_swd_glm)



prepare_data <- function(algorithm,
                         occ,
                         species = NULL,
                         x,
                         y,
                         raster_variables,
                         mask = NULL,
                         categorical_variables = NULL,
                         do_pca = FALSE,
                         variance_explained = 95,
                         min_explained = 5,
                         center = TRUE,
                         scale = TRUE,
                         write_pca = FALSE,
                         pca_directory = NULL,
                         exclude_from_pca = NULL,
                         n_background = 10000,
                         kfolds = 4,
                         weights = NULL,
                         min_number = 2,
                         min_continuous = NULL,
                         features = c("q", "lq", "lp", "qp", "lqp"),
                         reg_mult = c(0.1, 0.5, 1, 2, 3),
                         include_xy = TRUE,
                         write_file = FALSE,
                         file_name = NULL,
                         seed = 1) {
  #Check data
  if(length(algorithm) != 1){
    stop("'algorithm' must be a single value: either 'glm' or 'glmnet'")
  }

  if (!algorithm %in% c("glm", "glmnet")) {
    stop("'algorithm' must be 'glm' or 'glmnet'")
  }

  if (!inherits(occ, "data.frame")) {
    stop(paste0("Argument data must be a data.frame, not ",
                class(occ)))
  }

  if (!is.null(species) & !inherits(species, "character")) {
    stop(paste0("Argument species must be a character, not ",
                class(species)))
  }

  if (!inherits(x, "character")) {
    stop(paste0("Argument x must be a character, not ",
                class(x)))
  }

  if (!inherits(y, "character")) {
    stop(paste0("Argument y must be a character, not ",
                class(y)))
  }

  if(!(x %in% colnames(occ))){
    stop("The value provided for 'x' must match one of the column names in 'occ'.")
  }

  if(!(y %in% colnames(occ))){
    stop("The value provided for 'y' must match one of the column names in 'occ'.")
  }

  if (!inherits(raster_variables, "SpatRaster")) {
    stop(paste0("Argument raster_variables must be a SpatRaster, not ",
                class(raster_variables)))
  }

  if(!is.null(mask) & !inherits(mask, c("SpatRaster", "SpatVector",
                                        "SpatExtent"))){
    stop(paste0("Argument mask must be a SpatVector, SpatExtent or SpatRaster, not ",
                class(mask)))
  }

  if (!is.null(categorical_variables) & !inherits(categorical_variables, "character")) {
    stop(paste0("Argument categorical_variables must be a character, not ",
                class(categorical_variables)))
  }

  if (!is.logical(do_pca)) {
    stop(paste0("Argument do_pca must be logical, not ",
                class(do_pca)))
  }

  if (do_pca & !is.numeric(variance_explained)) {
    stop(paste0("Argument variance_explained must be numeric, not ",
                class(variance_explained)))
  }

  if (do_pca & !is.numeric(min_explained)) {
    stop(paste0("Argument min_explained must be numeric, not ",
                class(min_explained)))
  }

  if (do_pca & !is.logical(center)) {
    stop(paste0("Argument center must be logical, not ",
                class(center)))
  }

  if (do_pca & !is.logical(scale)) {
    stop(paste0("Argument scale must be logical, not ",
                class(scale)))
  }

  if (do_pca & !is.logical(write_pca)) {
    stop(paste0("Argument write_pca must be logical, not ",
                class(write_pca)))
  }

  if (!is.numeric(n_background)) {
    stop(paste0("Argument n_background must be numeric, not ",
                class(n_background)))
  }

  if (!is.numeric(kfolds)) {
    stop(paste0("Argument kfolds must be numeric, not ",
                class(kfolds)))
  }

  if (!is.null(weights) && nrow(occ) != length(weights)) {
    stop("Length of weights does not match the number of occurrences in occ")
  }

  if (!is.numeric(min_number)) {
    stop(paste0("Argument min_number must be numeric, not ",
                class(min_number)))
  }

  if (min_number > terra::nlyr(raster_variables)) {
    stop("'min_number' can't be greater than the number of variables in 'raster_variables'")
  }

  if (!is.null(min_continuous) & !is.numeric(min_continuous)) {
    stop(paste0("Argument min_continuous must be numeric, not ",
                class(min_continuous)))
  }

  if (!is.null(min_continuous)) {
    if(min_continuous > terra::nlyr(raster_variables)){
    stop("'min_continuous' can't be greater than the number of variables in 'raster_variables'")
  }}

  if (!is.character(features)) {
    stop(paste0("Argument features must be character, not ",
                class(features)))
  }

  check_features <- unique(unlist(strsplit(features, "")))
  check_features <- setdiff(check_features, c("l", "q", "p", "h", "t"))
  if(length(check_features) > 0){
    stop(paste(check_features, collapse = ", "), " is not a valid feature\n",
         "Features can be l, q, p, t, h, or combinations of these")
  }

  if (!is.numeric(reg_mult)) {
    stop(paste0("Argument reg_mult must be numeric, not ",
                class(reg_mult)))
  }

  if (!is.logical(include_xy)) {
    stop(paste0("Argument include_xy must be logical, not ",
                class(include_xy)))
  }

  if (!is.logical(write_file)) {
    stop(paste0("Argument write_file must be logical, not ",
                class(write_file)))
  }

  if (write_file & is.null(file_name)) {
    stop("If 'write_file = TRUE', you must specify a 'file_name'")
  }

  if (write_file & !is.character(file_name)) {
    stop("If 'write_file = TRUE', 'file_name' must be a character, not ",
         class(file_name))
  }

  if (!is.numeric(seed)) {
    stop(paste0("Argument seed must be numeric, not ",
                class(seed)))
  }

  if (!is.null(categorical_variables)){
    if (!categorical_variables %in% names(raster_variables)){
      stop("Categorical variables are not defined correctly. They must be a layer in 'raster_variables'")
    }
    continuous_variable_names <- setdiff(names(raster_variables), categorical_variables)
  } else {
    continuous_variable_names <- names(raster_variables)
  }

  sp_name <- species

  if (!is.null(mask)) {
    raster_variables <- terra::crop(raster_variables, mask, mask = TRUE)
  }

  if (do_pca) {
    if (!is.null(categorical_variables)){
      exclude_from_pca = c(categorical_variables, exclude_from_pca)
    } else {
      exclude_from_pca = exclude_from_pca
    }
    pca <- perform_pca(raster_variables, exclude_from_pca = exclude_from_pca,
                       project = FALSE, projection_data = NULL, out_dir = NULL,
                       overwrite = FALSE, progress_bar = FALSE,
                       center = center, scale = scale,
                       variance_explained = variance_explained,
                       min_explained = min_explained)
    pca$projection_directory <- NULL #Remove projection directory
    env <- pca$env
  } else {
    env <- raster_variables
    pca <- NULL
  }

  occ_var <- extract_occurrence_variables(occ, x, y, env)
  bg_var <- generate_background_variables(env, n_background)

  # combine occurrence and background data
  occ_bg <- rbind(occ_var, bg_var)
  occ_bg <- occ_bg[, c(1, 2, which(names(occ_bg) == "pr_bg"), (3:(ncol(occ_bg)-1))[-which(names(occ_bg) == "pr_bg")])]
  occ_bg <- handle_missing_data(occ_bg, weights)

  occ_bg_xy <- if (include_xy) occ_bg[, c(x, y)] else NULL
  occ_bg <- occ_bg[, -(1:2)]

  if (!is.null(categorical_variables)){
    occ_bg[categorical_variables] <- lapply(occ_bg[categorical_variables], factor)
  }


  k_f <- enmpa::kfold_partition(data = occ_bg, dependent = "pr_bg", k = kfolds,
                                seed = seed)

  #Formula grid
  formula_grid <- calibration_grid(occ_bg, min_number, min_continuous,
                                   categorical_var = categorical_variables,
                                   features, algorithm,
                                   reg_mult)

  data <- new_prepared_data(species, calibration_data = occ_bg, formula_grid,
                            kfolds = k_f, data_xy = occ_bg_xy,
                            continuous_variables = continuous_variable_names,
                            categorical_variables, weights, pca = pca$pca,
                            algorithm)

  if (write_file) {
    saveRDS(data, paste0(file_name, ".RDS"))
  }

  return(data)
}
