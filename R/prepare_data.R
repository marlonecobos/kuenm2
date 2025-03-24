#' Prepare data for model calibration
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
#'              n_background = 10000, bias_file = NULL, bias_effect = "direct",
#'              kfolds = 4, weights = NULL,
#'              min_number = 2, min_continuous = NULL,
#'              features = c("q", "lq", "lp", "qp", "lqp"),
#'              reg_mult = c(0.1, 0.5, 1, 2, 3), include_xy = TRUE,
#'              write_file = FALSE, file_name = NULL, seed = 1)
#'
#' @param algorithm (character) type algorithm, either "glm" or "maxnet".
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
#' @param bias_file (SpatRaster) a raster containing bias values (probability
#' weights) that influence the selection of background points. It must have the
#' same extent, resolution, and number of cells as the raster variables, unless
#' a mask is provided. Default is NULL.
#' @param bias_effect (character) a string specifying how the values in the
#' `bias_file` should be interpreted. If "direct", higher values in the bias
#' file increase the likelihood of selecting background points. If "inverse",
#' higher values decrease the likelihood. Default is "direct".
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
#' @param reg_mult (numeric) a vector of regularization parameters for maxnet.
#' Default is c(0.1, 1, 2, 3, 5).
#' @param include_xy (logical) whether to include the coordinates (longitude and
#' latitude) in the results from preparing data. Default is TRUE.
#' @param write_file (logical) whether to write the resulting prepared_data list
#' in a local directory. Default is FALSE.
#' @param file_name (character) the path or name of the folder where the
#' resulting list will be saved. This is only applicable if `write_file =
#' TRUE`. Default is NULL.
#' @param seed (numeric) integer value to specify an initial seed to split the
#' data and extract background. Default is 1.
#'
#' @return
#' An object of class `prepared_data` containing all elements to run a model
#' calibration routine. The elements include: species, calibration data,
#' a grid of model parameters, indices of k-folds for cross validation,
#' xy coordinates, names of continuous and categorical variables, weights,
#' results from PCA, and modeling algorithm.
#'
#' @importFrom enmpa kfold_partition aux_var_comb
#' @importFrom terra crop prcomp extract as.data.frame nlyr ncell ext res
#' @importFrom utils combn
#'
#' @export
#'
#' @examples
#' # Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' # Import occurrences
#' data(occ_data, package = "kuenm2")
#'
#' # Import a bias file
#' bias <- terra::rast(system.file("extdata", "bias_file.tif",
#'                                 package = "kuenm2"))
#'
#' # Prepare data for maxnet model
#' sp_swd <- prepare_data(algorithm = "maxnet", occ = occ_data,
#'                        species = occ_data[1, 1], x = "x", y = "y",
#'                        raster_variables = var,
#'                        categorical_variables = "SoilType",
#'                        n_background = 500, bias_file = bias,
#'                        features = c("l", "q", "p", "lq", "lqp"),
#'                        reg_mult = c(0.1, 1, 2, 3, 5))
#' print(sp_swd)
#'
#' # Prepare data for glm model
#' sp_swd_glm <- prepare_data(algorithm = "glm", occ = occ_data,
#'                            species = occ_data[1, 1], x = "x", y = "y",
#'                            raster_variables = var,
#'                            categorical_variables = "SoilType",
#'                            n_background = 500, bias_file = bias,
#'                            features = c("l", "q", "p", "lq", "lqp"))
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
                         bias_file = NULL,
                         bias_effect = "direct",
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
  if(length(algorithm) != 1) {
    stop("'algorithm' must be a single value: either 'glm' or 'maxnet'")
  } else {
    if (!algorithm %in% c("glm", "maxnet")) {
      stop("'algorithm' must be 'glm' or 'maxnet'")
    }
  }

  if (!inherits(occ, "data.frame")) {
    stop("Argument 'data' must be a 'data.frame'")
  }

  if (!is.null(species) & !inherits(species, "character")) {
    stop("Argument 'species' must be a character")
  }

  if (!inherits(x, "character") | !inherits(y, "character")) {
    stop("Arguments 'x' and 'y' must be of class character")
  }

  if(!x %in% colnames(occ) | !y %in% colnames(occ)) {
    stop("'x' and/or 'y' are not column names in 'occ'")
  }

  if (!inherits(raster_variables, "SpatRaster")) {
    stop("Argument 'raster_variables' must be a 'SpatRaster'")
  }

  if(!is.null(mask) &
     !inherits(mask, c("SpatRaster", "SpatVector", "SpatExtent"))) {
    stop("Argument 'mask' must be a 'SpatVector', 'SpatExtent' or 'SpatRaster'")
  }

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

  if (min_number > terra::nlyr(raster_variables)) {
    stop("'min_number' can't be greater than the number of layers in 'raster_variables'")
  }

  if (!is.null(min_continuous)) {
    if (!is.numeric(min_continuous)) {
      stop("Argument 'min_continuous' must be 'numeric'")
    } else {
      if (min_continuous > terra::nlyr(raster_variables)) {
        stop("'min_continuous' can't be greater than the number of layers in 'raster_variables'")
      }
    }
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

  if (!is.null(categorical_variables)) {
    if (any(!categorical_variables %in% names(raster_variables))) {
      stop("Defined 'categorical_variables' must be layers in 'raster_variables'")
    }
    continuous_variable_names <- setdiff(names(raster_variables),
                                         categorical_variables)
  } else {
    continuous_variable_names <- names(raster_variables)
  }

  if(!is.null(bias_file)){
    if(!inherits(bias_file, "SpatRaster")){
      stop("Argument bias_file must be a SpatRaster, not ", class(bias_file))
    }
    if(!inherits(bias_effect, "character")){
      stop("Argument bias_effect must be a character, not ", class(bias_file))
    }
    if(!(bias_effect %in% c("direct", "inverse")) | length(bias_effect) > 1) {
      stop("Argument bias_effect must be 'direct' or 'inverse'")
    }
    }

  sp_name <- species

  if (!is.null(mask)) {
    raster_variables <- terra::crop(raster_variables, mask, mask = TRUE)
    if(!is.null(bias_file)){
      bias_file <- terra::crop(bias_file, mask, mask = TRUE)
    }
  }

  if (do_pca) {
    if (!is.null(categorical_variables)){
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

  # combine occurrence and background data
  occ_bg <- rbind(occ_var, bg_var)
  occ_bg <- occ_bg[, c(1, 2, which(names(occ_bg) == "pr_bg"),
                       (3:(ncol(occ_bg)-1))[-which(names(occ_bg) == "pr_bg")])]
  occ_bg <- handle_missing_data(occ_bg, weights)

  occ_bg_xy <- if (include_xy) occ_bg[, c(x, y)] else NULL
  occ_bg <- occ_bg[, -(1:2)]

  if (!is.null(categorical_variables)) {
    occ_bg[categorical_variables] <- lapply(occ_bg[categorical_variables],
                                            factor)
  }


  k_f <- enmpa::kfold_partition(data = occ_bg, dependent = "pr_bg", k = kfolds,
                                seed = seed)

  #Formula grid
  formula_grid <- calibration_grid(occ_bg, min_number, min_continuous,
                                   categorical_var = categorical_variables,
                                   features, algorithm, reg_mult)

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
