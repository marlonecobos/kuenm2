#' prepare_data Class Constructor
#'
#' Creates a prepare_data object to hold the prepared data for model calibration.
#'
#' @param species Character name of the species.
#' @param calibration_data Data frame containing the calibration data.
#' @param kfolds List of k-fold partitions.
#' @param data_xy Data frame containing the XY coordinates.
#' @param continuous_variables Character vector of continuous variables.
#' @param categorical_variables Character vector of categorical variables.
#' @param weights Vector of weights for occurrences.
#' @param pca Principal component analysis results (if performed).
#' @return A prepare_data object.
#' @export

new_prepare_data <- function(species, calibration_data, kfolds, data_xy,
                                continuous_variables, categorical_variables,
                                weights, pca) {
  data <- list(
    species = species,
    calibration_data = calibration_data,
    kfolds = kfolds,
    data_xy = data_xy,
    continuous_variables = continuous_variables,
    categorical_variables = categorical_variables,
    weights = weights,
    pca = pca
  )
  class(data) <- "prepare_data"
  return(data)
}

#' Print Method for prepare_data Class
#'
#' Provides a summary of the prepare_data object.
#'
#' @param x A prepare_data object.
#' @param ... Additional arguments (currently not used).
#' @export

print.prepare_data <- function(x, ...) {
  cat("prepare_data object summary\n")
  cat("==========================\n")
  cat("Species:", x$species, "\n")
  cat("Number of occurrences:", nrow(x$calibration_data), "\n")
  cat("  - Number of presence points:", table(x$calibration_data$pr_bg)[2], "\n")
  cat("  - Number of background points:", table(x$calibration_data$pr_bg)[1], "\n")

  cat("k-Fold Cross-Validation:\n")
  cat("  - Number of folds:", length(x$kfolds), "\n")
  cat("Continuous Variables:\n")
  cat("  -", paste(x$continuous_variables, collapse = ", "), "\n")

  if (!is.null(x$categorical_variables) && length(x$categorical_variables) > 0) {
    cat("Categorical Variables:\n")
    cat("  -", paste(x$categorical_variables, collapse = ", "), "\n")
  } else {
    cat("Categorical Variables: None\n")
  }

  if (!is.null(x$pca)) {
    cat("PCA Information:\n")
    cat("  - Variables included:", paste(x$pca$pca$vars_in, collapse = ", "), "\n")
    cat("  - Number of PCA components:", length(x$pca$deviance_explained_cumsum), "\n")
  } else {
    cat("PCA Information: PCA not performed\n")
  }

  if (!is.null(x$weights)) {
    cat("Weights Information:\n")
    cat("  - Weights provided: Yes\n")
  } else {
    cat("Weights Information: No weights provided\n")
  }
}

#' Prepare Data for Model Calibration
#'
#' This function prepares data for model calibration, including optional PCA,
#' background point generation, and k-fold partitioning.
#'
#' @param occ A data frame containing occurrences.
#' @param species Character string specifying the species name.
#' @param x Character vector specifying the x-coordinate variable.
#' @param y Character vector specifying the y-coordinate variable.
#' @param spat_variables SpatRaster variables.
#' @param mask SpatRaster, SpatVector, SpatExtent, or any other object that has
#' a SpatExtent object to mask the spatial variables.
#' @param categorical_variables Character vector specifying categorical variables.
#' @param do_pca Logical, whether to perform PCA on the spatial variables.
#' @param deviance_explained Numeric, the percentage of deviance to be explained by PCA components.
#' @param min_explained Numeric, the minimum percentage of deviance explained by each PCA component.
#' @param center Logical, whether to center the variables before performing PCA.
#' @param scale Logical, whether to scale the variables before performing PCA.
#' @param write_pca Logical, whether to write the PCA output to a file.
#' @param nbg Numeric, number of background points to generate.
#' @param kfolds Numeric, number of folds for k-fold partitioning.
#' @param weights Numeric vector specifying weights for occurrences.
#' @param include_xy Logical, whether to include XY coordinates in the output.
#' @param write_files Logical, whether to write the results to a file.
#' @param file_name Character string specifying the output file name.
#' @param seed Numeric, the seed for random number generation.
#' @return A prepare_data object containing the prepared data.
#' @export

prepare_data <- function(occ,
                         species = NULL,
                         x,
                         y,
                         spat_variables,
                         mask = NULL,
                         categorical_variables = NULL,
                         do_pca = FALSE,
                         deviance_explained = 95,
                         min_explained = 5,
                         center = TRUE,
                         scale = TRUE,
                         write_pca = FALSE,
                         nbg = 10000,
                         kfolds = 4,
                         weights = NULL,
                         include_xy = TRUE,
                         write_files = FALSE,
                         file_name = NULL,
                         seed = 1) {

  if (write_files & is.null(file_name)) {
    file_name <- paste0("prepare_data_", species)
  }

  if (!is.null(weights) && nrow(occ) != length(weights)) {
    stop("Length of weights does not match the number of occurrences in occ")
  }

  if (!is.null(categorical_variables)){
    if (!categorical_variables %in% names(spat_variables)){
      stop("Categorical variables are not defined correctly.")
    }
    continuous_variable_names <- setdiff(names(spat_variables), categorical_variables)
  } else {
    continuous_variable_names <- names(spat_variables)
  }

  sp_name <- if (is.null(species)) as.character(occ[1, "species"]) else species

  if (!is.null(mask)) {
    spat_variables <- terra::crop(spat_variables, mask, mask = TRUE)
  }

  if (do_pca) {
    if (!is.null(categorical_variables)){
      exclude_from_pca = categorical_variables
    } else {
      exclude_from_pca = NULL
    }
    pca <- perform_pca(spat_variables, exclude_from_pca, center, scale,
                       deviance_explained, min_explained)
    env <- pca$env
  } else {
    env <- spat_variables
    pca <- NULL
  }


  occ_var <- extract_occurrence_variables(occ, x, y, env)
  bg_var <- generate_background_variables(env, nbg)

  # combine occurrence and background data
  occ_bg <- rbind(occ_var, bg_var)
  occ_bg <- occ_bg[, c(1, 2, which(names(occ_bg) == "pr_bg"), (3:(ncol(occ_bg)-1))[-which(names(occ_bg) == "pr_bg")])]
  occ_bg <- handle_missing_data(occ_bg, weights)

  occ_bg_xy <- if (include_xy) occ_bg[, 1:2] else NULL
  occ_bg <- occ_bg[, -(1:2)]

  if (!is.null(categorical_variables)){
    occ_bg[categorical_variables] <- lapply(occ_bg[categorical_variables], factor)
  }


  k_f <- enmpa::kfold_partition(data = occ_bg, dependent = "pr_bg", k = kfolds,
                                seed = seed)

  data <- new_prepare_data(sp_name, occ_bg, k_f, occ_bg_xy,
                              continuous_variable_names, categorical_variables,
                              weights, pca)

  if (write_files) {
    saveRDS(data, paste0(file_name, ".RDS"))
  }

  return(data)
}

# Helper function to perform PCA
perform_pca <- function(spat_variables, exclude_from_pca = NULL, center, scale,
                        deviance_explained, min_explained) {

  if (!is.null(exclude_from_pca)) {
    var_to_pca <-  setdiff(names(spat_variables), exclude_from_pca)
  } else {
    var_to_pca <- names(spat_variables)
  }


  pca <- terra::prcomp(spat_variables[[var_to_pca]], center = center, scale = scale)
  pca$x <- NULL
  pca$vars_in <- var_to_pca
  pca$vars_out <- exclude_from_pca

  d_exp <- cumsum(pca$sdev/sum(pca$sdev)) * 100
  d_exp <- d_exp[(pca$sdev/sum(pca$sdev) * 100) > min_explained]

  ind_exp <- if (max(d_exp) > deviance_explained) min(which(d_exp >= deviance_explained)) else length(d_exp)

  env <- predict(spat_variables[[var_to_pca]], pca, index = 1:ind_exp)

  if (!is.null(exclude_from_pca)) {
    env <- c(env, spat_variables[[exclude_from_pca]])
  }

  names(d_exp) <- paste0("PC", 1:ind_exp)

  return(list(env = env, pca = pca, deviance_explained_cumsum = d_exp))
  }

# Helper function to extract occurrence variables
extract_occurrence_variables <- function(occ, x, y, env) {
  xy <- as.matrix(occ[, c(x, y)])
  colnames(xy) <- c("x", "y")
  occ_var <- cbind(xy, terra::extract(x = env, y = xy))
  occ_var$pr_bg <- 1
  return(occ_var)
}

# Helper function to generate background variables
generate_background_variables <- function(env, nbg) {
  cell_samp <- terra::as.data.frame(env[[1]], na.rm = TRUE, cells = TRUE)[, "cell"]

  if (length(cell_samp) < nbg) {
    nbg <- length(cell_samp)
  }

  cell_samp <- sample(cell_samp, size = nbg, replace = FALSE)
  bg_var <- terra::extract(x = env, y = cell_samp, xy = TRUE)
  bg_var$pr_bg <- 0

  return(bg_var)
}

# Helper function to handle missing data
handle_missing_data <- function(occ_bg, weights) {
  na_rows <- which(!complete.cases(occ_bg))

  if (length(na_rows) > 0) {
    occ_bg <- occ_bg[-na_rows,]

    if (!is.null(weights)) {
      weights <- weights[-na_rows]
    }
    warning(length(na_rows), " rows were excluded from database because NAs were found")
  }
  return(occ_bg)
}

