#### HELPERS FOR PREPARE DATA ####
#' prepare_data Class Constructor
new_prepare_data <- function(species, calibration_data, formula_grid,
                             kfolds, data_xy, continuous_variables,
                             categorical_variables, weights, pca, model_type) {
  data <- list(
    species = species,
    calibration_data = calibration_data,
    formula_grid = formula_grid,
    kfolds = kfolds,
    data_xy = data_xy,
    continuous_variables = continuous_variables,
    categorical_variables = categorical_variables,
    weights = weights,
    pca = pca,
    model_type = model_type
  )
  class(data) <- "prepare_data"
  return(data)
}

#' Print Method for prepare_data Class
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

  model_type <- x$model_type

  #Print formula grid
  if(model_type == "glmnet"){
    cat("Calibration Grid (GLMNET)\n")
    cat("  - Features used:", paste(unique(x$formula_grid$Features), collapse = ", "), "\n")
    cat("  - Reg. multipliers used:", paste(unique(x$formula_grid$regm), collapse = ", "), "\n")
    cat("  - Number of combinations:", nrow(x$formula_grid), "\n")
    cat("  - Print formulas (n = 5):\n")
    print(head(x$formula_grid, 5))
  }

  if(model_type == "glm"){
    cat("Calibration Grid (GLM)\n")
    cat("  - Features used:", paste(unique(x$formula_grid$Features), collapse = ", "), "\n")
    cat("  - Number of combinations:", nrow(x$formula_grid), "\n")
    cat("  - Print (n = 5):\n")
    print(head(x$formula_grid, 5))
  }
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

#### Helper functions to create formula grid ####
calibration_grid <- function(occ_bg,
                             min_number = 2,
                             min_continuous = NULL,
                             categorical_var = NULL,
                             features = c("l", "q", "lq", "lqp", "p"),
                             model_type = c("glm", "glmnet"),
                             regm = c(0.1, 1, 2, 3, 5)) {

  # Validate the model_type input
  model_type <- match.arg(model_type)


  if (model_type == "glm") {
    # Call the GLM-specific function
    cal_grid_data <- calibration_grid_glm(occ_bg = occ_bg,
                                          min_number = min_number,
                                          min_continuous = min_continuous,
                                          categorical_var = categorical_var,
                                          features = features)
  } else if (model_type == "glmnet") {
    # Call the GLMNET-specific function
    cal_grid_data <- calibration_grid_glmnetmx(occ_bg = occ_bg,
                                               min_number = min_number,
                                               min_continuous = min_continuous,
                                               categorical_var = categorical_var,
                                               features = features,
                                               regm = regm)
  }

  return(cal_grid_data)
}

#' Calibration Grid Generation using GLMNET
calibration_grid_glmnetmx <- function(occ_bg,
                                      categorical_var = NULL,
                                      features = c("l", "q", "lq", "lqp", "p"),
                                      min_number = 2,
                                      min_continuous = NULL,
                                      regm = c(0.1, 1, 2, 3, 5)){

  #Get variable names
  var_names <- colnames(occ_bg[, -1, drop = FALSE])

  # Get variable combinations
  var_comb <- enmpa:::aux_var_comb(var_names = var_names, minvar = min_number)

  # Remove combinations according to minimum number of continuous variables
  if(!is.null(min_continuous) & !is.null(categorical_var)){
    n_cont <- sapply(var_comb, function(x) length(x[x != categorical_var]))
    var_comb <- var_comb[n_cont > min_continuous]
  }

  # Split features
  formula_x <- list()
  for(f_x in features) {
    if(grepl("p", f_x) & !is.null(categorical_var)) {
      var_comb_new <- var_comb[sapply(var_comb, function(x) sum(!x %in% categorical_var)) >= 2]
    } else {
      var_comb_new <- var_comb
    }
    for(vc in var_comb_new) {
      f_l <- prepare_formulas_glmnetmx(independent = vc, type = f_x, categorical_var = categorical_var)
      names(f_l) <- f_x
      formula_x <- c(formula_x, f_l)
    }
  }

  # Create a data frame with formulas and their corresponding features
  formula_d <- data.frame(Formulas = vapply(formula_x, function(f) paste(f, "-1"), character(1)),
                          Features = names(formula_x), stringsAsFactors = FALSE)

  # Expand the grid by combining formulas and regularization parameters
  f_grid <- expand.grid(Formulas = formula_d$Formulas, regm = regm, stringsAsFactors = FALSE)

  f_grid <- merge(f_grid, formula_d, by = "Formulas", sort = FALSE)
  f_grid$ID <- seq_len(nrow(f_grid))
  f_grid <- f_grid[, c("ID", "Formulas", "regm", "Features")]
  f_grid$Formulas <- as.character(f_grid$Formulas)

  return(f_grid)
}

#' Prepare Formulas for GLMNET
prepare_formulas_glmnetmx <- function(independent, type = "lqpht",
                                      categorical_var = NULL,
                                      minvar = 1, maxvar = NULL) {

  if (is.character(type)) {
    if (!all(unlist(strsplit(type, "")) %in% c("l", "p", "q", "h", "t"))) {
      stop("'type' must be: 'l', 'p', 'q', 'h', 't', or a combination of those.")
    }
  }  else {
    stop("'type' must be of class character.")
  }

  # Check if categorical variable is in independent variables. If not, set null
  if(!is.null(categorical_var)) {
    if(!(categorical_var %in% independent)) {
      categorical_var <- NULL
    }
  }

  if(!is.null(categorical_var)) {
    independent <- setdiff(independent, categorical_var)
  }

  predictors <- independent
  npred <- length(predictors)
  aux <- " "
  # Linear
  if (grepl("l", type)) {
    aux <- paste(aux, paste(predictors, collapse = " + "),
                 sep = " + ")
  }
  # Quadratic
  if (grepl("q", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("I(", predictors[v], "^2)"),
                   sep = " + ")
    }
  }
  # Product
  if (grepl("p", type)) {
    if (npred > 1) {
      inter_tab <- utils::combn(predictors, 2)
      aux_inter <- paste0(" + ", paste(apply(inter_tab,
                                             2, paste, collapse = ":"), collapse = " + "))
      if (length(aux_inter) > 0) {
        aux <- paste0(aux, aux_inter)
      }
    }
    else {
      if (grepl("l", type) | grepl("q", type)) {
        message("'p' is only possible with 2 or more independent variables.",
                "\nReturning other combinations.")
      }
      else {
        stop("'p' is only possible with 2 or more independent variables.",
             "\nTry other combinations of type.")
      }
    }
  }
  # Hinge
  if (grepl("h", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("hinge(", predictors[v], ")"),
                   sep = " + ")
    }
  }
  # Threshold
  if (grepl("t", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("thresholds(", predictors[v], ")"),
                   sep = " + ")
    }
  }
  out <- paste0("~",
                gsub(pattern = "  \\+ ", x = aux, replacement = ""))

  if(!is.null(categorical_var)) {
    out <- paste0(out, " + categorical(", categorical_var, ")")
  }
  return(out)
}

#' Calibration Grid Generation using GLM
calibration_grid_glm <- function(occ_bg,
                                 min_number = 2,
                                 min_continuous = NULL,
                                 categorical_var = NULL,
                                 features = c("l", "q", "lq", "lqp", "p")){

  #Get variable names
  var_names <- colnames(occ_bg[, -1, drop = FALSE])

  # Generate all combinations of variables
  var_combinations <- enmpa:::aux_var_comb(var_names = var_names, minvar = min_number)

  # Filter combinations based on the minimum number of continuous variables
  if (!is.null(min_continuous)) {
    continuous_counts <- sapply(var_combinations, function(vars) length(setdiff(vars, categorical_var)))
    var_combinations <- var_combinations[continuous_counts > min_continuous]
  }

  # Prepare formulas for each feature type
  formula_list <- list()

  for (feature_type in features) {
    if (grepl("p", feature_type) && !is.null(categorical_var)) {
      filtered_combinations <- var_combinations[sapply(var_combinations, function(vars) length(setdiff(vars, categorical_var)) >= 2)]
    } else {
      filtered_combinations <- var_combinations
    }

    for (vars in filtered_combinations) {
      formula <- prepare_formulas_glm(independent = vars, type = feature_type,
                                      categorical_var = categorical_var)
      formula_list <- c(formula_list, setNames(list(formula), feature_type))
    }
  }

  # Create a data frame with the formulas and their corresponding features
  formula_grid <- data.frame(
    ID = seq_along(formula_list),
    Formulas = unlist(formula_list),
    Features = names(formula_list),
    stringsAsFactors = FALSE
  )
  return(formula_grid)
}

#' Prepare Formulas for GLM
prepare_formulas_glm <- function(independent, type = "l", categorical_var = NULL) {

  # Validate the 'type' input
  valid_types <- c("l", "p", "q")
  if (!all(unlist(strsplit(type, "")) %in% valid_types)) {
    stop("'type' must be: 'l', 'p', 'q', or a combination of those.")
  }

  # Remove the categorical variable from independent variables if not present
  if (!is.null(categorical_var) && !(categorical_var %in% independent)) {
    categorical_var <- NULL
  }

  # Exclude the categorical variable from independent variables for GLM
  if (!is.null(categorical_var)) {
    independent <- setdiff(independent, categorical_var)
  }

  # Initialize the formula with an empty space
  formula_parts <- " "

  # Linear terms
  if (grepl("l", type)) {
    formula_parts <- paste(formula_parts, paste(independent, collapse = " + "), sep = " + ")
  }

  # Quadratic terms
  if (grepl("q", type)) {
    quadratic_terms <- paste0("I(", independent, "^2)")
    formula_parts <- paste(formula_parts, paste(quadratic_terms, collapse = " + "), sep = " + ")
  }

  # Product (interaction) terms
  if (grepl("p", type)) {
    if (length(independent) > 1) {
      interaction_terms <- apply(utils::combn(independent, 2), 2, paste, collapse = ":")
      formula_parts <- paste(formula_parts, paste(interaction_terms, collapse = " + "), sep = " + ")
    } else {
      message("'p' is only possible with 2 or more independent variables. Returning other combinations.")
    }
  }

  # Create the final formula
  formula <- paste0("~", gsub(pattern = "  \\+ ", replacement = "", x = formula_parts))

  # Add the categorical variable if provided
  if (!is.null(categorical_var)) {
    formula <- paste0(formula, " + categorical(", categorical_var, ")")
  }

  return(formula)
}
