#' Print Method for Calibration Grid GLM
#'
#' @param x An object of class `calibration_grid_glm`.
#' @export
print.calibration_grid_glm <- function(x, ...) {
  cat("Calibration Grid (GLM)\n")
  cat("Number of combinations:", nrow(x$data), "\n")
  cat("Print (n = 10):\n")
  print(head(x$data, 10))
}

#' Print Method for Calibration Grid GLMNET
#'
#' @param x An object of class `calibration_grid_glmnet`.
#' @export
print.calibration_grid_glmnet <- function(x, ...) {
  cat("Calibration Grid (GLMNET)\n")
  cat("Number of combinations:", nrow(x$data), "\n")
  cat("Print (n = 10):\n")
  print(head(x$data, 10))
}

#' Summary Method for Calibration Grid GLM
#'
#' @param object An object of class `calibration_grid_glm`.
#' @param n An integer specifying the number of formulas to display. Default = 10.
#' @export
summary.calibration_grid_glm <- function(object, n = 10, ...) {

  f <- unique(object$data$Formulas)
  if (length(f) < n){ n = length(f) }

  cat("Summary of Calibration Grid (GLM)\n")
  cat("Number of formulas used:", length(f), "\n")
  cat("Features used:", paste(unique(object$data$Features), collapse = ", "), "\n")
  cat(paste0("Print (n = ", n, "):\n"))
  print(head(unique(object$data$Formulas), n))
}

#' Summary method for calibration grid GLMNET
#'
#' @param object An object of class `calibration_grid_glmnet`.
#' @param n An integer specifying the number of formulas to display. Default = 10.
#' @export
summary.calibration_grid_glmnet <- function(object, n = 10,  ...) {

  f <- unique(object$data$Formulas)
  if (length(f) < n){ n = length(f) }

  cat("Summary of Calibration Grid (GLMNET)\n")
  cat("Number of combinations:", nrow(object$data), "\n")
  cat("Number of formulas used:", length(f), "\n")
  cat("Features used:", paste(unique(object$data$Features), collapse = ", "), "\n")
  cat("Regularization parameters:", paste(unique(object$data$regm), collapse = ", "), "\n")
  cat(paste0("Print (n = ", n, "):\n"))
  print(head(f, n))
}

#' Calibration grid generation
#'
#' Generates a grid of possible formulas and parameters for calibration using either GLM or GLMNET.
#'
#' @param swd A data frame containing species presence/absence and environmental variables.
#' @param var_names A character vector of variable names. Ignored if swd is provided.
#' @param min_number The minimum number of variables in a combination.
#' @param min_continuous The minimum number of continuous variables required in a combination.
#' @param categorical_var A character vector specifying categorical variables.
#' @param features A character vector of feature types (e.g., "l", "q", "p", "t", "h").
#' Note that threshold ("t") and hinges ("h") work only for glmnet.
#' @param model_type A character string specifying the model type to use, either "glm" or "glmnet".
#' @param regm A numeric vector of regularization parameters for glmnet.
#' @param write_file Logical, if TRUE the output grid is saved to a file.
#' @param file_name The name of the file to save the output grid if `write_file = TRUE`.
#' @return A data frame containing the calibration grid with possible formulas and parameters.
#' @export

calibration_grid <- function(swd = NULL,
                             var_names = NULL,
                             min_number = 2,
                             min_continuous = NULL,
                             categorical_var = NULL,
                             features = c("l", "q", "lq", "lqp", "p"),
                             model_type = c("glm", "glmnet"),
                             regm = c(0.1, 1, 2, 3, 5),
                             write_file = FALSE,
                             file_name = NULL) {

  # Validate the model_type input
  model_type <- match.arg(model_type)


  # Determine variable names if var_names is provided
  if(!is.null(swd)) {
    categorical_var = swd$categorical_variables
  }

  if (model_type == "glm") {
    # Call the GLM-specific function
    cal_grid_data <- calibration_grid_glm(swd = swd,
                                          var_names = var_names,
                                          min_number = min_number,
                                          min_continuous = min_continuous,
                                          categorical_var = categorical_var,
                                          features = features)
    cal_grid <- list(data = cal_grid_data)
    class(cal_grid) <- "calibration_grid_glm"

  } else if (model_type == "glmnet") {
    # Call the GLMNET-specific function
    cal_grid_data <- calibration_grid_glmnetmx(swd = swd,
                                               var_names = var_names,
                                               min_number = min_number,
                                               min_continuous = min_continuous,
                                               categorical_var = categorical_var,
                                               features = features,
                                               regm = regm)
    cal_grid <- list(data = cal_grid_data)
    class(cal_grid) <- "calibration_grid_glmnet"
  }

  # Optionally save the cal_grid to a file
  if (write_file) {
    if (is.null(file_name)) {
      file_name <- paste0("calibration_grid_", model_type)
    }
    saveRDS(cal_grid, paste0(file_name, ".RDS"))
  }

  return(cal_grid)
}

#' Calibration Grid Generation using GLMNET
#'
#' Generates a grid of possible formulas and parameters for calibration using GLMNET.
#'
#' @param swd A data frame containing species presence/absence and environmental variables.
#' @param var_names A character vector of variable names. Ignored if swd is provided.
#' @param categorical_var A character vector specifying categorical variables.
#' @param min_number The minimum number of variables in a combination.
#' @param min_continuous The minimum number of continuous variables required in a combination.
#' @param features A character vector of feature types (e.g., "l", "q", "p", "t", "h").
#' @param regm A numeric vector of regularization parameters for glmnet.
#' @return A data frame containing the calibration grid with possible formulas and parameters.
#' @export

calibration_grid_glmnetmx <- function(swd = NULL,
                                      var_names = NULL,
                                      categorical_var = NULL,
                                      features = c("l", "q", "lq", "lqp", "p"),
                                      min_number = 2,
                                      min_continuous = NULL,
                                      regm = c(0.1, 1, 2, 3, 5)){

  # Determine variable names if var_names is provided
  if(!is.null(swd) && is.null(var_names)) {
    var_names <- colnames(sp_swd$calibration_data[, -1, drop = FALSE])
  }

  # Get variable combinations
  var_comb <- enmpa:::aux_var_comb(var_names = var_names, minvar = min_number)

  # Remove combinations according to minimum number of continuous variables
  if(!is.null(min_continuous)){
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
#'
#' Auxiliary function to prepare GLMNET formulas based on the specified feature types.
#'
#' @param independent A character vector of independent variables.
#' @param type A character string specifying the type of formula to generate.
#' @param categorical_var A character vector specifying categorical variables.
#' @param minvar Minimum number of variables in a combination.
#' @param maxvar Maximum number of variables in a combination.
#' @return A character string representing the formula.
#' @export

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
#'
#' Generates a grid of possible formulas and parameters for calibration using GLM.
#'
#' @param swd A data frame containing species presence/absence and environmental variables.
#' @param var_names A character vector of variable names. Ignored if swd is provided.
#' @param min_number The minimum number of variables in a combination.
#' @param min_continuous The minimum number of continuous variables required in a combination.
#' @param categorical_var A character vector specifying categorical variables.
#' @param features A character vector of feature types (e.g., "l", "q", "p").
#' @return A data frame containing the calibration grid with possible formulas and parameters.
#' @export

calibration_grid_glm <- function(swd = NULL,
                                 var_names = NULL,
                                 min_number = 2,
                                 min_continuous = NULL,
                                 categorical_var = NULL,
                                 features = c("l", "q", "lq", "lqp", "p")){

  # Determine variable names if var_names is provided
  if(!is.null(swd) && is.null(var_names)) {
    var_names <- colnames(sp_swd$calibration_data[, -1, drop = FALSE])
  }

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
#'
#' Auxiliary function to prepare GLM formulas based on the specified feature types.
#'
#' @param independent A character vector of independent variables.
#' @param type A character string specifying the type of formula to generate.
#' @param categorical_var A character vector specifying categorical variables.
#' @return A character string representing the formula.
#' @export

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
