#' Print Method for prepared_data Class
#' @export
print.prepared_data <- function(x, ...) {
  cat("prepared_data object summary\n")
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
    cat("  - Variables included:", paste(x$pca$vars_in, collapse = ", "), "\n")
    cat("  - Number of PCA components:", length(x$pca$sdev), "\n")
  } else {
    cat("PCA Information: PCA not performed\n")
  }

  if (!is.null(x$weights)) {
    cat("Weights Information:\n")
    cat("  - Weights provided: Yes\n")
  } else {
    cat("Weights Information: No weights provided\n")
  }

  algorithm <- x$algorithm

  #Print formula grid
  if(algorithm == "glmnet"){
    cat("Calibration Grid (GLMNET)\n")
    cat("  - Features used:", paste(unique(x$formula_grid$Features), collapse = ", "), "\n")
    cat("  - Reg. multipliers used:", paste(unique(x$formula_grid$reg_mult), collapse = ", "), "\n")
    cat("  - Number of combinations:", nrow(x$formula_grid), "\n")
    cat("  - Print formulas (n = 5):\n")
    print(head(x$formula_grid, 5))
  }

  if(algorithm == "glm"){
    cat("Calibration Grid (GLM)\n")
    cat("  - Features used:", paste(unique(x$formula_grid$Features), collapse = ", "), "\n")
    cat("  - Number of combinations:", nrow(x$formula_grid), "\n")
    cat("  - Print (n = 5):\n")
    print(head(x$formula_grid, 5))
  }
}



#' Print Method for calibration_results Class
#' @export
print.calibration_results <- function(x, ...){
  cat(paste0("calibration_results object summary (", x$algorithm,")\n"))
  cat("=============================================================\n")
  cat("Species:", x$species, "\n")

  cat("Number of candidate models:", nrow(x$calibration_results$Summary), "\n")

  #Print summary
  cat("  - Models removed because they failed to fit:", length(x$summary$Errors), "\n")
  cat("  - Models removed with concave curves:", length(x$summary$Concave), "\n")
  cat("  - Models removed with non-significant values of pROC:",
      length(x$summary$Non_sig_pROC), "\n")
  cat("  - Models removed with omission rate >", paste0(x$omission_rate, "%:"),
      length(x$summary$High_omission_rate), "\n")
  cat("  - Models removed with delta AIC >", paste0(x$summary$delta_AIC, ":"),
      length(x$summary$High_AIC), "\n")
  cat("Selected models:", nrow(x$selected_models), "\n")

  cat("  - Print selected models (n = 5):\n")
  print(head(x$selected_models))
}



#' Print Method for fitted_models Class
#' @export
print.fitted_models <- function(x, ...){
  cat("fitted_models object summary\n")
  cat("==========================\n")
  cat("Species:", x$species, "\n")

  cat("Model type:", x$algorithm, "\n")
  cat("Number of fitted models:", length(x$Models), "\n")
  #Get number of replicates
  nr <- sum(grepl("Rep", names(x$Models[[1]])))
  cat("Models fitted with", nr, "replicates")
}



#' Print Method for prepared_proj Class
#' @export
print.projection_data <- function(x, ...){
  cat("projection_data object summary\n")
  cat("=============================\n")

  #Get times
  all_times <- x[c("Present", "Past", "Future")]
  all_times <- names(all_times[sapply(all_times, function(x) !is.null(x))])

  cat("Variables prepared to project models for",
      paste(all_times, collapse = " and "), "\n")

  #If Past...
  if("Past" %in% all_times){
    p <- names(x[["Past"]]) #Get periods
    g <- names(x[["Past"]][[1]]) #Get gcms
    cat("Past projections contain the following periods and GCMs:\n")
    cat("  - Periods:", paste(p, collapse = " | "), "\n")
    cat("  - GCMs:", paste(g, collapse = " | "), "\n")
  }

  #If Future...
  if("Future" %in% all_times){
    p <- names(x[["Future"]]) #Get periods
    s <- names(x[["Future"]][[1]]) #Get scenarios
    g <- names(x[["Future"]][[1]][[1]]) #Get gcms
    cat("Future projections contain the following periods, scenarios and GCMs:\n")
    cat("  - Periods:", paste(p, collapse = " | "), "\n")
    cat("  - Scenarios:", paste(s, collapse = " | "), "\n")
    cat("  - GCMs:", paste(g, collapse = " | "), "\n")
  }

  #Get folder
  f <- dirname(x[[1]][[1]])
  cat("All variables are located in the following root directory:\n")
  cat(f)
}



#' Print Method for model_projections Class
#' @export
print.model_projections <- function(x, ...){
  cat("model_projections object summary\n")
  cat("================================\n")

  #Get times
  all_times <- na.omit(unique(x$paths$Time))


  cat("Models projected for",
      paste(all_times, collapse = " and "), "\n")

  #If Past...
  if("Past" %in% all_times){
    p <- unique(x$paths$Period[x$paths$Time == "Past"])  #Get periods
    g <- unique(x$paths$GCM[x$paths$Time == "Past"]) #Get gcms
    cat("Past projections contain the following periods and GCMs:\n")
    cat("  - Periods:", paste(p, collapse = " | "), "\n")
    cat("  - GCMs:", paste(g, collapse = " | "), "\n")
  }

  #If Future...
  if("Future" %in% all_times){
    p <- unique(x$paths$Period[x$paths$Time == "Future"])  #Get periods
    s <- unique(x$paths$ssp[x$paths$Time == "Future"]) #Get scenarios
    g <- unique(x$paths$GCM[x$paths$Time == "Future"]) #Get gcms
    cat("Future projections contain the following periods, scenarios and GCMs:\n")
    cat("  - Periods:", paste(p, collapse = " | "), "\n")
    cat("  - Scenarios:", paste(s, collapse = " | "), "\n")
    cat("  - GCMs:", paste(g, collapse = " | "), "\n")
  }

  #Get folder
  f <- dirname(dirname(x$paths$output_path[1]))
  cat("All raster files containing the projection results are located in the following root directory:\n",
      f)
}
