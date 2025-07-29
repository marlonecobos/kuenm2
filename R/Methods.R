#' Predict method for glmnet_mx (maxnet) models
#'
#' @usage
#' predict.glmnet_mx(object, newdata, clamp = FALSE,
#'                   type = c("link", "exponential", "cloglog", "logistic",
#'                   "cumulative"))
#'
#' @name predict
#' @aliases predict,kuenm2_glmnet_mx-method
#'
#' @param object a glmnet_mx object.
#' @param newdata data to predict on.
#' @param clamp (logical) whether to clamp predictions. Default = FALSE.
#' @param type (character) type of prediction to be performed. Options are:
#' "link", "exponential", "cloglog", "logistic", and cumulative. Defaults to "link" if not defined.
#'
#' @export
#'
#' @returns
#' A glmnet_mx (maxnet) prediction.

predict.glmnet_mx <- function(object, newdata, clamp = FALSE,
                              type = c("link", "exponential", "cloglog",
                                       "logistic", "cumulative")) {
  if (clamp) {
    for (v in intersect(names(object$varmax), names(newdata))) {
      newdata[, v] <- pmin(pmax(newdata[, v], object$varmin[v]),
                           object$varmax[v])
    }
  }

  terms <- sub("hinge\\((.*)\\):(.*):(.*)$", "hingeval(\\1,\\2,\\3)",
               names(object$betas))
  terms <- sub("categorical\\((.*)\\):(.*)$", "categoricalval(\\1,\"\\2\")",
               terms)
  terms <- sub("thresholds\\((.*)\\):(.*)$", "thresholdval(\\1,\\2)",
               terms)
  f <- formula(paste("~", paste(terms, collapse = " + "), "-1"))

  mm <- model.matrix(f, data.frame(newdata))

  if (clamp) {
    mm <- t(pmin(pmax(t(mm), object$featuremins[names(object$betas)]),
                 object$featuremaxs[names(object$betas)]))
  }

  link <- (mm %*% object$betas) + object$alpha
  type <- match.arg(type)

  # return prediction
  if (type == "link") {
    return(link)
  }
  if (type == "exponential") {
    return(exp(link))
  }
  if (type == "cloglog") {
    return(1 - exp(0 - exp(object$entropy + link)))
  }
  if (type == "logistic") {
    return(1/(1 + exp(-object$entropy - link)))
  }
  if (type == "cumulative") {
   return(cumulative_predictions(predictions = exp(link)))
  }
}

#' Print method for kuenm2 objects
#'
#' @name print
#' @aliases print,kuenm2_prepared_data-method
#' @aliases print,kuenm2_calibration_results-method
#' @aliases print,kuenm2_fitted_models-method
#' @aliases print,kuenm2_projection_data-method
#' @aliases print,kuenm2_model_projections-method
#'
#' @param x an object of any of these classes: `prepared_data`,
#' `calibration_results`, `fitted_models`, `projection_data`, or
#' `model_projections`.
#' @param ... additional arguments affecting the summary produced. Ignored in
#' these functions.
#'
#' @export
#'
#' @rdname print
#'
#' @returns
#' A printed version of the object that summarizes the main elements contained.

# Do we need one for explore_back and explore_list????

print.prepared_data <- function(x, ...) {
  cat("prepared_data object summary\n")
  cat("============================\n")
  cat("Species:", x$species, "\n")
  cat("Number of Records:", nrow(x$calibration_data), "\n")
  cat("  - Presence:", table(x$calibration_data$pr_bg)[2], "\n")
  cat("  - Background:", table(x$calibration_data$pr_bg)[1], "\n")

  cat("Partition Method:", x$partition_method, "\n")
  if(x$partition_method %in% c("kfolds", "leave-one-out"))
    cat("  - Number of kfolds:", length(x$part_data), "\n")
  if(x$partition_method%in% c("subsample", "bootstrap")){
    cat("  - Number of replicates:", x$n_replicates, "\n")
    cat("  - Train proportion:", x$train_proportion, "\n")
  }

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
    cat("Weights:\n")
    cat("  - Weights provided: Yes\n")
  } else {
    cat("Weights: No weights provided\n")
  }

  #Print model calibration parameters
  cat("Calibration Parameters:\n")
  cat("  - Algorithm:", x$algorithm, "\n")
  cat("  - Number of candidate models:", nrow(x$formula_grid), "\n")
  cat("  - Features classes (responses):", paste(unique(x$formula_grid$Features), collapse = ", "), "\n")

  if (x$algorithm == "maxnet") {
    cat("  - Regularization multipliers:", paste(unique(x$formula_grid$R_multiplier), collapse = ", "), "\n")
  }
}



#' @rdname print
#' @export

print.calibration_results <- function(x, ...) {
  cat(paste0("calibration_results object summary (", x$algorithm,")\n"))
  cat("=============================================================\n")
  cat("Species:", x$species, "\n")

  cat("Number of candidate models:", nrow(x$calibration_results$Summary), "\n")

  #Print summary
  cat("  - Models removed because they failed to fit:", length(x$summary$Errors), "\n")
  cat("  - Models identified with concave curves:",
      sum(x$calibration_results$Summary$Is_concave), "\n")
  if(x$summary$Remove_concave){
    cat("  - Model with concave curves removed", "\n")
  } else {cat("  - Model with concave curves not removed", "\n")}
  cat("  - Models removed with non-significant values of pROC:",
      length(x$summary$Non_sig_pROC), "\n")
  cat("  - Models removed with omission error >", paste0(x$omission_rate, "%:"),
      length(x$summary$High_omission_rate), "\n")
  cat("  - Models removed with delta AIC >", paste0(x$summary$delta_AIC, ":"),
      length(x$summary$High_AIC), "\n")
  cat("Selected models:", nrow(x$selected_models), "\n")

  cat("  - Up to 5 printed here:\n")
  show_col <- c("ID", "Formulas", "Features", "R_multiplier",
                paste0("pval_pROC_at_", x$omission_rate, ".mean"),
                paste0("Omission_rate_at_", x$omission_rate, ".mean"),
                "dAIC", "Parameters")

  if (x$algorithm == "maxnet") {
    print(head(x$selected_models[, show_col]))
  } else {
    print(head(x$selected_models[, show_col[-4]]))
  }
}



#' @rdname print
#' @export

print.fitted_models <- function(x, ...) {
  cat("fitted_models object summary\n")
  cat("============================\n")
  cat("Species:", x$species, "\n")

  cat("Algortihm:", x$algorithm, "\n")
  cat("Number of fitted models:", length(x$Models), "\n")
  #Get number of replicates
  nr <- sum(grepl("Partition", names(x$Models[[1]])))
  if (nr > 0) {
    cat("Models fitted with", nr, "replicates")
  } else {
    cat("Only full models fitted, no replicates")
  }
}



#' @rdname print
#' @export

print.projection_data <- function(x, ...) {
  cat("projection_data object summary\n")
  cat("=============================\n")

  #Get times
  all_times <- x[c("Present", "Past", "Future")]
  all_times <- names(all_times[sapply(all_times, function(x) !is.null(x))])

  cat("Variables prepared to project models for",
      paste(all_times, collapse = " and "), "\n")

  #If Past...
  if ("Past" %in% all_times) {
    p <- names(x[["Past"]]) #Get periods
    g <- names(x[["Past"]][[1]]) #Get gcms
    cat("Past projections contain the following periods and GCMs:\n")
    cat("  - Periods:", paste(p, collapse = " | "), "\n")
    cat("  - GCMs:", paste(g, collapse = " | "), "\n")
  }

  #If Future...
  if ("Future" %in% all_times) {
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



#' @rdname print
#' @export

print.model_projections <- function(x, ...) {
  cat("model_projections object summary\n")
  cat("================================\n")

  #Get times
  all_times <- na.omit(unique(x$paths$Time))


  cat("Models projected for",
      paste(all_times, collapse = " and "), "\n")

  #If Past...
  if ("Past" %in% all_times) {
    p <- unique(x$paths$Period[x$paths$Time == "Past"])  #Get periods
    g <- unique(x$paths$GCM[x$paths$Time == "Past"]) #Get gcms
    cat("Past projections contain the following periods and GCMs:\n")
    cat("  - Periods:", paste(p, collapse = " | "), "\n")
    cat("  - GCMs:", paste(g, collapse = " | "), "\n")
  }

  #If Future...
  if ("Future" %in% all_times) {
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
