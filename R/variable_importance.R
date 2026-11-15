#' Variable importance
#'
#' @usage
#' variable_importance(models, modelID = NULL, by_terms = FALSE,
#'                     parallel = FALSE, ncores = NULL,
#'                     progress_bar = TRUE, progress_type = "timer",
#'                     use_lb = FALSE, verbose = TRUE)
#'
#' @param models an object of class `fitted_models` returned by the
#' \code{\link{fit_selected}}() function.
#' @param modelID (character). Default = NULL.
#' @param by_terms (logical) whether to calculate importance by model terms
#' (e.g., `bio1`, `I(bio1^2)`, `hinge(bio1)`) instead of aggregating by
#' variable. Default = FALSE.
#' @param parallel (logical) whether to calculate importance in parallel.
#' Default is FALSE.
#' @param ncores (numeric) number of cores to use for parallel processing.
#' Default is NULL and uses available cores - 1. This is only applicable if
#' `parallel = TRUE`.
#' @param progress_bar (logical) whether to display a progress bar during processing.
#' Default is TRUE.
#' @param progress_type (character) type of the progress bar. See
#' `pbapply::pbtypes()` for available progress bar types depending on operating
#' system. Default is "timer".
#' @param use_lb (logical) whether to use load balancing when running in
#' parallel. Default is FALSE.
#' @param verbose (logical) whether to display detailed messages during processing.
#' Default is TRUE.
#'
#' @return
#' A data.frame containing the relative contribution of each variable (or term
#' if `by_terms = TRUE`). An identification for distinct models is added if
#' `fitted` contains multiple models.
#'
#' @export
#'
#' @importFrom stats update as.formula deviance coef glm
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom pbapply pboptions pbsapply
#'
#' @seealso
#' [plot_importance()]
#'
#' @examples
#' # Example with maxnet
#' # Import example of fitted_models (output of fit_selected())
#' data(fitted_model_maxnet, package = "kuenm2")
#'
#' # Variable importance
#' imp_maxnet <- variable_importance(models = fitted_model_maxnet)
#'
#' # Plot
#' plot_importance(imp_maxnet)
#'
#' # Example with glm
#' # Import example of fitted_models (output of fit_selected())
#' data(fitted_model_glm, package = "kuenm2")
#'
#' # Variable importance
#' imp_glm <- variable_importance(models = fitted_model_glm)
#'
#' # Plot
#' plot_importance(imp_glm)

variable_importance <- function(models, modelID = NULL,
                                by_terms = FALSE,
                                parallel = FALSE,
                                ncores = NULL,
                                progress_bar = TRUE,
                                progress_type = "timer",
                                use_lb = FALSE,
                                verbose = TRUE) {
  # initial tests
  if (missing(models)) {
    stop("Argument 'models' must be defined.")
  }

  if (!is.null(modelID)) {
    if (!modelID %in% names(models[["Models"]])) {
      stop(paste0(
        "The 'ModelID' is not correct, check the following: [",
        paste(names(models[["Models"]]), collapse = ", ")),
        "]"
      )
    }
  }

  list_models <- models[["Models"]]
  model_info  <- models[["selected_models"]]
  data        <- models[["calibration_data"]]
  algorithm   <- models[["algorithm"]]

  if (is.null(modelID)) {
    models <- names(list_models)
    aux <- lapply(seq_along(models), function(y) {
      if (verbose) {
        message("\nCalculating variable contribution for model ", y, " of ",
                length(models))
      }
      var_importance_indmx(
        model = list_models[[y]][["Full_model"]],
        p = data$pr_bg,
        data = data,
        form = model_info[model_info$ID == gsub("Model_", "", models[y]), "Formulas"],
        rm = model_info[model_info$ID == gsub("Model_", "", models[y]), "R_multiplier"],
        algorithm = algorithm,
        by_terms = by_terms,
        parallel = parallel,
        ncores = ncores,
        progress_bar = progress_bar,
        progress_type = progress_type,
        use_lb = use_lb,
        verbose = verbose
      )
    })

    tab_contr <- do.call(rbind, aux)[, 1:2]
    tab_contr$Models <- rep(models, times = sapply(aux, nrow))
    rownames(tab_contr) <- NULL

  } else{
    tab_contr <- var_importance_indmx(
      model = list_models[[modelID]][["Full_model"]],
      p = data$pr_bg,
      data = data,
      form = model_info[model_info$ID == gsub("Model_", "", modelID), "Formulas"],
      rm = model_info[model_info$ID == gsub("Model_", "", modelID), "R_multiplier"],
      algorithm = algorithm,
      by_terms = by_terms,
      parallel = parallel,
      ncores = ncores,
      progress_bar = progress_bar,
      progress_type = progress_type,
      use_lb = use_lb,
      verbose = verbose
    )
  }

  return(tab_contr)
}


#
# Aux function to evaluate the Variable Contribution of the predictors
#

# get variable contribution for an individual model
var_importance_indmx <- function(model, p, data, form, rm, algorithm,
                                 by_terms,
                                 parallel,
                                 ncores,
                                 progress_bar,
                                 progress_type = "timer",
                                 use_lb = FALSE,
                                 verbose) {

  # initial tests
  if (missing(model)) {
    stop("Argument 'model' must be defined.")
  }

  # deviance of the full model
  dev_full <- if (algorithm == "maxnet") {
    deviance(model)[200]
  } else {
    deviance(model)
  }

  coefs <- if (algorithm == "maxnet") {
    names(model$betas)
  } else if (algorithm == "glm") {
    names(coef(model)[-1])
  }

  # Fix categorical, threshold and hinge terms
  coefs <- ifelse(
    grepl("^(thresholds|hinge|categorical)\\(", coefs),
    sub(":.*$", "", coefs),
    coefs
  )
  coefs <- unique(coefs)

  # Determine what to iterate: terms or variables
  if  (by_terms) {
    # By terms: each term is evaluated individually
    items_to_eval <- as.list(coefs)
    names(items_to_eval) <- coefs
  } else {
    # By variable: group all terms by their base variable
    var_names_list <- extract_variable_name(coefs)
    unique_vars <- unique(unlist(var_names_list))

    items_to_eval <- lapply(unique_vars, function(v) {
      belongs_to_var <- sapply(var_names_list, function(vars) v %in% vars)
      coefs[belongs_to_var]
    })
    names(items_to_eval) <- unique_vars
  }

  n_items <- length(items_to_eval)

  if (n_items == 1 && parallel) {
    parallel <- FALSE}

  # Store and set progress bar options safely for CRAN
  old_pbopt <- pbapply::pboptions(
    type = if (progress_bar) progress_type else "none",
    use_lb = use_lb
  )
  on.exit(pbapply::pboptions(old_pbopt), add = TRUE)

  if (parallel) {
    # Parallelization setup
    if (is.null(ncores)) {
      ncores <- max(1L, parallel::detectCores() - 1L)
      # Check if function is running in CRAN...
      chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
      if (nzchar(chk) && chk == "TRUE") {
        ncores <- min(ncores, 2L)
      }
    }
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else {
    cl <- NULL
  }


  # Fit the reduced models (either in parallel or sequentially)
  dev_reduction <- pbapply::pbsapply(
    X = items_to_eval,
    FUN = function(item, dev_full, p, data, form, rm, algorithm) {
      dev_full - get_red_devmx(item, p, data, form, rm, algorithm)
    },
    dev_full = dev_full,
    p = p,
    data = data,
    form = form,
    rm = rm,
    algorithm = algorithm,
    cl = cl
  )
  names(dev_reduction) <- names(items_to_eval)


  deviance_importance <- dev_reduction / sum(dev_reduction)

  # preparing results
  tab_contr <- data.frame(predictor = names(deviance_importance),
                          contribution = as.numeric(deviance_importance),
                          stringsAsFactors = FALSE)

  ord <- order(tab_contr$contribution, decreasing = TRUE)
  tab_contr <- tab_contr[ord, ]
  tab_contr$cum_contribution <- cumsum(tab_contr$contribution)
  rownames(tab_contr) <- NULL

  # returning results
  return(tab_contr)
}


#
# Aux function to extract the base variable name from model terms
#

extract_variable_name <- function(terms) {
  # Returns a list where each element contains the variable(s) for that term
  lapply(terms, function(term) {
    #  I(var^2) -> var
    if (grepl("^I\\(", term)) {
      var <- gsub("^I\\(([^\\^\\)]+).*\\)$", "\\1", term)
      return(var)
    }

    # Handle hinge, threshold, categorical: func(var) -> var
    if (grepl("^(hinge|thresholds|categorical)\\(", term)) {
      var <- gsub("^(hinge|thresholds|categorical)\\(([^\\)]+)\\)$", "\\2", term)
      return(var)
    }

    # Handle interactions: var1:var2 -> return BOTH variables
    if (grepl(":", term)) {
      vars <- strsplit(term, ":")[[1]]
      return(vars)
    }

    # Plain variable name
    return(term)
  })
}


# to get deviance of a model after excluding predictors
# reduce_var can be a single term or a vector of terms to remove
get_red_devmx <- function(reduce_var, p, data, f, rm, algorithm) {
  # Build the formula reduction string (remove all terms)
  reduce_terms <- paste(reduce_var, collapse = " - ")

  if (algorithm == "maxnet") {
    reduce_model <- glmnet_mx(
      p = p, data = data, regmult = rm,
      f = as.formula(paste(f, " - ", reduce_terms)),
      calculate_AIC = FALSE,
      addsamplestobackground = TRUE
    )
    return(deviance(reduce_model)[200])

  } else if (algorithm == "glm") {
    reduce_model <-
      glm_mx(formula = as.formula(paste("pr_bg", f, " - ", reduce_terms)),
             data = data)
    return(deviance(reduce_model))
  }
}



#' Summary plot for variable importance in models
#'
#' @description
#' See details in \code{\link[enmpa]{plot_importance}}
#'
#' @param x data.frame output from \code{\link{variable_importance}}().
#' @param xlab (character) a label for the x axis.
#' @param ylab (character) a label for the y axis.
#' @param main (character) main title for the plot.
#' @param extra_info (logical) when results are from more than one model, it adds information about the number of models using each predictor and the mean contribution found.
#' @param ... additional arguments passed to barplot or boxplot.
#'
#' @usage
#' plot_importance(x, xlab = NULL, ylab = "Relative contribution",
#'                 main = "Variable importance", extra_info = TRUE, ...)
#'
#' @return
#' No return value, called for side effects (a barplot or boxplot depending on
#' the number of models considered..
#'
#' @export

plot_importance <- enmpa::plot_importance
