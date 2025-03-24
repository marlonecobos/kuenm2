#' Variable importance
#'
#' @usage
#' var_importance(models, modelID = NULL, parallel = FALSE, ncores = 2,
#'                parallelType = "doSNOW", progress_bar = TRUE, verbose = TRUE)
#'
#' @param models an object of class `fitted_models` returned by the
#' \code{\link{fit_selected}}() function.
#' @param modelID (character). Default = NULL.
#' @param parallel (logical) whether to calculate importance in parallel.
#' Default is FALSE.
#' @param ncores (numeric) the number of cores to use for parallel processing.
#' Default is 2. This is only applicable if `parallel = TRUE`.
#' @param parallelType (character) the parallelization package to use:
#' either "doParallel" or "doSNOW". Default is "doSNOW". This is only applicable
#' if `parallel = TRUE`.
#' @param progress_bar (logical) whether to display a progress bar during processing.
#' Default is TRUE.
#' @param verbose (logical) whether to display detailed messages during processing.
#' Default is TRUE.
#'
#' @return
#' A data.frame containing the relative contribution of each variable. An
#' identification for distinct models is added if `fitted` contains multiple
#' models.
#'
#' @export
#'
#' @importFrom stats update as.formula deviance coef glm
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom foreach foreach `%dopar%`
#'
#' @examples
#' ##Example with maxnet
#' # Import example of fitted_models (output of fit_selected())
#' data("fitted_model_maxnet", package = "kuenm2")
#'
#' # Variable importance
#' imp_maxnet <- var_importance(models = fitted_model_maxnet)
#' # Plot using enmpa package
#' enmpa::plot_importance(imp_maxnet)
#'
#' ##Example with glm
#' # Import example of fitted_models (output of fit_selected())
#' data("fitted_model_glm", package = "kuenm2")
#'
#' # Variable importance
#' imp_glm <- var_importance(models = fitted_model_glm)
#' # Plot using enmpa package
#' enmpa::plot_importance(imp_glm)
#'
var_importance <- function(models, modelID = NULL,
                           parallel = FALSE,
                           ncores = 2,
                           parallelType = "doSNOW",
                           progress_bar = TRUE,
                           verbose = TRUE){
  # initial tests
  if (missing(models)) {
    stop("Argument 'model' must be defined.")
  }

  if (!is.null(modelID)){
    if (!modelID %in% names(models[["Models"]])){
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
  algorithm  <- models[["algorithm"]]

  if (is.null(modelID)){
    models <- names(list_models)
    aux <- lapply(seq_along(models), function(y) {
      if(verbose){
        message("\nCalculating variable contribution for model ", y, " of ",
                length(models))
      }
      var_importance_indmx(
        model = list_models[[y]][["Full_model"]],
        p = data$pr_bg,
        data = data,
        f = model_info[model_info$ID == gsub("Model_", "", models[y]), "Formulas"],
        rm = model_info[model_info$ID == gsub("Model_", "", models[y]), "reg_mult"],
        algorithm = algorithm,
        parallel = parallel,
        ncores = ncores,
        parallelType = parallelType,
        progress_bar = progress_bar,
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
      f = model_info[model_info$ID == gsub("Model_", "", modelID), "Formulas"],
      rm = model_info[model_info$ID == gsub("Model_", "", modelID), "reg_mult"],
      algorithm = algorithm,
      parallel = parallel,
      ncores = ncores,
      parallelType = parallelType,
      progress_bar = progress_bar,
      verbose = verbose
    )
  }

  return(tab_contr)
}


# Aux function to evaluate the Variable Contribution of the predictors
#

#to get deviance of a model after excluding predictors
get_red_devmx <- function(reduce_var, p, data, f, rm, algorithm) {
  if (algorithm == "maxnet") {
    reduce_model <- glmnet_mx(
      p = p, data = data, regmult = rm,
      f = as.formula(paste(f, " - ", reduce_var)),
      calculate_AIC = FALSE,
      addsamplestobackground = TRUE
    )
    return(deviance(reduce_model)[200])

  } else if (algorithm == "glm") {
    reduce_model <-
      glm_mx(formula = as.formula(paste("pr_bg", f, " - ", reduce_var)),
             data = data)
    return(deviance(reduce_model))
  }
}


# get variable contribution for an individual model
var_importance_indmx <- function(model, p, data, f, rm, algorithm,
                                 parallel,
                                 ncores,
                                 parallelType,
                                 progress_bar,
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

  #Fix categorical terms
  coefs <- unique(gsub("categorical\\(([^\\)]+)\\):[0-9]+", "categorical(\\1)",
                       coefs))

  # Progress bar setup
  if (progress_bar) {
    pb <- utils::txtProgressBar(0, length(coefs), style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
  }

  # Adjust parallelization based on number of tasks and cores
  if (length(coefs) == 1 & parallel) {
    parallel <- FALSE
  }
  if (length(coefs) < ncores & parallel) {
    ncores <- length(coefs)
  }

  # Setup parallel cluster
  if (parallel) {
    cl <- parallel::makeCluster(ncores)
    if (parallelType == "doParallel") {
      doParallel::registerDoParallel(cl)
      opts <- NULL
    } else if (parallelType == "doSNOW") {
      doSNOW::registerDoSNOW(cl)
      opts <- if (progress_bar) list(progress = progress) else NULL
    }
  } else {
    opts <- NULL
  }

  # Fit the best models (either in parallel or sequentially)
  if (parallel) {
    dev_reduction <- foreach::foreach(x = 1:length(coefs),
                                      .options.snow = opts,
                                      .combine = 'c'
    ) %dopar% {
      dev_full - get_red_devmx(coefs[x], p, data, f, rm, algorithm)
    }
  } else {
    dev_reduction <- c()
    for (x in 1:length(coefs)) {
      dev_reduction[x] <- dev_full - kuenm2:::get_red_devmx(coefs[x], p, data,
                                                            f, rm, algorithm)
      if (progress_bar) utils::setTxtProgressBar(pb, x)
    }
  }
  names(dev_reduction) <- coefs

  # Stop the cluster
  if (parallel) parallel::stopCluster(cl)

  # # deviance of the reduced models
  # dev_reduction <- sapply(coefs, function(variable) {
  #
  #   #abs(dev_full - get_red_devmx(variable, p, data, f, rm)) # negative values?
  #   dev_full - get_red_devmx(variable, p, data, f, rm, algorithm)
  # })

  deviance_importance <- dev_reduction / sum(dev_reduction)

  # preparing results
  tab_contr <- data.frame(predictor = names(deviance_importance),
                          stringsAsFactors = FALSE)
  tab_contr$contribution <- deviance_importance

  ord <- order(tab_contr$contribution, decreasing = TRUE)
  tab_contr <- tab_contr[ord,]
  tab_contr$cum_contribution <- cumsum(tab_contr$contribution)

  # returning results
  return(tab_contr)
}


#' Summary plot for variable importance in models
#'
#' @description
#' See details in \code{\link[enmpa]{plot_importance}}
#'
#' @usage
#' plot_importance(x, xlab = NULL, ylab = "Relative contribution",
#'                 main = "Variable importance", extra_info = TRUE, ...)
#' @export

plot_importance <- enmpa::plot_importance
