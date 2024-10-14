#' Variable importance
#'
#' @usage
#' var_importance(fitted, modelID = NULL)
#'
#'
#' @return
#' A data.frame containing the relative contribution of each variable. An
#' identification for distinct models is added if `fitted` contains multiple
#' models.
#'
#' @export
#'
#' @importFrom stats update as.formula deviance coef glm
#'

var_importance <- function(fitted, modelID = NULL){
  # initial tests
  if (missing(fitted)) {
    stop("Argument 'model' must be defined.")
  }

  list_models <- fitted[["Models"]]
  model_info  <- fitted[["selected_models"]]
  data        <- fitted[["calibration_data"]]
  model_type  <- fitted[["model_type"]]

  if (is.null(modelID)){
    models <- names(list_models)
    aux <- lapply(models, function(y) {
      var_importance_indmx(
        model = list_models[[y]][["Full_model"]],
        p = data$pr_bg,
        data = data,
        f = model_info[model_info$ID == gsub("Model_", "", y), "Formulas"],
        rm = model_info[model_info$ID == gsub("Model_", "", y), "regm"],
        model_type = model_type
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
      rm = model_info[model_info$ID == gsub("Model_", "", modelID), "regm"],
      model_type = model_type
    )
  }

  return(tab_contr)
}


# Aux function to evaluate the Variable Contribution of the predictors
#

#to get deviance of a model after excluding predictors
get_red_devmx <- function(reduce_var, p, data, f, rm, model_type) {
  if (model_type == "glmnet") {
    reduce_model <- glmnet_mx(
      p = p, data = data, regmult = rm,
      f = as.formula(paste(f, " - ", reduce_var)),
      calculate_AIC = FALSE,
      addsamplestobackground = TRUE
    )
    return(deviance(reduce_model)[200])

  } else if (model_type == "glm") {
    reduce_model <- suppressWarnings(
      glm(
        formula = as.formula(paste("pr_bg", f, " - ", reduce_var)),
        family = binomial(link = "cloglog"), data = data
      ))
    return(deviance(reduce_model))
  }
}


# get variable contribution for an individual model
var_importance_indmx <- function(model, p, data, f, rm, model_type) {

  # initial tests
  if (missing(model)) {
    stop("Argument 'model' must be defined.")
  }

  # deviance of the full model
  dev_full <- if (model_type == "glmnet") {
    deviance(model)[200]
  } else {
    deviance(model)
  }

  coefs <- if (model_type == "glmnet") {
    names(model$betas)
  } else if (model_type == "glm") {
    names(coef(model)[-1])
  }

  # deviance of the reduced models
  dev_reduction <- sapply(coefs, function(variable) {

    #abs(dev_full - get_red_devmx(variable, p, data, f, rm)) # negative values?
    dev_full - get_red_devmx(variable, p, data, f, rm, model_type)
  })

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


