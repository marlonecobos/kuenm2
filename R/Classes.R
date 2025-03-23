new_explore_back <- function(hist_m, hist_bg, hist_pr, mean_m, mean_bg, mean_pr,
                             cl_m, cl_bg, cl_pr, range_m, range_bg, range_pr) {
  # error checking

  val <- list(hist_m = hist_m, hist_bg = hist_bg, hist_pr = hist_pr,
              mean_m = mean_m, mean_bg = mean_bg, mean_pr = mean_pr,
              cl_m = cl_m, cl_bg = cl_bg, cl_pr = cl_pr, range_m = range_m,
              range_bg = range_bg, range_pr = range_pr)

  class(val) <- "explore_back"
  return(val)
}


new_explore_list <- function(summary, exploration_stats,
                             continuous_variables,
                             categorical_variables) {
  # error checking

  val <- list(summary = summary, exploration_stats = exploration_stats,
              continuous_variables = continuous_variables,
              categorical_variables = categorical_variables)

  class(val) <- "explore_calibration"
  return(val)
}


new_occ_cal <- function(summary, occurrences, block_limits = NULL,
                        calibration_regime = NULL) {
  # error checking

  val <- list(summary = summary, occurrences = occurrences,
              block_limits = block_limits,
              calibration_regime = calibration_regime)
  class(val) <- "occ_cal"
  return(val)
}



new_back_cal <- function(summary, background, background_bias = NULL,
                         training_regime = NULL) {
  # error checking

  val <- list(summary = summary, background = background,
              background_bias = background_bias,
              training_regime = training_regime)
  class(val) <- "back_cal"
  return(val)
}


new_lambdas <- function(lambdas_path){
  lambdas <- readLines(lambdas_path)
  maxcalcsV <- c("linearPredictorNormalizer",
                 "densityNormalizer",
                 "numBackgroundPoints",
                 "entropy")
  maxids <-  sapply(seq_along(maxcalcsV), function(x){
    ids <- stringr::str_detect(string = lambdas,maxcalcsV[x])
    return(which(ids))
  })
  if(length(maxids)==0L){
    stop("Please provide a valid lambdas file")
  }
  maxcalcs <- lapply(maxids, function(x){
    val <- stringr::str_split(lambdas[x],", ")[[1]]
    v <- as.numeric(val[2])
    return(v)
  })

  names(maxcalcs) <- maxcalcsV
  lambvals <- data.frame(stringr::str_split(string = lambdas[-maxids],
                                            pattern = "[,]",simplify = T))
  lambvals[,-1] <- as.numeric(as.matrix(lambvals[,-1]))
  names(lambvals) <- c("variable","lambda","min","max")
  feature <- sapply(lambvals$variable, function(x){
    qq <- which(stringr::str_detect(string =x ,"\\^"))
    pd <- which(stringr::str_detect(string =x ,"\\*"))
    fh <- which(stringr::str_detect(string =x ,"\\'"))
    rh <- which(stringr::str_detect(string =x ,"\\`"))
    th <- which(stringr::str_detect(string =x ,"\\<"))
    if(length(qq)){
      return("quadratic")
    } else if(length(pd)){
      return("product")
    } else if(length(fh)){
      return("forward_hinge")
    } else if(length(rh)){
      return("reverse_hinge")
    } else if(length(rh)){
      return("reverse_hinge")
    } else if(length(th)){
      return("threshold")
    } else{
      return("linear")
    }
  })
  lambvals$feature <- feature
  res <- list(lambdas_df=lambvals,maxmeta=maxcalcs)
  class(res) <- c("lambdas")
  return(res)
}



# prepared_data Class Constructor
new_prepared_data <- function(species, calibration_data, formula_grid,
                             kfolds, data_xy, continuous_variables,
                             categorical_variables, weights, pca, algorithm) {
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
    algorithm = algorithm
  )
  class(data) <- "prepared_data"
  return(data)
}


# calibration_results Class Constructor
new_calibration_results <- function(prepared_data, calibration_results,
                                    omission_rate, addsamplestobackground,
                                    weights, selected_models, summary) {

  fm <- c(prepared_data, calibration_results = calibration_results,
          omission_rate = omission_rate,
          addsamplestobackground = addsamplestobackground,
          weights = list(weights),
          selected_models = selected_models,
          summary = summary)

  class(fm) <- "calibration_results"
  return(fm)
}



# fitted_models Class Constructor
new_fitted_models <- function(species,
                              Models,
                              calibration_data,
                              continuous_variables,
                              categorical_variables,
                              selected_models,
                              weights,
                              pca,
                              addsamplestobackground,
                              omission_rate,
                              thresholds,
                              algorithm){
  data <- list(
    species = species,
    Models = Models,
    calibration_data = calibration_data,
    continuous_variables = continuous_variables,
    categorical_variables = categorical_variables,
    selected_models = selected_models,
    weights = weights,
    pca = pca,
    addsamplestobackground = addsamplestobackground,
    omission_rate = omission_rate,
    thresholds = thresholds,
    algorithm = algorithm
  )
  class(data) <- "fitted_models"
  return(data)
}



# prepared_proj Class Constructor
new_projection_data <- function(res_present, res_past, res_future, raster_pattern,
                                variables, pca){
  data <- list(Present = res_present,
               Past = res_past,
               Future = res_future,
               raster_pattern = raster_pattern,
               variables = variables,
               pca = pca)

  class(data) <- "projection_data"
  return(data)
}



# prepared_proj Class Constructor
new_model_projections <- function(paths, thresholds){
  data <- list(paths = paths,
               thresholds = thresholds)
  class(data) <- "model_projections"
  return(data)
}
