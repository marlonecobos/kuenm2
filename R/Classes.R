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


# prepared_data Class Constructor
new_prepared_data <- function(species, calibration_data, formula_grid,
                             part_data, partition_method, n_replicates,
                             train_proportion, data_xy,
                             continuous_variables,
                             categorical_variables, weights, pca, algorithm) {
  data <- list(
    species = species,
    calibration_data = calibration_data,
    formula_grid = formula_grid,
    part_data = part_data,
    partition_method = partition_method,
    n_replicates = n_replicates,
    train_proportion = train_proportion,
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
                              algorithm,
                              partition_method,
                              n_replicates,
                              train_proportion){
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
    algorithm = algorithm,
    partition_method = partition_method,
    n_replicates = n_replicates,
    train_proportion = train_proportion
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
