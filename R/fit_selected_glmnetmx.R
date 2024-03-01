# #Load packages
# library(terra)
# library(dplyr)
# library(pbapply)
# library(pbapply)
# library(foreach)
# library(parallel)
#
# #Load functions
# source("Functions/Metrics_Functions.R")
# source("Functions/eval_m.R")
# source("Functions/part_data.R")

#Prediction to dataframe or raster
#In df, each replicate in a column

# #Select 3 best models from candidate results
# cand_res <- read.csv("Models/Piper_fuligineum/candidate_results.csv")
# selected_models <- cand_res %>%
#   filter(is_concave == F, proc_pval.mean < 0.05,
#          Omission_rate_at_5.mean <= 0.05) %>% slice_min(AIC, n =3)
#
# #Import data to create the functions
# data <- read.csv("Models/Piper_fuligineum/occ_bg.csv")
# pr_bg <- "pr_bg"
# var_categorical = NULL
# replicates = TRUE
# n_replicates <- 5
# rep_type = "subsample"
# train_portion = 0.7
# write_models = TRUE
# out_dir = "Models/Piper_fuligineum/Best_models/" #Name of the folder to write final models
# parallel = TRUE
# ncores = 1
# progress_bar = TRUE
# parallelType = "doSNOW"
# verbose = TRUE
# to_export = c("aic_glmnetmx", "aic_maxnet", "eval_stats",
#               "get_formulas_maxnet",
#               "glmnet_mx", "kfold_part",
#               "maxnet.default.regularization",
#               "omrat_maxnet",
#               "predict.glmnet_mx", "empty_replicates",
#               "empty_summary",
#               "hinge", "hingeval", "thresholds",
#               "thresholdval", "categorical",
#               "categoricalval")




####Function to fit best models####
fit_selected_glmnetmx <- function(calibration_results,
                                  n_replicates = 5,
                                  rep_type = "kfold",
                                  train_portion = 0.7,
                                  write_models = FALSE, #Write files?
                                  file_name = NULL, #Name of the folder to write final models
                                  parallel = TRUE,
                                  ncores = 4,
                                  parallelType = "doSNOW",
                                  progress_bar = TRUE,
                                  verbose = TRUE) {
  #Args
  to_export <- c("aic_nk", "aic_ws",
                "eval_stats","glmnet_mx",
                "maxnet.default.regularization",
                "omrat","predict.glmnet_mx",
                "empty_replicates",
                "empty_summary", "hinge",
                "hingeval",
                "thresholds", "thresholdval",
                "categorical",
                "categoricalval", "fit_best_model", "n_replicates")

  #Extracts IDs from models
  m_ids <- calibration_results$selected_models$ID

  #Create grid of fitted models
  dfgrid <- expand.grid(models = m_ids, replicates = 1:n_replicates)
  n_tot <-  nrow(dfgrid)

  #Prepare data (index) to replicates
  if(n_replicates > 1) {
    #Partitioning data
    rep_data <- part_data(data = calibration_results$calibration_data,
                          pr_bg = "pr_bg",
                          train_portion = train_portion,
                          n_replicates = n_replicates,
                          method = rep_type)
  } else {
      rep_data <- NULL
    }


  #Show progress bar?
  if (isTRUE(progress_bar)) {
    pb <- txtProgressBar(0, n_tot, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n) }

  #Adjust parallelization according to the number of models and replicates
  if(n_tot == 1 & isTRUE(parallel)){
    parallel <- FALSE
  }
  if(n_tot < ncores & isTRUE(parallel)){
    ncores <- n_tot
  }

  #Set parallelization
  if(parallel) {
    #Make cluster
    cl <- parallel::makeCluster(ncores)

    if (parallelType == "doParallel") {
      doParallel::registerDoParallel(cl)
      opts <- NULL
    }

    if (parallelType == "doSNOW") {
      doSNOW::registerDoSNOW(cl)
      if (progress_bar)
        opts <- list(progress = progress)
      else opts <- NULL
    }
  } else {
    opts <- NULL}

  #In parallel (using %dopar%)
  if(parallel){
    best_models <- foreach(x = 1:n_tot,
                               .options.snow = opts,
                               .export = c(to_export)
    ) %dopar% {
      fit_best_model(x = x, dfgrid, calibration_results, data_x, n_replicates,
                     rep_data)
    }
  } else { #Not in parallel (using %do%)
    best_models <- vector("list", length = n_tot)
    # Loop for com barra de progresso manual
    for (x in 1:n_tot) {
      # Execute a função fit_eval_models
      best_models[[x]] <- fit_best_model(x = x, dfgrid, calibration_results,
                                         data_x, n_replicates,
                                         rep_data)

      # Sets the progress bar to the current state
      if(progress_bar){
        setTxtProgressBar(pb, x) }
    }  }

  #End cluster
  if(parallel){
  parallel::stopCluster(cl) }

  #Split list
  # Crie um vetor com o número de réplicas para cada modelo
  num_repl <- tapply(dfgrid$replicates, dfgrid$models, FUN = length)
  num_repl <- num_repl[match(dfgrid$models, names(num_repl))]

  # Split list
  best_models <- split(best_models, dfgrid$models)
  best_models <- lapply(best_models, function(sublist) {
    names(sublist) <- paste0("Rep_", seq_along(sublist))
    return(sublist)
  })

  #Rename models
  names(best_models) <- paste0("Model_", names(best_models))

  #Write models?
  if(write_models){
    saveRDS(best_models, file = paste0(file_name, ".RDS"))
  }

  return(best_models)
} #End of function


# #Test function
# #Load functions
# source("Functions/Metrics_Functions.R")
# source("Functions/eval_m.R")
# source("Functions/part_data.R")
#
# #Prediction to dataframe or raster
# #In df, each replicate in a column
#
# #Select 3 best models from candidate results
# cand_res <- read.csv("Models/Piper_fuligineum/candidate_results.csv")
# selected_models <- cand_res %>%
#   filter(is_concave == F, proc_pval.mean < 0.05,
#          Omission_rate_at_5.mean <= 0.05) %>% slice_min(AIC, n =3)
# occ_bg <- read.csv("Models/Piper_fuligineum/occ_bg.csv")
#
# #With kfold
# best_model <- fit_best(data = occ_bg,
#                        pr_bg = "pr_bg",
#                        selected_models = selected_models,
#                        var_categorical = NULL,
#                        replicates = TRUE,
#                        n_replicates = 4,
#                        train_portion = 0.7,
#                        rep_type = "kfold",
#                        parallel = TRUE,
#                        ncores = 5,
#                        progress_bar = TRUE,
#                        write_models = FALSE, #Write files?
#                        out_dir = NULL, #Name of the folder to write final models
#                        parallelType = "doSNOW",
#                        verbose = TRUE)
#
# #With subsample
# best_model2 <- fit_best(data = occ_bg,
#                        pr_bg = "pr_bg",
#                        selected_models = selected_models,
#                        var_categorical = NULL,
#                        replicates = TRUE,
#                        n_replicates = 5,
#                        train_portion = 0.7,
#                        rep_type = "subsample",
#                        parallel = TRUE,
#                        ncores = 5,
#                        progress_bar = TRUE,
#                        write_models = FALSE, #Write files?
#                        out_dir = NULL, #Name of the folder to write final models
#                        parallelType = "doSNOW",
#                        verbose = TRUE)
# #With bootstrap
# best_model3 <- fit_best(data = occ_bg,
#                         pr_bg = "pr_bg",
#                         selected_models = selected_models,
#                         var_categorical = NULL,
#                         replicates = TRUE,
#                         n_replicates = 10,
#                         train_portion = 0.7,
#                         rep_type = "bootstrap",
#                         parallel = TRUE,
#                         ncores = 5,
#                         progress_bar = TRUE,
#                         write_models = FALSE, #Write files?
#                         out_dir = NULL, #Name of the folder to write final models
#                         parallelType = "doSNOW",
#                         verbose = TRUE)
# #With kfold, no parallel and writing results
# best_model4 <- fit_best(data = occ_bg,
#                         pr_bg = "pr_bg",
#                         selected_models = selected_models,
#                         var_categorical = NULL,
#                         replicates = TRUE,
#                         n_replicates = 4,
#                         train_portion = 0.7,
#                         rep_type = "kfold",
#                         parallel = T,
#                         ncores = 5,
#                         progress_bar = TRUE,
#                         write_models = TRUE, #Write files?
#                         out_dir = "Models/Piper_fuligineum/Best_models/", #Name of the folder to write final models
#                         parallelType = "doSNOW",
#                         verbose = TRUE)
