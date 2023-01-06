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
                                  n_replicates = 1,
                                  rep_type = "kfold",
                                  train_portion = 0.7,
                                  write_models = FALSE, #Write files?
                                  file_name = NULL, #Name of the folder to write final models
                                  parallel = TRUE,
                                  ncores = 1,
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
                "categoricalval")

  #Extracts IDs from models
  m_ids <- calibration_results$selected_models$ID

  #Create grid of fitted models
  dfgrid <- expand.grid(models = m_ids, replicates = 1:n_replicates)
  n <-  nrow(dfgrid)

  #Prepare data (index) to replicates
  if(n_replicates > 1) {
    #Partitioning data
    rep_data <- part_data(data = calibration_results$calibration_data,
                          pr_bg = "pr_bg",
                          train_portion = train_portion,
                          n_replicates = n_replicates,
                          method = rep_type)
    }


  #Set parallelization
  if(parallel) {
    #Make cluster
    cl <- parallel::makeCluster(ncores)
    #Show progress bar?
    if (isTRUE(progress_bar)) {
      pb <- txtProgressBar(0, n, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n) }

    if (parallelType == "doParallel") {
      doParallel::registerDoParallel(cl)
      opts <- NULL
    }

    if (parallelType == "doSNOW") {
      doSNOW::registerDoSNOW(cl)
      if (isTRUE(progress_bar))
        opts <- list(progress = progress)
      else opts <- NULL
    }
  } else {
    opts <- NULL}

  #Fit models with replicates
  best_models <- foreach::foreach(x = 1:n, .options.snow = opts,
                                  .export = to_export) %dopar% {
                                    #Get grid
                                    grid_x <- dfgrid[x,]
                                    m_id <- grid_x$models
                                    rep_x <- grid_x$replicates

                                    #Get best model
                                    best_models_i <- calibration_results$selected_models[which(calibration_results$selected_models$ID == m_id),]
                                    #best_models_i <- selected_models[i,]
                                    best_formula <- best_models_i$Formulas
                                    best_regm <- best_models_i$regm

                                    #Get replicate, if necessary
                                    if(n_replicates > 1){
                                      rep_i <- rep_data[[rep_x]]
                                      data_x <- calibration_results$calibration_data[rep_i, ]
                                      } else { #Select i k-fold
                                        data_x <- data }
                                    #Run model
                                    mod_x <- glmnet_mx(p = data_x[,"pr_bg"], data = data_x,
                                                       f = as.formula(best_formula),
                                                       regmult = best_regm,
                                                       calculate_AIC = FALSE)

                                    #Only to make sure the IDs in list are correct
                                    mod_x$checkModel <- m_id
                                    mod_x$checkreplicate <- rep_x
                                    # #Predict
                                    # pred_x <- terra::predict(spat_var, mod_x, type = "cloglog",
                                    #                          na.rm = TRUE)
                                    return(mod_x) } #End of foreach

  #End cluster
  parallel::stopCluster(cl)

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
