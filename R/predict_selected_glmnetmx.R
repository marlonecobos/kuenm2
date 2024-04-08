#' Predict selected Maxent-like glmnet models
#'
#' @importFrom terra rast clamp predict median mean stdev diff writeRaster
#' @export

predict_selected_glmnetmx <- function(models,
                           spat_var,
                           write_files = FALSE,
                           write_replicates = FALSE,
                           out_dir = NULL,
                           consensus_per_model = TRUE,
                           consensus_general = TRUE,
                           consensus = c("median", "range", "mean", "stdev"), #weighted mean
                           clamping = FALSE,
                           var_to_clamp = NULL,
                           type = "cloglog",
                           overwrite = FALSE,
                           progress_bar = TRUE) {
  #Get models names
  models <- models[["Models"]]
  nm <- names(models)

  #Clamp
  if (clamping) {
    #Get minimum e maximum values from the models
    varmin <- models[[1]][[1]]$varmin[-1] #Get var min
    varmax <- models[[1]][[1]]$varmax[-1] #Get var min
    #Variables to clamp (all variables if var_to_clamp = NULL)
    if(is.null(var_to_clamp)) {
      var_to_clamp <- names(varmin)
    }

    #Clamp variables
    clamped_variables <- terra::rast(lapply(var_to_clamp, function(i){
      terra::clamp(spat_var[[i]], lower = varmin[i], upper = varmax[i],
                          values = TRUE)
    }))
    #Replace clamped variables in spat_var object
    spat_var[[names(clamped_variables)]] <- clamped_variables

  }

  #Get predictions for each replicate
  n_models <- length(models) #Number of models

  #Show progress bar?
  if(progress_bar) {
    pb <- txtProgressBar(0, n_models, style = 3)}

  #Create empty list
  p_models <- list()

  #Fill list with predictions
  for (i in seq_along(models)) {
    inner_list <- list()

    for (x in models[[i]]) {
      prediction <- terra::predict(spat_var, x, na.rm = TRUE, type = type)
      inner_list[[length(inner_list) + 1]] <- prediction
    }

    p_models[[length(p_models) + 1]] <- terra::rast(inner_list)
    #Set progress bar
    if(progress_bar){
      setTxtProgressBar(pb, i) }
  }

  #Rename models and replicates
  names(p_models) <- nm
  for (i in nm) {
    names(p_models[[i]]) <- names(models[[i]])
  }

  #Start to store results
  rep <- unlist(p_models)

  #Create empty list
  res <- list()
  #Get consensus by model
  if(consensus_per_model) {
  if("median" %in% consensus) {
    res$Consensus_per_model$median <- terra::rast(lapply(p_models, terra::median))
  }
  if("mean" %in% consensus) {
    res$Consensus_per_model$mean <- terra::rast(lapply(p_models, terra::mean))
  }
  if("stdev" %in% consensus) {
    res$Consensus_per_model$stdev <- terra::rast(lapply(p_models, terra::stdev))
  }
  if("range" %in% consensus) {
    res$Consensus_per_model$range <- terra::rast(lapply(p_models, function(r) {
      terra::diff(range(r))
    }))
    }
      }

  #Create empty list to general consensus
  gen_res <- list()
  #Get general consensus
  if(consensus_general & length(p_models) == 1 & consensus_per_model){

    if("median" %in% consensus) {
      gen_res$median <- res$Consensus_per_model$median
    }
    if("mean" %in% consensus) {
      gen_res$mean <- res$Consensus_per_model$mean
    }
    if("stdev" %in% consensus) {
      gen_res$stdev <- res$Consensus_per_model$stdev
    }
    if("range" %in% consensus) {
      gen_res$range <- res$Consensus_per_model$range
    }
  } else {
    if(consensus_general & length(p_models) == 1 & !consensus_per_model) {
      all_rep <- terra::rast(p_models)

      if("median" %in% consensus) {
        gen_res$median <- terra::median(all_rep)
      }
      if("mean" %in% consensus) {
        gen_res$mean <- terra::mean(all_rep)
      }
      if("stdev" %in% consensus) {
        gen_res$stdev <- terra::stdev(all_rep)
      }
      if("range" %in% consensus) {
        gen_res$range <- terra::diff(terra::range(all_rep))
      }
    }

    if(consensus_general  & length(p_models) > 1){
      all_rep <- terra::rast(p_models)

      if("median" %in% consensus) {
        gen_res$median <- terra::median(all_rep)
      }
      if("mean" %in% consensus) {
        gen_res$mean <- terra::mean(all_rep)
      }
      if("stdev" %in% consensus) {
        gen_res$stdev <- terra::stdev(all_rep)
      }
      if("range" %in% consensus) {
        gen_res$range <- terra::diff(terra::range(all_rep))
      }
    }
  }

  #Final list
  res <- lapply(1:length(nm), function(x) {
    mcs <- lapply(consensus, function(y) {
      res$Consensus_per_model[[y]][[x]]
    })
    mcs <- terra::rast(mcs)
    names(mcs) <- consensus

    list(Replicates = rep[[x]], Model_consensus = mcs)
  })

  names(res) <- nm

  res <- c(res, General_consensus = terra::rast(gen_res))

  #Write files?
  if(write_files) {
    if(!file.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }

  #Save files
  sapply(nm, function(i){
    #Write Replicates
    if(write_replicates) {
      terra::writeRaster(res[[i]]$Replicates,
                       file.path(out_dir,
                                 paste0(i, "_replicates.tiff")),
                overwrite = overwrite) }
    #Write consensus by model
    terra::writeRaster(res[[i]]$Model_consensus,
                file.path(out_dir,
                          paste0(i, "_consensus.tiff")),
                overwrite = overwrite)
    })
  #Write general consensus
  terra::writeRaster(res$General_consensus,
              file.path(out_dir, "General_consensus.tiff"),
              overwrite = overwrite)

    # sapply(r_i, function(x){
    #   writeRaster(x,
    #               paste0(out_dir, "/", i, "/", x, ".tiff"),
    #               overwrite = overwrite)
    #   })


  }
  return(res)
  }#End of function

# #Predict simple model
# source("Functions/Metrics_Functions.R")
# library(terra)
#
# models <- readRDS("Models/Piper_fuligineum/Best_models/Best_models.RDS")
# spat_var <- rast("Models/Piper_fuligineum/PCA_variables.tiff")
#
# #Predict with consensus per model and consensus general, without write files
# pm1 <- predict_models(models = models, spat_var = spat_var,
#                      write_files = FALSE,
#                      write_replicates = FALSE,
#                      out_dir = NULL,
#                      consensus_per_model = TRUE,
#                      consensus_general = TRUE,
#                      consensus = c("median", "range", "mean", "stdev"), #weighted mean
#                      type = "cloglog",
#                      overwrite = TRUE)
#
# #Predict with only consensus general, without write files
# pm2 <- predict_models(models = models, spat_var = spat_var,
#                       write_files = FALSE,
#                       write_replicates = FALSE,
#                       out_dir = NULL,
#                       consensus_per_model = FALSE,
#                       consensus_general = TRUE,
#                       consensus = c("median", "range", "mean", "stdev"), #weighted mean
#                       type = "cloglog",
#                       overwrite = TRUE)
#
# #Predict with only consensus per model, without write files
# pm3 <- predict_models(models = models, spat_var = spat_var,
#                       write_files = FALSE,
#                       write_replicates = FALSE,
#                       out_dir = NULL,
#                       consensus_per_model = TRUE,
#                       consensus_general = FALSE,
#                       consensus = c("median", "range", "mean", "stdev"), #weighted mean
#                       type = "cloglog",
#                       overwrite = TRUE)
# #Predict with consensus per model and consensus general, writing files, but not replicates
# pm4 <- predict_models(models = models, spat_var = spat_var,
#                       write_files = TRUE,
#                       write_replicates = FALSE,
#                       out_dir = "Models/Piper_fuligineum/Predictions",
#                       consensus_per_model = TRUE,
#                       consensus_general = TRUE,
#                       consensus = c("median", "range", "mean", "stdev"), #weighted mean
#                       type = "cloglog",
#                       overwrite = TRUE)
# #Predict with consensus per model and consensus general, writing consensus and replicates
# pm5 <- predict_models(models = models, spat_var = spat_var,
#                       write_files = TRUE,
#                       write_replicates = TRUE,
#                       out_dir = "Models/Piper_fuligineum/Predictions",
#                       consensus_per_model = TRUE,
#                       consensus_general = TRUE,
#                       consensus = c("median", "range", "mean", "stdev"), #weighted mean
#                       type = "cloglog",
#                       overwrite = TRUE)
