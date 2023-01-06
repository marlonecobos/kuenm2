# source("Functions/Metrics_Functions.R")
#
# models <- readRDS("Models/Myrcia_hatschbachii/Best_models.RDS")
# spat_var <- rast("Models/Myrcia_hatschbachii/PCA_var.tiff")
# consensus_per_model = TRUE
# consensus_general = TRUE
# consensus = c("median", "range", "mean", "stdev") #mean, weighted mean, variance
# write_files = TRUE
# write_replicates = FALSE
# type = "cloglog"
# out_dir = "Models/Myrcia_hatschbachii/Predictions"
# overwrite = TRUE

#If predict to other scenarios, keep replicates?

predict_selected_glmnetmx <- function(models,
                           spat_var,
                           write_files = FALSE,
                           write_replicates = FALSE,
                           out_dir = NULL,
                           consensus_per_model = TRUE,
                           consensus_general = TRUE,
                           consensus = c("median", "range", "mean", "stdev"), #weighted mean
                           type = "cloglog",
                           overwrite = FALSE) {
  #Get models names
  nm <- names(models)

  #Get predictions for each replicate
  p_models <- lapply(models, function(i){
    terra::rast(lapply(i, function(x) {
      terra::predict(spat_var, x, na.rm = TRUE,
                                 type = type)
    }))
  })
  names(p_models) <- nm

  #Start to store results
  rep <- unlist(p_models)

  #Create empty list
  res <- list()
  #Get consensus by model
  if(consensus_per_model) {
  if("median" %in% consensus) {
    res$Consensus_per_model$median <- rast(lapply(p_models, median))
  }
  if("mean" %in% consensus) {
    res$Consensus_per_model$mean <- rast(lapply(p_models, mean))
  }
  if("stdev" %in% consensus) {
    res$Consensus_per_model$stdev <- rast(lapply(p_models, stdev))
  }
  if("range" %in% consensus) {
    res$Consensus_per_model$range <- rast(lapply(p_models, function(r) {
      diff(range(r))
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
      all_rep <- rast(p_models)

      if("median" %in% consensus) {
        gen_res$median <- median(all_rep)
      }
      if("mean" %in% consensus) {
        gen_res$mean <- mean(all_rep)
      }
      if("stdev" %in% consensus) {
        gen_res$stdev <- stdev(all_rep)
      }
      if("range" %in% consensus) {
        gen_res$range <- diff(range(all_rep))
      }
    }

    if(consensus_general  & length(p_models) > 1){
      all_rep <- rast(p_models)

      if("median" %in% consensus) {
        gen_res$median <- median(all_rep)
      }
      if("mean" %in% consensus) {
        gen_res$mean <- mean(all_rep)
      }
      if("stdev" %in% consensus) {
        gen_res$stdev <- stdev(all_rep)
      }
      if("range" %in% consensus) {
        gen_res$range <- diff(range(all_rep))
      }
    }
  }

  #Final list
  res <- lapply(1:length(nm), function(x) {
    mcs <- lapply(consensus, function(y) {
      res$Consensus_per_model[[y]][[x]]
    })
    mcs <- rast(mcs)
    names(mcs) <- consensus

    list(Replicates = rep[[x]], Model_consensus = mcs)
  })

  names(res) <- nm

  res <- c(res, General_consensus = rast(gen_res))

  #Write files?
  if(write_files) {
    if(!file.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }

  #Save files
  sapply(nm, function(i){
    #Write Replicates
    if(write_replicates) {
    writeRaster(res[[i]]$Replicates,
                       file.path(out_dir,
                                 paste0(i, "_replicates.tiff")),
                overwrite = overwrite) }
    #Write consensus by model
    writeRaster(res[[i]]$Model_consensus,
                file.path(out_dir,
                          paste0(i, "_consensus.tiff")),
                overwrite = overwrite)
    })
  #Write general consensus
   writeRaster(res$General_consensus,
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
