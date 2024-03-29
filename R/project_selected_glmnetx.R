project_selected_glmnetx <- function(projection_file,
                                     models,
                                     out_dir = "Projection_results",
                                     consensus_per_model = TRUE,
                                     consensus_general = TRUE,
                                     consensus = c("median", "range", "mean", "stdev"), #weighted mean
                                     clamping = FALSE,
                                     var_to_clamp = NULL,
                                     type = "cloglog",
                                     overwrite = FALSE,
                                     progress_bar = TRUE,
                                     verbose = TRUE){
  #Extract variables from best models
  vars <- names(models[[1]][[1]]$samplemeans)[-1]

  if(!file.exists(out_dir)){
    dir.create(out_dir)
  }

  #Check scenarios to predict
  sc <- names(projection_file)

  #Project to present scenarios
  if("Present" %in% sc){
    if(verbose){
      message("Predicting models to Present scenarios...")
    }
    #Create folder
    present_dir <- file.path(out_dir, "Present/")
    present_sc <- names(projection_file[["Present"]])

    #Show progress bar?
    if(progress_bar) {
      pb <- txtProgressBar(min = 0, max = length(present_sc), style = 3)}

    for(i in 1:length(present_sc)) {
      p_i <- present_sc[i]
      present_sc_i <- projection_file[["Present"]][[p_i]]
      r <- rast(list.files(path = present_sc_i,
                           pattern = projection_file$Raster_pattern,
                           full.names = TRUE))
      #Create folder to save
      f_i <- file.path(present_dir, p_i)
      suppressWarnings(dir.create(f_i, recursive = TRUE))

      #Predict
      invisible(predict_selected_glmnetmx(models = models,
                     spat_var = r,
                     write_files = TRUE,
                     write_replicates = FALSE,
                     out_dir = f_i,
                     consensus_per_model = consensus_per_model,
                     consensus_general =   consensus_general,
                     consensus = consensus, #weighted mean
                     clamping = clamping,
                     var_to_clamp = var_to_clamp,
                     type = type,
                     overwrite = overwrite,
                     progress_bar = FALSE))
      #Set progress bar
      if(progress_bar){
        setTxtProgressBar(pb, i) }
    }
  } #End of present projections

  #Project to Past scenarios
  if("Past" %in% sc){
    if(verbose){
      message("\nPredicting models to Past scenarios...")
    }
    #Create folder
    past_dir <- file.path(out_dir, "Past/")
      #Get grid of projections
    df <- do.call(rbind, lapply(names(projection_file$Past), function(time) {
      time_data <- projection_file$Past[[time]]
      do.call(rbind, lapply(names(time_data), function(gcm) {
        data.frame(Time = time, GCM = gcm, Path = time_data[[gcm]], stringsAsFactors = FALSE)
      }))
    }))

    #Looping in the grid
    #Show progress bar?
    if(progress_bar) {
      pb <- txtProgressBar(0, nrow(df), style = 3)}

    for(i in 1:nrow(df)){
      time_i <- df$Time[i]
      gcm_i <- df$GCM[i]
      path_i <- df$Path[i]
      #Create folder
      f_i <- file.path(out_dir, "Past", time_i, gcm_i)
      suppressWarnings(dir.create(f_i, recursive = T))
      r <- rast(list.files(path_i, full.names = T,
                           pattern = projection_file$Raster_pattern))
      #Predict
      invisible(predict_selected_glmnetmx(models = models,
                                          spat_var = r,
                                          write_files = TRUE,
                                          write_replicates = FALSE,
                                          out_dir = f_i,
                                          consensus_per_model = consensus_per_model,
                                          consensus_general =   consensus_general,
                                          consensus = consensus, #weighted mean
                                          clamping = clamping,
                                          var_to_clamp = var_to_clamp,
                                          type = type,
                                          overwrite = overwrite,
                                          progress_bar = FALSE))
      #Set progress bar
      if(progress_bar){
        setTxtProgressBar(pb, i) }
    }
  } #End of past projections


  #Project to Future scenarios
  if("Future" %in% sc){
    if(verbose){
      message("\nPredicting models to Future scenarios...")
    }
    #Create grid of time-ssp-gcm
    df <- do.call(rbind, lapply(names(projection_file[["Future"]]), function(year_range) {
      year_range_data <- projection_file[["Future"]][[year_range]]
      do.call(rbind, lapply(names(year_range_data), function(ssp) {
        ssp_data <- year_range_data[[ssp]]
        do.call(rbind, lapply(names(ssp_data), function(gcm) {
          data.frame(Time = year_range, ssp = ssp, GCM = gcm, Path = ssp_data[[gcm]], stringsAsFactors = FALSE)
        }))      }))    }))

    #Looping in the grid
    #Show progress bar?
    if(progress_bar) {
      pb <- txtProgressBar(0, nrow(df), style = 3)}

    for(i in 1:nrow(df)){
      time_i <- df$Time[i]
      ssp_i <- df$ssp[i]
      gcm_i <- df$GCM[i]
      path_i <- df$Path[i]
      #Create folder
      f_i <- file.path(out_dir, "Future", time_i, ssp_i, gcm_i)
      suppressWarnings(dir.create(f_i, recursive = T))
      r <- rast(list.files(path_i, full.names = T,
                           pattern = projection_file$Raster_pattern))
      #Subset variables
      r <- r[[vars]]
      #Predict
      invisible(predict_selected_glmnetmx(models = models,
                                          spat_var = r,
                                          write_files = TRUE,
                                          write_replicates = FALSE,
                                          out_dir = f_i,
                                          consensus_per_model = consensus_per_model,
                                          consensus_general =   consensus_general,
                                          consensus = consensus, #weighted mean
                                          clamping = clamping,
                                          var_to_clamp = var_to_clamp,
                                          type = type,
                                          overwrite = overwrite,
                                          progress_bar = FALSE))
      #Set progress bar
      if(progress_bar){
        setTxtProgressBar(pb, i) }
    }
  }  #End of future projections
} #End of function

projection_file = readRDS("../test_kuenm2/Projection_file.RDS")
models <- readRDS("../test_kuenm2/Best_Models.RDS")
write_files = FALSE
write_replicates = FALSE
out_dir = "../test_kuenm2/Projection_results"
consensus_per_model = TRUE
consensus_general = TRUE
consensus = c("median", "range", "mean", "stdev") #weighted mean
clamping = FALSE
var_to_clamp = NULL
type = "cloglog"
overwrite = TRUE
verbose = TRUE
progress_bar = TRUE

# #Test function
# source("R/predict_selected_glmnetmx.R")
# source("R/helpers_glmnetmx.R")
# source("R/helpers_calibration_glmnetmx.R")
# pf = readRDS("../test_kuenm2/Projection_file.RDS")
# models <- readRDS("../test_kuenm2/Best_Models.RDS")
# out_dir = "../test_kuenm2/Projection_results"
#
# project_selected_glmnetx(projection_file = pf,
#                          models = models,
#                          out_dir = out_dir,
#                          consensus_per_model = TRUE,
#                          consensus_general = TRUE,
#                          consensus = c("median", "range", "mean", "stdev"), #weighted mean
#                          clamping = FALSE,
#                          var_to_clamp = NULL,
#                          type = "cloglog",
#                          overwrite = TRUE,
#                          progress_bar = TRUE,
#                          verbose = TRUE)
