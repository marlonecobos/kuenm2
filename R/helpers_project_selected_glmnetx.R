# ####Predict_to_past####
# predict_to_past <- function(i, df, models, write_replicates,
#                             consensus_per_model, consensus_general,
#                             consensus, clamping, var_to_clamp, type,
#                             overwrite, out_dir, projection_file) {
#   time_i <- df$Period[i]
#   gcm_i <- df$GCM[i]
#   path_i <- df$Path[i]
#   #Create folder
#   f_i <- normalizePath(file.path(out_dir, "Past", time_i, gcm_i))
#   suppressWarnings(dir.create(f_i, recursive = T))
#   r <- terra::rast(list.files(path_i, full.names = T,
#                               pattern = projection_file$Raster_pattern))
#   #Predict
#   invisible(predict_selected_glmnetmx(models = models,
#                                       spat_var = r,
#                                       write_files = TRUE,
#                                       write_replicates = write_replicates,
#                                       out_dir = f_i,
#                                       consensus_per_model = consensus_per_model,
#                                       consensus_general =   consensus_general,
#                                       consensus = consensus, #weighted mean
#                                       clamping = clamping,
#                                       var_to_clamp = var_to_clamp,
#                                       type = type,
#                                       overwrite = overwrite,
#                                       progress_bar = FALSE))
#   return(NULL)
# }
#
# ####Predict_to_future####
# predict_to_future <- function(i, df, dir_future) {
#   time_i <- df$Period[i]
#   ssp_i <- df$ssp[i]
#   gcm_i <- df$GCM[i]
#   path_i <- df$Path[i]
#   print("Até aqui funciona!!")
#   #Create folder
#   f_i <- normalizePath(file.path(dir_future, "Future", time_i, ssp_i, gcm_i))
#   print("Até aqui funciona 2222!!")
#   suppressWarnings(dir.create(f_i, recursive = T))
#   r <- terra::rast(list.files(path_i, full.names = T,
#                        pattern = projection_file$Raster_pattern))
#   #Predict
#   invisible(predict_selected_glmnetmx(models = models,
#                                       spat_var = r,
#                                       write_files = TRUE,
#                                       write_replicates = write_replicates,
#                                       out_dir = f_i,
#                                       consensus_per_model = consensus_per_model,
#                                       consensus_general =   consensus_general,
#                                       consensus = consensus, #weighted mean
#                                       clamping = clamping,
#                                       var_to_clamp = var_to_clamp,
#                                       type = type,
#                                       overwrite = overwrite,
#                                       progress_bar = FALSE))
#   return(NULL)
# }

#Predict to several scenarion
multiple_projections <- function(i, res_path, raster_pattern, par_list){
  res_i <- res_path[i,]
  input_i <- res_i$input_path
  output_i <- res_i$output_path
  #Import rasters
  r_i <- terra::rast(list.files(input_i, full.names = T,
                                pattern = raster_pattern))
  ####PCA ?##########
  #If pca is not NULL, predict pca
  if(!is.null(par_list$models$pca)){
    vars_in_pca <- names(par_list$models$pca$center)
    vars_out_pca <- setdiff(names(r_i), vars_in_pca)
    r_i_pca <- terra::predict(r_i[[vars_in_pca]], par_list$models$pca)
    if(length(vars_out_pca) == 0) {
      r_i <- r_i_pca} else {
        r_i <- c(r_i_pca, r_i[[ vars_out_pca]])
      }
  }
  ##################

  #Predict
  invisible(predict_selected_glmnetmx(models = par_list$models,
                                      spat_var = r_i,
                                      write_replicates = par_list$write_replicates,
                                      out_dir = output_i,
                                      consensus_per_model = par_list$consensus_per_model,
                                      consensus_general = par_list$consensus_general,
                                      consensus = par_list$consensus, #weighted mean
                                      clamping = par_list$clamping,
                                      var_to_clamp = par_list$var_to_clamp,
                                      type = par_list$type,
                                      overwrite = par_list$overwrite,
                                      progress_bar = FALSE,
                                      write_files = TRUE))
  return(invisible(NULL))
}

#Binarize model
binarize_values <- function(x, v, new_value = 1) {
  x <- (x >= v) * new_value
  return(x)
}

#Calculate threshold
calc_thr <- function(occ_suitability, thr = 0.1) {
  sort(occ_suitability)[round(length(occ_suitability) * thr) + 1]
}
# data <- readRDS("../test_kuenm2/Myrcia.RDS")
# r <- rast("../test_kuenm2/Projection_results/Present/Brazil/General_consensus.tiff")
# thr <- calc_thr(r = r, data = data, consensus = "median", thr = 0.1)
# thr

