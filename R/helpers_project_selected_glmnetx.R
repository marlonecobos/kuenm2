#Predict to several scenarion
multiple_projections <- function(i, res_path, raster_pattern, par_list){
  res_i <- res_path[i,]
  input_i <- res_i$input_path
  output_i <- res_i$output_path
  #Import rasters
  r_i <- terra::rast(list.files(input_i, full.names = T,
                                pattern = raster_pattern))

  #Predict
  invisible(predict_selected(models = par_list$models,
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


#Helper function to calculate variance coming from replicates by gcm
var_models_rep_by_gcm <- function(path){
  model_files <- list.files(path = path, pattern = "replicates",
                            full.names = TRUE)
  if(length(model_files) == 0){
    stop("Replicates not found.
Set by_replicate = FALSE or rerun project_selected_glmnetx() with write_replicates = TRUE")
  }

  if(length(model_files) > 1) {
    r_x <- lapply(model_files, terra::rast)
    #Take mean of replicates (1-1, 2-2, 3-3, etc...)
    n_replicates <- terra::nlyr(r_x[[1]])
    mean_replicates <- terra::rast(lapply(1:n_replicates, function(n){
      rep_n <- terra::mean(rast(lapply(r_x, function(x) x[[n]])))
    }))
    var_rep_x <- terra::app(mean_replicates, "var")} else {
      r_x <- terra::rast(model_files)
      var_rep_x <- terra::app(r_x, "var")}
  names(var_rep_x) <- basename(path)
  return(var_rep_x)
}

var_models_model_by_gcm <- function(path, consensus){
  r_x <- terra::rast(list.files(path = path, pattern = "Model_.*consensus",
                         full.names = TRUE))
  r_x <- r_x[[sapply(r_x, function(r) names(r) == consensus)]]
  var_x <- terra::app(r_x, "var")
  return(var_x)
}

var_models_across_gcm <- function(paths, consensus){
  # Read rasters
  r <- terra::rast(lapply(paths, function(x){
    terra::rast(file.path(x, "General_consensus.tiff"))[[consensus]]
  }))

  #plot(r)
  #Calculate variance
  v <- terra::app(r, "var")
  return(v)
}

var_models_across_ssp <- function(paths, consensus){
  # Read rasters
  r <- terra::rast(lapply(paths, function(x){
    terra::rast(file.path(x, "General_consensus.tiff"))[[consensus]]
  }))

  #plot(r)
  #Calculate variance
  v <- terra::app(r, "var")
  return(v)
}
