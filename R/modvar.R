#' Title
#'
#' @param projection_file
#' @param by_replicate
#' @param by_gcm
#' @param by_model
#' @param consensus
#' @param write_files
#' @param output_dir
#' @param return_rasters
#' @param progress_bar
#' @param verbose
#' @param overwrite
#'
#' @return
#' @export
#'
#' @examples
#'
modvar <- function(projection_file, #output of project_selected_glmnetx
                   by_replicate = TRUE,
                   by_gcm = TRUE,
                   by_model = TRUE,
                   consensus = "mean", #To calculate variance by_gcm and by_model
                   write_files = TRUE,
                   output_dir = NULL,
                   return_rasters = FALSE,
                   progress_bar = FALSE,
                   verbose = TRUE,
                   overwrite = FALSE){
  if(write_files){
    out_dir <- file.path(output_dir, "variance")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  #### Get data ####
  d <- projection_file[["paths"]]

  # Remove present
  #d <- d[d$Time != "Present",]

  # Get unique combinations of Time, Period, ssp and Scenario
  uc <- unique(d[, c("Time", "Period", "ssp", "Scenario")])
  #dplyr::tibble(uc)

  #Show progress bar?
  if (progress_bar) {
    pb <- txtProgressBar(0, nrow(uc), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n) }

  ####Iteration over combinations####
  res <- lapply(1:nrow(uc), function(p){
    #To test
    #p = 8
    ssp <- uc$ssp[p]
    period <- uc$Period[p]
    time <- uc$Time[p]
    scenario <- uc$Scenario[p]

    # Filter dataframe to current combination
    d_p <- d[
      (is.na(d$ssp) | d$ssp == ssp) &
        (is.na(d$Period) | d$Period == period) &
        (is.na(d$Scenario) | d$Scenario == scenario) &
        (is.na(d$Time) | d$Time == time), ]

    # Get all paths
    paths <- d_p$output_path

    #### By replicate ####
    if(by_replicate){
        if(verbose) {
        message("\nCalculating prediction variance coming from distinct replicates: scenario ", p, " of ", nrow(uc))
      }
      #Variance of means

      #### By replicates ####
      # Get variance of replicates in each gcm, than get the average across gcms
      var_rep_by_gcm <- terra::rast(lapply(paths, var_models_rep_by_gcm))
      var_rep <- terra::mean(var_rep_by_gcm)
      } else {#End of by_replicate
        var_rep <- NULL}

    #### By Model ####
    if(by_model){
      if(verbose) {
        message("Calculating prediction variance coming from distinct models: scenario ", p, " of ", nrow(uc))}
      # Get variance of models in each gcm, than get the average
      var_model_by_gcm <- terra::rast(lapply(paths, var_models_model_by_gcm, consensus))
      var_model <- terra::mean(var_model_by_gcm)
      names(var_model) <- "by_model"
    } else { #End of by model
      var_model <-  NULL}


    ####By GCM####
    if(by_gcm & period != "Present"){
      if(verbose) {
        message("Calculating prediction variance coming from distinct GCMs: scenario ", p, " of ", nrow(uc))}
      var_gcm <- var_models_across_gcm(paths = paths, consensus = consensus)
      names(var_gcm) <- "by_gcm"
    } else {
      var_gcm <- NULL}#End of by_gcm

    all_var <- terra::rast(c("by_rep" = var_rep,
                 "by_model" = var_model,
                 "by_gcm" = var_gcm))


    #Write results
    if(write_files){
      #Name of raster
      nr <- gsub("_NA", "",
                 paste(time, period, scenario, ssp, sep = "_"))
      terra::writeRaster(all_var,
                         filename = file.path(out_dir,
                                              paste0(nr, ".tiff")),
                         overwrite = overwrite)
    }

    #Progress bar
    if(progress_bar){
      setTxtProgressBar(pb, p) }

    if(return_rasters){
      return(all_var)} else {
        return(invisible(NULL))}
    }) #End of res

  if(return_rasters){
    names(res) <- gsub("Present_Present", "Present",
                       gsub("_NA", "",
                       paste(uc$Time, uc$Period, uc$ssp, uc$Scenario, sep = "_")))
    return(res) } else {
    return(invisible(NULL))
  }
} #End of function


# #Objects to test function internally
# library(terra)
# library(pbapply)
# by_replicate = TRUE
# by_gcm = TRUE
# by_model = TRUE
# projection_file = readRDS("../test_kuenm2/Projection_results/Projection_paths.RDS")
# consensus = "mean"
# write_files = TRUE
# output_dir = "../test_kuenm2/Projection_results/"
# return_rasters = FALSE
#
# #Test function
# v <- modvar(projection_file =projection_file, #output of project_selected_glmnetx
#             by_replicate = TRUE,
#             by_gcm = FALSE,
#             by_model = TRUE,
#             consensus = "mean", #To calculate variance by_gcm and by_model
#             write_files = FALSE,
#             output_dir = "../test_kuenm2/Projection_results/",
#             return_rasters = TRUE,
#             progress_bar = TRUE,
#             verbose = TRUE,
#             overwrite = TRUE)
#
