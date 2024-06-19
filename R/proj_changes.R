#Compute changes between scenarios
#' Title
#'
#' @param projection_paths
#' @param referece_id
#' @param consensus
#' @param include_id
#' @param user_thr
#' @param by_gcm
#' @param by_change
#' @param general_summary
#' @param force_resample
#' @param write_results
#' @param output_dir
#' @param overwrite
#' @param write_bin_models
#' @param return_rasters
#'
#' @return
#' @export
#'
#' @examples
proj_changes <- function(projection_paths, #Output from project_selected_glmnetx function
                         referece_id = 1, #In projection_paths, what is the id of the reference? (present time)
                         consensus = "median", #Consensus to use (median, mean, sd, etc)
                         include_id = NULL, #Which projections (by id) include (or exclude if negative)? Default "All"
                         user_thr = NULL, #User can specify his own threshold
                         by_gcm = TRUE, #Compute results by gcm?
                         by_change = TRUE, #Compute results by change (gain, loss and stability? )
                         general_summary = TRUE, #Get general summary?
                         force_resample = TRUE, #Force resample of projections to have same extension and resolution of reference raster (present)?
                         write_results = TRUE, #Write results?
                         output_dir = NULL,  #Output directory
                         overwrite = TRUE, #Overwrite rasters?
                         write_bin_models = FALSE, #Write binarized models?
                         return_rasters = FALSE){

  if(!return_rasters & !write_results) {
    warning("return_rasters and write_results cannot both be set to FALSE.
            Changing return_rasters to TRUE")
    return_rasters <- TRUE
  }

  #Create directory to save results, if necessary
  if(write_results){
    if(is.null(output_dir)){
      stop("If write_results = TRUE, you must define output_dir")
    }
    out_dir <- file.path(out_dir, "Projection_changes/")
    if(!file.exists(out_dir)){
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  }
  }

  #Extract threshold
  if(is.null(user_thr)) {
    thr <- projection_paths[["threshold"]][["consensus"]][[consensus]]
  } else {
    thr <- user_thr
  }

  #Remove threshold from projection paths
  projection_paths <- projection_paths[["paths"]]


  #Get only specified scenarios
  if(!is.null(include_id)) {
    to_include <- include_id[include_id > 0]
    if(length(to_include) > 0) {
      projection_paths <- projection_paths[projection_paths$id %in% to_include,]
    }
    to_exclude <- abs(include_id[include_id < 0])
    if(length(to_exclude) > 0) {
      projection_paths <- projection_paths[!(projection_paths$id %in%
                                               to_exclude),]
    }
  }


  #Get raster of reference
  dir_reference <- projection_paths$output_path[which(projection_paths$id ==
                                                        referece_id)]
  file_reference <- file.path(dir_reference, "General_consensus.tiff")
  r <- terra::rast(file_reference)[[consensus]]
  names(r) <- projection_paths$Time[which(projection_paths$id ==
                                                   referece_id)]


  ####Binarize models ####
  r_bin <- binarize_values(x = r, v = thr, new_value = 2)
  #Set levels
  levels(r_bin) <- data.frame(id = c(0, 2),
                          Result = c("Unsuitable",
                                     "Suitable"))
  names(r_bin) = "Present"

  #plot(r_bin)
  #Looping throughout projections to binarize
  pp <- projection_paths[projection_paths$Time %in% c("Past", "Future"),]

  #Binarize
  proj_bin <- rast(lapply(1:nrow(pp), function(i){
    #Get scenario projection
    d_i <- pp[i,]
    dir_change <- d_i$output_path
    #Get time
    proj_time <- d_i$Time
    file_change <- file.path(dir_change, "General_consensus.tiff")
    r_change <- rast(file_change)[[consensus]]

    #Force resample if necessary
    if(force_resample){
      if(any(res(r) != (res(r_change)) | ext(r) != ext(r_change))){
        r_change <- terra::resample(r_change, r_bin, method = "average")
      }}

    r_change <- binarize_values(x = r_change, v = thr, new_value = 1)
    levels(r_change) <- data.frame(id = c(0, 1),
                                   Result = c("Unsuitable",
                                              "Suitable"))
    return(r_change)
  }))
  #Rename binarizes scenarions
  names(proj_bin) <- paste(pp$Time, pp$Period, pp$ssp, pp$GCM, sep = "_")
  names(proj_bin) <- gsub("_NA_", "_", names(proj_bin))
  #plot(proj_bin[[1:6]])
  if(write_bin_models) {
    terra::writeRaster(x = c(r_bin, proj_bin),
                      filename = file.path(out_dir, "Binarized.tiff"),
                      overwrite = overwrite)
  }


  #Get single scenarios by Time and period
  sc <- unique(pp[,c("Time", "Period", "ssp")])

  ####Identify changes by scenario####
  if(by_gcm | by_change){
  # #Change values of r to compute gain and loss
  # r2 <- binarize_values(r, v = 1, new_value = 2)
  # r2[r2 == 0] <- 1
  # #plot(r2)

  #Table to set levels in raster
  cls <- data.frame(id = c(1, 2, 3, 0),
                    Result = c("Gain",
                               "Loss",
                               "Suitable-stable",
                               "Unsuitable-stable"))

  res_by_gcm <- lapply(proj_bin, function(i){
    #Compute changes
    r_result <- r_bin + i
    #Get legend
    levels(r_result) <- cls
    #plot(r_result)
    return(r_result)
  })

  #Rename res
  names(res_by_gcm) <- names(proj_bin)

  #Rasterize res
  res_by_gcm <- terra::rast(res_by_gcm)
  #plot(res_by_gcm[[1:8]])

  if(write_results & by_gcm){
    terra::writeRaster(x = res_by_gcm,
                       filename = file.path(out_dir, "Changes_by_GCM.tiff"),
                       overwrite = overwrite)
  }

  } #End of by_gcm


  ####Results by change####
  if(by_change) {
  res_by_change <- lapply(1:nrow(sc), function(i){
    sc_i <- sc[i,]
    scenario_i <- paste(sc_i$Time, sc_i$Period, sc_i$ssp, sep = "_")
    scenario_i <- gsub("_NA", "_", scenario_i)
    #Subset results
    res_i <- res_by_gcm[[grep(scenario_i, names(res_by_gcm))]]
    #Looping throught changes
    #Get available changes
    out <- unique(unlist(sapply(res_i, terra::unique, simplify = T)))
    sc_out <- lapply(out, function(x){
      #Get change value id
      out_id <- cls$id[cls$Result == x]
      res_i_bin <- (res_i == out_id) * 1

      #Sum results
      res_i_sum <- sum(res_i_bin) ####Import sum from terra
      #Levels
      l_sum <- data.frame(id = terra::unique(res_i_sum)$sum,
                          Result = paste0(x, " in ", terra::unique(res_i_sum)$sum, " GCMs"))
      levels(res_i_sum) <- l_sum
      return(res_i_sum)
    })
    names(sc_out) <- out
    return(rast(sc_out))
  })
  #Set names by scenarion
  names(res_by_change) <- gsub("_NA", "",
                             paste(sc$Time, sc$Period, sc$ssp, sep = "_"))
  #plot(res_by_change[[1]], main = names(res_by_change[1]))

  #Save results
  if(write_results){
    dir.create(file.path(out_dir, "Results_by_change"))
    sapply(names(res_by_change), function(z){
      terra::writeRaster(x = res_by_change[[z]],
                         filename = file.path(out_dir, "Results_by_change",
                                              paste0(z, ".tif")),
                         overwrite = overwrite)
      })}

  } else {res_by_change <- NULL} #End of by change

  #plot(res_by_change$Past_LGM)
  #plot(res_by_change$`Future_2041-2060_ssp245`)

  ####General summary####
  if(general_summary){
  res_summary <- lapply(1:nrow(sc), function(i){
    sc_i <- sc[i,]
    scenario_i <- paste(sc_i$Time, sc_i$Period, sc_i$ssp, sep = "_")
    scenario_i <- gsub("_NA", "_", scenario_i)
    #Subset results
    res_i <- proj_bin[[grep(scenario_i, names(proj_bin))]]
    #Presence value in present (number of gcms + 1)
    n_gcms <- terra::nlyr(res_i) + 1
    #Change values of present
    r_present <- (r_bin == 2) * n_gcms

    #Sum rasters to get results
    res_sum <- sum(c(r_present, res_i))

    # # preparing description table
    # vals <- sort(na.omit(unique(res_sum[])))
    # #Fix vals
    # if(length(vals) != max(vals)){
    #   vals <- seq(0, max(vals), 1)
    #   vals <- vals[vals != 2]
    # }

    vals <- seq(0, (n_gcms*2 - 1), 1)

    loss <- ceiling(max(vals)/2)
    l_val <- c(loss, vals[vals > loss & vals != max(vals)])
    g_val <- vals[vals < loss & vals != 0]
    gains <- paste0("gain in ", g_val, " GCMs")
    gains[gains == paste0("gain in ", n_gcms - 1, " GCMs")] <- "gain in all GCMs"
    losses <- paste0("loss in ", max(vals) - l_val, " GCMs")
    losses[losses == paste0("loss in ", n_gcms - 1, " GCMs")] <- "loss in all GCMs"
    descriptions <- c("stable, unsuitable in current period and all GCMs", gains,
                      losses, "stable, suitable in current period and all GCMs")
    res_table <- data.frame(Raster_value = vals, Description = descriptions)


    #Set levels
    levels(res_sum) <- res_table
    #plot(res_sum)
    return(res_sum)
  })
  #Set names by scenario
  names(res_summary) <- gsub("_NA", "",
                             paste(sc$Time, sc$Period, sc$ssp, sep = "_"))
  res_summary <- rast(res_summary) #Rasterize
  #plot(res_summary[[1:4]])

  if(write_results){
    terra::writeRaster(x = res_summary,
                       filename = file.path(out_dir, "Changes_summary.tiff"),
                       overwrite = overwrite)
  }



  } else {res_summary <- NULL} #End of general summary


  if(return_rasters){
  #Create final object
  res_final <- list(Binarized = c(r_bin, proj_bin),
                    Results_by_gcm = res_by_gcm,
                    Results_by_change = res_by_change,
                    Summary_changes = res_summary)
  return(res_final) } else {
    return(invisible(NULL))
  }

} #End of function


# #To test function internally
# projection_paths <- readRDS("../test_kuenm2/Projection_results/Projection_paths.RDS")
# referece_id = 2
# out_dir = "../test_kuenm2/Projection_results/"
# consensus = "median"
# force_resample = T
# output_dir <- "../test_kuenm2/Projection_results/"
# write_results = TRUE
# i = 12
# include_id = -1
# write_bin_models = FALSE
# thr = 0.1 #Threshold value (5%, 10%, etc...)
# data_occ = readRDS("../test_kuenm2/Myrcia.RDS") #Data created with prepare_data
# user_thr = NULL
# by_gcm = TRUE
# by_change = TRUE
# general_summary = TRUE
# overwrite = TRUE
# return_rasters = TRUE
