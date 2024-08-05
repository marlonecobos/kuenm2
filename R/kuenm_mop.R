kuenm_mop <- function(data,
                      variables = NULL,
                      projection_file,
                      out_dir,
                      type = "basic", calculate_distance = FALSE,
                      where_distance = "in_range", distance = "euclidean",
                      scale = FALSE, center = FALSE, fix_NA = TRUE, percentage = 1,
                      comp_each = 2000, tol = NULL, rescale_distance = FALSE,
                      parallel = FALSE, n_cores = NULL, progress_bar = FALSE,
                      overwrite = FALSE){
  #Get calibration data
  m <-  data$calibration_data[,-1]

  #Check scenarios to predict
  sc <- setdiff(names(projection_file), c("Raster_pattern", "do_pca"))

  #Get raster pattern to read
  raster_pattern <- projection_file$Raster_pattern

  #Get dataframe with path to predictions
  #Present
  if("Present" %in% sc){
    #Create folder
    present_dir <- file.path(out_dir, "Present/")
    present_sc <- names(projection_file[["Present"]])
    suppressWarnings({
      d_present <- data.frame(Time = "Present",
                              Period = "Present",
                              Scenario = present_sc,
                              input_path = unlist(projection_file[["Present"]]),
                              output_path = normalizePath(file.path(present_dir,
                                                                    present_sc)))
      })
  }
  #Past
  if("Past" %in% sc){
    #Create folder
    past_dir <- file.path(out_dir, "Past/")
    #Get grid of projections
    df_past <- do.call(rbind, lapply(names(projection_file$Past), function(time) {
      time_data <- projection_file$Past[[time]]
      do.call(rbind, lapply(names(time_data), function(gcm) {
        data.frame(Time = "Past", Period = time, GCM = gcm, Path = time_data[[gcm]], stringsAsFactors = FALSE)
      }))
    }))

    #Looping in the grid

    #Create dataframe with path to results
    suppressWarnings({
      d_past <- data.frame(Time = "Past",
                           Period = df_past$Period,
                           GCM = df_past$GCM,
                           input_path = df_past$Path,
                           output_path = normalizePath(file.path(past_dir, df_past$Period,
                                                                 df_past$GCM),
                                                       mustWork = FALSE))
      })
  }
  #Future
  ####Project to Future scenarios####
  if("Future" %in% sc){
    #Create folder
    future_dir <- file.path(out_dir, "Future/")

    #Create grid of time-ssp-gcm
    df_future <- do.call(rbind, lapply(names(projection_file[["Future"]]), function(year_range) {
      year_range_data <- projection_file[["Future"]][[year_range]]
      do.call(rbind, lapply(names(year_range_data), function(ssp) {
        ssp_data <- year_range_data[[ssp]]
        do.call(rbind, lapply(names(ssp_data), function(gcm) {
          data.frame(Time = "Future", Period = year_range, ssp = ssp,
                     GCM = gcm, Path = ssp_data[[gcm]],
                     stringsAsFactors = FALSE)
        }))      }))    }))


    #Create dataframe with path to results
    suppressWarnings({
      d_future <- data.frame(Time = df_future$Time,
                             Period = df_future$Period,
                             ssp = df_future$ssp,
                             GCM = df_future$GCM,
                             input_path = df_future$Path,
                             output_path = normalizePath(file.path(future_dir,
                                                                   df_future$Period,
                                                                   df_future$ssp,
                                                                   df_future$GCM),
                                                         mustWork = FALSE))
      })
    }

  #Get dataframe with path to each projection
  if(!("Present" %in% sc)){
    d_present <- NULL
  }
  if(!("Past" %in% sc)){
    d_past <- NULL
  }
  if(!("Future" %in% sc)){
    d_future <- NULL
  }

  #Return and write files with path
  res_path <- kuenm2:::bind_rows_projection(list(d_present, d_past, d_future))
  #Create ID
  res_path$id <- 1:nrow(res_path)

  #Show progress bar?
  if (progress_bar) {
    pb <- txtProgressBar(0, nrow(res_path), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n) }

  #Looping in projection file
  for(i in 1:nrow(res_path)){
    #Read rasters
    r <- terra::rast(
      list.files(path = res_path$input_path[i], pattern = raster_pattern,
                 full.names = TRUE))
    #Do PCA?
    if(!is.null(data$pca)){
      r <- terra::predict(r, data$pca)
    }

    #MOP
    mop_i <- mop::mop(m = m, g = r, type, calculate_distance,
                      where_distance, distance, scale,
                      center, fix_NA, percentage, comp_each,
                      tol, rescale_distance, parallel,
                      n_cores, progress_bar = FALSE)
    #Save results
    dir.create(dirname(res_path$output_path[i]), recursive = TRUE,
               showWarnings = FALSE)

    if(!is.null(mop_i$mop_basic))
    terra::writeRaster(mop_i$mop_basic,
                paste0(res_path$output_path[i], "_mopbasic.tif"),
                overwrite = overwrite)

    if(!is.null(mop_i$mop_simple)){
      terra::writeRaster(mop_i$mop_simple,
                  paste0(res_path$output_path[i], "_mopsimple.tif"),
                  overwrite = overwrite)}

    if(!is.null(mop_i$mop_detailed[[1]])){
      #Set levels
      lapply(names(mop_i$mop_detailed)[-1], function(x){
        if(x %in% c("towards_low_combined", "towards_high_combined")){
        levels(mop_i$mop_detailed[[x]]) <- list(data.frame(
          id = mop_i$mop_detailed$interpretation_combined$values,
          category = mop_i$mop_detailed$interpretation_combined$extrapolation_variables))}
      #Write raster
       terra::writeRaster(mop_i$mop_detailed[[x]],
                    paste0(res_path$output_path[i], "_mop_", x, ".tif"),
                    overwrite = overwrite)
        write.csv(mop_i$mop_detailed$interpretation_combined,
                  paste0(res_path$output_path[i], "_interpretation_combined.csv"))
      })
    }
    if(progress_bar){
      setTxtProgressBar(pb, i) }
  }

  #Return output path
  res_path$output_path <- dirname(res_path$output_path)
  return(unique(res_path))
}

# #Test function internally
# data = readRDS("../test_kuenm2/Myrcia.RDS")
# projection_file <- readRDS("../test_kuenm2/Projection_file.RDS")
# out_dir = "../test_kuenm2/mop_results"
# type = "basic";
# calculate_distance = FALSE;
# where_distance = "in_range"; distance = "euclidean";
# scale = FALSE; center = FALSE; fix_NA = TRUE; percentage = 1;
# comp_each = 2000; tol = NULL; rescale_distance = FALSE;
# parallel = FALSE; n_cores = NULL; progress_bar = FALSE
# overwrite = T
#
# #Test function
# kmop <- kuenm_mop(data = data, projection_file = projection_file,
#                   out_dir = out_dir, type = "detailed", calculate_distance = FALSE,
#                   where_distance = "in_range", distance = "euclidean",
#                   scale = FALSE, center = FALSE, fix_NA = TRUE, percentage = 1,
#                   comp_each = 2000, tol = NULL, rescale_distance = FALSE,
#                   parallel = FALSE, n_cores = NULL, progress_bar = TRUE,
#                   overwrite = TRUE)
