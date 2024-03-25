prepare_proj <- function(models, present_dir = NULL,
                         past_dir = NULL, past_time = NULL, past_gcm = NULL,
                         future_dir = NULL, future_time = NULL,
                         future_ssp = NULL, future_gcm = NULL,
                         filename,
                         raster_pattern = ".tif*") {

  #Extract variables from best models
  vars <- names(models[[1]][[1]]$samplemeans)[-1]

  ####Check directories with present projections####
  if(!is.null(present_dir)) {
    #Check folders
    if(!file.exists(present_dir)){
      stop(paste("present_dir", present_dir, "does not exist"))
    }

    #List internal directories
    internal_dir <- list.dirs(present_dir, recursive = T)[-1]

    if(length(internal_dir) == 0) {
      stop(paste("present_dir", present_dir, "is empty"))
    }

    #Check if there are several scenarios
    pdir <- list.dirs(present_dir)[-1]
    #Check if there is a file in the directory
    fdir <- list.files(present_dir, pattern = raster_pattern)
    if(length(pdir) == 0 & length(fdir) == 0) {
      stop("present_dir ", present_dir, " has no folders or files to project")
    }

    #To project using a single file
    if(length(pdir) == 0 & length(fdir) > 0) {
      r <- rast(file.path(present_dir, fdir))
      #Check absent vars
      abs_vars <- vars[!(vars %in% names(r))]
      if(length(abs_vars) > 0)
        stop("The following vars are absent from current projection folder", ":\n",
             paste(abs_vars, collapse = "\n"))
      #If everything is OK, create list with the path
      res_present <- list()
      res_present[["Present"]] <- present_dir
    }

    #To project using folders
    if(length(pdir) > 0) {
      #Check if variables are in the folders
      invisible(sapply(pdir, function(x){
        r <- list.files(x, pattern = raster_pattern, full.names = TRUE)
        if(length(r) == 0) {
          stop("The directory ", x, " has no ", raster_pattern, " files")
        }
        r <- rast(r)
        #Check absent vars
        abs_vars <- vars[!(vars %in% names(r))]
        if(length(abs_vars) > 0)
          stop("The following vars are absent from ", x, ":\n",
               paste(abs_vars, collapse = "\n"))
      }))

      #If everything is OK, create list with the path
      #Get scenarios
      sc <- gsub(present_dir, "", pdir)
      res_present <- list()
      for (scenario in sc) {
        res_present[[scenario]] <- file.path(present_dir, scenario)}
    }

  } #End of check present

  #####Check directories with future projections####
  if(!is.null(future_dir)) {

    #Check folders
    if(!file.exists(future_dir)){
      stop(paste("future_dir", future_dir, "does not exist"))
    }

    #List internal directories
    internal_dir <- list.dirs(future_dir, recursive = T)[-1]

    if(length(internal_dir) == 0) {
      stop(paste("future_dir", future_dir, "is empty"))
    }

  if(sum(!is.null(future_time),
         !is.null(future_ssp),
         !is.null(future_gcm)) != 3) {
    stop("To prepare projections for future time, you must set future_time, future_ssp and future_gcm arguments")
  }

  #Check folders
  expected_folders_future <- unlist(sapply(future_time, function(i){
    lapply(future_ssp, function(x){
      file.path(future_dir, i, x, future_gcm)
    })
  }, simplify = F, USE.NAMES = FALSE))
  future_exists <- unlist(sapply(expected_folders_future, file.exists,
                                 simplify = TRUE, USE.NAMES = FALSE))

  if(any(future_exists == FALSE)){
    eff <- expected_folders_future[future_exists == FALSE]
    stop("The following folders do not exist: ", "\n", paste(eff, collapse = "\n"))
  }

  #Check if variables are in the folders
  all_dir <- expected_folders_future
  invisible(sapply(all_dir, function(x){
    r <- list.files(x, pattern = raster_pattern, full.names = TRUE)
    if(length(r) == 0) {
      stop("The directory ", x, " has no ", raster_pattern, " files")
    }
    r <- rast(r)
    #Check absent vars
    abs_vars <- vars[!(vars %in% names(r))]
    if(length(abs_vars) > 0)
    stop("The following vars are absent from ", x, ":\n",
         paste(abs_vars, collapse = "\n"))
  }))

  #If everything is OK, create list with the path
  res_future <- list()
  for (year in future_time) {
    res_future[[year]] <- list()
  for (ssp in future_ssp) {
    res_future[[year]][[ssp]] <- list()
    for (gcm in future_gcm) {
      res_future[[year]][[ssp]][[gcm]] <- file.path(future_dir, year, ssp, gcm)
    }}}

  } #End of check future

  #####Check directories with past projections####
  if(!is.null(past_dir)) {

    #Check folders
    if(!file.exists(past_dir)){
      stop(paste("past_dir", past_dir, "does not exist"))
    }

    #List internal directories
    internal_dir <- list.dirs(past_dir, recursive = T)[-1]

    if(length(internal_dir) == 0) {
      stop(paste("past_dir", past_dir, "is empty"))
    }

    if(sum(!is.null(past_time),
           !is.null(past_gcm)) != 2) {
      stop("To prepare projections for past time, you must set past_time and past_gcm arguments")
    }


  #Check folders
  expected_folders_past <- unlist(sapply(past_time, function(i){
      file.path(past_dir, i, past_gcm)
  }, simplify = F, USE.NAMES = FALSE))
  past_exists <- unlist(sapply(expected_folders_past, file.exists,
                                 simplify = TRUE, USE.NAMES = FALSE))

  if(any(past_exists == FALSE)){
    eff <- expected_folders_past[past_exists == FALSE]
    stop("The following folders do not exist: ", "\n", paste(eff, collapse = "\n"))
  }

  #Check if variables are in the folders
  all_dir <- expected_folders_past
  invisible(sapply(all_dir, function(x){
    r <- list.files(x, pattern = raster_pattern, full.names = TRUE)
    if(length(r) == 0) {
      stop("The directory ", x, " has no ", raster_pattern, " files")
    }
    r <- rast(r)
    #Check absent vars
    abs_vars <- vars[!(vars %in% names(r))]
    if(length(abs_vars) > 0)
      stop("The following vars are absent from ", x, ":\n",
           paste(abs_vars, collapse = "\n"))
  }))

  #If everything is OK, create list with the path
  res_past <- list()
  for (year in past_time) {
    res_past[[year]] <- list()
      for (gcm in past_gcm) {
        res_past[[year]][[gcm]] <- file.path(past_dir, year, gcm)
      }}
  }#End of check past

  #Append results
  res <- list()
  if(!is.null(present_dir)){
    res[["Present"]] <- res_present
  }
  if(!is.null(past_dir)){
    res[["Past"]] <- res_past
  }
  if(!is.null(future_dir)){
    res[["Future"]] <- res_future
  }

  #Append raster pattern
  res[["Raster_pattern"]] <- raster_pattern

  #Save results as RDS
  saveRDS(res, paste0(filename, ".RDS"))

  return(res)
} #End of function


# #Check function
# bm <- readRDS("../test_kuenm2/Best_Models.RDS")
#
# pr <- prepare_proj(models = bm,
#              present_dir = "../test_kuenm2/Projections/Present/",
#              past_dir = "../test_kuenm2/Projections/Past/",
#              past_time = c("LGM", "MID"),
#              past_gcm = c("CCSM4", "MIROC-ESM", "MPI-ESM-P"),
#              future_dir = "../test_kuenm2/Projections/Future/",
#              future_time = c("2041-2060", "2081-2100"),
#              future_ssp = c("ssp245", "ssp585"),
#              future_gcm = c("BCC-CSM2-MR", "ACCESS-CM2", "CMCC-ESM2"),
#              filename = "../test_kuenm2/Projection_file",
#              raster_pattern = ".tif*")
#
#
#
# bm <- readRDS("../test_kuenm2/Best_Models.RDS")
# models <- bm
# future_dir <- "../test_kuenm2/Projections/Future/"
# future_time <- "2041-2060"
# future_ssp <- c("ssp245", "ssp585")
# future_gcm <- c("BCC-CSM2-MR", "ACCESS-CM2", "CMCC-ESM2")
# raster_pattern = ".tif*"
# past_dir <- "../test_kuenm2/Projections/Past/"
# past_time = c("LGM", "MID")
# past_gcm = c("CCSM4", "MIROC-ESM", "MPI-ESM-P")
# present_dir <- "../test_kuenm2/Projections/Present/"
# filename = "../test_kuenm2/Projection_file"
