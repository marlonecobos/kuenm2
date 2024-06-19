#' Projection of Maxent-like glmnet models to single or multiple scenarios
#'
#' @export

project_selected_glmnetx <- function(models,
                                     projection_file,
                                     out_dir = "Projection_results",
                                     #write_path = TRUE,
                                     consensus_per_model = TRUE,
                                     consensus_general = TRUE,
                                     consensus = c("median", "range", "mean", "stdev"), #weighted mean
                                     write_replicates = FALSE,
                                     clamping = FALSE,
                                     var_to_clamp = NULL,
                                     type = "cloglog",
                                     overwrite = FALSE,
                                     parallel = FALSE,
                                     ncores = 1,
                                     parallelType = "doSNOW",
                                     progress_bar = TRUE,
                                     verbose = TRUE){

  if(!any(c("median", "mean") %in% consensus)){
    stop("Consensus must have at least one of the options: 'median' or 'mean'")
  }

  #Save parameters in a list to send to foreach nodes#
  par_list <- list(models = models,
                   projection_file = projection_file,
                   consensus_per_model = consensus_per_model,
                   consensus_general = consensus_general,
                   consensus = consensus,
                   write_replicates = write_replicates,
                   clamping = clamping,
                   var_to_clamp = var_to_clamp,
                   type = type,
                   overwrite = overwrite)

  ####PREPARE DATAFRAME TO PREDICT####
  #Extract variables from best models
  vars <- names(models[["Models"]][[1]][[1]]$samplemeans)
  vars <- setdiff(vars, c("pr_bg", "fold"))

  if(!file.exists(out_dir)){
    dir.create(out_dir, showWarnings = FALSE)
  }
  #Normalize path
  out_dir <- normalizePath(out_dir)

  #Check scenarios to predict
  sc <- names(projection_file)

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
                                                                  present_sc)))})
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
                                                     mustWork = FALSE))})
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
                                                                 df_future$ssp, df_future$GCM),
                                                       mustWork = FALSE))})
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


  ####Configure parallelization####
  n_models <- nrow(res_path)
  if(n_models == 1 & isTRUE(parallel)){
    parallel <- FALSE
  } else {parallel <- TRUE}

  if(n_models < ncores & isTRUE(parallel)){
    ncores <- n_models
  } else {ncores = ncores}

  #Show progress bar?
  if (progress_bar) {
    pb <- txtProgressBar(0, n_models, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n) }


  #Set parallelization
  if(parallel) {
    #Make cluster
    cl <- parallel::makeCluster(ncores)

    if (parallelType == "doParallel") {
      doParallel::registerDoParallel(cl)
      opts <- NULL
    }

    if (parallelType == "doSNOW") {
      doSNOW::registerDoSNOW(cl)
      if (progress_bar)
        opts <- list(progress = progress)
      else opts <- NULL
    }
  } else {
    opts <- NULL}
  ###################################

  #Run predictions
  if(parallel){
    foreach(x = 1:n_models,
            .options.snow = opts
            ,
            .export = c("predict_selected_glmnetmx")
    ) %dopar% {
      multiple_projections(i = x, res_path, raster_pattern, par_list)
    }
  } else { #Not in parallel (using %do%)
    # Loop for com barra de progresso manual
    for (x in 1:n_models) {
      # Execute a função fit_eval_models
      multiple_projections(i = x, res_path, raster_pattern, par_list)

      # Sets the progress bar to the current state
      if(progress_bar){
        setTxtProgressBar(pb, x) }
    }
  }

  #Stop cluster
  if(parallel){parallel::stopCluster(cl)}

  #Append threshold to final results
  res_final <- list(paths = res_path,
                    threshold = models$thresholds)

  #Save
  saveRDS(res_final, file.path(out_dir, "Projection_paths.RDS"))
  return(res_final)
} #End of function

# #Test function internally
# models = readRDS("../test_kuenm2/Best_Models.RDS")
# projection_file = readRDS("../test_kuenm2/Projection_file.RDS")
# out_dir = "../test_kuenm2/Projection_results"
# #write_path = TRUE
# consensus_per_model = TRUE
# consensus_general = TRUE
# consensus = c("median", "range", "mean", "stdev") #weighted mean
# write_replicates = FALSE
# clamping = FALSE
# var_to_clamp = NULL
# type = "cloglog"
# overwrite = TRUE
# parallel = TRUE
# ncores = 8
# parallelType = "doSNOW"
# progress_bar = TRUE
# verbose = TRUE
