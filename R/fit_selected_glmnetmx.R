#' Fit selected Maxent-like glmnet models
#'
#' @export

fit_selected_glmnetmx <- function(calibration_results,
                                  n_replicates = 5,
                                  rep_type = "kfold",
                                  train_portion = 0.7,
                                  write_models = FALSE, #Write files?
                                  file_name = NULL, #Name of the folder to write final models
                                  parallel = FALSE,
                                  ncores = 2,
                                  parallelType = "doSNOW",
                                  progress_bar = TRUE,
                                  verbose = TRUE,
                                  seed = 42) {
  # #Args
  # to_export <- c("aic_nk", "aic_ws",
  #               "eval_stats","glmnet_mx",
  #               "maxnet.default.regularization",
  #               "omrat","predict.glmnet_mx",
  #               "empty_replicates",
  #               "empty_summary", "hinge",
  #               "hingeval",
  #               "thresholds", "thresholdval",
  #               "categorical",
  #               "categoricalval", "fit_best_model", "n_replicates")

  #Extracts IDs from models
  m_ids <- calibration_results$selected_models$ID

  ####Fit replicates####

  if(n_replicates > 1){
    if(verbose){
      message("Fitting replicates...")
    }
  #Create grid of fitted models
  dfgrid <- expand.grid(models = m_ids, replicates = 1:n_replicates)
  n_tot <-  nrow(dfgrid)

  #Prepare data (index) to replicates
  if(n_replicates > 1) {
    #Partitioning data
    rep_data <- part_data(data = calibration_results$calibration_data,
                          pr_bg = "pr_bg",
                          train_portion = train_portion,
                          n_replicates = n_replicates,
                          method = rep_type, seed = seed)
  } else {
      rep_data <- NULL
    }


  #Show progress bar?
  if (isTRUE(progress_bar)) {
    pb <- txtProgressBar(0, n_tot, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n) }

  #Adjust parallelization according to the number of models and replicates
  if(n_tot == 1 & isTRUE(parallel)){
    parallel <- FALSE
  }
  if(n_tot < ncores & isTRUE(parallel)){
    ncores <- n_tot
  }

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

  #In parallel (using %dopar%)
  if(parallel){
    best_models <- foreach(x = 1:n_tot,
                               .options.snow = opts
                           # ,
                           #     .export = c(to_export)
    ) %dopar% {
      kuenm2:::fit_best_model(x = x, dfgrid, calibration_results, data_x, n_replicates,
                     rep_data)
    }
  } else { #Not in parallel (using %do%)
    best_models <- vector("list", length = n_tot)
    # Loop for com barra de progresso manual
    for (x in 1:n_tot) {
      # Execute a função fit_eval_models
      best_models[[x]] <- kuenm2:::fit_best_model(x = x, dfgrid, calibration_results,
                                         data_x, n_replicates,
                                         rep_data)

      # Sets the progress bar to the current state
      if(progress_bar){
        setTxtProgressBar(pb, x) }
    }  }

  #End cluster
  if(parallel){
  parallel::stopCluster(cl) }

  #Split list
  # # Crie um vetor com o número de réplicas para cada modelo
  # num_repl <- tapply(dfgrid$replicates, dfgrid$models, FUN = length)
  # num_repl <- num_repl[match(dfgrid$models, names(num_repl))]

  # Split list
  best_models <- split(best_models, dfgrid$models)
  best_models <- lapply(best_models, function(sublist) {
    names(sublist) <- paste0("Rep_", seq_along(sublist))
    return(sublist)
  })

  #Rename models
  names(best_models) <- paste0("Model_", names(best_models))} else {
    best_models = list()
  }

  ####Fit full models####
  if(verbose){
    message("\nFitting full models...")
  }
  n_models <- length(m_ids)
  #Create grid of fitted models
  dfgrid <- expand.grid(models = m_ids, replicates = 1)
  #Adjust parallelization according to the number of models and replicates
  if(n_models == 1 & isTRUE(parallel)){
    parallel <- FALSE
  }
  if(n_models < ncores & isTRUE(parallel)){
    ncores <- n_models
  }

  #Show progress bar?
  if (isTRUE(progress_bar)) {
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


  #In parallel (using %dopar%)
  if(parallel){
    full_models <- foreach(x = 1:n_models,
                           .options.snow = opts
                           # ,
                           #    .export = c(to_export)
    ) %dopar% {
      kuenm2:::fit_best_model(x = x, dfgrid, calibration_results, data_x,
                     n_replicates = 1)
    }
  } else { #Not in parallel (using %do%)
    full_models <- vector("list", length = n_models)
    # Loop for com barra de progresso manual
    for (x in 1:n_models) {
      # Execute a função fit_eval_models
      full_models[[x]] <- kuenm2:::fit_best_model(x = x, dfgrid, calibration_results,
                                         data_x, n_replicates = 1)

      # Sets the progress bar to the current state
      if(progress_bar){
        setTxtProgressBar(pb, x) }
    }  }

  #Names models
  names(full_models) <- paste0("Model_", m_ids)

  #Append full models to replicates
  for(i in names(full_models)) {
    best_models[[i]]$Full_model <- full_models[[i]]
  }


  #End cluster
  if(parallel){
    parallel::stopCluster(cl) }

  #Compute threshold
  #Get models names
  nm <- names(best_models)

  #Get occurrences
  occ <- calibration_results$calibration_data[
    calibration_results$calibration_data$pr_bg == 1, -1]

  #Predict models to occurrences - Per replicate
  p_occ <- lapply(nm, function(x){
    m_x <- best_models[[x]]
    #Remove full model if replicates exist
    if(any(grepl("Rep", names(m_x)))){
      m_x$Full_model <- NULL
    }
    #Predict by replicate
    p_r <- sapply(m_x, function(i){
      predict.glmnet_mx(object = i, newdata = occ, type = "cloglog")
    })
    #Get mean and median
    p_mean <- apply(p_r, 1, mean, na.rm=T)
    p_median <- apply(p_r, 1, median, na.rm=T)
    #Save in a list
    p_res <- list(mean = p_mean,
                  median = p_median)
  })
  names(p_occ) <- nm
  #Get consensus
  mean_consensus <- apply(sapply(p_occ, function(x) x$mean),
                          1, mean, na.rm = TRUE)
  median_consensus <- apply(sapply(p_occ, function(x) x$median),
                          1, median, na.rm = TRUE)
  consensus <- list(mean = mean_consensus,
                    median = median_consensus)
  #Merge list
  p_occ <- c(p_occ, list(consensus = consensus))
  #Calculate threshold
  p_thr <- lapply(p_occ, function(model) {
      lapply(model, calc_thr,
             thr = calibration_results$omission_rate/100)
    })

  #Final results
  res <- list(species = calibration_results$species,
              Models = best_models,
              calibration_data = calibration_results$calibration_data,
              selected_models = calibration_results$selected_models,
              weights = calibration_results$weights,
              addsamplestobackground = calibration_results$addsamplestobackground,
              omission_rate = calibration_results$omission_rate,
              thresholds = p_thr)

  #Write models?
  if(write_models){
    saveRDS(res, file = paste0(file_name, ".RDS"))
  }

  return(res)
} #End of function


# # #Test function
# calibration_results = readRDS("../test_kuenm2/Candidate_and_Best_models.RDS")
# n_replicates = 5
# rep_type = "kfold"
# train_portion = 0.7
# write_models = FALSE #Write files?
# file_name = NULL #Name of the folder to write final models
# parallel = TRUE
# ncores = 5
# parallelType = "doSNOW"
# progress_bar = TRUE
# verbose = TRUE
