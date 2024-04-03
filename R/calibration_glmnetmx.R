#Looping through candidate models using cores
calibration_glmnetmx <- function(data, #Data in **CLASS??** format (includes weights)
                                 #pr_bg, #Column name with presence (1) or background (0)
                                 formula_grid, #Grid with formulas
                                 test_concave = TRUE, #Test concave curves in quadratic models?
                                 addsamplestobackground = TRUE,
                                 use_weights = FALSE,
                                 #weights = NULL,
                                 #folds = 4, #Columns name with k_folds or vector indicating k_folds
                                 parallel = TRUE,
                                 ncores = 1,
                                 progress_bar = TRUE, #Show progress bar? Only works if parallel_type = "doSNOW"
                                 write_summary = FALSE, #Write summary of each candidate evaluations?
                                 out_dir = NULL, #Name of the folder to write candidate evaluations
                                 parallel_type = "doSNOW",
                                 return_replicate = TRUE,
                                 omrat_thr = c(5, 10),
                                 omrat_threshold = 10,
                                 AIC = "ws",
                                 delta_aic = 2,
                                 allow_tolerance = TRUE, #If omission rate is higher than set, select the model with minimum omission rate
                                 tolerance = 0.01,
                                 skip_existing_models = FALSE, #Only works if write_summary = TRUE
                                 verbose = TRUE){

  #Args to parallel
  to_export <- c("aic_nk", "aic_ws", "eval_stats","glmnet_mx",
                 "maxnet.default.regularization", "omrat",
                 "predict.glmnet_mx", "empty_replicates",
                 "empty_summary", "hinge", "hingeval", "thresholds",
                 "thresholdval", "categorical", "categoricalval",
                 "fit_eval_concave", "fit_eval_models", "omrat_thr",
                 "omrat_threshold", "sel_best_models")

  #Check calibration data class and convert to dataframe if necessay
  if(is.matrix(data$calibration_data) | is.array(data$calibration_data)) {
    data$calibration_data <- as.data.frame(data$calibration_data)
  }

  #If write_summary = TRUE, create directory
  if(write_summary){
    if(!file.exists(out_dir))
      dir.create(out_dir)
  }

  #If skip_existing_models = TRUE, update grid
  if(skip_existing_models & write_summary) {
    ready_models <- list.files(path = out_dir, pattern = "summary",
                               full.names = T)
    ready_models <- do.call("rbind", lapply(seq_along(ready_models), function(i){
      read.csv( ready_models[i])
    }))
    run_models <- setdiff(formula_grid$ID, ready_models$ID)
    if(length(run_models) == 0) {
      stop(paste("All models completed. Check the folder:", out_dir))
    } else { #Update formula grid
      formula_grid <- formula_grid[formula_grid$ID %in% run_models, ]
    }
  }

  #Warning about samples added to background when weights is null
  if (verbose & addsamplestobackground & use_weights & is.null(data$weights)) {
  message("weights for samples added to background are the same as in samples.")
  }

  # Check weights
  if (use_weights & is.null(data$weights)) {
    message("'use_weights' = TRUE, but weights are not present in 'data'.\nSetting 'use_weights' = FALSE.")
    use_weights <- FALSE
  }

  # If weights is numeric, check if it has the same size of calibration data
    if(use_weights & length(data$weights) != nrow(data$calibration_data)) {
      stop("length of weights does not match number of rows in calibration_data")
    }


  if (parallel) {
  #Make cluster
  cl <- parallel::makeCluster(ncores) }

  #If test_concave = TRUE
  if (test_concave) {
    if(verbose){
      cat("\n
        Task 1 of 2: checking concave curves in quadratic models\n")
    }
    #Get only quadratic grids with higher regularization multiplier
    q_grids <- formula_grid[grepl("q", formula_grid$Features) &
                              formula_grid$regm == max(formula_grid$regm), ]
    #Set number of iteration
    n_tot <- nrow(q_grids)
    #If n = 0, do not run
    if(n_tot == 0) {
      warning("All quadratic models have been already tested")
    } else {
      #Show progress bar?
      if (isTRUE(progress_bar)) {
        pb <- txtProgressBar(min = 0, max = n_tot, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)}

      if (parallel & parallel_type == "doParallel") {
        doParallel::registerDoParallel(cl)
        opts <- NULL #Progress bar does not work with doParallel
      }

      if (parallel & parallel_type == "doSNOW") {
        doSNOW::registerDoSNOW(cl)
        if (isTRUE(progress_bar))
          opts <- list(progress = progress)
        else opts <- opts
      }

      #Test concave curves
      #In parallel (using %dopar%)
      if(parallel){
      results_concave <- foreach(
        x = 1:n_tot, .packages = c("glmnet", "enmpa"), .options.snow = opts,
        .export = c(to_export, "formula_grid", "q_grids","data",
                    "write_summary", "return_replicate")
      ) %dopar% {
        fit_eval_concave(x = x, q_grids, data, formula_grid,
                         omrat_thr,write_summary,
                         addsamplestobackground,
                         weights = data$weights,
                         return_replicate)
        }
      } else { #Not in parallel (using %do%)
        results_concave <- vector("list", length = n_tot)
        # Loop for com barra de progresso manual
        for (x in 1:n_tot) {
          # Execute a função fit_eval_models
          results_concave[[x]] <- fit_eval_concave(
            x = x, q_grids, data, formula_grid, omrat_thr = omrat_thr,
            write_summary = write_summary,
            addsamplestobackground = addsamplestobackground,
            weights = data$weights, return_replicate = return_replicate
          )

          # Sets the progress bar to the current state
          if(progress_bar){
            setTxtProgressBar(pb, x) }
        }
      }

  } #End of if(n > 0)
  }  #End of If test_concave = TRUE


  #Update grid
  if(!test_concave) {n_tot = 0}
  if(test_concave & n_tot > 0){

    #Convert results to dataframe
    #Replicate
    d_concave_rep <- do.call("rbind", lapply(results_concave,
                                            function(x) x$Replicates))
    row.names(d_concave_rep) <- NULL
    #Summary
    d_concave_sum <- do.call("rbind", lapply(results_concave,
                                             function(x) x$Summary))


    # d_concave <- do.call("rbind", results_concave)
    # d_concave <- d_concave[d_concave$is_concave == TRUE, ]

    #Identify formulas tested with concave curves
    formula_grid2 <- formula_grid[!(formula_grid$ID %in% d_concave_sum$ID), ]
  } else {
    formula_grid2 <- formula_grid
  }

  #Set number of iteration based on new grid
  n_tot <- nrow(formula_grid2)
  #If n == 0, skip non-quadratic models
  if(n_tot == 0) {
    warning("All non-quadratic models have been already tested")
  } else {

    #Show progress bar? - Update
    if (progress_bar) {
      pb <- txtProgressBar(0, n_tot, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n) }

    if (parallel_type == "doParallel") {
      doParallel::registerDoParallel(cl)
      opts <- NULL
    }

    if (parallel_type == "doSNOW") {
      doSNOW::registerDoSNOW(cl)
      if (isTRUE(progress_bar))
        opts <- list(progress = progress)
      else opts <- NULL
    }

    if(verbose) {
      if(test_concave) {
        cat("\nTask 2 of 2: calibrating non-quadratic models and quadratic models
        without concave curves\n
        ")
      } else {
        cat("
        Task 1 of 1: calibrating models\n")
      } }

    ####Results non-concave####
    #Test concave curves
    #In parallel (using %dopar%)
    if(parallel){
    results <- foreach(
      x = 1:n_tot, .packages = c("glmnet", "enmpa"), .options.snow = opts,
      .export = c(to_export, "formula_grid2", "q_grids", "data",
                  "write_summary", "return_replicate")
    ) %dopar% {
      fit_eval_models(x, formula_grid2, data,
                      formula_grid, omrat_thr,
                      write_summary,
                      addsamplestobackground = addsamplestobackground,
                      weights = data$weights,
                      return_replicate)
      }
    } else { #Not in parallel (using %do%)
      results <- vector("list", length = n_tot)
      # Loop for com barra de progresso manual
      for (x in 1:n_tot) {
        # Execute a função fit_eval_models
        results[[x]] <- fit_eval_models(
          x, formula_grid2 = g, data = data, formula_grid = g, omrat_thr,
          addsamplestobackground =  addsamplestobackground,
          weights = data$weights, write_summary, return_replicate
        )

        # Sets the progress bar to the current state
        if(progress_bar){
        setTxtProgressBar(pb, x) }
      }
    }

    #Stop cluster
    if(parallel){
    parallel::stopCluster(cl) }

    #Convert to dataframe
    #Replicate
    d_res_rep <- do.call("rbind", lapply(results, function(x) x$Replicates))
    row.names(d_res_rep) <- NULL
    #Summary
    d_res_sum <- do.call("rbind", lapply(results, function(x) x$Summary))

    # Join results with results concave, if it exists
    if(test_concave) {
      replicates_final <- rbind(d_concave_rep, d_res_rep)
      summary_final <- rbind(d_concave_sum, d_res_sum)
      res_final <- list(All_results = replicates_final,
                        Summary = summary_final)
    } else {
      res_final <- list(All_results =  d_res_rep,
                        Summary = d_res_sum)
    }

  } #End of if(n == 0)

  #Select best models
  bm <- sel_best_models(cand_models = res_final$Summary, #dataframe with Candidate results
                        test_concave = test_concave, #Remove models with concave curves?
                        omrat_threshold = omrat_threshold, #Omission rate (test points) used to select the best models
                        allow_tolerance = allow_tolerance, #If omission rate is higher than set, select the model with minimum omission rate
                        tolerance = tolerance, #If allow tollerance, select the model with minimum omission rate + tolerance
                        AIC = AIC, #Which AIC? japones (nk) or Warrien (ws?
                        significance = 0.05, #Significante to select models based on pROC
                        verbose = verbose, #Show messages?
                        delta_aic = delta_aic, #Delta AIC to select best models
                        save_file = F, #Save file with best models?
                        file_name = NULL)
  #Concatenate final results

  return(c(data, calibration_results = list(res_final),
           addsamplestobackground = addsamplestobackground,
           weights = list(data$weights), selected_models = list(bm)))

} #End of function

