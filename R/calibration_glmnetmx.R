#Looping through candidate models using cores
calibration_glmnetmx <- function(data, #Data in **CLASS??** format
          #pr_bg, #Column name with presence (1) or background (0)
          formula_grid, #Grid with formulas
          test_convex = TRUE, #Test convex curves in quadratic models?
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
          skip_existing_models = FALSE, #Only works if write_summary = TRUE
          verbose = TRUE){

  #Args to parallel
  to_export <- c("aic_nk", "aic_ws", "eval_stats","glmnet_mx",
                "maxnet.default.regularization", "omrat",
                "predict.glmnet_mx", "empty_replicates",
                "empty_summary", "hinge", "hingeval", "thresholds",
                "thresholdval", "categorical", "categoricalval")

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

  #Make cluster
  cl <- parallel::makeCluster(ncores)

  #If test_convex = TRUE
  if(test_convex){
    if(verbose){
      cat("\n
        Task 1 of 2: checking convex curves in quadratic models\n")
    }
    #Get only quadratic grids with higher regularization multiplier
    q_grids <- formula_grid[grepl("q", formula_grid$Features) &
                              formula_grid$regm == max(formula_grid$regm), ]

    #Set number of iteration
    n <- nrow(q_grids)
    #If n = 0, do not run
    if(n == 0) {
      warning("All quadratic models have been already tested")
    } else {
      #Show progress bar?
      if (isTRUE(progress_bar)) {
        pb <- txtProgressBar(0, n, style = 3)
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
  #Start results convex
 results_convex <- foreach::foreach(
        x = 1:n, .options.snow = opts,
        .export = to_export) %dopar% {
 #Get grid x
 grid_x <- q_grids[x,] #Get i candidate model
 formula_x <- as.formula(grid_x$Formulas) #Get formula from grid x
 reg_x <- grid_x$regm #Get regularization multiplier from grid x

 #Complete model with AIC
 m_aic <- try(glmnet_mx(p = data$calibration_data$pr_bg,
                        data = data$calibration_data,
           f = formula_x, regmult = reg_x, calculate_AIC = T),
     silent = TRUE)
 if(any(class(m_aic) == "try-error")) {
      npar <- NA
      AICc <- NA
      is_c <- NA
      mods <- NA
      class(mods) <- "try-error"
      } else {

 #Get number of parameters
 npar <- length(m_aic$betas)

 #Calculate AIC from Warren
 vals <- predict.glmnet_mx(m_aic,
                           data$calibration_data[data$calibration_data$pr_bg == 1,],
                           type = "exponential")

 AICc <- aic_ws(pred_occs = vals, ncoefs = npar)


 #Check if model has convex curves
 m_betas <- m_aic$betas
 # #Check plot
 # p_aic <- predict.glmnet_mx(m_aic, data, type = "cloglog")
 # plot_aic <- cbind(data, pred = p_aic)
 # plot_aic[,c(5,14)] %>% plot()

 #Select only quadratic betas
 q_betas <- m_betas[grepl("\\^2", names(m_betas))]
 #Check if is convex
 if(length(q_betas) == 0) {
   is_c <- FALSE} else {
     is_c <- any(q_betas > 0)
   }
      }

 #If is convex, write results and check another combination
 if(isTRUE(is_c) | is.na(is_c)){
   ####Save metrics in a dataframe
   all_reg <- unique(formula_grid$regm)
   grid_q <- do.call("rbind",
         lapply(seq_along(all_reg), function(k){
         grid_x_i <- grid_x
         grid_x_i$regm <- all_reg[k]
         #Check id
         grid_x_i$ID<- formula_grid[formula_grid$Formulas ==
                                      grid_x$Formulas &
                                      formula_grid$regm == all_reg[k], "ID"]
        return(grid_x_i)
  }))

   df_eval_q <- empty_replicates(omrat_thr = omrat_thr,
    n_row = nrow(grid_q)*length(data$kfolds),
    replicates = names(data$kfolds), is_c = is_c)

   df_eval_q2 <- cbind(grid_q, df_eval_q)

   #Summarize results
   eval_final_q <- empty_summary(omrat_thr = omrat_thr,
      is_c = is_c)
   eval_final_q <- cbind(grid_q, eval_final_q)

 } else {#If is not convex, keep calculating metrics
   bgind <- which(data$calibration_data == 0)
   mods <- lapply(1:length(data$kfolds), function(i) {
     notrain <- -data$kfolds[[i]]
     data_i <- data$calibration_data[notrain,]
    #Run model
     mod_i <- glmnet_mx(p = data_i$pr_bg, data = data_i,
     f = formula_x, regmult = reg_x,
     calculate_AIC = FALSE)

     #Predict model only to background
     pred_i <- as.numeric(predict(object = mod_i,
                                  newdata = data$calibration_data,
                                  clamp = FALSE, type = "cloglog"))

     #Extract suitability in train and test points
     suit_val_cal <- pred_i[unique(c(notrain, -bgind))]
     suit_val_eval <- pred_i[which(!-notrain %in% bgind)]

     ####Calculate omission rate following kuenm####
     om_rate <- omrat(threshold = omrat_thr,
         pred_train = suit_val_cal,
         pred_test = suit_val_eval)
     #Calculate pROC following enmpa
     proc_i <- enmpa::proc_enm(test_prediction = suit_val_eval,
  prediction = pred_i)


    ####Save metrics in a dataframe
     df_eval_q <- data.frame(Replicate = i,
         t(data.frame(om_rate)),
         proc_auc_ratio = proc_i$pROC_summary[1],
         proc_pval = proc_i$pROC_summary[2],
         AIC_nk = m_aic$AIC,
         AIC_ws = AICc,
         npar = npar,
         is_convex = is_c,
         row.names = NULL)
     df_eval_q2 <- cbind(grid_x, df_eval_q)
     return(df_eval_q2)
   })
   names(mods) <- names(data$kfolds)
   #Return evaluation final
   eval_final_q <- do.call("rbind", mods)

   #Summarize results?
   eval_final_q_summary <- eval_stats(eval_final_q)

 } #End of else

 #If write_summary = T...
 if(write_summary){
   write.csv(eval_final_q_summary, file.path(out_dir,
        paste0("Summary_cand_model_", grid_x$ID, ".csv")),
    row.names = F)
 }

 #Return replicates?
 if(!return_replicate)
   eval_final_q <- NULL

 return(list(All_results = eval_final_q,
             Summary = eval_final_q_summary))


        } #End of first foreach
    } #End of if(n > 0)
  }  #End of If test_convex = TRUE

  #Update grid
  if(!test_convex) {n = 0}
  if(test_convex & n > 0){

      #Convert results to dataframe
    #Replicate
    d_convex_rep <- do.call("rbind", lapply(results_convex,
                                          function(x) x$Replicates))
    row.names(d_convex_rep) <- NULL
    #Summary
    d_convex_sum <- do.call("rbind", lapply(results_convex, function(x) x$Summary))


    # d_convex <- do.call("rbind", results_convex)
    # d_convex <- d_convex[d_convex$is_convex == TRUE, ]

    #Identify formulas tested with convex curves
    formula_grid2 <- formula_grid[!(formula_grid$ID %in% d_convex_sum$ID), ]
  } else {
    formula_grid2 <- formula_grid
  }

  #Set number of iteration based on new grid
  n <- nrow(formula_grid2)
  #If n == 0, skip non-quadratic models
  if(n == 0) {
    warning("All non-quadratic models have been already tested")
  } else {

    #Show progress bar? - Update
    if (isTRUE(progress_bar)) {
      pb <- txtProgressBar(0, n, style = 3)
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
      if(test_convex) {
        cat("\nTask 2 of 2: calibrating non-quadratic models and quadratic models
        without convex curves\n
        ")
      } else {
        cat("
        Task 1 of 1: calibrating models\n")
      } }
    ####Results non-convex####
    results <- foreach::foreach(x = 1:n, .options.snow = opts,
   .export = to_export) %dopar% {
     #Get grid x
     grid_x <- formula_grid2[x,] #Get i candidate model
     formula_x <- as.formula(grid_x$Formulas) #Get formula from grid x
     reg_x <- grid_x$regm #Get regularization multiplier from grid x

     #Complete model with AIC
     #Complete model with AIC
     m_aic <- try(glmnet_mx(p = data$calibration_data$pr_bg,
                            data = data$calibration_data,
                            f = formula_x, regmult = reg_x, calculate_AIC = T),
                  silent = TRUE)

     if(any(class(m_aic) == "try-error")) {
       npar <- NA
       AICc <- NA
       is_c <- NA
       mods <- NA
       class(mods) <- "try-error"
     } else
       {

       #Get number of parameters
       npar <- length(m_aic$betas)

       #Calculate AIC from Warren
       vals <- predict.glmnet_mx(m_aic,
                                 data$calibration_data[data$calibration_data$pr_bg == 1,],
                                 type = "exponential")

       AICc <- aic_ws(pred_occs = vals, ncoefs = npar)

       #Check if model has convex curves
       m_betas <- m_aic$betas
       #Select only quadratic betas
       q_betas <- m_betas[grepl("\\^2", names(m_betas))]
       #Check if is convex
       if(length(q_betas) == 0) {
         is_c <- FALSE} else {
  is_c <- any(q_betas > 0)
         }
       #Get background index (again, if necessary)
       bgind <- which(data$calibration_data == 0)
       mods <- try(lapply(1:length(data$kfolds), function(i) {
         notrain <- -data$kfolds[[i]]
         data_i <- data$calibration_data[notrain,]
         #Run model
         mod_i <- glmnet_mx(p = data_i$pr_bg, data = data_i,
                            f = formula_x, regmult = reg_x,
                            calculate_AIC = FALSE)

         #Predict model only to background
         pred_i <- as.numeric(predict(object = mod_i,
                                      newdata = data$calibration_data,
                                      clamp = FALSE, type = "cloglog"))

         # #Predict model to spatraster, only to check
         # vars_r <- terra::rast("Models/Araucaria_angustifolia/PCA_variables.tiff")
         # pred_r <- terra::predict(vars_r, mod_i, type = "cloglog", na.rm = TRUE)
         # plot(pred_r)

         #Extract suitability in train and test points
         suit_val_cal <- pred_i[unique(c(notrain, -bgind))]
         suit_val_eval <- pred_i[which(!-notrain %in% bgind)]

         ####Calculate omission rate following kuenm####
         om_rate <- omrat(threshold = omrat_thr,
                          pred_train = suit_val_cal,
                          pred_test = suit_val_eval)
         #Calculate pROC following enmpa
         proc_i <- enmpa::proc_enm(test_prediction = suit_val_eval,
                                   prediction = pred_i)

         ####Save metrics in a dataframe
         df_eval <- data.frame(Replicate = i, t(data.frame(om_rate)),
                               proc_auc_ratio = proc_i$pROC_summary[1],
                               proc_pval = proc_i$pROC_summary[2],
                               AIC_nk = m_aic$AIC,
                               AIC_ws = AICc,
                               npar = npar,
                               is_convex = is_c,
                               row.names = NULL)
         df_eval2 <- cbind(grid_x, df_eval)
         return(df_eval2)
       }), silent = TRUE)
       names(mods) <- names(data$kfolds)
       }

     #####Create empty dataframe if mods is an error####
     if(class(mods) == "try-error") {
       eval_final <- cbind(grid_x, empty_replicates(omrat_thr = omrat_thr,
     n_row = length(data$kfolds),
     replicates = names(data$kfolds), is_c = is_c))
     } else{
       #Return evaluation final
       eval_final <- do.call("rbind", mods) }

     #Summarize results
     if(class(mods) == "try-error") {
         eval_final <- cbind(grid_x,
           empty_summary(omrat_thr = omrat_thr, is_c = is_c))
       } else {
         eval_final_summary <- eval_stats(eval_final) }


     #If write_summary = T...
     if(write_summary){
       write.csv(eval_final_summary, file.path(out_dir,
 paste0("Summary_cand_model_", grid_x$ID, ".csv")),
        row.names = F)
     }


     if(!return_replicate)
       eval_final <- NULL

     return(list(All_results = eval_final,
                 Summary = eval_final_summary)) } #End of foreach

    #Stop cluster
    parallel::stopCluster(cl)

    #Convert to dataframe
      #Replicate
    d_res_rep <- do.call("rbind", lapply(results, function(x) x$Replicates))
    row.names(d_res_rep) <- NULL
      #Summary
    d_res_sum <- do.call("rbind", lapply(results, function(x) x$Summary))

    # Join results with results convex, if it exists
    if(test_convex) {
      replicates_final <- rbind(d_convex_rep, d_res_rep)
      summary_final <- rbind(d_convex_sum, d_res_sum)
      res_final <- list(All_results = replicates_final,
                        Summary = summary_final)
    } else {
      res_final <- list(All_results =  d_res_rep,
                        Summary = d_res_sum)
    }

  } #End of if(n == 0)

  #Select best models
  bm <- sel_best_models(cand_models = res_final$Summary, #dataframe with Candidate results
                        test_convex = test_convex, #Remove models with concave curves?
                        omrat_threshold = omrat_threshold, #Omission rate (test points) used to select the best models
                        allow_tolerance = T, #If omission rate is higher than set, select the model with minimum omission rate
                        tolerance = 0.01, #If allow tollerance, select the model with minimum omission rate + tolerance
                        AIC = "nk", #Which AIC? japones (nk) or Warrien (ws?
                        significance = 0.05, #Significante to select models based on pROC
                        verbose = verbose, #Show messages?
                        delta_aic = 2, #Delta AIC to select best models
                        save_file = F, #Save file with best models?
                        file_name = NULL)
  #Concatenate final results
  return(c(data, calibration_results = list(res_final),
             selected_models = list(bm)))

} #End of function

