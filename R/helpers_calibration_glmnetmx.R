#' Summarize evaluation results
#'
#' @importFrom stats aggregate
#' @importFrom enmpa proc_enm
#' @export

eval_stats <- function(calib_results, omrat_thr){
  omission_rates <- paste0("Omission_rate_at_", omrat_thr)
  toagg <- c(omission_rates,
               "proc_auc_ratio", "proc_pval")
  agg_by <- c("Formulas", "regm", "Features")
  to_keep <- c("ID", "Formulas", "regm", "Features", "AIC_nk",
                 "AIC_ws", "npar", "is_concave")
  agg_formula <- paste("~", paste(agg_by, collapse = " + "))

  #Get summary
  xy <- lapply(toagg, function(x) {
    do.call(
      data.frame,
      stats::aggregate(as.formula(paste(x, agg_formula)),
                       data = calib_results, FUN = function(y) {
                         c(mean = round(mean(y), 4), sd = round(sd(y), 4))
                       }, na.action = NULL)
    )
  })

  #Summarize stats
  stats <- Reduce(function(x, y) merge(x, y, by = agg_by), xy)

  stats_AICS <- calib_results[!duplicated(calib_results[,to_keep]),][,to_keep]
  stats_final <- merge(stats, stats_AICS, by = agg_by)
  return(stats_final)
}


####Create empty dataframes####
empty_replicates <- function(omrat_thr, n_row = 4, replicates = 1:4,
                             is_c = NA) {
  column_names <- c("Replicate", paste0("Omission_rate_at_", omrat_thr),
                    "proc_auc_ratio", "proc_pval", "AIC_nk", "AIC_ws", "npar",
                    "is_concave")
  df_eval_q <- data.frame(matrix(NA, nrow = n_row, ncol = length(column_names)))
  colnames(df_eval_q) <- column_names
  df_eval_q$Replicate <- replicates
  df_eval_q$is_concave = is_c

  return(df_eval_q)
}

empty_summary <- function(omrat_thr, is_c) {
  om_means <- paste0("Omission_rate_at_", omrat_thr, ".mean")
  om_sd <- paste0("Omission_rate_at_", omrat_thr, ".sd")
  column_names <- c(om_means, om_sd,
                    "proc_auc_ratio.mean", "proc_auc_ratio.sd", "proc_pval.mean",
                    "proc_pval.sd", "AIC_nk", "AIC_ws", "npar", "is_concave")
  eval_final_q  <- data.frame(matrix(NA, nrow = 1, ncol = length(column_names)))
  colnames(eval_final_q) <- column_names
  eval_final_q$is_concave = is_c
  return(eval_final_q)
}

#Function to test concave curves
fit_eval_concave <- function(x, q_grids, data, formula_grid, omrat_thr,
                             write_summary, addsamplestobackground, weights,
                             return_replicate, AIC) {
  grid_x <- q_grids[x, ]
  formula_x <- as.formula(grid_x$Formulas)
  reg_x <- grid_x$regm

  #Complete model with AIC
  m_aic <- try(glmnet_mx(p = data$calibration_data$pr_bg,
                         data = data$calibration_data,
                         f = formula_x, regmult = reg_x,
                         addsamplestobackground = addsamplestobackground,
                         weights = weights,
                         calculate_AIC = T,
                         AIC = AIC),
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
    vals <- predict.glmnet_mx(
      m_aic, data$calibration_data[data$calibration_data$pr_bg == 1,],
      type = "exponential"
    )

    AICc <- aic_ws(pred_occs = vals, ncoefs = npar)


    #Check if model has concave curves
    #Select only quadratic betas
    q_betas <- m_aic$betas[grepl("\\^2", names(m_aic$betas))]
    #Check if is concave
    if(length(q_betas) == 0) {
      is_c <- FALSE} else {
        is_c <- any(q_betas > 0)
      }
  }

  #If is concave, write results and check another combination
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
    eval_final_q_summary <- cbind(grid_q, eval_final_q)

  } else {#If is not concave, keep calculating metrics
    bgind <- which(data$calibration_data == 0)
    mods <- lapply(1:length(data$kfolds), function(i) {
      notrain <- -data$kfolds[[i]]
      data_i <- data$calibration_data[notrain,]
      #Set weights per k-fold
      if(!is.null(weights)) {
        weights_i <- weights[notrain]
      } else {
        weights_i <- weights
      }
      #Run model
      mod_i <- glmnet_mx(p = data_i$pr_bg, data = data_i,
                         f = formula_x, regmult = reg_x,
                         addsamplestobackground = addsamplestobackground,
                         weights = weights_i,
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
                              is_concave = is_c,
                              row.names = NULL)
      return(cbind(grid_x, df_eval_q))
    })
    names(mods) <- names(data$kfolds)
    #Return evaluation final
    eval_final_q <- do.call("rbind", mods)

    #Summarize results?
    eval_final_q_summary <- eval_stats(eval_final_q, omrat_thr)

  } #End of else

  #If write_summary = T...
  if(write_summary){
    write.csv(eval_final_q_summary,
              file.path(out_dir,
                        paste0("Summary_cand_model_", grid_x$ID, ".csv")),
              row.names = F)
  }

  #Return replicates?
  if(!return_replicate)
    eval_final_q <- NULL

  return(list(All_results = eval_final_q,
              Summary = eval_final_q_summary))
}

#Function to test all models curves (except quadratic models when test_concave = TRUE)
fit_eval_models <- function(x, #Each line of the formula grid
                            formula_grid, #Formula grid (output of calibration_grid_glmnetmx)
                            data, #Data to fit models (output of prepare_data)
                            omrat_thr, #Omission rate threshold
                            write_summary, #Write individual summary for each line in formula grid?
                            addsamplestobackground, #Add samples to background?
                            weights, #Add weights?
                            return_replicate, #Return results of replicates
                            AIC) { #Which AIC: ws (Warren) or nk (Nikimura)
  #Get grid x
  grid_x <- formula_grid[x,] #Get i candidate model
  formula_x <- as.formula(grid_x$Formulas) #Get formula from grid x
  reg_x <- grid_x$regm #Get regularization multiplier from grid x

  #Complete model with AIC
  m_aic <- try(glmnet_mx(p = data$calibration_data$pr_bg,
                         data = data$calibration_data,
                         f = formula_x, regmult = reg_x,
                         addsamplestobackground = addsamplestobackground,
                         weights = weights,
                         calculate_AIC = T,
                         AIC = AIC),
               silent = TRUE)

  #If it's not possible fit the model, return evaluation metric as NA
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
    vals <- predict.glmnet_mx(
      m_aic, data$calibration_data[data$calibration_data$pr_bg == 1,],
      type = "exponential"
    )

    AICc <- aic_ws(pred_occs = vals, ncoefs = npar)

    #Check if model has concave curves
    m_betas <- m_aic$betas
    #Select only quadratic betas
    q_betas <- m_betas[grepl("\\^2", names(m_betas))]
    #Check if is concave
    if(length(q_betas) == 0) {
      is_c <- FALSE} else {
        is_c <- any(q_betas > 0)
      }
    #Get background index (again, if necessary)
    bgind <- which(data$calibration_data == 0)
    #Fit model using k-folds
    mods <- try(lapply(1:length(data$kfolds), function(i) {
      notrain <- -data$kfolds[[i]]
      data_i <- data$calibration_data[notrain,]
      #Set weights per k-fold
      if(!is.null(weights)) {
        weights_i <- weights[notrain]
      } else {
        weights_i <- weights
      }
      #Run model
      mod_i <- glmnet_mx(p = data_i$pr_bg, data = data_i,
                         f = formula_x, regmult = reg_x,
                         addsamplestobackground = addsamplestobackground,
                         weights = weights_i,
                         calculate_AIC = FALSE)

      #Predict model only to background
      pred_i <- as.numeric(predict.glmnet_mx(object = mod_i,
                                   newdata = data$calibration_data,
                                   clamp = FALSE, type = "cloglog"))

      #Extract suitability in train and test points
      suit_val_cal <- pred_i[unique(c(notrain, -bgind))]
      suit_val_eval <- pred_i[which(!-notrain %in% bgind)]

      ####Calculate omission rate following kuenm####
      #This function is in evaluation_functions.R
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
                            is_concave = is_c,
                            row.names = NULL)
      return(cbind(grid_x, df_eval))
    }), silent = TRUE)
    }

  #####Create empty dataframe if mods is an error####
  if(class(mods) == "try-error") {
    eval_final <- cbind(grid_x,
                        empty_replicates(omrat_thr = omrat_thr,
                                         n_row = length(data$kfolds),
                                         replicates = names(data$kfolds),
                                         is_c = is_c))
  } else{
    #Return evaluation final
    names(mods) <- names(data$kfolds)
    eval_final <- do.call("rbind", mods) }

  #Summarize results using eval_stats function
  if(class(mods) == "try-error") {
    eval_final_summary <- cbind(grid_x,
                                empty_summary(omrat_thr = omrat_thr,
                                              is_c = is_c))
  } else {
    eval_final_summary <- eval_stats(eval_final, omrat_thr) }


  #If write_summary = T...
  if(write_summary){
    write.csv(eval_final_summary,
              file.path(out_dir,
                        paste0("Summary_cand_model_", grid_x$ID, ".csv")),
              row.names = F)
  }


  if(!return_replicate)
    eval_final <- NULL

  return(list(All_results = eval_final,
              Summary = eval_final_summary))
}

#Fit best models
fit_best_model <- function(x, dfgrid, calibration_results, n_replicates,
                           rep_data){

  #Get grid
  grid_x <- dfgrid[x,]
  m_id <- grid_x$models
  rep_x <- grid_x$replicates

  #Get best model
  best_models_i <- calibration_results$selected_models[which(
    calibration_results$selected_models$ID == m_id),]
  #best_models_i <- selected_models[i,]
  best_formula <- best_models_i$Formulas
  best_regm <- best_models_i$regm

  #Get replicate, if necessary
  if(n_replicates > 1){
    rep_i <- rep_data[[rep_x]]
    data_x <- calibration_results$calibration_data[rep_i, ]
  } else { #Select i k-fold
    data_x <- calibration_results$calibration_data }
  #Run model
  mod_x <- glmnet_mx(p = data_x[,"pr_bg"], data = data_x,
                     f = as.formula(best_formula),
                     regmult = best_regm,
                     addsamplestobackground = calibration_results[["addsamplestobackground"]],
                     weights = calibration_results[["weights"]],
                     calculate_AIC = FALSE)

  #Only to make sure the IDs in list are correct
  mod_x$checkModel <- m_id
  mod_x$checkreplicate <- rep_x
  # #Predict
  # pred_x <- terra::predict(spat_var, mod_x, type = "cloglog",
  #                          na.rm = TRUE)
  return(mod_x)
}

#Bind rows to get path for each projection in project_selected_glmnetx
bind_rows_projection <- function(data_frames) {
  all_columns <- c("Time", "Period", "Scenario", "ssp", "GCM",
                   "input_path", "output_path")
  result <- NULL

  for (df in data_frames) {
    if (exists("df") && !is.null(df)) {  # Check if df exists and is not NULL
      missing_columns <- setdiff(all_columns, colnames(df))
      for (col in missing_columns) {
        df[[col]] <- NA
      }
      result <- rbind(result, df[, all_columns]) # Reorder columns as specified
    }
  }

  # Remove columns with only NA values
  result <- result[, colSums(is.na(result)) < nrow(result)]

  return(result)
}

####Get formulas####
prepare_formulas_glmnetmx <- function(independent, type = "lqpht",
                                      categorical_var = NULL,
                                      minvar = 1, maxvar = NULL)
{
  if (is.character(type)) {
    if (!all(unlist(strsplit(type, "")) %in% c("l", "p",
                                               "q", "h", "t"))) {
      stop("'type' must be: 'l', 'p', 'q', 'h', 't', or a combination of those")
    }
  }  else {
    stop("'type' must be of class character.")
  }

  #Check if categorical variable is in independente variables. If not, set null
  if(!is.null(categorical_var)) {
    if(!(categorical_var %in% independent)) {
      categorical_var <- NULL
    }  }

  if(!is.null(categorical_var)) {
    independent <- setdiff(independent, categorical_var)
  }

  predictors <- independent
  npred <- length(predictors)
  aux <- " "
  #Linear
  if (grepl("l", type)) {
    aux <- paste(aux, paste(predictors, collapse = " + "),
                 sep = " + ")
  }
  #Quadratic
  if (grepl("q", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("I(", predictors[v], "^2)"),
                   sep = " + ")
    }
  }
  #Product
  if (grepl("p", type)) {
    if (npred > 1) {
      inter_tab <- utils::combn(predictors, 2)
      aux_inter <- paste0(" + ", paste(apply(inter_tab,
                                             2, paste, collapse = ":"), collapse = " + "))
      if (length(aux_inter) > 0) {
        aux <- paste0(aux, aux_inter)
      }
    }
    else {
      if (grepl("l", type) | grepl("q", type)) {
        message("'p' is is only possible with 2 or more independent variables.",
                "\nReturning other combinations.")
      }
      else {
        stop("'p' is is only possible with 2 or more independent variables.",
             "\nTry other combinations of type.")
      }
    }
  }
  #Hinge
  if (grepl("h", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("hinge(", predictors[v], ")"),
                   sep = " + ")
    }
  }
  #Threshold
  if (grepl("t", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("thresholds(", predictors[v], ")"),
                   sep = " + ")
    }
  }
  out <- paste0("~",
                gsub(pattern = "  \\+ ", x = aux, replacement = ""))

  if(!is.null(categorical_var)) {
    out <- paste0(out, " + categorical(", categorical_var, ")")
  }
  return(out)
}
