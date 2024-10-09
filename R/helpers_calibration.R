#' Summarize evaluation results
#'
#' @importFrom stats aggregate glm as.formula
#' @importFrom enmpa proc_enm
#' @export

eval_stats <- function(calib_results, omrat_thr, model_type) {

  # Define omission rates and columns to aggregate
  omission_rates <- paste0("Omission_rate_at_", omrat_thr)
  toagg <- c(omission_rates, "proc_auc_ratio", "proc_pval")

  # Aggregation groups depending on the model type
  if (model_type == "glmnet") {
    agg_by <- c("Formulas", "regm", "Features")
    to_keep <- c("ID", "Formulas", "regm", "Features", "AIC_nk", "AIC_ws",
                 "npar", "is_concave")

  } else if (model_type == "glm") {
    agg_by <- c("Formulas", "Features")
    to_keep <- c("ID", "Formulas", "Features", "AIC", "npar", "is_concave")

  } else {
    stop("Unsupported model type. Please use 'glmnet' or 'glm'.")
  }

  # Aggregation formula
  agg_formula <- paste("~", paste(agg_by, collapse = " + "))

  # Get summary statistics
  xy <- lapply(toagg, function(x) {
    do.call(
      data.frame,
      stats::aggregate(as.formula(paste(x, agg_formula)),
                       data = calib_results, FUN = function(y) {
                         c(mean = round(mean(y), 4), sd = round(sd(y), 4))
                       }, na.action = NULL)
    )
  })

  # Summarize statistics and merge
  stats <- Reduce(function(x, y) merge(x, y, by = agg_by), xy)

  # Keep AIC and related information
  stats_AICS <- calib_results[!duplicated(calib_results[, to_keep]), ][, to_keep]
  stats_final <- merge(stats, stats_AICS, by = agg_by)

  return(stats_final)
}

empty_replicates <- function(omrat_thr, n_row = 4, replicates = 1:4,
                             is_c = NA, model_type) {

  # Define column names based on model type
  if (model_type == "glmnet") {
    column_names <- c("Replicate", paste0("Omission_rate_at_", omrat_thr),
                      "proc_auc_ratio", "proc_pval", "AIC_nk", "AIC_ws", "npar",
                      "is_concave")
  } else if (model_type == "glm") {
    column_names <- c("Replicate", paste0("Omission_rate_at_", omrat_thr),
                      "proc_auc_ratio", "proc_pval", "AIC", "npar", "is_concave")
  } else {
    stop("Unsupported model type. Please use 'glmnet' or 'glm'.")
  }

  # Create an empty dataframe
  df_eval_q <- data.frame(matrix(NA, nrow = n_row, ncol = length(column_names)))
  colnames(df_eval_q) <- column_names

  # Assign replicate values and concavity status
  df_eval_q$Replicate <- replicates
  df_eval_q$is_concave <- is_c

  return(df_eval_q)
}

empty_summary <- function(omrat_thr, is_c, model_type) {

  # Omission rates column names
  om_means <- paste0("Omission_rate_at_", omrat_thr, ".mean")
  om_sd <- paste0("Omission_rate_at_", omrat_thr, ".sd")

  # Base column names depending on the model type
  if (model_type == "glmnet") {
    column_names <- c(om_means, om_sd,
                      "proc_auc_ratio.mean", "proc_auc_ratio.sd",
                      "proc_pval.mean", "proc_pval.sd",
                      "AIC_nk", "AIC_ws", "npar", "is_concave")
  } else if (model_type == "glm") {
    column_names <- c(om_means, om_sd,
                      "proc_auc_ratio.mean", "proc_auc_ratio.sd",
                      "proc_pval.mean", "proc_pval.sd",
                      "AIC", "npar", "is_concave")
  } else {
    stop("Unsupported model type. Please use 'glmnet' or 'glm'.")
  }

  # Create the empty dataframe
  eval_final_q <- data.frame(matrix(NA, nrow = 1, ncol = length(column_names)))
  colnames(eval_final_q) <- column_names
  eval_final_q$is_concave <- is_c

  return(eval_final_q)
}

fit_eval_concave <- function(x, q_grids, data, formula_grid, omrat_thr,
                             write_summary, addsamplestobackground, weights,
                             return_replicate, model_type,
                             AIC = NULL) {

  grid_x <- q_grids[x, ]

  if (model_type == "glmnet") {
    # For glmnet model
    formula_x <- as.formula(grid_x$Formulas)
    reg_x <- grid_x$regm
    m_aic <- try(glmnet_mx(p = data$calibration_data$pr_bg,
                           data = data$calibration_data,
                           f = formula_x, regmult = reg_x,
                           addsamplestobackground = addsamplestobackground,
                           weights = weights,
                           calculate_AIC = TRUE,
                           AIC = AIC),
                 silent = TRUE)

  } else if (model_type == "glm") {
    # For glm model
    formula_x <- as.formula(paste("pr_bg ", grid_x$Formulas))
    m_aic <- glm_mx(formula = formula_x, family = binomial(link = "cloglog"),
                    data = data$calibration_data, weights = weights)

  } else {
    stop("Unsupported model type. Please use 'glmnet' or 'glm'.")
  }

  if (any(class(m_aic) == "try-error")) {
    npar <- NA
    AICc <- NA
    is_c <- NA
    mods <- NA
    class(mods) <- "try-error"
  } else {
    # Get number of parameters
    npar <- if (model_type == "glmnet") {
      length(m_aic$betas)
    } else{
      length(m_aic$coefficients[-1])
    }

    if (model_type == "glmnet") {
      vals <- predict.glmnet_mx(object = m_aic,
                                newdata = data$calibration_data[
                                  data$calibration_data$pr_bg == 1, ],
                                type = "exponential")
      AICc <- aic_ws(pred_occs = vals, ncoefs = npar)
    }

    # Check for concave curves (quadratic terms)
    q_betas <- if (model_type == "glmnet") {
      m_aic$betas[grepl("\\^2", names(m_aic$betas))]
    } else {
      m_aic$coefficients[-1][grepl("\\^2", names(m_aic$coefficients[-1]))]
    }

    is_c <- if (length(q_betas) == 0) FALSE else any(q_betas > 0)
  }

  # Handle concave results
  if (isTRUE(is_c) | is.na(is_c)) {
    # If concave, return grid
    grid_q <- if (model_type == "glmnet") {
      all_reg <- unique(formula_grid$regm)
      do.call("rbind", lapply(seq_along(all_reg), function(k) {
        grid_x_i <- grid_x
        grid_x_i$regm <- all_reg[k]
        grid_x_i$ID <- formula_grid[formula_grid$Formulas == grid_x$Formulas &
                                      formula_grid$regm == all_reg[k], "ID"]
        return(grid_x_i)
      }))
    } else {
      formula_grid[x, ]
    }

    df_eval_q <- empty_replicates(omrat_thr = omrat_thr,
                                  n_row = nrow(grid_q) * length(data$kfolds),
                                  replicates = names(data$kfolds),
                                  is_c = is_c,
                                  model_type = model_type)
    df_eval_q2 <- cbind(grid_q, df_eval_q)
    eval_final_q <- empty_summary(omrat_thr = omrat_thr, is_c = is_c,
                                  model_type = model_type)
    eval_final_q_summary <- cbind(grid_q, eval_final_q)

  } else {
    # If not concave, calculate metrics
    bgind <- which(data$calibration_data == 0)
    mods <- lapply(1:length(data$kfolds), function(i) {
      notrain <- -data$kfolds[[i]]
      data_i <- data$calibration_data[notrain, ]

      if (!is.null(data$weights)){
        weights_i <- data$weights[notrain]
      } else {
        weights_i <- NULL
      }

      if (model_type == "glmnet") {
        mod_i <- glmnet_mx(p = data_i$pr_bg, data = data_i,
                           f = formula_x, regmult = reg_x,
                           addsamplestobackground = addsamplestobackground,
                           weights = weights_i, calculate_AIC = FALSE)
      } else if (model_type == "glm") {
        mod_i <- glm_mx(formula = formula_x, family = binomial(link = "cloglog"),
                        data = data_i, weights = weights_i)
      }

      pred_i <- if (model_type == "glmnet") {
        predict.glmnet_mx(object = mod_i, newdata = data$calibration_data,
                          clamp = FALSE, type = "cloglog")
      } else if (model_type == "glm") {
        enmpa::predict_glm(model = mod_i, newdata = data$calibration_data,
                           type = "response")
      }

      # Calculate metrics (omission rate, pROC)
      suit_val_cal <- pred_i[unique(c(notrain, -bgind))]
      suit_val_eval <- pred_i[which(!-notrain %in% bgind)]
      om_rate <- omrat(threshold = omrat_thr, pred_train = suit_val_cal,
                       pred_test = suit_val_eval)
      proc_i <- enmpa::proc_enm(test_prediction = suit_val_eval,
                                prediction = pred_i)

      df_eval_q <-  if (model_type == "glmnet") {
        data.frame(Replicate = i,
                   t(data.frame(om_rate)),
                   proc_auc_ratio = proc_i$pROC_summary[1],
                   proc_pval = proc_i$pROC_summary[2],
                   AIC_nk = m_aic$AIC,
                   AIC_ws = AICc,
                   npar = npar,
                   is_concave = is_c,
                   row.names = NULL)
      } else {
        data.frame(Replicate = i,
                   t(data.frame(om_rate)),
                   proc_auc_ratio = proc_i$pROC_summary[1],
                   proc_pval = proc_i$pROC_summary[2],
                   AIC = m_aic$aic,
                   npar = npar,
                   is_concave = is_c,
                   row.names = NULL)
      }
      return(cbind(grid_x, df_eval_q))
    })
    names(mods) <- names(data$kfolds)
    eval_final_q <- do.call("rbind", mods)
    eval_final_q_summary <- eval_stats(eval_final_q, omrat_thr, model_type)
  }

  # Write summary if requested
  if (write_summary) {
    write.csv(eval_final_q_summary,
              file.path(out_dir, paste0("Summary_cand_model_", grid_x$ID, ".csv")),
              row.names = FALSE)
  }

  # Return final results
  if (!return_replicate)
    eval_final_q <- NULL

  return(list(All_results = eval_final_q, Summary = eval_final_q_summary))
}

fit_eval_models <- function(x, formula_grid, data, omrat_thr, write_summary,
                            addsamplestobackground, weights, return_replicate,
                            model_type, AIC) {

  # Arguments:
  # x: Each line of the formula grid
  # formula_grid: Formula grid (output of calibration_grid_glmnetmx)
  # data: Data to fit models (output of prepare_data)
  # omrat_thr: Omission rate threshold
  # write_summary: Write individual summary for each line in formula grid?
  # addsamplestobackground: Add samples to background?
  # weights: Add weights
  # return_replicate: Return results of replicates
  # model_type: Type of model, either glmnet or glm
  # AIC: AIC for glmnet (can be ws or nk)

  grid_x <- formula_grid[x,] # Get i candidate model

  if (model_type == "glmnet") {
    # Fit glmnet model
    reg_x <- grid_x$regm # Get regularization multiplier for glmnet
    formula_x <- as.formula(grid_x$Formulas) # Get formula from grid x

    m_aic <- try(glmnet_mx(p = data$calibration_data$pr_bg,
                           data = data$calibration_data,
                           f = formula_x, regmult = reg_x,
                           addsamplestobackground = addsamplestobackground,
                           weights = weights, calculate_AIC = TRUE, AIC = AIC
    ),
    silent = TRUE)
  } else if (model_type == "glm") {
    # Fit glm model
    formula_x <- as.formula(paste("pr_bg", grid_x$Formulas))
    m_aic <- glm_mx(formula = formula_x,
                    family = binomial(link = "cloglog"),
                    data = data$calibration_data,
                    weights = weights)
  } else {
    stop("Unsupported model type. Please use 'glmnet' or 'glm'.")
  }

  # Handle errors during model fitting
  if (any(class(m_aic) == "try-error")) {
    npar <- NA
    AICc <- NA
    is_c <- NA
    mods <- NA
    class(mods) <- "try-error"
  } else {
    # Get number of parameters
    npar <- if (model_type == "glmnet") {
      length(m_aic$betas)
    } else{
      length(m_aic$coefficients[-1])
    }

    # Calculate AIC for glmnet or glm
    AICc <- if (model_type == "glmnet") {
      vals <- predict.glmnet_mx(object = m_aic,
                                newdata = data$calibration_data[
                                  data$calibration_data$pr_bg == 1, ],
                                type = "exponential")
      aic_ws(pred_occs = vals, ncoefs = npar)

    }

    # Check for concave curves (quadratic terms)
    q_betas <- if (model_type == "glmnet") {
      m_aic$betas[grepl("\\^2", names(m_aic$betas))]
    } else {
      m_aic$coefficients[-1][grepl("\\^2", names(m_aic$coefficients[-1]))]
    }

    is_c <- if (length(q_betas) == 0) FALSE else any(q_betas > 0)

    # Get background index
    bgind <- which(data$calibration_data == 0)

    # Fit models using k-fold cross-validation
    mods <- try(lapply(1:length(data$kfolds), function(i) {
      notrain <- -data$kfolds[[i]]
      data_i <- data$calibration_data[notrain,]

      # Set weights per k-fold
      if (!is.null(weights)){
        weights_i <- weights[notrain]
      } else {
        weights_i <- NULL
      }

      if (model_type == "glmnet") {
        # Run glmnet model
        cat(weights_i)
        mod_i <- glmnet_mx(p = data_i$pr_bg, data = data_i,
                           f = formula_x, regmult = reg_x,
                           addsamplestobackground = addsamplestobackground,
                           weights = weights_i, calculate_AIC = FALSE)
      } else {
        # Run glm model
        mod_i <- glm_mx(formula = formula_x,
                        family = binomial(link = "cloglog"),
                        data = data_i, weights = weights_i)
      }

      # Predict model
      pred_i <- if (model_type == "glmnet") {
        as.numeric(predict.glmnet_mx(object = mod_i,
                                     newdata = data$calibration_data,
                                     clamp = FALSE, type = "cloglog"))
      } else if (model_type == "glm") {
        enmpa::predict_glm(model = mod_i,
                           newdata = data$calibration_data,
                           type = "response")
      }

      # Extract suitability in train and test points
      suit_val_cal <- pred_i[unique(c(notrain, -bgind))]
      suit_val_eval <- pred_i[which(!-notrain %in% bgind)]

      # Calculate omission rate and pROC
      om_rate <- omrat(threshold = omrat_thr, pred_train = suit_val_cal,
                       pred_test = suit_val_eval)
      proc_i <- enmpa::proc_enm(test_prediction = suit_val_eval,
                                prediction = pred_i)

      # Save metrics in a dataframe
      df_eval <-  if (model_type == "glmnet") {
        data.frame(Replicate = i,
                   t(data.frame(om_rate)),
                   proc_auc_ratio = proc_i$pROC_summary[1],
                   proc_pval = proc_i$pROC_summary[2],
                   AIC_nk = m_aic$AIC,
                   AIC_ws = AICc,
                   npar = npar,
                   is_concave = is_c,
                   row.names = NULL)
      } else if (model_type == "glm") {
        data.frame(Replicate = i,
                   t(data.frame(om_rate)),
                   proc_auc_ratio = proc_i$pROC_summary[1],
                   proc_pval = proc_i$pROC_summary[2],
                   AIC = m_aic$aic,
                   npar = npar,
                   is_concave = is_c,
                   row.names = NULL)
      }
      return(cbind(grid_x, df_eval))
    }), silent = TRUE)
  }

  ##### Handle errors and summarize results #####
  if (class(mods) == "try-error") {
    eval_final <- cbind(grid_x,
                        empty_replicates(omrat_thr = omrat_thr,
                                         n_row = length(data$kfolds),
                                         replicates = names(data$kfolds),
                                         is_c = is_c, model_type = model_type))
  } else {
    # Combine evaluation results
    names(mods) <- names(data$kfolds)
    eval_final <- do.call("rbind", mods)
  }

  # Summarize results using eval_stats
  eval_final_summary <- if (class(mods) == "try-error") {
    cbind(grid_x, empty_summary(omrat_thr = omrat_thr, is_c = is_c,
                                model_type = model_type))
  } else {
    eval_stats(eval_final, omrat_thr, model_type = model_type)
  }

  # Write summary if requested
  if (write_summary) {
    write.csv(eval_final_summary,
              file.path(out_dir, paste0("Summary_cand_model_", grid_x$ID, ".csv")),
              row.names = FALSE)
  }

  # Return replicates?
  if (!return_replicate) eval_final <- NULL
  return(list(All_results = eval_final, Summary = eval_final_summary))
}

#Fit best models (pending to modify)
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


#' Print Method for prepare_data Class
#' @export
print.calibration_results <- function(x, ...){
  cat(paste0("calibration_results object summary (", x$model_type,")\n"))
  cat("=============================================================\n")
  cat("Species:", x$species, "\n")

  cat("Number of candidate models:", nrow(x$calibration_results$Summary), "\n")

  #Print summary
  cat("  - Models removed because they failed to fit:", length(x$summary$Errors), "\n")
  cat("  - Models removed with concave curves:", length(x$summary$Concave), "\n")
  cat("  - Models removed with non-significant values of pROC:",
      length(x$summary$Non_sig_pROC), "\n")
  cat("  - Models removed with omission rate >", paste0(x$omission_rate, "%:"),
      length(x$summary$High_omission_rate), "\n")
  cat("  - Models removed with delta AIC >", paste0(x$summary$delta_AIC, ":"),
      length(x$summary$High_AIC), "\n")
  cat("Selected models:", nrow(x$selected_models), "\n")

  cat("  - Print selected models (n = 5):\n")
  print(head(x$selected_models))
}
