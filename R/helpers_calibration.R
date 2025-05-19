
eval_stats <- function(cal_res, error_considered, algorithm) {

  # Arguments:
  # cal_res: Calibration results
  # error_considered: Omission rate threshold
  # algorithm: Type of model, either maxnet or glm

  # Check arguments:
  if (any(missing(cal_res), missing(error_considered), missing(algorithm))) {
    stop("Arguments 'cal_res', 'error_considered', 'algorithm' must be defined.")
  }
  if (!inherits(cal_res, "data.frame")) {
    stop("Argument 'cal_res' must be a 'dataframe'.")
  }

  if (!inherits(error_considered, "numeric")) {
    stop("Argument 'error_considered' must be 'numeric'.")
  }

  if (!inherits(algorithm, "character")) {
    stop("Argument 'algorithm' must be a 'character'.")
  }

  if (!(algorithm %in% c("glm", "maxnet"))) {
    stop("Argument algorithm must be 'glm' or 'maxnet'")
  }
  ####

  # Define omission rates and proc to aggregate
  omission_rates <- paste0("Omission_rate_at_", error_considered)
  proc_values <- paste0("Mean_AUC_ratio_at_", error_considered)
  pval_values <- paste0("pval_pROC_at_", error_considered)

  toagg <- c(omission_rates, proc_values, pval_values)

  # Aggregation groups depending on the model type
  if (algorithm == "maxnet") {
    agg_by <- c("Formulas", "R_multiplier", "Features")
    to_keep <- c("ID", "Formulas", "R_multiplier", "Features", "AICc",
                 "Parameters", "Is_concave")

  } else if (algorithm == "glm") {
    agg_by <- c("Formulas", "Features")
    to_keep <- c("ID", "Formulas", "Features", "AIC", "Parameters", "Is_concave")

  } else {
    stop("Unsupported model type. Please use 'maxnet' or 'glm'.")
  }

  # Aggregation formula
  agg_formula <- paste("~", paste(agg_by, collapse = " + "))

  # Get summary statistics
  xy <- lapply(toagg, function(x) {
    do.call(
      data.frame,
      stats::aggregate(as.formula(paste(x, agg_formula)),
                       data = cal_res, FUN = function(y) {
                         c(mean = round(mean(y), 4), sd = round(sd(y), 4))
                       }, na.action = NULL)
    )
  })

  # Summarize statistics and merge
  stats <- Reduce(function(x, y) merge(x, y, by = agg_by), xy)

  # Keep AIC and related information
  stats_AICS <- cal_res[!duplicated(cal_res[, to_keep]), ][, to_keep]
  stats_final <- merge(stats, stats_AICS, by = agg_by)

  #colnames(stats_final) <- gsub("\\.", "_", colnames(stats_final))

  return(stats_final)
}



empty_replicates <- function(error_considered, n_row, replicates,
                             is_c, algorithm) {

  # Arguments:
  # omission_rate: Omission rate threshold
  # n_row: Number of rows
  # replicates: Replicates
  # is_c: Concavity status
  # algorithm: Type of model, either maxnet or glm

  # Check arguments:
  if (any(missing(error_considered), missing(n_row), missing(replicates),
          missing(is_c), missing(algorithm))) {
    stop("Arguments 'error_considered', 'n_row', 'replicates', 'is_c', and 'algorithm' must be defined.")
  }
  if (!inherits(error_considered, "numeric")) {
    stop("Argument 'error_considered' must be 'numeric'.")
  }
  if (!is.numeric(n_row)) {
    stop("Argument 'n_row' must be 'numeric'.")
  }
  if (all(!inherits(replicates, "numeric") & !inherits(replicates, "character"))) {
    stop("Argument 'replicates' must be 'numeric' or 'character'.")
  }

  if (!is.na(is_c)) {
    if (!inherits(is_c, "logical")) {
      stop("Argument 'is_c' must be NA or 'logical'.")
    }
  }

  if (!inherits(algorithm, "character")) {
    stop("Argument 'algorithm' must be a 'character'.")
  }

  if (!(algorithm %in% c("glm", "maxnet"))) {
    stop("Argument algorithm must be 'glm' or 'maxnet'")
  }
  ####

  # Define column names based on model type
  if (algorithm == "maxnet") {
    column_names <- c("Fold",
                      paste0("Omission_rate_at_", error_considered),
                      paste0("Mean_AUC_ratio_at_", error_considered),
                      paste0("pval_pROC_at_", error_considered),
                      "AICc", "Parameters", "Is_concave")
  } else if (algorithm == "glm") {
    column_names <- c("Fold",
                      paste0("Omission_rate_at_", error_considered),
                      paste0("Mean_AUC_ratio_at_", error_considered),
                      paste0("pval_pROC_at_", error_considered),
                      "AIC", "Parameters", "Is_concave")
  } else {
    stop("Unsupported model type. Please use 'maxnet' or 'glm'.")
  }

  # Create an empty dataframe
  df_eval_q <- data.frame(matrix(NA, nrow = n_row, ncol = length(column_names)))
  colnames(df_eval_q) <- column_names

  # Assign Fold values and concavity status
  df_eval_q$Fold <- replicates
  df_eval_q$Is_concave <- is_c

  return(df_eval_q)
}



empty_summary <- function(error_considered, is_c, algorithm) {


  # Arguments:
  # omission_rate: Omission rate threshold
  # is_c: Concavity status
  # algorithm: Type of model, either maxnet or glm

  # Check arguments:
  if (any(missing(error_considered), missing(is_c), missing(algorithm))) {
    stop("Arguments 'error_considered', 'is_c', and 'algorithm' must be defined.")
  }
  if (!inherits(error_considered, "numeric")) {
    stop("Argument 'error_considered' must be 'numeric'.")
  }

  if (!is.na(is_c)) {
    if (!inherits(is_c, "logical")) {
      stop("Argument 'is_c' must be NA or 'logical'.")
    }
  }

  if (!inherits(algorithm, "character")) {
    stop("Argument 'algorithm' must be a 'character'.")
  }

  if (!(algorithm %in% c("glm", "maxnet"))) {
    stop("Argument algorithm must be 'glm' or 'maxnet'")
  }
  ####

  # Omission rates column names
  om_means <- paste0("Omission_rate_at_", error_considered, ".mean")
  om_sd <- paste0("Omission_rate_at_", error_considered, ".sd")

  #Proc columns names
  auc_means <- paste0("Mean_AUC_ratio_at_", error_considered, ".mean")
  auc_sd <- paste0("Mean_AUC_ratio_at_", error_considered, ".sd")
  pval_means <- paste0("pval_pROC_at_", error_considered, ".mean")
  pval_sd <- paste0("pval_pROC_at_", error_considered, ".sd")

  # Base column names depending on the model type
  if (algorithm == "maxnet") {
    column_names <- c(om_means, om_sd,
                      auc_means, auc_sd,
                      pval_means, pval_sd,
                      "AICc", "Parameters", "Is_concave")
  } else if (algorithm == "glm") {
    column_names <- c(om_means, om_sd,
                      auc_means, auc_sd,
                      pval_means, pval_sd,
                      "AIC", "Parameters", "Is_concave")
  } else {
    stop("Unsupported model type. Please use 'maxnet' or 'glm'.")
  }

  # Create the empty dataframe
  eval_final_q <- data.frame(matrix(NA, nrow = 1, ncol = length(column_names)))
  colnames(eval_final_q) <- column_names
  eval_final_q$Is_concave <- is_c

  return(eval_final_q)
}



fit_eval_concave <- function(x, q_grids, data, formula_grid, error_considered, omission_rate,
                             write_summary, addsamplestobackground, weights = NULL,
                             return_all_results, algorithm,
                             proc_for_all, out_dir) {

  # Arguments:
  # x: Each line of the formula grid
  # q_grids: Formula grid for quadratic terms
  # data: Data to fit models (output of prepare_data)
  # formula_grid: Formula grid
  # omission_rate: Omission rate threshold
  # write_summary: Write individual summary for each line in formula grid?
  # addsamplestobackground: Add samples to background?
  # weights: Add weights
  # return_all_results: Return results of replicates
  # algorithm: Type of model, either maxnet or glm

  # Check arguments

  if (!inherits(q_grids, "data.frame")) {
    stop("Argument 'q_grids' must be a data.frame.")
  }

  if (!inherits(data, "prepared_data")) {
    stop("Argument 'data' must be a prepared_data object.")
  }

  if (!inherits(formula_grid, "data.frame")) {
    stop("Argument 'formula_grid' must be a data.frame.")
  }

  if (!inherits(error_considered, "numeric")) {
    stop("Argument error_considered must be numeric.")
  }

  if (!inherits(omission_rate, "numeric")) {
    stop("Argument omission_rate must be numeric.")
  }

  if (!inherits(write_summary, "logical")) {
    stop("Argument write_summary must be logical.")
  }

  if (!inherits(addsamplestobackground, "logical")) {
    stop("Argument 'addsamplestobackground' must be 'logical'.")
  }

  if (!is.null(weights)) {
    if (!inherits(weights, "numeric")) {
      stop("Argument 'weights' must be NULL or 'numeric'.")
    }}

  if (!inherits(return_all_results, "logical")) {
    stop("Argument 'return_all_results' must be 'logical'.")
  }

  if (!inherits(algorithm, "character")) {
    stop("Argument 'algorithm' must be a 'character'.")
  }

  if (!(algorithm %in% c("glm", "maxnet"))) {
    stop("Argument 'algorithm' must be 'glm' or 'maxnet'.")
  }

  if (!inherits(proc_for_all, "logical")) {
    stop("Argument 'proc_for_all' must be 'logical'.")
  }

  ####

  grid_x <- q_grids[x, ]

  if (algorithm == "maxnet") {
    # For maxnet model
    formula_x <- as.formula(grid_x$Formulas)
    reg_x <- grid_x$R_multiplier
    m_aic <- try(glmnet_mx(p = data$calibration_data$pr_bg,
                           data = data$calibration_data,
                           f = formula_x, regmult = reg_x,
                           addsamplestobackground = addsamplestobackground,
                           weights = weights),
                 silent = TRUE)

  } else if (algorithm == "glm") {
    # For glm model
    formula_x <- as.formula(paste("pr_bg ", grid_x$Formulas))
    m_aic <- glm_mx(formula = formula_x, family = binomial(link = "cloglog"),
                    data = data$calibration_data, weights = weights)

  } else {
    stop("Unsupported model type. Please use 'maxnet' or 'glm'.")
  }

  if (any(class(m_aic) == "try-error")) {
    npar <- NA
    AICc <- NA
    is_c <- NA
    mods <- NA
    class(mods) <- "try-error"
  } else {
    # Get number of parameters
    npar <- if (algorithm == "maxnet") {
      length(m_aic$betas)
    } else{
      length(m_aic$coefficients[-1])
    }

    if (algorithm == "maxnet") {
      vals <- predict.glmnet_mx(object = m_aic,
                                newdata = data$calibration_data[
                                  data$calibration_data$pr_bg == 1, ],
                                type = "exponential")
      AICc <- aic_ws(pred_occs = vals, ncoefs = npar)
    }

    # Check for concave curves (quadratic terms)
    q_betas <- if (algorithm == "maxnet") {
      m_aic$betas[grepl("\\^2", names(m_aic$betas))]
    } else {
      m_aic$coefficients[-1][grepl("\\^2", names(m_aic$coefficients[-1]))]
    }

    is_c <- if (length(q_betas) == 0) FALSE else any(q_betas > 0)
  }

  # Handle concave results
  if (isTRUE(is_c) | is.na(is_c)) {
    # If concave, return grid
    grid_q <- if (algorithm == "maxnet") {
      all_reg <- unique(formula_grid$R_multiplier)
      do.call("rbind", lapply(seq_along(all_reg), function(k) {
        grid_x_i <- grid_x
        grid_x_i$R_multiplier <- all_reg[k]
        grid_x_i$ID <- formula_grid[formula_grid$Formulas == grid_x$Formulas &
                                      formula_grid$R_multiplier == all_reg[k], "ID"]
        return(grid_x_i)
      }))
    } else {
      formula_grid[x, ]
    }

    df_eval_q <- empty_replicates(error_considered = error_considered,
                                  n_row = nrow(grid_q) * length(data$kfolds),
                                  replicates = names(data$kfolds),
                                  is_c = is_c,
                                  algorithm = algorithm)
    df_eval_q2 <- cbind(grid_q, df_eval_q)
    eval_final_q <- empty_summary(error_considered = error_considered, is_c = is_c,
                                  algorithm = algorithm)
    eval_final_q_summary <- reorder_stats_columns(cbind(grid_q, eval_final_q),
                                                  error_considered)

  } else {
    # If not concave, calculate metrics
    bgind <- which(data$calibration_data == 0)
    mods <- lapply(1:length(data$kfolds), function(i) {
      notrain <- -data$kfolds[[i]]
      data_i <- data$calibration_data[notrain, ]

      if (!is.null(data$weights)) {
        weights_i <- data$weights[notrain]
      } else {
        weights_i <- NULL
      }

      if (algorithm == "maxnet") {
        mod_i <- glmnet_mx(p = data_i$pr_bg, data = data_i,
                           f = formula_x, regmult = reg_x,
                           addsamplestobackground = addsamplestobackground,
                           weights = weights_i, calculate_AIC = FALSE)
      } else if (algorithm == "glm") {
        mod_i <- glm_mx(formula = formula_x, family = binomial(link = "cloglog"),
                        data = data_i, weights = weights_i)
      }

      pred_i <- if (algorithm == "maxnet") {
        as.numeric(predict.glmnet_mx(object = mod_i,
                                     newdata = data$calibration_data,
                                     clamp = FALSE, type = "cloglog"))
      } else if (algorithm == "glm") {
        enmpa::predict_glm(model = mod_i, newdata = data$calibration_data,
                           type = "response")
      }

      # Calculate metrics (omission rate, pROC)
      suit_val_cal <- pred_i[unique(c(notrain, -bgind))]
      suit_val_eval <- pred_i[which(!-notrain %in% bgind)]
      om_rate <- omrat(threshold = error_considered, pred_train = suit_val_cal,
                       pred_test = suit_val_eval)

      #Calculate PROC? ...
      if (proc_for_all) {
        proc_i <- lapply(error_considered, function(omr) {
          proc_omr <- enmpa::proc_enm(test_prediction = suit_val_eval,
                                      prediction = pred_i,
                                      threshold = omr)$pROC_summary
          names(proc_omr) <- c(paste0("Mean_AUC_ratio_at_", omr),
                               paste0("pval_pROC_at_", omr))
          return(proc_omr)
        })
        proc_i <- unlist(proc_i)} else {
          #Or fill PROC with NA
          proc_i <- rep(NA, length(error_considered) * 2)
          names(proc_i) <- c(paste0("Mean_AUC_ratio_at_", error_considered),
                             paste0("pval_pROC_at_", error_considered))
        }



      df_eval_q <-  if (algorithm == "maxnet") {
        data.frame(Fold = i,
                   t(om_rate),
                   t(proc_i),
                   AICc = AICc,
                   Parameters = npar,
                   Is_concave = is_c,
                   row.names = NULL)
      } else {
        data.frame(Fold = i,
                   t(om_rate),
                   t(proc_i),
                   AIC = m_aic$aic,
                   Parameters = npar,
                   Is_concave = is_c,
                   row.names = NULL)
      }
      return(cbind(grid_x, df_eval_q))
    })
    names(mods) <- names(data$kfolds)
    eval_final_q <- do.call("rbind", mods)
    eval_final_q_summary <- reorder_stats_columns(eval_stats(eval_final_q,
                                                             error_considered,
                                                             algorithm),
                                                  error_considered = error_considered)
  }

  # Write summary if requested
  if (write_summary) {
    write.csv(eval_final_q_summary,
              file.path(out_dir, paste0("Summary_cand_model_", grid_x$ID, ".csv")),
              row.names = FALSE)
  }

  # Return final results
  if (!return_all_results)
    eval_final_q <- NULL

  return(list(All_results = eval_final_q, Summary = eval_final_q_summary))
}



fit_eval_models <- function(x, formula_grid, data, error_considered, omission_rate,
                            write_summary, addsamplestobackground, weights = NULL,
                            return_all_results, algorithm,
                            proc_for_all, out_dir) {

  # Arguments:
  # x: Each line of the formula grid
  # formula_grid: Formula grid (output of calibration_grid)
  # data: Data to fit models (output of prepare_data)
  # omission_rate: Omission rate threshold
  # write_summary: Write individual summary for each line in formula grid?
  # addsamplestobackground: Add samples to background?
  # weights: Add weights
  # return_all_results: Return results of replicates
  # algorithm: Type of model, either maxnet or glm

  # Check rguments

  if (!inherits(data, "prepared_data")) {
    stop("Argument 'data' must be a 'prepared_data' object.")
  }

  if (!inherits(formula_grid, "data.frame")) {
    stop("Argument 'formula_grid' must be a 'data.frame'.")
  }

  if (!inherits(error_considered, "numeric")) {
    stop("Argument 'error_considered' must be 'numeric'.")
  }

  if (!inherits(omission_rate, "numeric")) {
    stop("Argument 'omission_rate' must be 'numeric'.")
  }

  if (!inherits(write_summary, "logical")) {
    stop("Argument 'write_summary' must be 'logical'.")
  }

  if (!inherits(addsamplestobackground, "logical")) {
    stop("Argument 'addsamplestobackground' must be 'logical'.")
  }

  if (!is.null(weights)) {
    if (!inherits(weights, "numeric")) {
      stop("Argument 'weights' must be NULL or 'numeric'.")
    }}

  if (!inherits(return_all_results, "logical")) {
    stop("Argument 'return_all_results' must be 'logical'.")
  }

  if (!inherits(algorithm, "character")) {
    stop("Argument 'algorithm' must be a 'character'.")
  }

  if (!(algorithm %in% c("glm", "maxnet"))) {
    stop("Argument 'algorithm' must be 'glm' or 'maxnet'.")
  }

  if (!inherits(proc_for_all, "logical")) {
    stop("Argument 'proc_for_all' must be 'logical'.")
  }
  ####

  grid_x <- formula_grid[x,] # Get i candidate model

  if (algorithm == "maxnet") {
    # Fit maxnet model
    reg_x <- grid_x$R_multiplier # Get regularization multiplier for maxnet
    formula_x <- as.formula(grid_x$Formulas) # Get formula from grid x

    m_aic <- try(glmnet_mx(p = data$calibration_data$pr_bg,
                           data = data$calibration_data,
                           f = formula_x, regmult = reg_x,
                           addsamplestobackground = addsamplestobackground,
                           weights = weights),
    silent = TRUE)
  } else if (algorithm == "glm") {
    # Fit glm model
    formula_x <- as.formula(paste("pr_bg", grid_x$Formulas))
    m_aic <- glm_mx(formula = formula_x,
                    family = binomial(link = "cloglog"),
                    data = data$calibration_data,
                    weights = weights)
  } else {
    stop("Unsupported model type. Please use 'maxnet' or 'glm'.")
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
    npar <- if (algorithm == "maxnet") {
      length(m_aic$betas)
    } else{
      length(m_aic$coefficients[-1])
    }

    # Calculate AIC for maxnet or glm
    AICc <- if (algorithm == "maxnet") {
      vals <- predict.glmnet_mx(object = m_aic,
                                newdata = data$calibration_data[
                                  data$calibration_data$pr_bg == 1, ],
                                type = "exponential")
      aic_ws(pred_occs = vals, ncoefs = npar)

    }

    # Check for concave curves (quadratic terms)
    q_betas <- if (algorithm == "maxnet") {
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
      if (!is.null(weights)) {
        weights_i <- weights[notrain]
      } else {
        weights_i <- NULL
      }

      if (algorithm == "maxnet") {
        # Run maxnet model
        mod_i <- try(glmnet_mx(p = data_i$pr_bg, data = data_i,
                           f = formula_x, regmult = reg_x,
                           addsamplestobackground = addsamplestobackground,
                           weights = weights_i, calculate_AIC = FALSE))
      } else {
        # Run glm model
        mod_i <- try(glm_mx(formula = formula_x,
                        family = binomial(link = "cloglog"),
                        data = data_i, weights = weights_i))
      }

      # Predict model
      pred_i <- if (algorithm == "maxnet" & inherits(mod_i, "glmnet_mx")) {
        as.numeric(predict.glmnet_mx(object = mod_i,
                                     newdata = data$calibration_data, clamp = FALSE,
                                     type = "cloglog"))
      } else if (algorithm == "glm") {
        enmpa::predict_glm(model = mod_i, newdata = data$calibration_data,
                           type = "response")
      } else if (inherits(mod_i, "try-error")){
        rep(NA, nrow(data$calibration_data))
      }

      # Calculate metrics (omission rate, pROC)
      if(inherits(mod_i, "try-error")){
        om_rate <- rep(NA, length(omission_rate))
        names(om_rate) <- paste0("Omission_rate_at_", omission_rate)
      } else {
        suit_val_cal <- pred_i[unique(c(notrain, -bgind))]
        suit_val_eval <- pred_i[which(!-notrain %in% bgind)]
        om_rate <- omrat(threshold = omission_rate, pred_train = suit_val_cal,
                         pred_test = suit_val_eval)}

      #Calculate PROC? ...
      if(proc_for_all & !inherits(mod_i, "try-error")) {
        proc_i <- lapply(omission_rate, function(omr) {
          proc_omr <- enmpa::proc_enm(test_prediction = suit_val_eval,
                                      prediction = pred_i,
                                      threshold = omr)$pROC_summary
          names(proc_omr) <- c(paste0("Mean_AUC_ratio_at_", omr),
                               paste0("pval_pROC_at_", omr))
          return(proc_omr)
        })
        proc_i <- unlist(proc_i)} else {
          #Or fill PROC with NA
          proc_i <- rep(NA, length(error_considered) * 2)
          names(proc_i) <- c(paste0("Mean_AUC_ratio_at_", error_considered),
                             paste0("pval_pROC_at_", error_considered))
        }


      # Save metrics in a dataframe
      df_eval <-  if (algorithm == "maxnet") {
        data.frame(Fold = i,
                   t(om_rate),
                   t(proc_i),
                   AICc = AICc,
                   Parameters = npar,
                   Is_concave = is_c,
                   row.names = NULL)
      } else if (algorithm == "glm") {
        data.frame(Fold = i,
                   t(om_rate),
                   t(proc_i),
                   AIC = m_aic$aic,
                   Parameters = npar,
                   Is_concave = is_c,
                   row.names = NULL)
      }
      return(cbind(grid_x, df_eval))
    }), silent = TRUE)
  }

  ##### Handle errors and summarize results #####
  if (inherits(mods, "try-error")) {
    eval_final <- cbind(grid_x,
                        empty_replicates(error_considered = error_considered,
                                         n_row = length(data$kfolds),
                                         replicates = names(data$kfolds),
                                         is_c = is_c, algorithm = algorithm))
  } else {
    # Combine evaluation results
    names(mods) <- names(data$kfolds)
    eval_final <- do.call("rbind", mods)
  }

  # Summarize results using eval_stats
  eval_final_summary <- if (inherits(mods, "try-error")) {
    reorder_stats_columns(cbind(grid_x, empty_summary(error_considered, is_c,
                                                      algorithm)),
                          error_considered = error_considered)
  } else {
    reorder_stats_columns(eval_stats(eval_final, error_considered, algorithm),
                          error_considered = error_considered)
  }

  # Write summary if requested
  if (write_summary) {
    write.csv(eval_final_summary,
              file.path(out_dir, paste0("Summary_cand_model_", grid_x$ID, ".csv")),
              row.names = FALSE)
  }

  # Return replicates?
  if (!return_all_results) {
    eval_final <- NULL
  }

  return(list(All_results = eval_final, Summary = eval_final_summary))
}


fit_best_model <- function(x, dfgrid, cal_res, n_replicates = 1,
                           rep_data = NULL, algorithm = "maxnet") {

  # Arguments:
  # x: index of the grid
  # dfgrid: dataframe with the grid
  # cal_res: output of the calibration function
  # n_replicates: number of replicates
  # rep_data: data splitting (replicated data)
  # algorithm: Type of model, either "maxnet" or "glm"

  # Check arguments
  if (missing(x)) {
    stop("Argumen 'x' must be defined.")
  }
  if (missing(dfgrid)) {
    stop("Argumen 'dfgrid' must be defined.")
  }
  if (missing(cal_res)) {
    stop("Argumen 'cal_res' must be defined.")
  }

  if (!inherits(dfgrid, "data.frame")) {
    stop("Argument 'dfgrid' must be a 'data.frame.'")
  }

  if (!inherits(cal_res, "calibration_results")) {
    stop("Argument 'cal_res' must be a 'calibration_results' object.")
  }

  if (!inherits(n_replicates, "numeric")) {
    stop("Argument 'n_replicates' must be 'numeric'.")
  }
  if (n_replicates > 1) {
    if (!inherits(rep_data, "list")) {
      stop("Argument 'rep_data' must be a 'list'.")
    }
  }

  if (!inherits(algorithm, "character")) {
    stop("Argument 'algorithm' must be a 'character'.")
  }

  if (!(algorithm %in% c("glm", "maxnet"))) {
    stop("Argument algorithm must be 'glm' or 'maxnet'.")
  }
  ####

  # Get the grid information
  grid_x <- dfgrid[x, ]
  m_id <- grid_x$models
  rep_x <- grid_x$replicates

  # Get the best model's formula and parameters from calibration results
  best_model <- cal_res$selected_models[cal_res$selected_models$ID == m_id, ]
  best_formula <- best_model$Formulas

  if (algorithm == "maxnet") {
    best_regm <- best_model$R_multiplier  # Regularization multiplier for maxnet
  }

  # Select data for the Fold, or use the entire calibration data
  # if n_replicates == 1
  if (n_replicates > 1) {
    rep_i <- rep_data[[rep_x]]
    data_x <- cal_res$calibration_data[rep_i, ]
  } else {
    data_x <- cal_res$calibration_data
  }

  # Fit the model based on the algorithm
  if (algorithm == "maxnet") {
    # Run maxnet model
    mod_x <- glmnet_mx(p = data_x$pr_bg, data = data_x,
                       f = as.formula(best_formula),
                       regmult = best_regm,
                       addsamplestobackground = cal_res$addsamplestobackground,
                       weights = cal_res$weights,
                       calculate_AIC = FALSE)

  } else if (algorithm == "glm") {
    # Run glm model
    mod_x <- glm_mx(formula = as.formula(paste("pr_bg ", best_formula)),
                    family = binomial(link = "cloglog"),
                    data = data_x,
                    weights = cal_res$weights)

    #mod_x$data <- NULL # avoid store redundant info
  }
  # Assign model ID and Fold number for tracking
  mod_x$checkModel <- m_id
  mod_x$checkReplicate <- rep_x

  return(mod_x)
}



# Bind rows to get path for each projection in project_selected_glmnetx
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



# Reorder columns in stats final
reorder_stats_columns <- function(stats_final, error_considered) {

  if (!inherits(stats_final, "data.frame")) {
    stop("Argument 'stats_final' must be a 'data.frame'.")
  }

  if (!inherits(error_considered, "numeric")) {
    stop("Argument 'error_considered' must be 'numeric'.")
  }

  first_cols<- intersect(c("ID", "Formulas", "R_multiplier", "Features"),
                         colnames(stats_final))
  metric_cols <- c(paste0("Omission_rate_at_", error_considered, ".mean"),
                   paste0("Omission_rate_at_", error_considered, ".sd"),
                   #Proc columns names
                   paste0("Mean_AUC_ratio_at_", error_considered, ".mean"),
                   paste0("Mean_AUC_ratio_at_", error_considered, ".sd"),
                   paste0("pval_pROC_at_", error_considered, ".mean"),
                   paste0("pval_pROC_at_", error_considered, ".sd"))
  ordered_metric_cols <- unlist(lapply(error_considered, function(rate) {
    grep(paste0("_", rate, "."), metric_cols, value = TRUE)
  }))
  last_cols <- setdiff(colnames(stats_final), c(first_cols, metric_cols))
  orders_cols <- c(first_cols, ordered_metric_cols, last_cols)
  return(stats_final[,orders_cols])
}

#### PROC ####
proc <- function(x, formula_grid, data, error_considered = 10,
                 addsamplestobackground = TRUE, weights = NULL,
                 algorithm) {

  #Check arguments
  if (missing(x)) {
    stop("Argumen 'x' must be defined.")
  }
  if (missing(formula_grid)) {
    stop("Argumen 'formula_grid' must be defined.")
  }
  if (missing(data)) {
    stop("Argumen 'data' must be defined.")
  }
  if (!inherits(formula_grid, "data.frame")) {
    stop("Argument 'formula_grid' must be a 'data.frame'.")
  }

  if (!inherits(error_considered, "numeric")) {
    stop("Argument 'error_considered' must be 'numeric'.")
  }

  if (!inherits(addsamplestobackground, "logical")) {
    stop("Argument 'addsamplestobackground' must be 'logical'.")
  }

  if (!is.null(weights)) {
    if (!inherits(weights, "numeric")) {
      stop("Argument 'weights' must be NULL or 'numeric'.")
    }}

  if (!inherits(algorithm, "character")) {
    stop("Argument 'algorithm' must be a 'character'.")
  }

  if (!(algorithm %in% c("glm", "maxnet"))) {
    stop("Argument algorithm must be 'glm' or 'maxnet'.")
  }
  ####

  grid_x <- formula_grid[x,] # Get i candidate model
  #Get formula and reg
  formula_x <- as.formula(grid_x$Formulas)
  reg_x <- grid_x$R_multiplier

  if (algorithm == "glm") {
    formula_x <- as.formula(paste("pr_bg ", grid_x$Formulas))
  }

  # Get background index
  bgind <- which(data$calibration_data$pr_bg == 0)

  # Fit models using k-fold cross-validation
  mods <- try(lapply(1:length(data$kfolds), function(i) {
    notrain <- -data$kfolds[[i]]
    data_i <- data$calibration_data[notrain,]

    # Set weights per k-fold
    if (!is.null(weights)) {
      weights_i <- weights[notrain]
    } else {
      weights_i <- NULL
    }

    if (algorithm == "maxnet") {
      # Run glmnet model
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
    pred_i <- if (algorithm == "maxnet") {
      as.numeric(predict.glmnet_mx(object = mod_i,
                                   newdata = data$calibration_data,
                                   clamp = FALSE, type = "cloglog"))
    } else if (algorithm == "glm") {
      enmpa::predict_glm(model = mod_i,
                         newdata = data$calibration_data,
                         type = "response")
    }

    # Extract suitability in train and test points
    suit_val_cal <- pred_i[unique(c(notrain, -bgind))]
    suit_val_eval <- pred_i[which(!-notrain %in% bgind)]

    #Proc
    proc_i <- lapply(error_considered, function(omr) {
      proc_omr <- enmpa::proc_enm(test_prediction = suit_val_eval,
                                  prediction = pred_i,
                                  threshold = omr)$pROC_summary
      names(proc_omr) <- c(paste0("Mean_AUC_ratio_at_", omr),
                           paste0("pval_pROC_at_", omr))
      return(proc_omr)
    })
    proc_i <- unlist(proc_i)

    # Save metrics in a dataframe
    df_proc <-  if (algorithm == "maxnet") {
      data.frame(Fold = i,
                 t(proc_i),
                 row.names = NULL)
    } else if (algorithm == "glm") {
      data.frame(Fold = i,
                 t(proc_i),
                 row.names = NULL)
    }
    return(df_proc)
  }), silent = TRUE)
  #Get summary
  proc_df <- do.call("rbind", mods)

  #Get means and sd
  means <- sapply(proc_df[, -1], mean)  # Excluindo a coluna Fold
  sds <- sapply(proc_df[, -1], sd)

  #Create new dataframe
  proc_df <- data.frame(
    t(c(means, sds))
  )
  #Rename
  names(proc_df) <- c(paste0(names(means), ".mean"), paste0(names(sds), ".sd"))

  #Append model ID
  proc_df$ID <- grid_x$ID

  return(proc_df)
}
