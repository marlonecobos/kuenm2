#' Summarize evaluation results
#'
#' @importFrom stats aggregate glm as.formula
#' @importFrom enmpa proc_enm
#' @export

eval_stats <- function(cal_res, omission_rate, algorithm) {

  # Arguments:
  # cal_res: Calibration results
  # omission_rate: Omission rate threshold
  # algorithm: Type of model, either maxnet or glm

  # Check arguments:
  if(!inherits(cal_res, "data.frame")){
    stop("Argument cal_res must be a dataframe, not ", class(cal_res))
  }

  if(!inherits(omission_rate, "numeric")){
    stop("Argument omission_rate must be numeric, not ", class(omission_rate))
  }

  if(!inherits(algorithm, "character")){
    stop("Argument algorithm must be a character, not ", class(algorithm))
  }

  if(!(algorithm %in% c("glm", "maxnet"))){
    stop("Argument algorithm must be 'glm' or 'maxnet'")
  }
  ####

  # Define omission rates and proc to aggregate
  omission_rates <- paste0("Omission_rate_at_", omission_rate)
  proc_values <- paste0("Mean_AUC_ratio_at_", omission_rate)
  pval_values <- paste0("pval_pROC_at_", omission_rate)

  toagg <- c(omission_rates, proc_values, pval_values)

  # Aggregation groups depending on the model type
  if (algorithm == "maxnet") {
    agg_by <- c("Formulas", "reg_mult", "Features")
    to_keep <- c("ID", "Formulas", "reg_mult", "Features", "AIC_nk", "AIC_ws",
                 "npar", "is_concave")

  } else if (algorithm == "glm") {
    agg_by <- c("Formulas", "Features")
    to_keep <- c("ID", "Formulas", "Features", "AIC", "npar", "is_concave")

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



empty_replicates <- function(omission_rate, n_row, replicates,
                             is_c, algorithm) {

  # Arguments:
  # omrat_thr: Omission rate threshold
  # n_row: Number of rows
  # replicates: Replicates
  # is_c: Concavity status
  # algorithm: Type of model, either maxnet or glm

  # Check arguments:
  if(!inherits(omission_rate, "numeric")){
    stop("Argument omission_rate must be numeric, not ", class(omission_rate))
  }

  if(all(!inherits(replicates, "numeric") & !inherits(replicates, "character"))){
    stop("Argument replicates must be numeric or character, not ", class(replicates))
  }

  if(!is.na(is_c)){
    if(!inherits(is_c, "logical")){
    stop("Argument is_c must be NA or logical, not ", class(is_c))}
  }

  if(!inherits(algorithm, "character")){
    stop("Argument algorithm must be a character, not ", class(algorithm))
  }

  if(!(algorithm %in% c("glm", "maxnet"))){
    stop("Argument algorithm must be 'glm' or 'maxnet'")
  }
  ####

  # Define column names based on model type
  if (algorithm == "maxnet") {
    column_names <- c("Replicate",
                      paste0("Omission_rate_at_", omission_rate),
                      paste0("Mean_AUC_ratio_at_", omission_rate),
                      paste0("pval_pROC_at_", omission_rate),
                      "AIC_nk", "AIC_ws", "npar",
                      "is_concave")
  } else if (algorithm == "glm") {
    column_names <- c("Replicate",
                      paste0("Omission_rate_at_", omission_rate),
                      paste0("Mean_AUC_ratio_at_", omission_rate),
                      paste0("pval_pROC_at_", omission_rate),
                      "AIC", "npar", "is_concave")
  } else {
    stop("Unsupported model type. Please use 'maxnet' or 'glm'.")
  }

  # Create an empty dataframe
  df_eval_q <- data.frame(matrix(NA, nrow = n_row, ncol = length(column_names)))
  colnames(df_eval_q) <- column_names

  # Assign replicate values and concavity status
  df_eval_q$Replicate <- replicates
  df_eval_q$is_concave <- is_c

  return(df_eval_q)
}



empty_summary <- function(omission_rate, is_c, algorithm) {

  # Arguments:
  # omrat_thr: Omission rate threshold
  # is_c: Concavity status
  # algorithm: Type of model, either maxnet or glm

  # Check arguments:
  if(!inherits(omission_rate, "numeric")){
    stop("Argument omission_rate must be numeric, not ", class(omission_rate))
  }

  if(!is.na(is_c)){
    if(!inherits(is_c, "logical")){
      stop("Argument is_c must be NA or logical, not ", class(is_c))}
  }

  if(!inherits(algorithm, "character")){
    stop("Argument algorithm must be a character, not ", class(algorithm))
  }

  if(!(algorithm %in% c("glm", "maxnet"))){
    stop("Argument algorithm must be 'glm' or 'maxnet'")
  }
  ####

  # Omission rates column names
  om_means <- paste0("Omission_rate_at_", omission_rate, ".mean")
  om_sd <- paste0("Omission_rate_at_", omission_rate, ".sd")

  #Proc columns names
  auc_means <- paste0("Mean_AUC_ratio_at_", omission_rate, ".mean")
  auc_sd <- paste0("Mean_AUC_ratio_at_", omission_rate, ".sd")
  pval_means <- paste0("pval_pROC_at_", omission_rate, ".mean")
  pval_sd <- paste0("pval_pROC_at_", omission_rate, ".sd")

  # Base column names depending on the model type
  if (algorithm == "maxnet") {
    column_names <- c(om_means, om_sd,
                      auc_means, auc_sd,
                      pval_means, pval_sd,
                      "AIC_nk", "AIC_ws", "npar", "is_concave")
  } else if (algorithm == "glm") {
    column_names <- c(om_means, om_sd,
                      auc_means, auc_sd,
                      pval_means, pval_sd,
                      "AIC", "npar", "is_concave")
  } else {
    stop("Unsupported model type. Please use 'maxnet' or 'glm'.")
  }

  # Create the empty dataframe
  eval_final_q <- data.frame(matrix(NA, nrow = 1, ncol = length(column_names)))
  colnames(eval_final_q) <- column_names
  eval_final_q$is_concave <- is_c

  return(eval_final_q)
}



fit_eval_concave <- function(x, q_grids, data, formula_grid, omission_rate, omrat_thr,
                             write_summary, addsamplestobackground, weights = NULL,
                             return_replicate, algorithm, AIC_option,
                             proc_for_all) {

  # Arguments:
  # x: Each line of the formula grid
  # q_grids: Formula grid for quadratic terms
  # data: Data to fit models (output of prepare_data)
  # formula_grid: Formula grid
  # omrat_thr: Omission rate threshold
  # write_summary: Write individual summary for each line in formula grid?
  # addsamplestobackground: Add samples to background?
  # weights: Add weights
  # return_replicate: Return results of replicates
  # algorithm: Type of model, either maxnet or glm
  # AIC_option: AIC for maxnet (can be ws or nk)

  # Check arguments

  if(!inherits(q_grids, "data.frame")){
    stop("Argument q_grids must be a data.frame, not ", class(q_grids))
  }

  if(!inherits(data, "prepared_data")){
    stop("Argument data must be a prepared_data object, not ", class(prepared_data))
  }

  if(!inherits(formula_grid, "data.frame")){
    stop("Argument formula_grid must be a data.frame, not ", class(formula_grid))
  }

  if(!inherits(omission_rate, "numeric")){
    stop("Argument omission_rate must be numeric, not ", class(omission_rate))
  }

  if(!inherits(omrat_thr, "numeric")){
    stop("Argument omrat_thr must be numeric, not ", class(omrat_thr))
  }

  if(!inherits(write_summary, "logical")){
    stop("Argument write_summary must be logical, not ", class(write_summary))
  }

  if(!inherits(addsamplestobackground, "logical")){
    stop("Argument addsamplestobackground must be logical, not ", class(addsamplestobackground))
  }

  if(!is.null(weights)){
    if(!inherits(weights, "numeric")){
    stop("Argument weights must be NULL or numeric, not ", class(weights))
  }}

  if(!inherits(return_replicate, "logical")){
    stop("Argument return_replicate must be logical, not ", class(return_replicate))
  }

  if(!inherits(algorithm, "character")){
    stop("Argument algorithm must be a character, not ", class(algorithm))
  }

  if(!(algorithm %in% c("glm", "maxnet"))){
    stop("Argument algorithm must be 'glm' or 'maxnet'")
  }

  if(!inherits(AIC_option, "character")){
    stop("Argument AIC_option must be a character, not ", class(AIC_option))
  }

  if(!(AIC_option %in% c("ws", "nk"))){
    stop("Argument AIC_option must be 'ws' or 'nk'")
  }

  if(!inherits(proc_for_all, "logical")){
    stop("Argument proc_for_all must be logical, not ", class(proc_for_all))
  }

  ####

  grid_x <- q_grids[x, ]

  if (algorithm == "maxnet") {
    # For maxnet model
    formula_x <- as.formula(grid_x$Formulas)
    reg_x <- grid_x$reg_mult
    m_aic <- try(glmnet_mx(p = data$calibration_data$pr_bg,
                           data = data$calibration_data,
                           f = formula_x, regmult = reg_x,
                           addsamplestobackground = addsamplestobackground,
                           weights = weights,
                           calculate_AIC = TRUE,
                           AIC_option = AIC_option),
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
      all_reg <- unique(formula_grid$reg_mult)
      do.call("rbind", lapply(seq_along(all_reg), function(k) {
        grid_x_i <- grid_x
        grid_x_i$reg_mult <- all_reg[k]
        grid_x_i$ID <- formula_grid[formula_grid$Formulas == grid_x$Formulas &
                                      formula_grid$reg_mult == all_reg[k], "ID"]
        return(grid_x_i)
      }))
    } else {
      formula_grid[x, ]
    }

    df_eval_q <- empty_replicates(omission_rate = omission_rate,
                                  n_row = nrow(grid_q) * length(data$kfolds),
                                  replicates = names(data$kfolds),
                                  is_c = is_c,
                                  algorithm = algorithm)
    df_eval_q2 <- cbind(grid_q, df_eval_q)
    eval_final_q <- empty_summary(omission_rate = omission_rate, is_c = is_c,
                                  algorithm = algorithm)
    eval_final_q_summary <- reorder_stats_columns(cbind(grid_q, eval_final_q),
                                                  omission_rate)

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
      om_rate <- kuenm2:::omrat(threshold = omission_rate, pred_train = suit_val_cal,
                       pred_test = suit_val_eval)

      #Calculate PROC? ...
      if(proc_for_all){
        proc_i <- lapply(omission_rate, function(omr){
          proc_omr <- enmpa::proc_enm(test_prediction = suit_val_eval,
                                      prediction = pred_i,
                                      threshold = omr)$pROC_summary
          names(proc_omr) <- c(paste0("Mean_AUC_ratio_at_", omr),
                               paste0("pval_pROC_at_", omr))
          return(proc_omr)
        })
        proc_i <- unlist(proc_i)} else {
          #Or fill PROC with NA
          proc_i <- rep(NA, length(omission_rate) * 2)
          names(proc_i) <- c(paste0("Mean_AUC_ratio_at_", omission_rate),
                             paste0("pval_pROC_at_", omission_rate))
        }


      # #Proc
      # proc_i <- lapply(omission_rate, function(omr){
      #   proc_omr <- enmpa::proc_enm(test_prediction = suit_val_eval,
      #                   prediction = pred_i, threshold = omr)$pROC_summary
      #   names(proc_omr) <- c(paste0("Mean_AUC_ratio_at_", omr),
      #                        paste0("pval_pROC_at_", omr))
      #   return(proc_omr)
      # })
      # proc_i <- unlist(proc_i)


      df_eval_q <-  if (algorithm == "maxnet") {
        data.frame(Replicate = i,
                   t(om_rate),
                   t(proc_i),
                   AIC_nk = m_aic$AIC,
                   AIC_ws = AICc,
                   npar = npar,
                   is_concave = is_c,
                   row.names = NULL)
      } else {
        data.frame(Replicate = i,
                   t(om_rate),
                   t(proc_i),
                   AIC = m_aic$aic,
                   npar = npar,
                   is_concave = is_c,
                   row.names = NULL)
      }
      return(cbind(grid_x, df_eval_q))
    })
    names(mods) <- names(data$kfolds)
    eval_final_q <- do.call("rbind", mods)
    eval_final_q_summary <- reorder_stats_columns(eval_stats(eval_final_q,
                                                             omission_rate,
                                                             algorithm),
                                                  omission_rate = omission_rate)
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



fit_eval_models <- function(x, formula_grid, data, omission_rate, omrat_thr,
                            write_summary, addsamplestobackground, weights = NULL,
                            return_replicate, algorithm, AIC_option,
                            proc_for_all) {

  # Arguments:
  # x: Each line of the formula grid
  # formula_grid: Formula grid (output of calibration_grid)
  # data: Data to fit models (output of prepare_data)
  # omrat_thr: Omission rate threshold
  # write_summary: Write individual summary for each line in formula grid?
  # addsamplestobackground: Add samples to background?
  # weights: Add weights
  # return_replicate: Return results of replicates
  # algorithm: Type of model, either maxnet or glm
  # AIC_option: AIC_option for maxnet (can be ws or nk)

  # Check rguments

  if(!inherits(data, "prepared_data")){
    stop("Argument data must be a prepared_data object, not ", class(prepared_data))
  }

  if(!inherits(formula_grid, "data.frame")){
    stop("Argument formula_grid must be a data.frame, not ", class(formula_grid))
  }

  if(!inherits(omission_rate, "numeric")){
    stop("Argument omission_rate must be numeric, not ", class(omission_rate))
  }

  if(!inherits(omrat_thr, "numeric")){
    stop("Argument omrat_thr must be numeric, not ", class(omrat_thr))
  }

  if(!inherits(write_summary, "logical")){
    stop("Argument write_summary must be logical, not ", class(write_summary))
  }

  if(!inherits(addsamplestobackground, "logical")){
    stop("Argument addsamplestobackground must be logical, not ", class(addsamplestobackground))
  }

  if(!is.null(weights)){
    if(!inherits(weights, "numeric")){
      stop("Argument weights must be NULL or numeric, not ", class(weights))
    }}

  if(!inherits(return_replicate, "logical")){
    stop("Argument return_replicate must be logical, not ", class(return_replicate))
  }

  if(!inherits(algorithm, "character")){
    stop("Argument algorithm must be a character, not ", class(algorithm))
  }

  if(!(algorithm %in% c("glm", "maxnet"))){
    stop("Argument algorithm must be 'glm' or 'maxnet'")
  }

  if(!inherits(AIC_option, "character")){
    stop("Argument AIC_option must be a character, not ", class(AIC_option))
  }

  if(!(AIC_option %in% c("ws", "nk"))){
    stop("Argument AIC_option must be 'ws' or 'nk'")
  }

  if(!inherits(proc_for_all, "logical")){
    stop("Argument proc_for_all must be logical, not ", class(proc_for_all))
  }
  ####

  grid_x <- formula_grid[x,] # Get i candidate model

  if (algorithm == "maxnet") {
    # Fit maxnet model
    reg_x <- grid_x$reg_mult # Get regularization multiplier for maxnet
    formula_x <- as.formula(grid_x$Formulas) # Get formula from grid x

    m_aic <- try(glmnet_mx(p = data$calibration_data$pr_bg,
                           data = data$calibration_data,
                           f = formula_x, regmult = reg_x,
                           addsamplestobackground = addsamplestobackground,
                           weights = weights, calculate_AIC = TRUE,
                           AIC_option = AIC_option),
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
      if (!is.null(weights)){
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
      if(proc_for_all & !inherits(mod_i, "try-error")){
        proc_i <- lapply(omission_rate, function(omr){
          proc_omr <- enmpa::proc_enm(test_prediction = suit_val_eval,
                                      prediction = pred_i,
                                      threshold = omr, ...)$pROC_summary
          names(proc_omr) <- c(paste0("Mean_AUC_ratio_at_", omr),
                               paste0("pval_pROC_at_", omr))
          return(proc_omr)
        })
        proc_i <- unlist(proc_i)} else {
          #Or fill PROC with NA
          proc_i <- rep(NA, length(omission_rate) * 2)
          names(proc_i) <- c(paste0("Mean_AUC_ratio_at_", omission_rate),
                             paste0("pval_pROC_at_", omission_rate))
        }

      # #Proc
      # proc_i <- lapply(omission_rate, function(omr){
      #   proc_omr <- enmpa::proc_enm(test_prediction = suit_val_eval,
      #                               prediction = pred_i, threshold = omr)$pROC_summary
      #   names(proc_omr) <- c(paste0("Mean_AUC_ratio_at_", omr),
      #                        paste0("pval_pROC_at_", omr))
      #   return(proc_omr)
      # })
      # proc_i <- unlist(proc_i)


      # Save metrics in a dataframe
      df_eval <-  if (algorithm == "maxnet") {
        data.frame(Replicate = i,
                   t(om_rate),
                   t(proc_i),
                   AIC_nk = m_aic$AIC,
                   AIC_ws = AICc,
                   npar = npar,
                   is_concave = is_c,
                   row.names = NULL)
      } else if (algorithm == "glm") {
        data.frame(Replicate = i,
                   t(om_rate),
                   t(proc_i),
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
                        empty_replicates(omission_rate = omission_rate,
                                         n_row = length(data$kfolds),
                                         replicates = names(data$kfolds),
                                         is_c = is_c, algorithm = algorithm))
  } else {
    # Combine evaluation results
    names(mods) <- names(data$kfolds)
    eval_final <- do.call("rbind", mods)
  }

  # Summarize results using eval_stats
  eval_final_summary <- if (class(mods) == "try-error") {
    reorder_stats_columns(cbind(grid_x, empty_summary(omission_rate, is_c,
                                                      algorithm)),
                          omission_rate = omission_rate)
  } else {
    reorder_stats_columns(eval_stats(eval_final, omission_rate, algorithm),
                          omission_rate = omission_rate)
  }

  # Write summary if requested
  if (write_summary) {
    write.csv(eval_final_summary,
              file.path(out_dir, paste0("Summary_cand_model_", grid_x$ID, ".csv")),
              row.names = FALSE)
  }

  # Return replicates?
  if (!return_replicate) {
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
    stop("Argumen 'x' must be defined")
  }
  if (missing(dfgrid)) {
    stop("Argumen 'dfgrid' must be defined")
  }
  if (missing(cal_res)) {
    stop("Argumen 'cal_res' must be defined")
  }

  if(!inherits(dfgrid, "data.frame")){
    stop("Argument dfgrid must be a data.frame, not ", class(dfgrid))
  }

  if(!inherits(cal_res, "calibration_results")){
    stop("Argument cal_res must be a data.frame, not ", class(calibration_results))
  }

  if(!inherits(n_replicates, "numeric")){
    stop("Argument n_replicates must be numeric, not ", class(n_replicates))
  }
  if (n_replicates > 1) {
    if(!inherits(rep_data, "list")){
      stop("Argument rep_data must be a list, not ", class(rep_data))
    }
  }

  if(!inherits(algorithm, "character")){
    stop("Argument algorithm must be a character, not ", class(algorithm))
  }

  if(!(algorithm %in% c("glm", "maxnet"))){
    stop("Argument algorithm must be 'glm' or 'maxnet'")
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
    best_regm <- best_model$reg_mult  # Regularization multiplier for maxnet
  }

  # Select data for the replicate, or use the entire calibration data
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
  # Assign model ID and replicate number for tracking
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
reorder_stats_columns <- function(stats_final, omission_rate){

  if(!inherits(stats_final, "data.frame")){
    stop("Argument stats_final must be a data.frame, not ", class(stats_final))
  }

  if(!inherits(omission_rate, "numeric")){
    stop("Argument omission_rate must be numeric, not ", class(omission_rate))
  }

  first_cols<- intersect(c("ID", "Formulas", "reg_mult", "Features"),
                         colnames(stats_final))
  metric_cols <- c(paste0("Omission_rate_at_", omission_rate, ".mean"),
                   paste0("Omission_rate_at_", omission_rate, ".sd"),
                   #Proc columns names
                   paste0("Mean_AUC_ratio_at_", omission_rate, ".mean"),
                   paste0("Mean_AUC_ratio_at_", omission_rate, ".sd"),
                   paste0("pval_pROC_at_", omission_rate, ".mean"),
                   paste0("pval_pROC_at_", omission_rate, ".sd"))
  ordered_metric_cols <- unlist(lapply(omission_rate, function(rate) {
    grep(paste0("_", rate, "."), metric_cols, value = TRUE)
  }))
  last_cols <- setdiff(colnames(stats_final), c(first_cols, metric_cols))
  orders_cols <- c(first_cols, ordered_metric_cols, last_cols)
  return(stats_final[,orders_cols])
}

#### PROC ####
proc <- function(x, formula_grid, data, omission_rate = 10,
                 addsamplestobackground = TRUE, weights = NULL,
                 algorithm){

  #Check arguments
  if(!inherits(formula_grid, "data.frame")){
    stop("Argument formula_grid must be a data.frame, not ", class(formula_grid))
  }

  if(!inherits(omission_rate, "numeric")){
    stop("Argument omission_rate must be numeric, not ", class(omission_rate))
  }

  if(!inherits(addsamplestobackground, "logical")){
    stop("Argument addsamplestobackground must be logical, not ", class(addsamplestobackground))
  }

  if(!is.null(weights)){
    if(!inherits(weights, "numeric")){
      stop("Argument weights must be NULL or numeric, not ", class(weights))
    }}

  if(!inherits(algorithm, "character")){
    stop("Argument algorithm must be a character, not ", class(algorithm))
  }

  if(!(algorithm %in% c("glm", "maxnet"))){
    stop("Argument algorithm must be 'glm' or 'maxnet'")
  }
  ####

  grid_x <- formula_grid[x,] # Get i candidate model
  #Get formula and reg
  formula_x <- as.formula(grid_x$Formulas)
  reg_x <- grid_x$reg_mult

  if(algorithm == "glm"){
    formula_x <- as.formula(paste("pr_bg ", grid_x$Formulas))
  }

  # Get background index
  bgind <- which(data$calibration_data$pr_bg == 0)

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
    proc_i <- lapply(omission_rate, function(omr){
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
      data.frame(Replicate = i,
                 t(proc_i),
                 row.names = NULL)
    } else if (algorithm == "glm") {
      data.frame(Replicate = i,
                 t(proc_i),
                 row.names = NULL)
    }
    return(df_proc)
  }), silent = TRUE)
  #Get summary
  proc_df <- do.call("rbind", mods)

  #Get means and sd
  means <- sapply(proc_df[, -1], mean)  # Excluindo a coluna Replicate
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
