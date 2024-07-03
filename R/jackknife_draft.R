#' Jackknife
#'
#' @param fitted object.
#' @param modelID (character) name of  ModelID from the fitted object.
#' @param cv (numeric) number of k-folds froms cross validation. Default = 5.
#'
#' @return list.
#' @export
#'
#'

jackknife <- function(fitted, modelID, cv = 5){

  # initial tests
  if (missing(fitted) | missing(modelID)) {
    stop("Arguments 'fitted', and 'modelID' must be defined.")
  }

  if (!modelID %in% names(fitted[["Models"]])){
    stop(paste0(
      "The 'ModelID' is not correct, check the following: [",
      paste(names(fitted[["Models"]]), collapse = ", ")),
      "]"
    )
  }

  model_info  <- fitted[["selected_models"]]
  model <- fitted[["Models"]][[modelID]][["Full_model"]]
  data <- fitted[["calibration_data"]]

  # get the model predictors
  ff <- model_info[model_info$ID == gsub("Model_", "",modelID), ]
  features <- names(model$betas)

  # get the formula combinations
  comb <- utils::combn(features, length(features) - 1)

  withon <- sapply(1:length(features), function(x) {
    enmpa::get_formulas(dependent = "", rev(features)[x], type = "l", mode = "complex")
  })

  without <- sapply(1:length(features), function(x) {

    if (length(comb[, x]) > 1) {
      paste(enmpa::get_formulas(dependent = "", comb[, x], type = "l", mode = "complex"), "-1")
    } else {
      enmpa::get_formulas(dependent = "", comb[, x], type = "l", mode = "complex")
    }
  })

  metrics <- c("Omission_rate_at_5", "AIC_ws", "Residual_deviance")

  full <- eval_model_jk(formula = ff$Formulas,
                        reg = ff$regm,
                        data = fitted$calibration_data,
                        cv = cv,
                        omrat_thr = fitted$omission_rate,
                        addsamplestobackground = fitted$addsamplestobackground,
                        weights = fitted$weights)
  full$AIC_nk <- NULL

  withon_r <- lapply(withon, function(x) {
    out <- eval_model_jk(formula = x,
                         reg = ff$regm,
                         data = fitted$calibration_data,
                         cv = cv,
                         omrat_thr = fitted$omission_rate,
                         addsamplestobackground = fitted$addsamplestobackground,
                         weights = fitted$weights)[, metrics]

    return(apply(out, 2, mean))

  })

  without_r <- lapply(without, function(x) {
    out <- eval_model_jk(formula = x,
                         reg = ff$regm,
                         data = fitted$calibration_data,
                         cv = cv,
                         omrat_thr = fitted$omission_rate,
                         addsamplestobackground = fitted$addsamplestobackground,
                         weights = fitted$weights)[, metrics]

    return(apply(out, 2, mean))
  })

  withon_r <- do.call(rbind, withon_r)
  without_r <- do.call(rbind, without_r)

  rownames(withon_r) <- features
  rownames(without_r) <- features


  return(list(Full_model_stats = full,
              Formula = ff$Formulas,
              Without = without_r,
              With_only = withon_r)
  )
}

# Function to evaluate model
eval_model_jk <- function(formula, reg, data, cv = 3,  omrat_thr = 5,
                          addsamplestobackground, weights) {

  #Complete model with AIC
  m_aic <- suppressWarnings(
    glmnet_mx(p = data$pr_bg,
              data = data,
              f = as.formula(formula),
              regmult = reg,
              addsamplestobackground = addsamplestobackground,
              weights = weights,
              calculate_AIC = FALSE)
  )


  #Get number of parameters
  npar <- length(m_aic$betas)

  #Calculate AIC from Warren
  vals <- predict.glmnet_mx(
    m_aic, data[data$pr_bg == 1,],
    type = "exponential"
  )

  AICc <- kuenm2:::aic_ws(pred_occs = vals, ncoefs = npar)
  Residual_deviance <- deviance(m_aic)[200]

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
  bgind <- which(data$pr_bg == 0)

  data_partition <- enmpa::kfold_partition(data = data, k = cv,
                                           dependent = "pr_bg")

  mods <- try(lapply(1:cv, function(i) {
    notrain <- -data_partition[[i]]
    data_i <- data[notrain,]
    #Set weights per k-fold
    if(!is.null(weights)) {
      weights_i <- weights[notrain]
    } else {
      weights_i <- weights
    }

    #Run model
    mod_i <-suppressWarnings(
      glmnet_mx(p = data_i$pr_bg,
                        data = data_i,
                        f = as.formula(formula),
                        regmult = reg,
                        addsamplestobackground = addsamplestobackground,
                        weights = weights_i,
                        calculate_AIC = FALSE)
    )

    #Predict model only to background
    pred_i <- as.numeric(predict.glmnet_mx(object = mod_i,
                                           newdata = data,
                                           clamp = FALSE,
                                           type = "cloglog")
    )

    #Extract suitability in train and test points
    suit_val_cal <- pred_i[unique(c(notrain, -bgind))]
    suit_val_eval <- pred_i[which(!-notrain %in% bgind)]

    #### Calculate omission rate following kuenm ####
    om_rate <- kuenm2:::omrat(threshold = omrat_thr,
                              pred_train = suit_val_cal,
                              pred_test = suit_val_eval)

    #Calculate pROC following enmpa
    proc_i <- enmpa::proc_enm(test_prediction = suit_val_eval,
                              prediction = pred_i)

    ####Save metrics in a dataframe
    df_eval <- data.frame(Replicate = i,
                          t(data.frame(om_rate)),
                          proc_auc_ratio = proc_i$pROC_summary[1],
                          proc_pval = proc_i$pROC_summary[2],
                          AIC_nk = m_aic$AIC,
                          AIC_ws = AICc,
                          npar = npar,
                          is_concave = is_c,
                          Residual_deviance = Residual_deviance,
                          row.names = NULL)

    df_eval2 <- cbind(Formula = formula, reg = reg, df_eval)
    return(df_eval2)
  }), silent = TRUE)


  #Return evaluation final
  eval_final <- do.call("rbind", mods)

  return(eval_final)
}

# Function to plot jackknife
plot_jk <- function(x,
                    metric = "Omission_rate_at_5",
                    legend = TRUE,
                    colors = c("cyan", "blue", "red"),
                    xlab = NULL, xlim = NULL, main = NULL){

  if (missing(x)) {
    stop("Arguments 'x' must be defined.")
  }

  if (!metric %in% c("Omission_rate_at_5", "AIC_ws", "Residual_deviance")){

  }


  # get matrix
  m <- as.matrix(cbind(x$Without[, metric], x$With_only[, metric]))
  colnames(m) <- c("without", "with_only")

  fm <- mean(x[["Full_model_stats"]][[metric]])

  if(is.null(xlab)){xlab = metric}
  if(is.null(main)){main = paste("Jackknife of", metric)}

  if (is.null(xlim)){
    xlim = c(0, max(c(max(m), fm)) + 0.10*max(c(max(m), fm)))
    }

  # Adjust margins to accommodate variable names and legend
  original_par <- par(no.readonly = TRUE)  # Save the original par settings

  if (legend){
    par(mar = c(5, 8, 4, 8) + 0.1)
  } else {
    par(mar = c(5, 8, 4, 6) + 0.1)
  }

  barplot(t(m),
          beside = TRUE,
          horiz = TRUE,
          xlim = xlim,
          las=1,
          col =  colors[1:2],
          xlab = xlab,
          cex.names = 1,
          main = main)
  abline(v = fm, col = colors[3], lty = "dotted", lwd = 2)

  if (legend){
    legend("topright",
           inset=c(-0.285, -0.075), xpd = TRUE, horiz = FALSE, bty = "n",
           cex = 0.95,
           legend=c("Without variable", "With only variable", "Full model"),
           fill = colors)
  }
  # Restore the original par settings
  par(original_par)
}




