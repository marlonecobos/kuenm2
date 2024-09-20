#' Select best among candidate models
#'
#' @export

sel_best_models <- function(cand_models,
                     test_concave = TRUE,
                     # omrat = 5,
                     omrat_threshold = 5,
                     allow_tolerance = T,
                     tolerance = 0.01,
                     AIC = "nk",
                     significance = 0.05,
                     verbose = TRUE,
                     delta_aic = 2,
                     save_file = T,
                     file_name = NULL){

  #Remove column with the other AIC
  if(AIC == "nk") {
    AIC <- "AIC_nk"
    cand_models$AIC_ws <- NULL}

  if(AIC == "ws") {
    AIC <- "AIC_ws"
    cand_models$AIC_nk <- NULL}

  #Which omission rate?
  om_thr <- paste0("Omission_rate_at_", omrat_threshold, ".mean")

  #How many models are being filtered?
  if(verbose){
    message("\nFiltering ", nrow(cand_models), " models")
  }

  #Remove models with errors
  na_models <- cand_models[is.na(cand_models$is_concave), "ID"]
  if(verbose){
    message("Removing ", length(na_models),
            " model(s) because they failed to fit")
  }
  cand_models <- cand_models[!is.na(cand_models$is_concave),]


  #If test concave = TRUE, remove concave curves
  if(test_concave) {
    concave_models <- cand_models[cand_models$is_concave, "ID"]
    if(verbose){
      message("Removing ", length(concave_models),
              " model(s) with concave curves")
    }
    cand_models <- cand_models[!cand_models$is_concave, ]
  } else {concave_models <- 0}


  #Subset models with significativa pROC
  insig_proc <- cand_models[cand_models$proc_pval.mean >= significance |
                            is.na(cand_models$proc_pval.mean), "ID"]
  if(verbose){
    message("Removing ", length(insig_proc),
            " model(s) with non-significant values of pROC")
  }
  cand_models <- cand_models[cand_models$proc_pval.mean < significance &
                                !is.na(cand_models$proc_pval.mean),]

  #Subset models by omission rate
  high_omr <- cand_models[cand_models[,om_thr] > omrat_threshold/100, "ID"]

  cand_om <- cand_models[cand_models[,om_thr] <= omrat_threshold/100,]
  if(verbose){
    message(nrow(cand_om), " models were selected with omission rate below ",
            omrat_threshold, "%")
  }

  if(nrow(cand_om) == 0 & !allow_tolerance) {
  stop("There is no models with values of omission rate below than ",
       omrat_threshold, ".\nTry with allow_tolerance = T")
  }

  #If 0 models were selected and allow tolerance
  if(nrow(cand_om) == 0 & allow_tolerance) {
    min_thr <- min(cand_models[,om_thr])
    cand_om <- subset(cand_models, cand_models[,om_thr] <= min_thr + tolerance)
    high_omr <- cand_models[cand_models[,om_thr] > min_thr + tolerance, "ID"]
    if(verbose){
      message("Minimum value of omission rate (",  round(min_thr*100, 1), "%) is above the selected theshold (", (omrat_threshold),"%).\nApplying tolerance and selecting ", nrow(cand_om), " models with omission rate <", round(min_thr*100 + tolerance, 1), "%")
    }
  }

  #Calculate delta AIC
  cand_om$dAIC <- cand_om[, AIC] - min(cand_om[, AIC])
  #Select delta AIC
  high_aic <- cand_om[cand_om$dAIC > delta_aic, "ID"]
  cand_final <- cand_om[cand_om$dAIC <= delta_aic,]

  if(verbose){
    message("Selecting ", nrow(cand_final), " final model(s) with delta AIC <",
            delta_aic)
  }

  if(save_file == T){
    if(is.null(file_name)){
      stop("file name need to be defined")
      }
  write.csv(cand_final, paste0(file_name, ".csv"),
              row.names = F)
  }

  #Final resuls
  sel_res <- list(cand_final = cand_final,
                  summary = list("delta_AIC" = delta_aic,
                                 "Errors" = na_models,
                                 "Concave" = concave_models,
                                 "Non_sig_pROC" = insig_proc,
                                 "High_omission_rate"= high_omr,
                                 "High_AIC" = high_aic,
                                 "Selected" = cand_final$ID))

  return(sel_res)
  }

# # #Test function
# #With minimum omission rate below the selected threshold
# bm <- sel_best_models(cand_models = cr,
#                        test_concave = TRUE,
#                        omrat = 5,
#                        omrat_threshold = 5, #5%
#                        allow_tolerance = T,
#                        tolerance = 0.01,
#                        AIC = "nk",
#                        significance = 0.05,
#                        verbose = TRUE,
#                        save_file = T,
#                        output_dir = dir,
#                        delta_aic = 2)
# #Save best model
# write.csv(bm, "Models/Piper_fuligineum/selected_models.csv", row.names = F)
#
# #With minimum omission rate above the selected threshold, allowing tolerance
# bm2 <- sel_best_models(cand_models = cr,
#                       test_concave = TRUE,
#                       omrat = 5,
#                       omrat_threshold = 1, #1%
#                       allow_tolerance = T,
#                       tolerance = 0.01,
#                       AIC = "nk",
#                       significance = 0.05,
#                       verbose = TRUE,
#                       delta_aic = 2,
#                       save_file = F,
#                       output_dir = NUL)
#
# #With minimum omission rate above the selected threshold, allowing tolerance
# bm3 <- sel_best_models(cand_models = cr,
#                       test_concave = TRUE,
#                       omrat = 5,
#                       omrat_threshold = 1, #1%
#                       allow_tolerance = F,
#                       tolerance = 0.01,
#                       AIC = "nk",
#                       significance = 0.05,
#                       verbose = TRUE,
#                       delta_aic = 2,
#                       save_file = F,
#                       output_dir = NUL)
#
#
#
#
