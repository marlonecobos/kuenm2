# # # ####Import calibration results####
# # cr <- read.csv("Models/Piper_fuligineum/calibration_results_v2.csv")
# # dir <- "Models/Piper_fuligineum/"
#
# #Test function
# cand_models <- cr

####Function to select best models####
sel_best_models <- function(cand_models,
                     test_convex = TRUE,
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

  if(AIC == "nk") {
    AIC <- "AIC_nk"}

  if(AIC == "ws") {
    AIC <- "AIC_ws"}


  #Which omission rate?
  om_thr <- paste0("Omission_rate_at_", omrat_threshold, ".mean")

  #How many models are being filtered?
  if(verbose){
    message("\nFiltering ", nrow(cand_models), " models")
  }

  #If test convex = TRUE, remove convex curves
  if(test_convex) {
    if(verbose){
      message("Removing ", nrow(subset(cand_models, is_convex == TRUE)),
              " models with convex curves")
    }
    cand_models <- subset(cand_models, is_convex == FALSE)
  }

  #Remove NAs from results
  if(verbose){
    message("Removing ", nrow(cand_models) - nrow(na.omit(cand_models)),
            " models because NAs were found")
  }
  cand_models <- na.omit(cand_models)

  #Subset models with significativa pROC
  if(verbose){
    message("Removing ", nrow(subset(cand_models,
                                     proc_pval.mean >= significance)),
            " models with non-significant values of pROC")
  }
  cand_models <- subset(cand_models, proc_pval.mean <= significance)

  #Subset models by omission rate
  cand_om <- subset(cand_models, cand_models[,om_thr] <= omrat_threshold/100)
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
    if(verbose){
      message("Minimum value of omission rate (",  round(min_thr*100, 1), "%) is above the selected theshold (", (omrat_threshold),"%).\nApplying tolerance and selecting ", nrow(cand_om), " models with omission rate <", round(min_thr*100 + tolerance, 1), "%")
    }
  }

  #Calculate delta AIC
  cand_om$dAIC <- cand_om[, AIC] - min(cand_om[, AIC])
  #Select delta AIC
  cand_final <- subset(cand_om, cand_om$dAIC <= delta_aic)

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

  return(cand_final)
  }

# # #Test function
# #With minimum omission rate below the selected threshold
# bm <- sel_best_models(cand_models = cr,
#                        test_convex = TRUE,
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
#                       test_convex = TRUE,
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
#                       test_convex = TRUE,
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
