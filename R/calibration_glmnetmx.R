#' #' Calibration, evaluation and selection of candidate models
#' #'
#' #' @description
#' #' This function fits and validates candidate models using the data and grid of
#' #' formulas prepared with `prepare_data()`. Then, it selects the best models
#' #' based on concave curves (optional), omission rate, and AIC values.
#' #'
#' #' @usage
#' #' calibration_glmnetmx(data, test_concave = TRUE,
#' #'                      addsamplestobackground = TRUE, use_weights = FALSE,
#' #'                      parallel = TRUE, ncores = 1, parallel_type = "doSNOW",
#' #'                      progress_bar = TRUE, write_summary = FALSE,
#' #'                      out_dir = NULL, skip_existing_models = FALSE,
#' #'                      return_replicate = TRUE, omrat_threshold = 10,
#' #'                      AIC = "ws", delta_aic = 2, allow_tolerance = TRUE,
#' #'                      tolerance = 0.01, verbose = TRUE)
#' #'
#' #' @param data an object of class `prepare_data` returned by the prepare_data()
#' #' function.
#' #' @param formula_grid an object of class `formula_grid` containing the grid of
#' #' formulas returned by the calibration_grid_glmnetmx() function.
#' #' @param test_concave (logical) whether to test for and remove candidate models
#' #' presenting concave curves. Default is TRUE.
#' #' @param addsamplestobackground (logical) whether to add to the background any
#' #' presence sample that is not already there. Default is TRUE.
#' #' @param use_weights (logical) whether to apply the weights present in the
#' #' data. Default is FALSE.
#' #' @param parallel (logical) whether to fit the candidate models in parallel.
#' #' Default is FALSE.
#' #' @param ncores (numeric) number of cores to use for parallel processing.
#' #' Default is 1. This is only applicable if `parallel = TRUE`.
#' #' @param parallel_type (character) the package to use for parallel processing:
#' #' "doParallel" or "doSNOW". Default is "doSNOW". This is only applicable if
#' #' `parallel = TRUE`.
#' #' @param progress_bar (logical) whether to display a progress bar during
#' #' processing. Default is TRUE.
#' #' @param write_summary (logical) whether to save the evaluation results for
#' #' each candidate model to disk. Default is FALSE.
#' #' @param out_dir (character) the file name, with or without a path, for saving
#' #' the evaluation results for each candidate model. This is only applicable if
#' #' `write_summary = TRUE`.
#' #' @param skip_existing_models (logical) whether to check for and skip candidate
#' #' models that have already been fitted and saved in `out_dir`. This is only
#' #' applicable if `write_summary = TRUE`. Default is FALSE.
#' #' @param return_replicate (logical) whether to return the evaluation results
#' #' for each replicate. Default if FALSE, meaning only the summary (mean and
#' #' standard deviation) of the evaluation results will be returned.
#' #' @param omrat_threshold (numeric) a value from 0 to 100 representing the
#' #' percentage of potential error (E) that the data due to any source of
#' #' uncertainty. Default = 10.
#' #' @param AIC (character) the type of AIC to be calculated: "ws" for AIC
#' #' proposed by Warren and Seifert (2011), or "nk" for AIC proposed by Ninomiya
#' #' and Kawano (2016). Default is "ws". See References for details.
#' #' @param delta_aic the value of delta AIC used as a threshold to select models.
#' #' Default is 2.
#' #' @param allow_tolerance (logical) whether to allow selection of models with
#' #' minimum values of omission rates even if their omission rate surpasses the
#' #' `omrat_threshold`. This is only applicable if  all candidate models have
#' #' omission rates higher than the `omrat_threshold`. Default is TRUE.
#' #' @param tolerance (numeric) The value added to the minimum omission rate if it
#' #' exceeds the `omrat_threshold`. If `allow_tolerance = TRUE`, selected models
#' #' will have an omission rate equal to or less than the minimum rate plus this
#' #' tolerance. Default is 0.01.
#' #' @param verbose (logical) whether to display messages during processing.
#' #' Default is TRUE.
#' #'
#' #' @importFrom parallel makeCluster stopCluster
#' #' @importFrom doParallel registerDoParallel
#' #' @importFrom doSNOW registerDoSNOW
#' #' @importFrom foreach foreach `%dopar%`
#' #' @importFrom utils txtProgressBar setTxtProgressBar
#' #' @importFrom stats aggregate
#' #' @importFrom glmnet glmnet.control glmnet
#' #'
#' #' @export
#' #'
#' #' @return
#' #' An object of class 'calibration_results' containing the following elements:
#' #' - species: a character string with the name of the species.
#' #' - calibration data: a data.frame containing a column (`pr_bg`) that
#' #' identifies occurrence points (1) and background points (0), along with the
#' #' corresponding values of predictor variables for each point.
#' #' - formula_grid: data frame containing the calibration grid with possible
#' #' formulas and parameters.
#' #' - kfolds: a list of vectors with row indices corresponding to each fold.
#' #' - data_xy: a data.frame with occurrence and background coordinates.
#' #' - continuous_variables: a character indicating the continuous variables.
#' #' - categorical_variables: a character, categorical variable names (if used).
#' #' - weights: a numeric vector specifying weights for data_xy (if used).
#' #' - pca: if a principal component analysis was performed with variables, a list
#' #' of class "prcomp". See ?stats::prcomp() for details.
#' #' - model_type: the model type (glm or glmnet)
#' #' - calibration_results: a list containing a data frame with all evaluation
#' #' metrics for all replicates (if `return_replicate = TRUE`) and a summary of
#' #' the evaluation metrics for each candidate model.
#' #' - omission_rate: The omission rate determined by `omrat_threshold`.
#' #' - addsampletobackground: a logical value indicating whether any presence
#' #' sample not already in the background was added. Default is TRUE.
#' #' - selected_models:  data frame with the ID and the summary of evaluation
#' #' metrics for the selected models.
#' #' - summary: A list containing the delta AIC values for model selection, and
#' #' the ID values of models that failed to fit, had concave curves,
#' #' non-significant pROC values, omission rates above the threshold, delta AIC
#' #' values above the threshold, and the selected models.
#' #'
#' #' @references
#' #' Ninomiya, Yoshiyuki, and Shuichi Kawano. "AIC for the Lasso in generalized
#' #' linear models." (2016): 2537-2560.
#' #'
#' #' Warren, D. L., & Seifert, S. N. (2011). Ecological niche modeling in Maxent:
#' #' the importance of model complexity and the performance of model selection
#' #' criteria. Ecological applications, 21(2), 335-342.
#' #'
#' #' @examples
#' #' # Import occurrences
#' #' data(occ_data, package = "kuenm2")
#' #'
#' #' # Import variables
#' #' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#' #'                                package = "kuenm2"))
#' #'
#' #' # Use only variables 1, 2 and 3
#' #' var <- var[[1:3]]
#' #'
#' #' # Prepare data
#' #' sp_swd <- prepare_data(model_type = "glmnet", occ = occ_data,
#' #'                        species = occ_data[1, 1], x = "x", y = "y",
#' #'                        spat_variables = var, mask = NULL,
#' #'                        categorical_variables = NULL,
#' #'                        do_pca = FALSE, deviance_explained = 95,
#' #'                        min_explained = 5, center = TRUE, scale = TRUE,
#' #'                        write_pca = FALSE, output_pca = NULL, nbg = 100,
#' #'                        kfolds = 4, weights = NULL, min_number = 2,
#' #'                        min_continuous = NULL,
#' #'                        features = c("l", "lq"),
#' #'                        regm = 1,
#' #'                        include_xy = TRUE,
#' #'                        write_file = FALSE, file_name = NULL,
#' #'                        seed = 1)
#' #'
#' #' # Calibrate models
#' #' m <- calibration_glmnetmx(data = sp_swd,
#' #'                           test_concave = TRUE,
#' #'                           parallel = FALSE,
#' #'                           ncores = 1,
#' #'                           progress_bar = TRUE,
#' #'                           write_summary = FALSE,
#' #'                           out_dir = NULL,
#' #'                           parallel_type = "doSNOW",
#' #'                           return_replicate = TRUE,
#' #'                           omrat_threshold = 10,
#' #'                           allow_tolerance = TRUE,
#' #'                           tolerance = 0.01,
#' #'                           AIC = "ws",
#' #'                           delta_aic = 2,
#' #'                           skip_existing_models = FALSE,
#' #'                           verbose = TRUE)
#' #' m
#'
#' calibration_glmnetmx <- function(data,
#'                                 test_concave = TRUE,
#'                                 addsamplestobackground = TRUE,
#'                                 use_weights = FALSE,
#'                                 parallel = TRUE,
#'                                 ncores = 1,
#'                                 parallel_type = "doSNOW",
#'                                 progress_bar = TRUE,
#'                                 write_summary = FALSE,
#'                                 out_dir = NULL,
#'                                 skip_existing_models = FALSE,
#'                                 return_replicate = TRUE,
#'                                 omrat_threshold = 10,
#'                                 AIC = "ws",
#'                                 delta_aic = 2,
#'                                 allow_tolerance = TRUE,
#'                                 tolerance = 0.01,
#'                                 verbose = TRUE) {
#'
#'   #Check calibration data class and convert to dataframe if necessay
#'   if(is.matrix(data$calibration_data) | is.array(data$calibration_data)) {
#'     data$calibration_data <- as.data.frame(data$calibration_data)
#'   }
#'
#'   #If write_summary = TRUE, create directory
#'   if(write_summary){
#'     if(!file.exists(out_dir))
#'       dir.create(out_dir)
#'   }
#'
#'   model_type <- data$model_type
#'   formula_grid <- data$formula_grid
#'
#'   #If skip_existing_models = TRUE, update grid
#'   if(skip_existing_models & write_summary) {
#'     ready_models <- list.files(path = out_dir, pattern = "summary", full.names = T)
#'     ready_models <- do.call("rbind", lapply(seq_along(ready_models), function(i){
#'       read.csv( ready_models[i])
#'     }))
#'     run_models <- setdiff(formula_grid$ID, ready_models$ID)
#'     if(length(run_models) == 0) {
#'       stop(paste("All models completed. Check the folder:", out_dir))
#'     } else { #Update formula grid
#'       formula_grid <- formula_grid[formula_grid$ID %in% run_models, ]
#'     }
#'   }
#'
#'   # Warning about samples added to background when weights are null
#'   if (verbose & addsamplestobackground & use_weights & is.null(data$weights)) {
#'   message("Weights for samples added to background are the same as in samples.")
#'   }
#'
#'   # Check weights
#'   if (use_weights & is.null(data$weights)) {
#'     message("'use_weights' = TRUE, but weights are not present in 'data'.\nSetting 'use_weights' = FALSE.")
#'     use_weights <- FALSE
#'   }
#'
#'   # If weights is numeric, check if it has the same size of calibration data
#'     if(use_weights & length(data$weights) != nrow(data$calibration_data)) {
#'       stop("length of weights does not match number of rows in calibration_data")
#'     }
#'
#'   #Make cluster
#'   if (parallel) {
#'     cl <- parallel::makeCluster(ncores)
#'     }
#'
#'   # Task 1: Checking concave curves in quadratic models
#'   if (test_concave) {
#'     if (verbose) {
#'       cat("\n Task 1 of 2: checking concave curves in quadratic models\n")
#'     }
#'
#'     q_grids <- formula_grid[grepl("q", formula_grid$Features), ]
#'     n_tot <- nrow(q_grids)
#'
#'     if(n_tot == 0) {
#'       warning("All quadratic models have been already tested")
#'     } else {
#'       if (isTRUE(progress_bar)) {
#'         pb <- txtProgressBar(min = 0, max = n_tot, style = 3)
#'         progress <- function(n) setTxtProgressBar(pb, n)
#'         opts <- list(progress = progress)}
#'
#'       if (parallel & parallel_type == "doParallel") {
#'         doParallel::registerDoParallel(cl)
#'         opts <- NULL #Progress bar does not work with doParallel
#'       }
#'
#'       if (parallel & parallel_type == "doSNOW") {
#'         doSNOW::registerDoSNOW(cl)
#'         if (isTRUE(progress_bar))
#'           opts <- list(progress = progress)
#'         else opts <- opts
#'       }
#'
#'       #Test concave curves
#'       #In parallel (using %dopar%)
#'       if(parallel){
#'       results_concave <- foreach::foreach(
#'         x = 1:n_tot, .packages = c("glmnet", "enmpa"), .options.snow = opts
#'         # ,
#'         #
#'         # .export = c(to_export, "formula_grid", "q_grids","data",
#'         #             "write_summary", "return_replicate")
#'       ) %dopar% {
#'         fit_eval_concave(x = x, q_grids, data, formula_grid,
#'                          omrat_thr = omrat_threshold,
#'                          write_summary = write_summary,
#'                          addsamplestobackground = addsamplestobackground,
#'                          weights = data$weights,
#'                          return_replicate = return_replicate,
#'                          AIC = AIC,
#'                          model_type = model_type)
#'         }
#'       } else { #Not in parallel (using %do%)
#'         results_concave <- vector("list", length = n_tot)
#'         # Loop for com barra de progresso manual
#'         for (x in 1:n_tot) {
#'           # Execute a função fit_eval_models
#'           results_concave[[x]] <- fit_eval_concave(
#'             x = x, q_grids, data, formula_grid, omrat_thr = omrat_threshold,
#'             write_summary = write_summary,
#'             addsamplestobackground = addsamplestobackground,
#'             weights = data$weights, return_replicate = return_replicate,
#'             AIC = AIC, model_type = model_type
#'           )
#'
#'           # Sets the progress bar to the current state
#'           if(progress_bar){
#'             setTxtProgressBar(pb, x) }
#'         }
#'       }
#'
#'   } #End of if(n > 0)
#'   }  #End of If test_concave = TRUE
#'
#'
#'   #Update grid
#'   if(!test_concave) {n_tot = 0}
#'   if(test_concave & n_tot > 0){
#'
#'     #Convert results to dataframe
#'     #Replicate
#'     d_concave_rep <- do.call("rbind", lapply(results_concave, function(x) x$Replicates))
#'     row.names(d_concave_rep) <- NULL
#'     #Summary
#'     d_concave_sum <- do.call("rbind", lapply(results_concave, function(x) x$Summary))
#'
#'
#'     # d_concave <- do.call("rbind", results_concave)
#'     # d_concave <- d_concave[d_concave$is_concave == TRUE, ]
#'
#'     #Identify formulas tested with concave curves
#'     formula_grid <- formula_grid[!(formula_grid$ID %in% d_concave_sum$ID), ]
#'   }
#'
#'   #Set number of iteration based on new grid
#'   n_tot <- nrow(formula_grid)
#'   #If n == 0, skip non-quadratic models
#'   if(n_tot == 0) {
#'     warning("All non-quadratic models have been already tested")
#'   } else {
#'
#'     #Show progress bar? - Update
#'     if (progress_bar) {
#'       pb <- txtProgressBar(0, n_tot, style = 3)
#'       progress <- function(n) setTxtProgressBar(pb, n) }
#'
#'     if (parallel_type == "doParallel") {
#'       doParallel::registerDoParallel(cl)
#'       opts <- NULL
#'     }
#'
#'     if (parallel_type == "doSNOW") {
#'       doSNOW::registerDoSNOW(cl)
#'       if (isTRUE(progress_bar))
#'         opts <- list(progress = progress)
#'       else opts <- NULL
#'     }
#'
#'     if(verbose) {
#'       if(test_concave) {
#'         cat("\nTask 2 of 2: calibrating non-quadratic models and quadratic models
#'         without concave curves\n
#'         ")
#'       } else {
#'         cat("
#'         Task 1 of 1: calibrating models\n")
#'       } }
#'
#'     ####Results non-concave####
#'     #Test concave curves
#'     #In parallel (using %dopar%)
#'     if(parallel){
#'     results <- foreach(
#'       x = 1:n_tot, .packages = c("glmnet", "enmpa"), .options.snow = opts
#'       # ,
#'       # .export = c(to_export, "formula_grid", "q_grids", "data",
#'       #             "write_summary", "return_replicate")
#'     ) %dopar% {
#'       fit_eval_models(x, formula_grid, data,
#'                       omrat_thr = omrat_threshold,
#'                       write_summary = write_summary,
#'                       addsamplestobackground = addsamplestobackground,
#'                       weights = data$weights,
#'                       return_replicate = return_replicate, AIC = AIC)
#'       }
#'     } else { #Not in parallel (using %do%)
#'       results <- vector("list", length = n_tot)
#'       # Loop for com barra de progresso manual
#'       for (x in 1:n_tot) {
#'         # Execute a função fit_eval_models
#'         results[[x]] <- fit_eval_models(
#'           x, formula_grid = formula_grid, data = data, omrat_thr = omrat_threshold,
#'           addsamplestobackground =  addsamplestobackground,
#'           weights = data$weights, write_summary, return_replicate,
#'           AIC = AIC
#'         )
#'
#'         # Sets the progress bar to the current state
#'         if(progress_bar){
#'         setTxtProgressBar(pb, x) }
#'       }
#'     }
#'
#'     #Stop cluster
#'     if(parallel){
#'     parallel::stopCluster(cl) }
#'
#'     #Convert to dataframe
#'     #Replicate
#'     d_res_rep <- do.call("rbind", lapply(results, function(x) x$Replicates))
#'     row.names(d_res_rep) <- NULL
#'     #Summary
#'     d_res_sum <- do.call("rbind", lapply(results, function(x) x$Summary))
#'
#'     # Join results with results concave, if it exists
#'     if(test_concave) {
#'       replicates_final <- rbind(d_concave_rep, d_res_rep)
#'       summary_final <- rbind(d_concave_sum, d_res_sum)
#'       res_final <- list(All_results = replicates_final,
#'                         Summary = summary_final)
#'     } else {
#'       res_final <- list(All_results =  d_res_rep,
#'                         Summary = d_res_sum)
#'     }
#'
#'   } #End of if(n == 0)
#'
#'   #Select best models
#'   bm <- sel_best_models(cand_models = res_final$Summary, #dataframe with Candidate results
#'                         test_concave = test_concave, #Remove models with concave curves?
#'                         omrat_threshold = omrat_threshold, #Omission rate (test points) used to select the best models
#'                         allow_tolerance = allow_tolerance, #If omission rate is higher than set, select the model with minimum omission rate
#'                         tolerance = tolerance, #If allow tollerance, select the model with minimum omission rate + tolerance
#'                         AIC = AIC, #Which AIC? japones (nk) or Warrien (ws?
#'                         significance = 0.05, #Significante to select models based on pROC
#'                         verbose = verbose, #Show messages?
#'                         delta_aic = delta_aic, #Delta AIC to select best models
#'                         save_file = F, #Save file with best models?
#'                         file_name = NULL)
#'   #Concatenate final results
#'   fm <- c(data, calibration_results = list(res_final),
#'           omission_rate = omrat_threshold,
#'           addsamplestobackground = addsamplestobackground,
#'           weights = list(data$weights), selected_models = list(bm$cand_final),
#'           summary = list(bm$summary))
#'   class(fm) <- "calibration_results"
#'
#'   return(fm)
#'
#' } #End of function
#'
