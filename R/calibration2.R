#' Calibration, evaluation, and selection of candidate models
#'
#' @description
#' This function fits and validates candidate models using the data and grid of
#' formulas prepared with \code{\link{prepare_data}}(). It supports both `glm` and `maxnet`
#' model types. The function then selects the best models based on concave
#' curves (optional), partial ROC, omission rate, and AIC values.
#'
#' @usage
#' calibration2(data, test_concave = TRUE, proc_for_all = FALSE,
#'              extrapolation_factor = 0.1,
#'              var_limits = NULL, averages_from = "pr" ,
#'              addsamplestobackground = TRUE,
#'              use_weights = FALSE, parallel = TRUE, ncores = 4,
#'              parallel_type = "doSNOW", progress_bar = TRUE,
#'              write_summary = FALSE,
#'              out_dir = NULL, skip_existing_models = FALSE,
#'              return_replicate = TRUE, omission_rate= 10, omrat_threshold = 10,
#'              AIC = "ws", delta_aic = 2, allow_tolerance = TRUE,
#'              tolerance = 0.01, verbose = TRUE)
#'
#' @param data an object of class `prepared_data` returned by the
#' \code{\link{prepare_data()}} function. It contains the calibration data,
#' formulas grid, kfolds, and model type.
#' @param test_concave (logical) whether to test for and remove candidate models
#' presenting concave curves. Default is TRUE.
#' @param proc_for_all (logical) whether to apply partial ROC tests to all
#' candidate models or only to the selected models. Default is FALSE, meaning
#' that tests are applied only to the selected models.
#' @param extrapolation_factor (numeric) a multiplier used to calculate the
#' extrapolation range for each variable when detecting concave curves. Larger
#' values allow broader extrapolation beyond the observed data range, while
#' smaller values restrict the range. Default is 0.1. See details.
#' @param var_limits (list) A named list specifying the lower and/or upper limits
#' for some variables. The first value represents the lower limit, and the
#' second value represents the upper limit. Default is \code{NULL}, meaning no
#' specific limits are applied, and the range will be calculated using the
#' \code{extrapolation_factor}. See details.
#' @param averages_from (character) specifies how the averages or modes of the
#' variables are calculated. Available options are "pr" (to calculate averages
#' from the presence localities) or "pr_bg" (to use the combined set of presence
#' and background localities). Default is "pr". See details.
#' @param addsamplestobackground (logical) whether to add to the background any
#' presence sample that is not already there. Default is TRUE.
#' @param use_weights (logical) whether to apply the weights present in the
#' data. Default is FALSE.
#' @param parallel (logical) whether to fit the candidate models in parallel.
#' Default is FALSE.
#' @param ncores (numeric) number of cores to use for parallel processing.
#' Default is 1. This is only applicable if `parallel = TRUE`.
#' @param parallel_type (character) the package to use for parallel processing:
#' "doParallel" or "doSNOW". Default is "doSNOW". This is only applicable if
#' `parallel = TRUE`.
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is TRUE.
#' @param write_summary (logical) whether to save the evaluation results for
#' each candidate model to disk. Default is FALSE.
#' @param out_dir (character) the file name, with or without a path, for saving
#' the evaluation results for each candidate model. This is only applicable if
#' `write_summary = TRUE`.
#' @param skip_existing_models (logical) whether to check for and skip candidate
#' models that have already been fitted and saved in `out_dir`. This is only
#' applicable if `write_summary = TRUE`. Default is FALSE.
#' @param return_replicate (logical) whether to return the evaluation results
#' for each replicate. Default is TRUE, meaning evaluation results for each
#' replicate will be returned.
#' @param omission_rate (numeric) values from 0 to 100 representing the
#' percentage of potential error due to any source of uncertainty. This value is
#' used to calculate the omission rate. Default is 10. See details.
#' @param omrat_threshold (numeric) the maximum omission rate a candidate model
#' can have to be considered a best model. This value must match one of the
#' values specified in omrat. Defaut is 10.
#' @param AIC (character) the type of AIC to be calculated: "ws" for AIC
#' proposed by Warren and Seifert (2011), or "nk" for AIC proposed by Ninomiya
#' and Kawano (2016). Default is "ws". See References for details.
#' @param delta_aic (numeric) the value of delta AIC used as a threshold to
#' select models. Default is 2.
#' @param allow_tolerance (logical) whether to allow selection of models with
#' minimum values of omission rates even if their omission rate surpasses the
#' `omrat_threshold`. This is only applicable if all candidate models have
#' omission rates higher than the `omrat_threshold`. Default is TRUE.
#' @param tolerance (numeric) The value added to the minimum omission rate if it
#' exceeds the `omrat_threshold`. If `allow_tolerance = TRUE`, selected models
#' will have an omission rate equal to or less than the minimum rate plus this
#' tolerance. Default is 0.01.
#' @param verbose (logical) whether to display messages during processing.
#' Default is TRUE.
#' @param ... Additional arguments passed to \code{\link[enmpa]{proc_enm}} for
#' calculating partial ROC. See ?enmpa::proc_enm
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach `%dopar%`
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats aggregate
#' @importFrom glmnet glmnet.control glmnet
#'
#' @export
#'
#' @return
#' An object of class 'calibration_results' containing the following elements:
#' - species: a character string with the name of the species.
#' - calibration data: a data.frame containing a column (`pr_bg`) that
#' identifies occurrence points (1) and background points (0), along with the
#' corresponding values of predictor variables for each point.
#' - formula_grid: data frame containing the calibration grid with possible
#' formulas and parameters.
#' - kfolds: a list of vectors with row indices corresponding to each fold.
#' - data_xy: a data.frame with occurrence and background coordinates.
#' - continuous_variables: a character indicating the continuous variables.
#' - categorical_variables: a character, categorical variable names (if used).
#' - weights: a numeric vector specifying weights for data_xy (if used).
#' - pca: if a principal component analysis was performed with variables, a list
#' of class "prcomp". See ?stats::prcomp() for details.
#' - algorithm: the model type (glm or maxnet)
#' - calibration_results: a list containing a data frame with all evaluation
#' metrics for all replicates (if `return_replicate = TRUE`) and a summary of
#' the evaluation metrics for each candidate model.
#' - omission_rate: The omission rate determined by `omrat_threshold` for
#' selecting best models.
#' - addsamplestobackground: a logical value indicating whether any presence
#' sample not already in the background was added.
#' - selected_models: data frame with the ID and the summary of evaluation
#' metrics for the selected models.
#' - summary: A list containing the delta AIC values for model selection, and
#' the ID values of models that failed to fit, had concave curves,
#' non-significant pROC values, omission rates above the threshold, delta AIC
#' values above the threshold, and the selected models.
#'
#' @references
#' Ninomiya, Yoshiyuki, and Shuichi Kawano. "AIC for the Lasso in generalized
#' linear models." (2016): 2537-2560.
#'
#' Warren, D. L., & Seifert, S. N. (2011). Ecological niche modeling in Maxent:
#' the importance of model complexity and the performance of model selection
#' criteria. Ecological applications, 21(2), 335-342.
#'
#' @details
#' Partial ROC is calculated following Peterson et al.
#' (2008; http://dx.doi.org/10.1016/j.ecolmodel.2007.11.008).
#'
#' Omission rates are calculated using models trained with separate testing data
#' subsets. Users can specify multiple omission rates to be calculated
#' (e.g., c(5, 10)), though only one can be used as the threshold for selecting
#' the best models.
#'
#' Model complexity (AIC) is assessed using models generated with the complete
#' set of occurrences.
#'
#' Concave curves are identified by analyzing the beta coefficients of quadratic
#' terms within the variable's range. The range for extrapolation is calculated
#' as the difference between the variable's maximum and minimum values in the
#' model, multiplied by the extrapolation factor. A concave curve is detected
#' when the beta coefficient is positive, and the vertex — where the curve
#' changes direction — lies between the lower and upper limits of the variable.
#'
#' Users can specify the lower and upper limits for certain variables using
#' \code{var_limits}. For example, if \code{var_limits = list("bio12" = c(0, NA),
#' "bio15" = c(0, 100))}, the lower limit for \code{bio12} will be 0, and the
#' upper limit will be calculated using the extrapolation factor. Similarly,
#' the lower and upper limits for \code{bio15} will be 0 and 100, respectively.
#'
#' For calculating the vertex position, a response curve for a given variable is
#' generated with all other variables set to their mean values (or mode for
#' categorical variables). These values are calculated either from the presence
#' localities (if \code{averages_from = "pr"}) or from the combined set of
#' presence and background localities (if \code{averages_from = "pr_bg"}).
#'
#' @examples
#' # Import occurrences
#' data(occ_data, package = "kuenm2")
#'
#' # Import variables
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#' #Remove categorical variable
#' var <- var[[1:4]]
#'
#' #### maxnet ####
#' # Prepare data for maxnet model
#' sp_swd <- prepare_data(algorithm = "maxnet", occ = occ_data,
#'                        species = occ_data[1, 1], x = "x", y = "y",
#'                        raster_variables = var,
#'                        n_background = 100,
#'                        features = c("l", "lq"),
#'                        reg_mult = 1)
#'
#'
#' # Calibrate maxnet models
#' m <- calibration2(data = sp_swd,
#'                   test_concave = TRUE,
#'                   proc_for_all = FALSE,
#'                   extrapolation_factor = 0.1,
#'                   var_limits = list("bio_7" = c(0, NA),
#'                                     "bio_12" = c(0, NA),
#'                                     "bio_15" = c(0, 100)),
#'                   averages_from = "pr",
#'                   parallel = FALSE,
#'                   ncores = 1,
#'                   progress_bar = TRUE,
#'                   write_summary = FALSE,
#'                   out_dir = NULL,
#'                   parallel_type = "doSNOW",
#'                   return_replicate = TRUE,
#'                   omission_rate = c(5, 10),
#'                   omrat_threshold = 10,
#'                   allow_tolerance = TRUE,
#'                   tolerance = 0.01,
#'                   AIC = "ws",
#'                   delta_aic = 2,
#'                   skip_existing_models = FALSE,
#'                   verbose = TRUE)
#' m
#'
#' #### GLM ####
#' # Prepare data for glm model
#' sp_swd_glm <- prepare_data(algorithm = "glm", occ = occ_data,
#'                            species = occ_data[1, 1], x = "x", y = "y",
#'                            raster_variables = var,
#'                            n_background = 100,
#'                            features = c("l", "lq"))
#'
#' # Calibrate glm models
#' m_glm <- calibration2(data = sp_swd_glm,
#'                       test_concave = TRUE,
#'                       proc_for_all = FALSE,
#'                       extrapolation_factor = 0.1,
#'                       var_limits = list("bio_7" = c(0, NA),
#'                                         "bio_12" = c(0, NA),
#'                                         "bio_15" = c(0, 100)),
#'                       averages_from = "pr",
#'                       parallel = FALSE,
#'                       ncores = 1,
#'                       progress_bar = TRUE,
#'                       write_summary = FALSE,
#'                       out_dir = NULL,
#'                       parallel_type = "doSNOW",
#'                       return_replicate = TRUE,
#'                       omission_rate = c(5, 10),
#'                       omrat_threshold = 10,
#'                       allow_tolerance = TRUE,
#'                       tolerance = 0.01,
#'                       AIC = "ws",
#'                       delta_aic = 2,
#'                       skip_existing_models = FALSE,
#'                       verbose = TRUE)
#' m_glm
#'
calibration2 <- function(data,
                        test_concave = TRUE,
                        proc_for_all = FALSE,
                        extrapolation_factor = 0.1,
                        var_limits = NULL,
                        averages_from = "pr",
                        addsamplestobackground = TRUE,
                        use_weights = FALSE,
                        parallel = TRUE,
                        ncores = 4,
                        parallel_type = "doSNOW",
                        progress_bar = TRUE,
                        write_summary = FALSE,
                        out_dir = NULL,
                        skip_existing_models = FALSE,
                        return_replicate = TRUE,
                        omission_rate = 10,
                        omrat_threshold = 10,
                        AIC = "ws",
                        delta_aic = 2,
                        allow_tolerance = TRUE,
                        tolerance = 0.01,
                        verbose = TRUE, ...) {

  # Convert calibration data to dataframe if necessary
  if (is.matrix(data$calibration_data) || is.array(data$calibration_data)) {
    data$calibration_data <- as.data.frame(data$calibration_data)
  }

  # If write_summary = TRUE, create directory
  if (write_summary && !file.exists(out_dir)) {
    dir.create(out_dir)
  }

  # Define global vars
  algorithm <- data$algorithm
  formula_grid <- data$formula_grid
  weights <- data$weights

  # Skip existing models if requested
  if (skip_existing_models && write_summary) {
    ready_models <- list.files(path = out_dir, pattern = "summary", full.names = TRUE)
    ready_models <- do.call("rbind", lapply(seq_along(ready_models), function(i) {
      read.csv(ready_models[i])
    }))
    run_models <- setdiff(formula_grid$ID, ready_models$ID)

    if (length(run_models) == 0) {
      stop(paste("All models completed. Check the folder:", out_dir))
    } else {
      formula_grid <- formula_grid[formula_grid$ID %in% run_models, ]
    }
  }

  # Warning about samples added to background when weights are null
  if (verbose & addsamplestobackground & use_weights & is.null(weights)) {
    message("Weights for samples added to background are the same as in samples.")
  }

  # Validate weights
  if (use_weights && is.null(weights)) {
    message("'use_weights' = TRUE, but weights are not present in 'data'.\n",
            "Setting 'use_weights' = FALSE.")
    use_weights <- FALSE
  }

  if (use_weights && length(weights) != nrow(data$calibration_data)) {
    stop("Length of weights does not match number of rows in calibration_data.")
  }

  # Parallelization setup
  if (parallel) {
    if (!(parallel_type %in% c("doSNOW", "doParallel"))) {
      stop("Invalid parallel_type. Use 'doSNOW' or 'doParallel'.")
    }
    cl <- parallel::makeCluster(ncores)
  }

  # Task 1: Checking concave curves in quadratic models_________________________
  # ____________________________________________________________________________

  if (test_concave) {
    if (verbose) {
      cat("\n Task 1/2: checking concave curves in quadratic models\n")
    }

    q_grids <- formula_grid[grepl("q", formula_grid$Features), ]
    n_tot <- as.numeric(nrow(q_grids))

    if(n_tot == 0) {
      warning("All quadratic models have been already tested")
    } else {
      if (isTRUE(progress_bar)) {
        pb <- txtProgressBar(min = 0, max = n_tot, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)} else {opts <- NULL}

      if (parallel & parallel_type == "doParallel") {
        doParallel::registerDoParallel(cl)
        opts <- NULL # Progress bar does not work with doParallel
      }

      if (parallel & parallel_type == "doSNOW") {
        doSNOW::registerDoSNOW(cl)
        if (isTRUE(progress_bar))
          opts <- list(progress = progress)
        else opts <- opts
      }

      # Execute fit_eval_concave in parallel or sequentially
      if(parallel){
        results_concave <- foreach::foreach(
          x = 1:n_tot,
          .packages = c("glmnet", "enmpa"),
          .options.snow = opts) %dopar% {
            kuenm2:::fit_eval_concave2(x = x, q_grids, data, formula_grid,
                             omission_rate = omission_rate,
                             omrat_thr = omrat_threshold,
                             write_summary = write_summary,
                             addsamplestobackground = addsamplestobackground,
                             weights = weights,
                             return_replicate = return_replicate,
                             algorithm = algorithm, AIC = AIC,
                             extrapolation_factor = extrapolation_factor,
                             var_limits = var_limits,
                             averages_from = averages_from, proc_for_all)
          }
      } else {
        results_concave <- vector("list", length = n_tot)
        for (x in 1:n_tot) {
          results_concave[[x]] <-
            kuenm2:::fit_eval_concave2(x = x, q_grids, data, formula_grid,
                             omission_rate = omission_rate,
                             omrat_thr = omrat_threshold,
                             write_summary = write_summary,
                             addsamplestobackground = addsamplestobackground,
                             weights = weights,
                             return_replicate = return_replicate,
                             algorithm = algorithm, AIC = AIC,
                             extrapolation_factor = extrapolation_factor,
                             var_limits = var_limits,
                             averages_from = averages_from,
                             proc_for_all)
          # Sets the progress bar to the current state
          if(progress_bar) setTxtProgressBar(pb, x)
        }
      }
    } # End of if(n > 0)
  } # End of If test_concave = TRUE

  # Update formula grid after concave test
  if(!test_concave) {n_tot = 0}
  if(test_concave & n_tot > 0){

    # Convert results to dataframe
    d_concave_rep <- do.call("rbind", lapply(results_concave,
                                             function(x) x$Replicates))
    row.names(d_concave_rep) <- NULL
    d_concave_sum <- do.call("rbind", lapply(results_concave,
                                             function(x) x$Summary))
    formula_grid <- formula_grid[!(formula_grid$ID %in% d_concave_sum$ID), ]
  }

  # Task 2: Fitting remaining models ___________________________________________
  # ____________________________________________________________________________

  n_tot <- nrow(formula_grid)

  if(n_tot == 0) {
    warning("All non-quadratic models have been already tested")
  } else {

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
        cat("\nTask 2/2: calibrating non-quadratic models and quadratic models",
            "without concave curves\n")
      } else {
        cat("Task 1/1: calibrating models\n")
      }
    }

    if(parallel){
      results <- foreach(
        x = 1:n_tot,
        .packages = c("glmnet", "enmpa"),
        .options.snow = opts ) %dopar% {
          fit_eval_models2(x, formula_grid, data,
                          omission_rate = omission_rate,
                          omrat_thr = omrat_threshold,
                          write_summary = write_summary,
                          addsamplestobackground = addsamplestobackground,
                          weights = weights,
                          return_replicate = return_replicate, AIC = AIC,
                          algorithm = algorithm,
                          extrapolation_factor = extrapolation_factor,
                          var_limits = var_limits,
                          averages_from = averages_from, proc_for_all)
        }
    } else {
      results <- vector("list", length = n_tot)
      for (x in 1:n_tot) {
        results[[x]] <-
          fit_eval_models2(x, formula_grid = formula_grid, data = data,
                          omission_rate = omission_rate,
                          omrat_thr = omrat_threshold,
                          addsamplestobackground =  addsamplestobackground,
                          weights = weights, write_summary,
                          return_replicate, AIC = AIC, algorithm = algorithm,
                          extrapolation_factor = extrapolation_factor,
                          var_limits = var_limits,
                          averages_from = averages_from, proc_for_all)

        if(progress_bar) setTxtProgressBar(pb, x)
      }
    }

    # Stop cluster
    if(parallel) parallel::stopCluster(cl)

    d_res_rep <- do.call("rbind", lapply(results, function(x) x$Replicates))
    row.names(d_res_rep) <- NULL
    d_res_sum <- do.call("rbind", lapply(results, function(x) x$Summary))

    # Join results with results concave, if it exists
    if(test_concave) {
      replicates_final <- rbind(d_concave_rep, d_res_rep)
      summary_final <- rbind(d_concave_sum, d_res_sum)
      res_final <- list(All_results = replicates_final,
                        Summary = summary_final)
    } else {
      res_final <- list(All_results =  d_res_rep, Summary = d_res_sum)
    }
  }# End of if(n == 0)

  if(n_tot == 0){
    res_final <- list(All_results = d_concave_rep,
                      Summary = d_concave_sum)}

  # Select the best models
  bm <- sel_best_models(cand_models = res_final$Summary,
                        test_concave = test_concave,
                        calc_proc = !proc_for_all,
                        data = data,
                        omrat_threshold = omrat_threshold,
                        allow_tolerance = allow_tolerance,
                        tolerance = tolerance, AIC_option = AIC,
                        significance = 0.05, verbose = verbose,
                        delta_aic = delta_aic,
                        algorithm = algorithm,
                        parallel = parallel, ncores = ncores,
                        parallel_type = parallel_type,
                        progress_bar = progress_bar)

  # Concatenate final results
  fm <- c(data, calibration_results = list(res_final),
          omission_rate = omrat_threshold,
          addsamplestobackground = addsamplestobackground,
          weights = list(weights),
          selected_models = list(bm$cand_final),
          summary = list(bm$summary))

  class(fm) <- "calibration_results"
  return(fm)
}

fit_eval_concave2 <- function(x, q_grids, data, formula_grid, omission_rate, omrat_thr,
                             write_summary, addsamplestobackground, weights,
                             return_replicate, algorithm,
                             AIC, extrapolation_factor, var_limits,
                             averages_from = "pr", proc_for_all, ...) {

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
  # AIC: AIC for maxnet (can be ws or nk)

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
                           AIC = AIC),
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
    is_c <- detect_concave(model = m_aic, calib_data = data$calibration_data,
                           extrapolation_factor = extrapolation_factor,
                           var_limits = var_limits,
                           averages_from = averages_from,
                           plot = FALSE)
    is_c <- any(sapply(is_c, function(x) x$is_concave))
  }

  # Handle concave results
  if (isTRUE(is_c) | is.na(is_c)) {
    # If concave, return grid
    grid_q <- if (algorithm == "maxnet") {
      all_reg <- formula_grid$reg_mult[formula_grid$Formulas == grid_x$Formulas &
                                     formula_grid$Features == grid_x$Features]
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
                                  n_row = as.numeric(nrow(grid_q) * length(data$kfolds)),
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
    bgind <- which(data$calibration_data$pr_bg == 0)
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
      om_rate <- omrat(threshold = omission_rate, pred_train = suit_val_cal,
                       pred_test = suit_val_eval)


      #Calculate PROC? ...
      if(proc_for_all){
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

fit_eval_models2 <- function(x, formula_grid, data, omission_rate, omrat_thr,
                            write_summary, addsamplestobackground, weights,
                            return_replicate, algorithm, AIC,
                            extrapolation_factor,
                            var_limits, averages_from, proc_for_all, ...) {
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
  # AIC: AIC for maxnet (can be ws or nk)

  grid_x <- formula_grid[x,] # Get i candidate model

  if (algorithm == "maxnet") {
    # Fit maxnet model
    reg_x <- grid_x$reg_mult # Get regularization multiplier for maxnet
    formula_x <- as.formula(grid_x$Formulas) # Get formula from grid x

    m_aic <- try(glmnet_mx(p = data$calibration_data$pr_bg,
                           data = data$calibration_data,
                           f = formula_x, regmult = reg_x,
                           addsamplestobackground = addsamplestobackground,
                           weights = weights, calculate_AIC = TRUE, AIC = AIC
    ),
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
    is_c <- detect_concave(model = m_aic, calib_data = data$calibration_data,
                           extrapolation_factor = extrapolation_factor,
                           var_limits = var_limits,
                           averages_from = averages_from, plot = FALSE)
    is_c <- any(sapply(is_c, function(x) x$is_concave))

    # # Check for concave curves (quadratic terms)
    # q_betas <- if (algorithm == "maxnet") {
    #   m_aic$betas[grepl("\\^2", names(m_aic$betas))]
    # } else {
    #   m_aic$coefficients[-1][grepl("\\^2", names(m_aic$coefficients[-1]))]
    # }
    #
    # is_c <- if (length(q_betas) == 0) FALSE else any(q_betas > 0)

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
        # Run maxnet model
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

      # Calculate omission rate and pROC
      om_rate <- omrat(threshold = omission_rate, pred_train = suit_val_cal,
                       pred_test = suit_val_eval)
      #Proc
      # proc_i <- lapply(omission_rate, function(omr){
      #   proc_omr <- enmpa::proc_enm(test_prediction = suit_val_eval,
      #                               prediction = pred_i,
      #                               threshold = omr, ...)$pROC_summary
      #   names(proc_omr) <- c(paste0("Mean_AUC_ratio_at_", omr),
      #                        paste0("pval_pROC_at_", omr))
      #   return(proc_omr)
      # })
      # proc_i <- unlist(proc_i)

      #Calculate PROC? ...
      if(proc_for_all){
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
                                         n_row = as.numeric(length(data$kfolds)),
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
  if (!return_replicate) eval_final <- NULL
  return(list(All_results = eval_final, Summary = eval_final_summary))
}
