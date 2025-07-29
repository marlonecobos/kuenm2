#' Evaluate models with independent data
#'
#' @description
#' This function evaluates the selected models using independent data (i.e.,
#' data not used during model calibration). The function computes omission rate
#' and pROC, and optionally assesses whether the environmental conditions in the
#' independent data are analogous (i.e., within the range) to those in the
#' calibration data.
#'
#' @usage independent_evaluation(fitted_models, new_data,
#'                               consensus = c("mean", "median"),
#'                               type = "cloglog", extrapolation_type = "E",
#'                               var_to_clamp = NULL, perform_mop = TRUE,
#'                               mop_type = "detailed",
#'                               calculate_distance = TRUE,
#'                               where_distance = "all",
#'                               return_predictions = TRUE,
#'                               return_binary = TRUE,
#'                               progress_bar = FALSE, ...)
#'
#' @param fitted_models an object of class `fitted_models` returned by the
#' \code{\link{fit_selected}}() function.
#' @param new_data a `data.frame` containing environmental variables for
#' independent test records. The column names must correspond exactly to the
#' environmental variables used to fit the selected models, and each row to an
#' individual test record.
#' @param consensus (character) vector specifying the types of consensus to
#' use. Available options are `"median"` and `"mean"`. Default is
#' `c("median", "mean")`.
#' @param type (character) the format of prediction values. For `maxnet` models,
#' valid options are `"raw"`, `"cumulative"`, `"logistic"`, and `"cloglog"`. For
#' `glm` models, valid options are `"response"` and `"cloglog"`. Default is
#' `"cloglog"`.
#' @param extrapolation_type (character) extrapolation type of model. Models can
#' be transferred with three options: free extrapolation ('E'), extrapolation
#' with clamping ('EC'), and no extrapolation ('NE'). Default = 'E'. See details.
#' @param var_to_clamp (character) vector specifying which variables to clamp or
#' not extrapolate. Only applicable if extrapolation_type is "EC" or "NE".
#' Default is `NULL`, meaning all variables will be clamped or not extrapolated.
#' @param perform_mop (logical) whether to execute a Mobility-Oriented Parity
#' (MOP) analysis. This analysis assesses if the environmental conditions in the
#' `new_data` are analogous (within ranges) to those in the calibration data.
#' Defaults to `TRUE`.
#' @param mop_type (character) type of MOP analysis to be performed. Options
#' available are "basic", "simple" and "detailed". Default is 'simples'. See
#' \code{\link{projection_mop}}() for more details.
#' @param calculate_distance (logical) whether to calculate distances
#' (dissimilarities) between new_data and calibration data. Default is TRUE.
#' @param where_distance (character) specifies which values in `new_data` should
#' be used to calculate distances. Options are: "in_range" (only conditions
#' within the calibration range), "out_range" (only conditions outside the
#' calibration range), and "all" (all conditions). Default is "all".
#' @param return_predictions (logical) whether to return continuous predictions
#' at the locations of independent records in `new_data`. Default is TRUE.
#' @param return_binary (logical) whether to return binary predictions
#' at the locations of independent records in `new_data`. The predictions are
#' binarized using the respective thresholds stores in `fitted_models`. Default
#' is TRUE.
#' @param progress_bar (logical) whether to display a progress bar during
#' mop processing. Default is FALSE.
#' @param ... additional arguments passed to \code{\link[mop]{mop}()}.
#'
#' @importFrom mop mop
#' @importFrom fpROC auc_metrics
#' @importFrom stats setNames
#' @return
#' A list containing the following elements:
#'
#' - **evaluation**: A `data.frame` with omission rate and pROC values for each
#' selected model and for the consensus.
#' - **mop_results**: (Only if `perform_mop = TRUE`) An object of class
#' `mop_results`, with metrics of environmental similarity between calibration
#' and independent data.
#' - **predictions**: (Only if `return_predictions = TRUE`) A `list` of
#' `data.frames` containing continuous and binary predictions at the independent
#' record locations, along with MOP distances, an indicator of whether
#' environmental conditions at each location fall within the calibration range,
#' and the identity of the variables that have lower and higher values than the
#' calibration range. If the `fitted_models` object includes categorical
#' variables, the returned `data.frame` will also contain columns indicating
#' which values in `new_data` were not present in the calibration data.
#' @examples
#' # Example with maxnet
#' # Import example of fitted_models (output of fit_selected())
#' data("fitted_model_maxnet", package = "kuenm2")
#'
#' # Import independent records to evaluate the models
#' data("new_occ", package = "kuenm2")
#'
#' # Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' #Extract variables to occurrences
#' new_data <- extract_occurrence_variables(occ = new_occ, x = "x", y = "y",
#'                                          raster_variables = var)
#'
#' #Add some fake data beyond the limits of calibration ranges
#' fake_data <- data.frame("pr_bg" = c(1, 1, 1),
#'                         "x" = c(NA, NA, NA),
#'                         "y" = c(NA, NA, NA),
#'                         "bio_1" = c(10, 15, 23),
#'                         "bio_7" = c(12, 16, 20),
#'                         "bio_12" = c(2300, 2000, 1000),
#'                         "bio_15" = c(30, 40, 50),
#'                         "SoilType" = c(1, 1, 1))
#' new_data <- rbind(new_data, fake_data)
#'
#'
#' # Evaluate models with independent data
#' res_ind <- independent_evaluation(fitted_models = fitted_model_maxnet,
#'                                   new_data = new_data)
#'
#' @export
independent_evaluation <- function(fitted_models, new_data,
                                   consensus = c("mean", "median"),
                                   type = "cloglog",
                                   extrapolation_type = "E",
                                   var_to_clamp = NULL,
                                   perform_mop = TRUE,
                                   mop_type = "detailed",
                                   calculate_distance = TRUE,
                                   where_distance = "all",
                                   return_predictions = TRUE,
                                   return_binary = TRUE,
                                   progress_bar = FALSE,
                                   ...){
  #### Check data ####
  if (missing(fitted_models)) {
    stop("Argument 'fitted_models' must be defined.")
  }
  if (!inherits(fitted_models, "fitted_models")) {
    stop("Argument 'fitted_models' must be a 'fitted_models' object.")
  }
  if (missing(new_data)) {
    stop("Argument 'new_data' must be defined.")
  }
  if (!inherits(new_data, "data.frame")) {
    stop("Argument 'new_data' must be a 'data.frame' object.")
  }

  if (!inherits(consensus, "character")) {
    stop("Argument 'consensus' must be a 'character'.")
  }
  consensus_out <- setdiff(consensus, c("median", "mean"))
  if (length(consensus_out) > 0) {
    stop("Invalid 'consensus' provided.",
         "\nAvailable options are 'median' and 'mean'.")
  }

  if(fitted_models$algorithm == "maxnet"){
    if (!any(c("raw", "cumulative", "logistic", "cloglog") %in% type)) {
      stop("Invalid 'type' provided.",
           "\nAvailable options for maxnet fitted_models are: 'raw', 'cumulative',
           'logistic', or 'cloglog'.")
    }
    if(type == "raw")
      type <-  "exponential"
  }

  if(fitted_models$algorithm == "glm"){
    if (!any(c("response", "cloglog") %in% type)) {
      stop("Invalid 'type' provided.",
           "\nAvailable options for glm fitted_models are 'response' or 'cloglog'.")
    }
    if(type == "cloglog")
      type = "link"
  }

  if(length(extrapolation_type) > 1){
    stop("Extrapolation type accepts only one of these values: 'E', 'EC', or
         'NE'")
  }

  extrapolation_out <- setdiff(extrapolation_type, c("E", "EC", "NE"))
  if (length(extrapolation_out) > 0) {
    stop("Invalid 'extrapolation type' provided.",
         "\nAvailable options are: 'E', 'EC', and 'NE'.")
  }

  if (extrapolation_type %in% c("E", "EC") & !is.null(var_to_clamp) &
      !inherits(var_to_clamp, "character")) {
    stop("Argument 'var_to_clamp' must be NULL or 'character'.")
  }


  if(!inherits(perform_mop, "logical")){
    stop("Argument 'perform_mop' must be 'logical'")
  }

  if (!inherits(mop_type, "character")) {
    stop("Argument 'mop_type' must be a 'character'.")
  }

  mop_type_out <- setdiff(mop_type, c("basic", "simple", "detailed"))
  if (length(mop_type_out) > 0) {
    stop("Invalid 'mop_type' provided.",
         "\nAvailable options are: 'basic', 'simple', or 'detailed'.")
  }

  if(!inherits(calculate_distance, "logical")){
    stop("Argument 'calculate_distance' must be 'logical'.")
  }

  distance_out <- setdiff(where_distance, c("in_range", "out_range", "all"))
  if (length(distance_out) > 0) {
    stop("Invalid 'where_distance' provided.",
         "\nAvailable options are: 'in_range', 'out_range', and 'all'.")
  }

  if(!inherits(return_predictions, "logical")){
    stop("Argument 'return_predictions' must be 'logical'.")
  }

  if(!inherits(return_binary, "logical")){
    stop("Argument 'return_binary' must be 'logical'.")
  }

  if(!inherits(progress_bar, "logical")){
    stop("Argument 'progress_bar' must be 'logical'.")
  }

  #Check variables
  v <- unique(unlist(sapply(fitted_models$Models, function(x)
    names(x$Full_model$betas)[-1],
    simplify = F)))
  v <- gsub("I\\((.*?)\\^2\\)", "\\1", v) #Remove quadratic pattern
  v <- v[!grepl("categorical", v)] #Remove categorical pattern
  v <- unique(unlist(strsplit(v, ":"))) #Remove product pattern

  #Check variables absent from new_data
  v_out <- setdiff(v, colnames(new_data))
  if(length(v_out) > 0){
    stop("The following variables used to fit the models are absent from the 'new_data:\n'", paste(v_out, collapse = "; "))
  }

  #Predict to independent records
  pred_test <- predict_selected(models = fitted_models,
                                raster_variables = new_data[,v],
                                consensus = consensus,
                                extrapolation_type = extrapolation_type,
                                var_to_clamp = var_to_clamp,
                                type = type,
                                progress_bar = FALSE)
  #Save names
  nm <- names(pred_test)

  #Predict to background
  bg_data <- fitted_models$calibration_data
  pred_bg <- predict_selected(models = fitted_models,
                              raster_variables = bg_data[,v],
                              consensus = consensus,
                              extrapolation_type = extrapolation_type,
                              var_to_clamp = var_to_clamp,
                              type = type,
                              progress_bar = FALSE)

  #Get only consensus predictions
  pred_test <- lapply(names(pred_test), function(i){
    if(i == "General_consensus"){
      return(pred_test[[i]])
    } else {
      return(pred_test[[i]]$Model_consensus)
    }
  })
  names(pred_test) <- nm

  #Get thresholds
  thr <- lapply(names(fitted_models$thresholds), function(i){
    if(i != "General_consensus"){
      return(fitted_models$thresholds[[i]])
    } else {
      return(fitted_models$thresholds$consensus)
    }
  })
  names(thr) <- names(fitted_models$thresholds)
  names(thr)[names(thr) == "consensus"] <- "General_consensus"

  res <- lapply(names(pred_test), function(i){
    #print(i)
    #Get pred test i
    p_i <- pred_test[[i]]

    #Get consensus predictions
    if(i != "General_consensus"){
        #p_i <- p_i$Model_consensus
        bg_i <- pred_bg[[i]]$Model_consensus
    } else {
      bg_i <- pred_bg[[i]]
    }

    #Calculate omission rate
    omr_i <- sapply(names(p_i), function(x){
      sum(p_i[[x]] < thr[[i]][[x]])/length(p_i[[x]])
    })

    #Calculate proc
    proc_i <- lapply(names(p_i), function(x){
      fpROC::auc_metrics(test_prediction = p_i[[x]],
                         prediction = bg_i[[x]],
                         threshold = fitted_models$omission_rate)$summary[, 4:5]
    })
    names(proc_i) <- names(p_i$Model_consensus)
    proc_i <- as.data.frame(do.call(rbind, proc_i))

    #Save results
    df_i <- data.frame(Model = i,
                       consensus = names(p_i),
                       omr = omr_i)
    df_i <- cbind(df_i, proc_i)
    colnames(df_i)[3] <- paste0("Omission_rate_at_", fitted_models$omission_rate)
    row.names(df_i) <- NULL
    return(df_i)
  })
  res <- do.call(rbind, res)

  if(perform_mop){
    mop_res <- mop_with_records(train_data = bg_data,
                                test_data = new_data,
                                variables = v,
                                categorical_variables = fitted_models$categorical_variables,
                                mop_type = mop_type,
                                calculate_distance = calculate_distance,
                                where_distance = where_distance,
                                progress_bar = progress_bar, ...)
    } else {
    mop_res <- NULL}

  if(return_predictions){
    predictions <- list()
    predictions[["continuous"]] <- do.call(cbind, pred_test)

    if(return_binary){
      # Initializing a new list to store the binarized data
      pred_test_binarized <- list()

      # Iterating over each model/consensus and applying binarization
      for (model_name in names(pred_test)) {
        # Initialize a sublist for the current model within pred_test_binarized
        pred_test_binarized[[model_name]] <- list()
        # Get the pred_test dataframe for the current model
        current_pred_data <- pred_test[[model_name]]
        # Get the threshold values for the current model
        current_thr_data <- thr[[model_name]]
        # Iterate over the 'mean' and 'median' metrics
        for (metric in consensus) {
          binarized_values <- as.numeric(current_pred_data[[metric]] >= current_thr_data[[metric]])
          # Add the binarized values to the model's sublist
          pred_test_binarized[[model_name]][[metric]] <- binarized_values
        }
        # Convert the model's sublist to a dataframe to maintain the original structure
        pred_test_binarized[[model_name]] <- as.data.frame(pred_test_binarized[[model_name]])
      }

      predictions$binary <- do.call(cbind, pred_test_binarized)

    }

  } else { #End of return_predictions
    predictions <- NULL
  }


  #Append mop results to predictions
  if(perform_mop && return_predictions){
    for(pred_type in names(predictions)){
      predictions[[pred_type]] <- cbind(predictions[[pred_type]], mop_res$mop_records)
    }
  } else if (perform_mop && !return_predictions){
    predictions <- mop_res$mop_records
  }


  #Final results
  final_res <- list("evaluation" = res,
                    "mop_results" = mop_res$mop_results,
                    "predictions" = predictions)

  return(final_res)
}


#Consensus = mean
#Model 192: Higher omission rate (0.4, > 0.1) and significant pROC value
#Model 189: Higher omission rate (0.4, > 0.1) and significant pROC value
#General consensus: Higher omission rate (0.4, > 0.1) and significant pROC value

#Mop summary
#All independent data are within the ranges of calibration data.
#The mean distance was, varying between and




