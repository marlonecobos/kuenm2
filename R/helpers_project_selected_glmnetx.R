#Predict to several scenarion
multiple_projections <- function(i, res_path, raster_pattern, par_list) {
  res_i <- res_path[i,]
  input_i <- res_i$input_path
  output_i <- res_i$output_path
  #Import rasters
  r_i <- terra::rast(list.files(input_i, full.names = T,
                                pattern = raster_pattern))

  #Mask?
  if (!is.null(par_list$mask)) {
    #If in parallel, unwrap mask
    if (inherits(par_list$mask, "PackedSpatVector")) {
      par_list$mask <- terra::unwrap(par_list$mask)
    }
    r_i <- terra::crop(x = r_i, y = par_list$mask, mask = TRUE)
  }

  #Predict
  invisible(
    predict_selected(models = par_list$models,
                     new_variables = r_i,
                     write_partitions = par_list$write_partitions,
                     out_dir = output_i,
                     consensus_per_model = par_list$consensus_per_model,
                     consensus_general = par_list$consensus_general,
                     consensus = par_list$consensus,  # weighted mean
                     extrapolation_type = par_list$extrapolation_type,
                     var_to_clamp = par_list$var_to_clamp,
                     type = par_list$type,
                     overwrite = par_list$overwrite,
                     progress_bar = FALSE,
                     write_files = TRUE)
  )
  return(invisible(NULL))
}


#Binarize model
binarize_values <- function(x, v, new_value = 1) {
  x <- (x >= v) * new_value
  return(x)
}


#Calculate threshold
calc_thr <- function(occ_suitability, thr = 0.1) {
  sort(occ_suitability)[round(length(occ_suitability) * thr) + 1]
}
# data <- readRDS("../test_kuenm2/Myrcia.RDS")
# r <- rast("../test_kuenm2/Projection_results/Present/Brazil/General_consensus.tif")
# thr <- calc_thr(r = r, data = data, consensus = "median", thr = 0.1)
# thr


#Helper function to calculate variance coming from partitions by gcm
var_models_rep_by_gcm <- function(path) {
  model_files <- list.files(path = path, pattern = "partitions",
                            full.names = TRUE)
  if (length(model_files) == 0) {
    stop("Partitions not found.",
         "\nSet by_partition = FALSE or rerun project_selected() with write_partitions = TRUE.")
  }

  if (length(model_files) > 1) {
    r_x <- lapply(model_files, terra::rast)
    #Take mean of replicates (1-1, 2-2, 3-3, etc...)
    n_replicates <- terra::nlyr(r_x[[1]])
    mean_replicates <- terra::rast(lapply(1:n_replicates, function(n) {
      rep_n <- terra::mean(rast(lapply(r_x, function(x) x[[n]])))
    }))
    var_rep_x <- terra::app(mean_replicates, "var")} else {
      r_x <- terra::rast(model_files)
      var_rep_x <- terra::app(r_x, "var")}
  names(var_rep_x) <- basename(path)
  return(var_rep_x)
}

var_models_model_by_gcm <- function(path, consensus) {
  r_x <- terra::rast(list.files(path = path, pattern = "Model_.*consensus",
                         full.names = TRUE))
  r_x <- r_x[[sapply(r_x, function(r) names(r) == consensus)]]
  if (terra::nlyr(r_x) > 1) {
  var_x <- terra::app(r_x, "var") } else {
    var_x <- r_x * 0
  }
  return(var_x)
}


var_models_across_gcm <- function(paths, consensus) {
  # Read rasters
  r <- terra::rast(lapply(paths, function(x) {
    terra::rast(file.path(x, "General_consensus.tif"))[[consensus]]
  }))

  #plot(r)
  #Calculate variance
  v <- terra::app(r, "var")
  return(v)
}

var_models_across_ssp <- function(paths, consensus) {
  # Read rasters
  r <- terra::rast(lapply(paths, function(x) {
    terra::rast(file.path(x, "General_consensus.tif"))[[consensus]]
  }))

  #plot(r)
  #Calculate variance
  v <- terra::app(r, "var")
  return(v)
}


#### Check scenarios to predict and get dataframe with paths ####
check_pred_scenarios <- function(projection_data, out_dir) {
  #Check scenarios to predict
  sc <- names(projection_data[sapply(projection_data, function(x) !is.null(x))])

  #Get raster pattern to read
  raster_pattern <- projection_data$Raster_pattern

  #Get dataframe with path to predictions
  #Present
  if ("Present" %in% sc) {
    #Create folder
    present_dir <- file.path(out_dir, "Present/")
    present_sc <- names(projection_data[["Present"]])
    suppressWarnings({
      d_present <- data.frame(
        Time = "Present",
        Period = "Present",
        Scenario = present_sc,
        input_path = unlist(projection_data[["Present"]]),
        output_path = normalizePath(file.path(present_dir,
                                              present_sc))
      )
    })
  }
  #Past
  if ("Past" %in% sc) {
    #Create folder
    past_dir <- file.path(out_dir, "Past/")
    #Get grid of projections
    df_past <- do.call(
      rbind,
      lapply(names(projection_data$Past), function(time) {
        time_data <- projection_data$Past[[time]]
        do.call(rbind, lapply(names(time_data), function(gcm) {
          data.frame(Time = "Past", Period = time, GCM = gcm,
                     Path = time_data[[gcm]], stringsAsFactors = FALSE)
        }))
      })
    )

    #Looping in the grid

    #Create dataframe with path to results
    suppressWarnings({
      d_past <- data.frame(
        Time = "Past",
        Period = df_past$Period,
        GCM = df_past$GCM,
        input_path = df_past$Path,
        output_path = normalizePath(file.path(past_dir, df_past$Period,
                                              df_past$GCM),
                                    mustWork = FALSE)
      )
    })
  }
  #Future
  ####Project to Future scenarios####
  if ("Future" %in% sc) {
    #Create folder
    future_dir <- file.path(out_dir, "Future/")

    #Create grid of time-ssp-gcm
    df_future <- do.call(
      rbind,
      lapply(names(projection_data[["Future"]]), function(year_range) {
        year_range_data <- projection_data[["Future"]][[year_range]]
        do.call(rbind, lapply(names(year_range_data), function(ssp) {
          ssp_data <- year_range_data[[ssp]]
          do.call(rbind, lapply(names(ssp_data), function(gcm) {
            data.frame(Time = "Future", Period = year_range, ssp = ssp,
                       GCM = gcm, Path = ssp_data[[gcm]],
                       stringsAsFactors = FALSE)
          }))
        }))
      })
    )


    #Create dataframe with path to results
    suppressWarnings({
      d_future <- data.frame(
        Time = df_future$Time,
        Period = df_future$Period,
        Scenario = df_future$ssp,
        GCM = df_future$GCM,
        input_path = df_future$Path,
        output_path = normalizePath(file.path(future_dir, df_future$Period,
                                              df_future$ssp, df_future$GCM),
                                    mustWork = FALSE)
      )
    })
  }

  #Get dataframe with path to each projection
  if (!("Present" %in% sc)) {
    d_present <- NULL
  }
  if (!("Past" %in% sc)) {
    d_past <- NULL
  }
  if (!("Future" %in% sc)) {
    d_future <- NULL
  }

  #Return and write files with path
  res_path <- bind_rows_projection(list(d_present, d_past, d_future))
  #Create ID
  res_path$id <- 1:nrow(res_path)
  return(res_path)
}

#### Function to get cumulative predictions ####
cumulative_predictions <- function(predictions){
  original_order <- order(order(predictions)) #Get original order
  sorted_predictions <- sort(predictions) #Sort raw predictions
  cumulative_sum <- cumsum(sorted_predictions) #Cumulative sum
  max_cumulative_sum <- max(cumulative_sum, na.rm = TRUE) #Maximum cumulative value
  if (max_cumulative_sum == 0) {
    # Se todos os valores de entrada forem zero, o resultado também deve ser zero
    normalized_cumulative <- rep(0, length(predictions))
  } else {
    # 5. Normalizar a soma cumulativa para o intervalo de 0 a 100
    # Formula: (valor_atual / valor_maximo) * 100
    normalized_cumulative <- (cumulative_sum / max_cumulative_sum) * 100
  }
  return(normalized_cumulative[original_order])
}

helper_organize_proj <- function(r, mask, variable_names, fixed_variables,
                                 check_extent, resample_to_present,
                                 r_present =NULL, categorical_variables = NULL,
                                 file_name){

  #Mask variable, if necessary
  if(!is.null(mask)){
    r <- terra::crop(r, mask, mask = TRUE)
  }

  #Append fixed variables, if necessary
  if(!is.null(fixed_variables)){
    #Check if they have the same resolution
    if(any(terra::res(fixed_variables) != terra::res(r))){
      stop("Resolution of fixed_variables are different from the resolution of ",
file_name)
    }

    if(!is.null(mask)){
      fixed_variables <- terra::crop(fixed_variables, mask, mask = TRUE)
    }
    #Check extent
    if(check_extent){
      if(terra::ext(fixed_variables) != terra::ext(r)){
        terra::ext(fixed_variables) <- terra::ext(r)
      }
    }
    #Append variables
    r <- c(r, fixed_variables)
  }

  var_out <- setdiff(variable_names, names(r))
  if(length(var_out) > 0){
    stop("The following variable(s) is/are absent from the ", file_name, ":\n",
         paste(var_out, collapse = ", "))    }

  if(resample_to_present){
    if(is.null(categorical_variables)){
      r <- terra::resample(r, r_present, method = "bilinear")
    } else {
      r_cont <- terra::resample(r[[setdiff(names(r), categorical_variables)]],
                                r_present, method = "bilinear")
      r_cat <- terra::as.factor(r[[categorical_variables]])
      r_cat <- terra::resample(r_cat,
                               r_present, method = "near")
      r <- c(r_cont, r_cat)
    }
  }

  return(r)
}
