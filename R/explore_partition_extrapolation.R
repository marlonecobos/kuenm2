#' Analysis of extrapolation risks in partitions using the MOP metric
#'
#' @description
#' This function calculates environmental dissimilarities and identifies
#' non-analogous conditions by comparing the training data against the test data
#' for each partition, using the MOP (Mobility-Oriented Parity) metric.
#'
#' @usage explore_partition_extrapolation(data, include_train_background = TRUE,
#'                                        include_test_background = FALSE,
#'                                        variables = NULL,
#'                                        mop_type = "detailed",
#'                                        calculate_distance = TRUE,
#'                                        where_distance = "all",
#'                                        return_spatial = TRUE,
#'                                        crs = "+init=epsg:4326",
#'                                        progress_bar = FALSE, ...)
#'
#' @param data an object of class `prepared_data` returned by the
#' `prepare_data()` function.
#' @param include_train_background (logical) whether to include the background
#' points used in training to define the environmental range of the training
#' data. If set to FALSE, only the environmental conditions of the training
#' presence records will be considered. Default is TRUE, meaning both presence
#' and background points are used.
#' @param include_test_background (logical) whether to compute MOP for both the
#' test presence records and the background points not used during training.
#' Default is FALSE, meaning MOP will be calculated only for the test presences.
#' @param variables (character) names of the variables to be used in the MOP
#' calculation. Default is NULL, meaning all variables in `data` will be used.
#' @param mop_type (character) type of MOP analysis to be performed. Options
#' available are "basic", "simple" and "detailed". Default is 'simples'. See
#' \code{\link{projection_mop}}() for more details.
#' @param calculate_distance (logical) whether to calculate distances
#' (dissimilarities) between train and test data. Default is TRUE.
#' @param where_distance (character) specifies which values in train data should
#' be used to calculate distances. Options are: "in_range" (only conditions
#' within the train range), "out_range" (only conditions outside the
#' train range), and "all" (all conditions). Default is "all".
#' @param return_spatial (logical) whether to return a `SpatVector` showing the
#' dissimilarities and the spatial distribution of test data that falls within
#' and outside the range of the training data. Default is TRUE.
#' @param crs The coordinate reference system to spatialize the results. Only
#' applicable if `return_spatial = TRUE`. Default is "+init=epsg:4326".
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is FALSE.
#' @param ... additional arguments passed to \code{\link[mop]{mop}()}.
#'
#' @return
#' A `data.frame` containing:
#' - MOP distances (if `calculate_distance = TRUE`);
#' - an indicator of whether environmental conditions at each test record fall
#' within the training range;
#' - the number of variables outside the training range;
#' - the names of variables with values lower or higher than the training range;
#' - if the `prepared_data` object includes categorical variables, it will also
#' contain columns indicating which values in the testing data were not present
#' in the training data.
#'
#' If `return_raster_result = TRUE`, it also returns a `SpatRaster` showing the
#' spatial distribution of test data that falls within and outside the range of
#' the training data.
#' @export
#' @importFrom terra vect crs rasterize mask coltab
#' @importFrom mop mop
#' @importFrom stats setNames
#' @examples
#' #Prepare data
#' # Import occurrences
#' data(occ_data, package = "kuenm2")
#'
#' # Import raster layers
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' # Prepare data for maxnet model
#' sp_swd <- prepare_data(algorithm = "maxnet", occ = occ_data,
#'                        x = "x", y = "y",
#'                        raster_variables = var,
#'                        species = occ_data[1, 1],
#'                        n_background = 100,
#'                        categorical_variables = "SoilType",
#'                        features = c("l", "lq"),
#'                        r_multiplier = 1,
#'                        partition_method = "kfolds")
#'
#' # Analysis of extrapolation risks in partitions
#' res <- explore_partition_extrapolation(data = sp_swd,
#'                                        raster_variables = var,
#'                                        include_test_background = TRUE)
#' #Plot spatial spatial distribution of test data
#' terra::plot(res$Spatial_results)
#'
explore_partition_extrapolation <- function(data,
                                            include_train_background = TRUE,
                                            include_test_background = FALSE,
                                            variables = NULL,
                                            mop_type = "detailed",
                                            calculate_distance = TRUE,
                                            where_distance = "all",
                                            return_spatial = TRUE,
                                            crs = "+init=epsg:4326",
                                            progress_bar = FALSE, ...){
  #Check data
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (!inherits(data, "prepared_data")) {
    stop("'data' must be a 'prepared_data' object.")
  }
  if (!inherits(include_train_background, "logical")) {
    stop("'include_train_background' must be 'logical'.")
  }
  if (!inherits(include_test_background, "logical")) {
    stop("'include_test_background' must be 'logical'.")
  }
  if (!is.null(variables)) {
    if(!inherits(variables, "logical")) {
    stop("'variables' must be 'character' or NULL.")}
    var_out <- setdiff(variables, colnames(data$calibration_data))
    if(length(var_out) > 0){
      stop("Some 'variables' provided are absent from data")
    }
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

  if(!inherits(return_spatial, "logical")){
    stop("Argument 'return_spatial' must be 'logical'.")
  }


  if(!inherits(progress_bar, "logical")){
    stop("Argument 'progress_bar' must be 'logical'.")
  }

  #Remove categorical if necessary
  if(!is.null(data$categorical_variables)){
    v <- setdiff(colnames(data$calibration_data), c("pr_bg",
                                                    data$categorical_variables))
  } else {
    v <- setdiff(colnames(data$calibration_data),
                 data$categorical_variables)
  }

  #Remove variables, if necessary
  if(!is.null(variables)){
    v <- intersect(variables, v)
  }

  #Create list to save partition results
  res <- list()

  #Get number of partitions
  n_partitions <- length(data$part_data)

  if (progress_bar) {
    pb <- txtProgressBar(min = 0, max = n_partitions, style = 3)
    progress <- function(n) {
      setTxtProgressBar(pb, n)
    }
    opts <- list(progress = progress)
  } else {
    opts <- NULL
  }

  for(i in 1:n_partitions){
    #Get partition name
    i_name <- names(data$part_data[i])

    #Get test and train
    partition_i <- data$part_data[[i_name]]
    test_i <- data$calibration_data[partition_i, ]
    train_i <- data$calibration_data[-partition_i, ]

    #Exclude background, if necessary
    if(!include_train_background){
      train_i <- train_i[train_i$pr_bg == 1, ]
    }
    if(!include_test_background){
      test_i <- test_i[test_i$pr_bg == 1, ]
    }

    #Run analysis
    res_i <- mop_with_records(train_data = train_i, test_data = test_i,
                              variables = v,
                              categorical_variables = data$categorical_variables,
                              mop_type = mop_type,
                              calculate_distance = calculate_distance,
                              where_distance = where_distance,
                              progress_bar = FALSE, ...)

    res[[i_name]][["test_data"]] <- res_i$mop_records
    #Append pr_bg information
    if(include_test_background){
      res[[i_name]][["test_data"]]$pr_bg <- data$calibration_data[partition_i,
                                                                  "pr_bg"]
    } else {
      res[[i_name]][["test_data"]]$pr_bg <- 1
    }


    #Append xy information
    if(return_spatial){
      if(!include_test_background){
        xy <- cbind(data$data_xy[partition_i, ],
                    "pr_bg" = data$calibration_data[partition_i, "pr_bg"])
        xy <- xy[xy$pr_bg == 1, c("x", "y")]
      } else {
        xy <- data$data_xy[partition_i, ]
      }

    res[[i_name]][["test_data"]]<- cbind(xy,
                                         res[[i_name]][["test_data"]])}

    # #Append index, for test
    # res[[i_name]][["test_data"]]$index <- partition_i

    # Sets the progress bar to the current state
    if (progress_bar) setTxtProgressBar(pb, i)

  } #End of for i

  #Summarize data
  all_res <- lapply(names(res), function(i) {
    cbind("Partition" = i, res[[i]]$test_data)
  })
  all_res <- do.call("rbind", all_res)

  #Organize columns
  c_order <- intersect(c("Partition", "pr_bg", "x", "y"), colnames(all_res))
  c_order <- c(c_order, setdiff(colnames(all_res), c_order))
  all_res <- all_res[, c_order]

  # #Spatialize results
  # if(return_spatial){
  #
  # res_spatial <- list()

  # #Categories and colors
  # l <- data.frame(id = 0:5,
  #                 category = c("Calibration area",
  #                              "Test data within range",
  #                              "Test data out of range",
  #                              "Background test within range",
  #                              "Background test out of range",
  #                              "Train data"))
  # cores <- data.frame(value = 0:5, col = c("grey90",
  #                                          "#009E73",
  #                                          "#D55E00",
  #                                          "#0072B2",
  #                                          "#F0E442",
  #                                          "black"))


  # for (i in names(res)) {
  #   #print(i)
  #   #Filter results
  #   rep_i <- all_res[all_res$Partition == i, ]
  #
  #   #Add attributed
  #   rep_i$category <- NA # initialize column
  #
  #   # Set categories
  #   rep_i$category[rep_i$pr_bg == 1 & rep_i$inside_range] <- "Test presence within range"
  #   rep_i$category[rep_i$pr_bg == 1 & !rep_i$inside_range] <- "Test presence out of range"
  #   rep_i$category[rep_i$pr_bg == 0 & rep_i$inside_range] <- "Test background within range"
  #   rep_i$category[rep_i$pr_bg == 0 & !rep_i$inside_range] <- "Test background out of range"
  #
  #   #Get only columns of interest
  #   rep_i <- rep_i[,c("x", "y", "category", "mop_distance")]
  #
  #   # #Append train points
  #   # train_index <- -data$part_data[[i]]
  #   # rep_train_i <- cbind(data$data_xy[train_index,],
  #   #                      "pr_bg" = data$calibration_data[train_index, "pr_bg"])
  #   # rep_train_i <- rep_train_i[rep_train_i$pr_bg == 1,]
  #   # rep_train_i$category <- 5
  #   # #Merge results
  #   # rep_i<- rbind(rep_i, rep_train_i[,colnames(rep_i)])
  #
  #   #Spatialize
  #   rep_vect <- terra::vect(rep_i, geom = c(x = "x", y = "y"),
  #                           crs = crs)
  #
  #
  #   # #Get template raster
  #   # template_raster <- raster_variables[[1]]
  #   # template_raster[!is.na(template_raster)] <- 0
  #   #
  #   # #Vectorize raster
  #   # spatial_points <- terra::vect(rep_i[,c("x", "y", "category")],
  #   #                               geom=c("x", "y"),
  #   #                               crs= terra::crs(raster_variables))
  #   # #Rasterize data
  #   # rasterized_data <- suppressWarnings(terra::rasterize(spatial_points,
  #   #                                                      template_raster,
  #   #                                                      field="category",
  #   #                                                      fun="min",
  #   #                                                      background = 0,
  #   #                                                      na.rm = TRUE))
  #   # #Identify calibration area
  #   # rasterized_data <- terra::mask(rasterized_data, template_raster)
  #   #
  #   # #Set levels
  #   # levels(rasterized_data) <- l
  #   #
  #   # #Set colors
  #   # terra::coltab(rasterized_data) <- cores
  #   res_spatial[[i]] <- rep_vect
  # }
  # } #End of return_spatial
  #

  #Returnin final results
  res_final <- list("Mop_results" = all_res,
                    "calibration_data" = data$calibration_data,
                    "categorical_variables" = data$categorical_variables,
                    "pca" = data$pca)
  class(res_final) <- "explore_partition"
  return(res_final)
  } #End of function





