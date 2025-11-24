#' Analysis of extrapolation risks in projections using the MOP metric
#'
#' @description
#' Calculates the mobility-oriented parity metric and other sub-products to
#' represent dissimilarities and non-analogous conditions when comparing a set
#' of reference conditions (M) against model projection conditions (G).
#'
#' @usage
#' projection_mop(data, projection_data, out_dir,
#'                subset_variables = FALSE, mask = NULL, type = "basic",
#'                na_in_range = TRUE, calculate_distance = FALSE,
#'                where_distance = "in_range", distance = "euclidean",
#'                scale = FALSE, center = FALSE, fix_NA = TRUE, percentage = 1,
#'                comp_each = 2000, tol = NULL, rescale_distance = FALSE,
#'                parallel = FALSE, ncores = NULL, progress_bar = TRUE,
#'                overwrite = FALSE)
#'
#' @param data an object of class `fitted_models` returned by the
#' [fit_selected()] function or an object of class `prepared_data` returned by
#' the [prepare_data()] function.
#' @param projection_data an object of class `projection_data` returned by the
#' [prepare_projection()]function. This file contains the paths to
#' the rasters representing each scenario.
#' @param out_dir (character) a path to a root directory for saving the raster
#' file of each projection.
#' @param subset_variables (logical) whether to include in the analysis only the
#' variables present in the selected models. Default is FALSE
#' @param mask (SpatRaster, SpatVector, or SpatExtent) spatial object used to
#' mask the variables (optional). Default is NULL.
#' @param type (character) type of MOP analysis to be performed. Options
#' available are "basic", "simple" and "detailed". See Details for further
#' information.
#' @param na_in_range (logical) whether to assign `NA` to regions within the
#' projected area (G) where environmental conditions fall within the range of
#' the calibration data (M). If `TRUE` (default), these regions are assigned
#' `NA`. If `FALSE`, they are assigned `0` in the simple and basic MOP outputs,
#' and `"within ranges"` in the detailed MOP output.
#' @param calculate_distance (logical) whether to calculate distances
#' (dissimilarities) between m and g. The default, FALSE, runs rapidly and does
#' not assess dissimilarity levels.
#' @param where_distance (character) where to calculate distances, considering
#' how conditions in g are positioned in comparison to the range of conditions
#' in m. Options available are "in_range", "out_range" and "all". Default is
#' "in_range".
#' @param distance (character) which distances are calculated, euclidean or
#' mahalanobis. Only applicable if calculate_distance = TRUE.
#' @param scale (logical or numeric) whether to scale as in
#' \code{\link[base]{scale}}. Default is FALSE.
#' @param center (logical or numeric) whether to center as in
#' \code{\link[base]{scale}}. Default is FALSE.
#' @param fix_NA (logical) whether to fix layers so cells with NA values are
#' the same in all layers. Setting to FALSE may save time if the rasters are
#' big and have no NA matching problems. Default is TRUE.
#' @param percentage (numeric) percentage of \code{m} closest conditions used
#' to derive mean environmental distances to each combination of conditions in
#' \code{g}.
#' @param comp_each (numeric) number of combinations in \code{g} to be used for
#' distance calculations at a time. Increasing this number requires more RAM
#' @param tol (numeric) tolerance to detect linear dependencies when calculating
#' Mahalanobis distances. The default, NULL, uses `.Machine$double.eps`.
#' @param rescale_distance (logical) whether to re-scale distances 0-1.
#' Re-scaling prevents comparisons of dissimilarity values obtained from runs
#' with different values of \code{percentage}.
#' @param parallel (logical) whether to fit the candidate models in parallel.
#' Default is FALSE.
#' @param ncores (numeric) number of cores to use for parallel processing.
#' Default is NULL and uses available cores - 1. This is only applicable if
#' `parallel = TRUE`.
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is TRUE.
#' @param overwrite (logical) whether to overwrite SpatRaster if they already
#' exists. Only applicable if `write_files` is set to TRUE. Default is FALSE.
#'
#' @details
#' `type` options return results that differ in the detail of how non-analogous
#' conditions are identified.
#' - **basic** - makes calculation as proposed by Owens et al. (2013)
#' <doi:10.1016/j.ecolmodel.2013.04.011>.
#' - **simple** - calculates how many variables in the set of interest are
#' non-analogous to those in the reference set.
#' - **detailed** - calculates five additional extrapolation metrics. See
#' `mop_detailed` under `Value` below for full details.
#'
#' `where_distance` options determine what values should be used to calculate
#' dissimilarity
#' - **in_range** - only conditions inside `m` ranges
#' - **out_range** - only conditions outside `m` ranges
#' - **all** - all conditions
#'
#' When the variables used to represent conditions have different units,
#' scaling and centering are recommended. This step is only valid when Euclidean
#' distances are used.
#'
#' @return
#' An object of class `mop_projections`, with the root directory and the dataframe
#' containing the file paths where the results were stored for each scenario.
#' The paths contain the following files:
#' - **summary** - a data.frame with details of the data used in the analysis:
#'     - *variables* - names of variables considered.
#'     - *type* - type of MOP analysis performed.
#'     - *scale* - value according to the argument \code{scale}.
#'     - *center* - value according to the argument \code{center}.
#'     - *calculate_distance* - value according to the argument
#'     \code{calculate_distance}.
#'     - *distance* - option regarding distance used.
#'     - *percentage* - percentage of \code{m} used as reference for
#'     distance calculation.
#'     - *rescale_distance* - value according to the argument
#'     \code{rescale_distance}.
#'     - *fix_NA* - value according to the argument \code{fix_NA}.
#'     - *N_m* - total number of elements (cells with values or valid
#'     rows) in \code{m}.
#'     - *N_g* - total number of elements (cells with values or valid
#'     rows) in \code{g}.
#'     - *m_min* - minimum values (lower limit) of the variables in reference
#'     conditions (\code{m}).
#'     - *m_max* - maximum values (upper limit) of the variables in reference
#'     conditions (\code{m}).
#' - **mop_distances** - if \code{calculate_distance} = TRUE, a SpatRaster or
#' vector with distance values for the set of interest (\code{g}). Higher values
#' represent greater dissimilarity compared to the set of reference (\code{m}).
#' - **mop_basic** - a SpatRaster or vector, for the set of interest,
#' representing conditions in which at least one of the variables is
#' non-analogous to the set of reference. Values should be: 1 for non-analogous
#' conditions, and NA for conditions inside the ranges of the reference set.
#' - **mop_simple** - a SpatRaster or vector, for the set of interest,
#' representing how many variables in the set of interest are non-analogous to
#' those in the reference set. NA is used for conditions inside the ranges of
#' the reference set.
#' - **mop_detailed** - a list containing:
#'     - *interpretation_combined* - a data.frame to help identify combinations
#'     of variables in *towards_low_combined* and *towards_high_combined* that
#'     are non-analogous to \code{m}.
#'     - *towards_low_end* - a SpatRaster or matrix for all variables
#'     representing where non-analogous conditions were found towards low values
#'     of each variable.
#'     - *towards_high_end* - a SpatRaster or matrix for all variables
#'     representing where non-analogous conditions were found towards high
#'     values of each variable.
#'     - *towards_low_combined* - a SpatRaster or vector with values
#'     representing the identity of the variables found to have non-analogous
#'     conditions towards low values. If vector, interpretation requires the use
#'     of the data.frame *interpretation_combined*.
#'     - *towards_high_combined* - a SpatRaster or vector with values
#'     representing the identity of the variables found to have non-analogous
#'     conditions towards high values. If vector, interpretation requires the
#'     use of the data.frame *interpretation_combined*.
#'
#' @export
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom terra rast predict writeRaster mask
#' @importFrom mop mop
#'
#' @seealso
#' [organize_future_worldclim()], [prepare_projection()]
#'
#' @examples
#' # Step 1: Organize variables for current projection
#' ## Import current variables (used to fit models)
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' ## Create a folder in a temporary directory to copy the variables
#' out_dir_current <- file.path(tempdir(), "Current_raw4")
#' dir.create(out_dir_current, recursive = TRUE)
#'
#' ## Save current variables in temporary directory
#' terra::writeRaster(var, file.path(out_dir_current, "Variables.tif"))
#'
#'
#' # Step 2: Organize future climate variables (example with WorldClim)
#' ## Directory containing the downloaded future climate variables (example)
#' in_dir <- system.file("extdata", package = "kuenm2")
#'
#' ## Create a folder in a temporary directory to copy the future variables
#' out_dir_future <- file.path(tempdir(), "Future_raw4")
#'
#' ## Organize and rename the future climate data (structured by year and GCM)
#' ### 'SoilType' will be appended as a static variable in each scenario
#' organize_future_worldclim(input_dir = in_dir, output_dir = out_dir_future,
#'                           name_format = "bio_", fixed_variables = var$SoilType)
#'
#' # Step 3: Prepare data to run multiple projections
#' ## An example with maxnet models
#' ## Import example of fitted_models (output of fit_selected())
#' data(fitted_model_maxnet, package = "kuenm2")
#'
#' ## Prepare projection data using fitted models to check variables
#' pr <- prepare_projection(models = fitted_model_maxnet,
#'                          present_dir = out_dir_current,
#'                          future_dir = out_dir_future,
#'                          future_period = c("2041-2060", "2081-2100"),
#'                          future_pscen = c("ssp126", "ssp585"),
#'                          future_gcm = c("ACCESS-CM2", "MIROC6"),
#'                          raster_pattern = ".tif*")
#'
#' # Step 4: Perform MOP for all projection scenarios
#' ## Create a folder to save MOP results
#' out_dir <- file.path(tempdir(), "MOPresults")
#' dir.create(out_dir, recursive = TRUE)
#'
#' ## Run MOP
#' kmop <- projection_mop(data = fitted_model_maxnet, projection_data = pr,
#'                        out_dir = out_dir, type = "detailed")

projection_mop <- function(data,
                           projection_data,
                           out_dir,
                           subset_variables = FALSE,
                           mask = NULL,
                           type = "basic",
                           na_in_range = TRUE,
                           calculate_distance = FALSE,
                           where_distance = "all",
                           distance = "euclidean",
                           scale = FALSE,
                           center = FALSE,
                           fix_NA = TRUE,
                           percentage = 1,
                           comp_each = 2000,
                           tol = NULL,
                           rescale_distance = FALSE,
                           parallel = FALSE,
                           ncores = NULL,
                           progress_bar = TRUE,
                           overwrite = FALSE) {
  #Check data
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (missing(projection_data)) {
    stop("Argument 'projection_data' must be defined.")
  }
  if (missing(out_dir)) {
    stop("Argument 'out_dir' must be defined.")
  }

  if (!inherits(data, "prepared_data") && !inherits(data, "fitted_models")) {
    stop("Argument 'data' must be a 'prepared_data' or 'fitted_models' object.")
  }
  if (!inherits(projection_data, "projection_data")) {
    stop("Argument 'projection_data' must be a 'projection_data' object.")
  }
  if (!inherits(out_dir, "character")) {
    stop("Argument 'out_dir' must be a 'character'.")
  }

  if (subset_variables & !inherits(data, "fitted_models")) {
    stop("If 'subset_variables' = TRUE, 'data' must be 'fitted_models' object.")
  }
  if (!is.null(mask) & !inherits(mask, c("SpatVector", "SpatVector", "SpatExtent"))) {
    stop("Argument 'mask' must be a 'SpatVector', 'SpatVector', or 'SpatExtent'.")
  }

  if (!inherits(type, "character")) {
    stop("Argument 'type' must be a 'character'.")
  }

  type_out <- setdiff(type, c("basic", "simple", "detailed"))
  if (length(type_out) > 0) {
    stop("Invalid 'type' provided.",
         "\nAvailable options are: 'basic', 'simple', or 'detailed'.")
  }

  if (!inherits(where_distance, "character")) {
    stop("Argument 'where_distance' must be a 'character'.")
  }

  where_distance_out <- setdiff(where_distance, c("in_range", "out_range", "all"))
  if (length(where_distance_out) > 0) {
    stop("Invalid 'where_distance' provided.",
         "\nThe available options are: 'in_range', 'out_range', or 'all'.")
  }
  if (!inherits(distance, "character")) {
    stop("Argument 'distance' must be a 'character'.")
  }

  distance_out <- setdiff(distance, c("euclidean", "mahalanobis"))
  if (length(distance_out) > 0) {
    stop("Invalid distance provided.",
         "\nThe available options are: 'euclidean' or 'mahalanobis'")
  }
  if (!inherits(scale, c("logical", "numeric"))) {
    stop("Argument 'scale' must be 'logical' or 'numeric'.")
  }
  if (!inherits(center, c("logical", "numeric"))) {
    stop("Argument 'center' must be 'logical' or 'numeric'.")
  }
  if (!inherits(percentage, "numeric")) {
    stop("Argument 'percentage' must be 'numeric'.")
  }

  if(!inherits(calculate_distance, "logical")){
    stop("Argument 'calculate_distance' must be 'logical'.")
  }

  #Get calibration data
  m <-  data$calibration_data[, -1]

  # exclude categorical variables from M
  if (!is.null(data$categorical_variables)) {
    m <- m[, !colnames(m) %in% data$categorical_variables]
  }

  #Check scenarios to predict
  sc <- setdiff(names(projection_data), c("raster_pattern", "pca"))
  sc <- names(projection_data[sc][sapply(projection_data[sc],
                                       function(x) !is.null(x))])

  #Get raster pattern to read
  raster_pattern <- projection_data$raster_pattern

  #Check scenarios for predicting
  res_path <- check_pred_scenarios(projection_data = projection_data,
                                   out_dir = out_dir)

  #Show progress bar?
  if (progress_bar) {
    pb <- utils::txtProgressBar(0, nrow(res_path), style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
  }

  #Looping in projection file
  for(i in 1:nrow(res_path)) {
    #Read rasters
    r <- terra::rast(
      list.files(path = res_path$input_path[i], pattern = raster_pattern,
                 full.names = TRUE))

    if (!is.null(mask)) {
      r <- terra::crop(x = r, y = mask, mask = TRUE)
    }

    # exclude categorical variables from G
    if (!is.null(data$categorical_variables)) {
      nc <- setdiff(names(r), data$categorical_variables)
      r <- r[[nc]]
    }

    #Do PCA?
    if (!is.null(data$pca)) {
      r <- terra::predict(r[[data$pca$vars_in]], data$pca)
    }

    #Subset variables?
    if (subset_variables) {
      v <- unique(unlist(sapply(data$Models, function(x)
        names(x$Full_model$betas)[-1],
                  simplify = F)))
      v <- gsub("I\\((.*?)\\^2\\)", "\\1", v) #Remove quadratic pattern
      v <- v[!grepl("categorical", v)] #Remove categorical pattern
      v <- unique(unlist(strsplit(v, ":"))) #Remove product pattern

      # remove categorical variables again just in case
      if (!is.null(data$categorical_variables)) {
        v <- setdiff(v, data$categorical_variables)
      }

      r <- r[[v]]
      m <- m[, v]
    }

   #See if there are variables in G absent in M
    g_out_m <- setdiff(names(r), names(m))
    if(length(g_out_m) > 0) {
      r <- r[[names(m)]]
      warning(length(g_out_m), " variable(s) in 'G' were removed because they are not present in 'M'.")
    }

    #MOP
    if (parallel) {
      if (is.null(ncores)) {
        ncores <- max(1, parallel::detectCores() - 1)
      }
    }

    mop_i <- mop::mop(m = m, g = r, type = type,
                      calculate_distance = calculate_distance,
                      where_distance = where_distance, distance = distance,
                      scale = scale, center = center, fix_NA = fix_NA,
                      percentage = percentage, comp_each = comp_each,
                      tol = tol, rescale_distance = rescale_distance,
                      parallel = parallel, n_cores = ncores,
                      progress_bar = FALSE)
    #Save results
    dir.create(dirname(res_path$output_path[i]), recursive = TRUE,
               showWarnings = FALSE)
    #Get summary
    s <- data.frame(variables = paste(mop_i$summary$variables, collapse = ", "),
                    type, scale, center, calculate_distance, distance,
                    percentage,
                    rescale_distance, fix_NA, N_m = mop_i$summary$N_m,
                    N_g = mop_i$summary$N_g,
                    m_min = paste(round(mop_i$summary$m_ranges[1, ], 1),
                                  collapse = ", "),
                    m_max = paste(round(mop_i$summary$m_ranges[2, ], 1),
                                  collapse = ", "))

    write.csv(s, paste0(res_path$output_path[i], "summary.csv"))

    if (!is.null(mop_i$mop_distances)) {
      if (all(!is.na(terra::minmax(mop_i$mop_distances)))) {
        if(!na_in_range){
          #Mask mop_distance
          mop_i$mop_distances[is.na(mop_i$mop_distances)] <- 0
          mop_i$mop_distances <- terra::mask(mop_i$mop_distances, r[[1]])}
        #Write
        terra::writeRaster(mop_i$mop_distances,
                           paste0(res_path$output_path[i], "_mop_distances.tif"),
                           overwrite = overwrite)
      }
    }

    if (!is.null(mop_i$mop_basic)) {
      if (all(!is.na(terra::minmax(mop_i$mop_basic)))) {
        if(!na_in_range){
          #Mask mop_basic
          mop_i$mop_basic[is.na(mop_i$mop_basic)] <- 0
          mop_i$mop_basic <- terra::mask(mop_i$mop_basic, r[[1]])}
        #Save
        terra::writeRaster(mop_i$mop_basic,
                           paste0(res_path$output_path[i], "_mopbasic.tif"),
                           overwrite = overwrite)
      }
    }

    if (!is.null(mop_i$mop_simple)) {
      if (all(!is.na(terra::minmax(mop_i$mop_basic)))) {
        if(!na_in_range){
          #Mask mop_simple
          mop_i$mop_simple[is.na(mop_i$mop_simple)] <- 0
          mop_i$mop_simple <- terra::mask(mop_i$mop_simple, r[[1]])
          levels(mop_i$mop_simple) <- NULL}

        terra::writeRaster(mop_i$mop_simple,
                           paste0(res_path$output_path[i], "_mopsimple.tif"),
                           overwrite = overwrite)
      }
    }

    if (!is.null(mop_i$mop_detailed[[1]])) {
      #Set levels
      lapply(names(mop_i$mop_detailed)[-1], function(x) {
        res_x <- mop_i$mop_detailed[[x]]
        #Remove variables with all values in range
        all_in_range <- sapply(res_x, function(z){
          !any(is.na(minmax(z)))
        })
        if(sum(all_in_range) > 0){
          res_x <- res_x[[all_in_range]]

          if(!na_in_range){
            #Mask mop_detailed
            res_x <- mop_i$mop_detailed[[x]]
            res_x[is.na(res_x)] <- 0
            res_x <- terra::mask(res_x, r[[1]], )

            if(all(terra::is.factor(res_x))){
            #Add level no risk
            levels(res_x) <- rbind(data.frame(id = 0, category = "Within ranges"),
                                   levels(res_x)[[1]])}
          }
          terra::writeRaster(res_x,
                             paste0(res_path$output_path[i], "_mop_", x, ".tif"),
                             overwrite = overwrite)
        }
      })
      #Save interpretation_combined
      if(!na_in_range){
        interpretation_combined <- rbind(
          data.frame(values = 0, extrapolation_variables = "Within ranges"),
          mop_i$mop_detailed$interpretation_combined)
        write.csv(interpretation_combined,
                  paste0(res_path$output_path[i],
                         "_interpretation_combined.csv"))
        } else {
          write.csv(mop_i$mop_detailed$interpretation_combined,
                    paste0(res_path$output_path[i],
                           "_interpretation_combined.csv"))
        }

    }
    if (progress_bar) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  #Return root directory and output path
  res_path$output_path <- dirname(res_path$output_path)
  res_final <- list("paths" = unique(res_path),
                    "root_directory" = out_dir)
  class(res_final) <- "mop_projections"

  #Save mop_projections object
  saveRDS(res_final, paste0(out_dir, "/mop_projections.rds"))

  return(res_final)
}

