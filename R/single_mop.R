#' Analysis of extrapolation risks using the MOP metric (for single scenario)
#'
#' @description
#' Calculates the mobility-oriented parity metric and other sub-products to
#' represent dissimilarities and non-analogous conditions when comparing a set
#' of reference conditions (M) against another set of scenario conditions (G).
#'
#' @usage
#' single_mop(data, new_variables, subset_variables = FALSE,
#'            mask = NULL, type = "basic", na_in_range = TRUE,
#'            calculate_distance = FALSE, where_distance = "in_range",
#'            distance = "euclidean", scale = FALSE, center = FALSE,
#'            fix_NA = TRUE, percentage = 1, comp_each = 2000, tol = NULL,
#'            rescale_distance = FALSE, parallel = FALSE, ncores = NULL,
#'            progress_bar = TRUE, write_files = FALSE, out_dir = NULL,
#'            overwrite = FALSE)
#'
#' @param data an object of class `fitted_models` returned by the
#' [fit_selected()] function, an object of class `prepared_data` returned by
#' the [prepare_data()] function, or an object of class `calibration_results`
#' returned by the [calibration()] function.
#' @param new_variables a SpatRaster or data.frame of predictor variables.
#' The names of these variables must match those used to prepare the date or
#' calibrate the models provided in `data`.
#' @param subset_variables (logical) whether to include in the analysis only the
#' variables present in the selected models. Only applicable if `data` is a
#' `fitted_models` object. Default is FALSE
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
#' @param parallel (logical) whether to compute MOP in parallel. Default is
#' FALSE.
#' @param ncores (numeric) number of cores to use for parallel processing.
#' Default is NULL and uses available cores - 1. This is only applicable if
#' `parallel = TRUE`.
#' @param progress_bar (logical) whether to display a progress bar during
#' distance calculation. Only applicable if calculate_distance is TRUE.
#' Default is TRUE.
#' @param write_files (logical) whether to save the MOP results (SpatRasters and
#' data.frames) to disk. Default is FALSE.
#' @param out_dir (character) directory path where results will be saved.
#' Only relevant if `write_files = TRUE`.
#' @param overwrite (logical) whether to overwrite SpatRasters if they already
#' exist. Only applicable if `write_files = TRUE`. Default is FALSE.
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
#' @returns
#' An object of class `mop_results` containing:
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
#'     - *m_ranges* - the range (minimum and maximum values) of the variable in
#'     reference conditions (\code{m})
#' - **mop_distances** - if \code{calculate_distance} = TRUE, a SpatRaster
#' with distance values for the set of interest (\code{g}). Higher values
#' represent greater dissimilarity compared to the set of reference (\code{m}).
#' - **mop_basic** - a SpatRaster for the set of interest representing
#' conditions in which at least one of the variables is non-analogous to the set
#' of reference. Values should be: 1 for non-analogous conditions, and NA for
#' conditions inside the ranges of the reference set.
#' - **mop_simple** - a SpatRaster, for the set of interest, representing how
#' many variables in the set of interest are non-analogous to those in the
#' reference set. NA is used for conditions inside the ranges of the reference
#' set.
#' - **mop_detailed** - a list containing:
#'     - *interpretation_combined* - a data.frame to help identify combinations
#'     of variables in *towards_low_combined* and *towards_high_combined* that
#'     are non-analogous to \code{m}.
#'     - *towards_low_end* - a SpatRaster for all variables representing where
#'     non-analogous conditions were found towards low values of each variable.
#'     - *towards_high_end* - a SpatRaster for all variables representing where
#'     non-analogous conditions were found towards high values of each variable.
#'     - *towards_low_combined* - a SpatRaster with values representing the
#'     identity of the variables found to have non-analogous conditions towards
#'     low values.
#'     - *towards_high_combined* - a SpatRaster with values representing the
#'     identity of the variables found to have non-analogous conditions towards
#'     high values.
#'
#' @export
#'
#' @importFrom terra crop predict minmax mask writeRaster is.factor
#' @importFrom parallel detectCores
#' @importFrom mop mop
#'
#' @seealso
#' [projection_mop()]
#' @examples
#' # Import an example of fitted models (output of fit_selected())
#' data("fitted_model_maxnet", package = "kuenm2")
#'
#' # Import variables under a new set of conditions
#' # Here, future climate data
#' future_scenario <- terra::rast(system.file("extdata",
#'                                           "wc2.1_10m_bioc_ACCESS-CM2_ssp585_2081-2100.tif",
#'                                           package = "kuenm2"))
#'
#'
#'
#' # Rename variables to match the variable names in the fitted models
#' names(future_scenario) <- sub("bio0", "bio", names(future_scenario))
#' names(future_scenario) <- sub("bio", "bio_", names(future_scenario))
#'
#' # Run MOP
#' sm <- single_mop(data = fitted_model_maxnet, new_variables = future_scenario,
#'                  type = "detailed")

single_mop <- function(data,
                       new_variables,
                       subset_variables = FALSE,
                       mask = NULL,
                       type = "basic",
                       na_in_range = TRUE,
                       calculate_distance = FALSE,
                       where_distance = "in_range",
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
                       write_files = FALSE,
                       out_dir = NULL,
                       overwrite = FALSE) {

  #Check data
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }

  if (!inherits(data, "prepared_data") && !inherits(data, "fitted_models") &&
      !inherits(data, "calibration_results")) {
    stop("Argument 'data' must be a 'prepared_data', 'fitted_models' or 'calibration_results' object.")
  }

  if (missing(new_variables)) {
    stop("Argument 'new_variables' must be defined.")
  }

  if (!inherits(new_variables, "SpatRaster")) {
    stop("Argument 'new_variables' must be a 'SpatRaster' object.")
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

  if(!inherits(fix_NA, "logical")){
    stop("Argument 'fix_NA' must be 'logical'.")
  }

  if(!inherits(parallel, "logical")){
    stop("Argument 'parallel' must be 'logical'.")
  }

  if(!inherits(progress_bar, "logical")){
    stop("Argument 'progress_bar' must be 'logical'.")
  }

  if(!inherits(write_files, "logical")){
    stop("Argument 'write_files' must be 'logical'.")
  }

  if(write_files){
    if(is.null(out_dir)){
      stop("If 'write_files = TRUE', 'out_dir' must be specified")
    }

    if (!inherits(out_dir, "character")) {
      stop("Argument 'out_dir' must be a 'character'.")
    }

    if(!file.exists(out_dir)){
      stop("Directory specified in 'out_dir' does not exist.")
    }

    if(!inherits(overwrite, "logical")){
      stop("Argument 'overwrite' must be 'logical'.")
    }
  }


  #### Function starts here ####

  #Get calibration data
  m <-  data$calibration_data[, -1]

  # exclude categorical variables from M
  if (!is.null(data$categorical_variables)) {
    m <- m[, !colnames(m) %in% data$categorical_variables]
  }


  # Mask
  if (!is.null(mask)) {
    new_variables <- terra::crop(x = new_variables, y = mask, mask = TRUE)
  }

  # exclude categorical variables from G
  if (!is.null(data$categorical_variables)) {
    nc <- setdiff(names(new_variables), data$categorical_variables)
    new_variables <- new_variables[[nc]]
  }

  #Do PCA?
  if (!is.null(data$pca)) {
    new_variables <- terra::predict(new_variables[[data$pca$vars_in]], data$pca)
  }

  #Subset variables?
  if (subset_variables && inherits(data, "fitted_models")) {
    v <- unique(unlist(sapply(data$Models, function(x)
      names(x$Full_model$betas),
      simplify = F)))
    v <- gsub("I\\((.*?)\\^2\\)", "\\1", v) #Remove quadratic pattern
    v <- v[!grepl("categorical", v)] #Remove categorical pattern
    v <- unique(unlist(strsplit(v, ":"))) #Remove product pattern

    # remove categorical variables again just in case
    if (!is.null(data$categorical_variables)) {
      v <- setdiff(v, data$categorical_variables)
    }

    new_variables <- new_variables[[v]]
    m <- m[, v, drop = FALSE]
  }


  #See if there are variables in G absent in M
  g_out_m <- setdiff(names(new_variables), names(m))
  if(length(g_out_m) > 0) {
    new_variables <- new_variables[[names(m)]]
    warning(length(g_out_m), " variable(s) in 'new_variables' were removed because they are not present in 'M'.")
  }

  #MOP
  if (parallel) {
    if (is.null(ncores)) {
      ncores <- max(1, parallel::detectCores() - 1)
    }
  }

  mop_i <- mop::mop(m = m, g = new_variables, type = type,
                    calculate_distance = calculate_distance,
                    where_distance = where_distance, distance = distance,
                    scale = scale, center = center, fix_NA = fix_NA,
                    percentage = percentage, comp_each = comp_each,
                    tol = tol, rescale_distance = rescale_distance,
                    parallel = parallel, n_cores = ncores,
                    progress_bar = progress_bar)

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
  if(write_files){
    write.csv(s, paste0(out_dir, "/summary.csv"))
  }


  if (!is.null(mop_i$mop_distances)) {
    if (all(!is.na(terra::minmax(mop_i$mop_distances)))) {
      if(!na_in_range){
        #Mask mop_distance
        mop_i$mop_distances[is.na(mop_i$mop_distances)] <- 0
        mop_i$mop_distances <- terra::mask(mop_i$mop_distances,
                                           new_variables[[1]])}
      # Write
      if(write_files){
        terra::writeRaster(mop_i$mop_distances,
                           paste0(out_dir, "/mop_distance.tif"),
                           overwrite = overwrite)
      }
    }
  }

  if (!is.null(mop_i$mop_basic)) {
    if (all(!is.na(terra::minmax(mop_i$mop_basic)))) {
      if(!na_in_range){
        #Mask mop_basic
        mop_i$mop_basic[is.na(mop_i$mop_basic)] <- 0
        mop_i$mop_basic <- terra::mask(mop_i$mop_basic, new_variables[[1]])}
      if(write_files){
        terra::writeRaster(mop_i$mop_basic,
                           paste0(out_dir, "/mop_basic.tif"),
                           overwrite = overwrite)
      }
    }
  }

  if (!is.null(mop_i$mop_simple)) {
    if (all(!is.na(terra::minmax(mop_i$mop_basic)))) {
      if(!na_in_range){
        #Mask mop_simple
        mop_i$mop_simple[is.na(mop_i$mop_simple)] <- 0
        mop_i$mop_simple <- terra::mask(mop_i$mop_simple, new_variables[[1]])
        levels(mop_i$mop_simple) <- NULL}
      #
      if(write_files){
        terra::writeRaster(mop_i$mop_simple,
                           paste0(out_dir, "/mop_simple.tif"),
                           overwrite = overwrite)
      }
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
          res_x <- terra::mask(res_x, new_variables[[1]], )

          if(all(terra::is.factor(res_x))){
            #Add level no risk
            levels(res_x) <- rbind(data.frame(id = 0,
                                              category = "Within ranges"),
                                   levels(res_x)[[1]])}
        }
        if(write_files){
          terra::writeRaster(res_x,
                             paste0(out_dir, "/mop_", x, ".tif"),
                             overwrite = overwrite)
        }
      }
    })

    #Save interpretation_combined
    if(write_files){
      if(!na_in_range){
        interpretation_combined <- rbind(
          data.frame(values = 0, extrapolation_variables = "Within ranges"),
          mop_i$mop_detailed$interpretation_combined)
        write.csv(interpretation_combined,
                  paste0(out_dir,
                         "/interpretation_combined.csv"))
      } else {
        write.csv(mop_i$mop_detailed$interpretation_combined,
                  paste0(out_dir,
                         "/interpretation_combined.csv"))
      }

    }

  }
  return(mop_i)
}
