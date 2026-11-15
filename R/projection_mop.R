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
#'                progress_type = "timer", overwrite = FALSE)
#'
#' @param data an object of class `fitted_models` returned by the
#' [fit_selected()] function or an object of class `prepared_data` returned by
#' the [prepare_data()] function.
#' @param projection_data an object of class `projection_data` returned by the
#' [prepare_projection()] function. This object contains the paths to
#' the rasters representing each scenario.
#' @param out_dir (character) a path to a root directory for saving the raster
#' file of each projection.
#' @param subset_variables (logical) whether to include in the analysis only the
#' variables present in the selected models. Only applicable if `data` is a
#' `fitted_models` object. Default is FALSE.
#' @param mask (SpatRaster, SpatVector, or SpatExtent) spatial object used to
#' mask the variables (optional). Default is NULL.
#' @param type (character) type of MOP analysis to be performed. Options
#' available are `"basic"`, `"simple"`, and `"detailed"`. See Details for
#' further information. Default is `"basic"`.
#' @param na_in_range (logical) whether to assign `NA` to regions within the
#' projected area (G) where environmental conditions fall within the range of
#' the calibration data (M). If `TRUE` (default), these regions are assigned
#' `NA`. If `FALSE`, they are assigned `0` in the simple and basic MOP outputs,
#' and `"Within ranges"` in the detailed MOP output.
#' @param calculate_distance (logical) whether to calculate distances
#' (dissimilarities) between m and g. The default, FALSE, runs rapidly and does
#' not assess dissimilarity levels.
#' @param where_distance (character) where to calculate distances, considering
#' how conditions in g are positioned in comparison to the range of conditions
#' in m. Options available are `"in_range"`, `"out_range"`, and `"all"`.
#' Default is `"in_range"`.
#' @param distance (character) which distances are calculated, `"euclidean"` or
#' `"mahalanobis"`. Only applicable if `calculate_distance = TRUE`. Default is
#' `"euclidean"`.
#' @param scale (logical or numeric) whether to scale as in
#' \code{\link[base]{scale}}. Default is FALSE.
#' @param center (logical or numeric) whether to center as in
#' \code{\link[base]{scale}}. Default is FALSE.
#' @param fix_NA (logical) whether to fix layers so cells with NA values are
#' the same in all layers. Setting to FALSE may save time if the rasters are
#' big and have no NA matching problems. Default is TRUE.
#' @param percentage (numeric) percentage of \code{m} closest conditions used
#' to derive mean environmental distances to each combination of conditions in
#' \code{g}. Default is 1.
#' @param comp_each (numeric) number of combinations in \code{g} to be used for
#' distance calculations at a time. Increasing this number requires more RAM.
#' Default is 2000.
#' @param tol (numeric) tolerance to detect linear dependencies when calculating
#' Mahalanobis distances. The default, NULL, uses `.Machine$double.eps`.
#' @param rescale_distance (logical) whether to re-scale distances 0-1.
#' Re-scaling prevents comparisons of dissimilarity values obtained from runs
#' with different values of \code{percentage}. Default is FALSE.
#' @param parallel (logical) whether to run distance calculations in parallel
#' within `mop::mop()`. Default is FALSE.
#' @param ncores (numeric) number of cores to use for parallel processing.
#' Default is NULL and uses available cores - 1. This is only applicable if
#' `parallel = TRUE`.
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is TRUE.
#' @param progress_type (character) type of the progress bar. See
#' `pbapply::pbtypes()` for available progress bar types depending on operating
#' system. Default is `"timer"`.
#' @param overwrite (logical) whether to overwrite SpatRaster files if they
#' already exist. Default is FALSE.
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
#' dissimilarity:
#' - **in_range** - only conditions inside `m` ranges
#' - **out_range** - only conditions outside `m` ranges
#' - **all** - all conditions
#'
#' When the variables used to represent conditions have different units,
#' scaling and centering are recommended. This step is only valid when Euclidean
#' distances are used.
#'
#' @return
#' An object of class `mop_projections`, with the root directory and the data.frame
#' containing the file paths where the results were stored for each scenario.
#'
#' @export
#'
#' @importFrom utils write.csv
#' @importFrom stats coef
#' @importFrom terra rast predict writeRaster mask crop minmax is.factor
#' @importFrom parallel detectCores
#' @importFrom pbapply pboptions pblapply
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
#'                           name_format = "bio_", static_variables = var$SoilType)
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
                           progress_type = "timer",
                           overwrite = FALSE) {
  # Check data
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

  if (subset_variables && !inherits(data, "fitted_models")) {
    stop("If 'subset_variables' = TRUE, 'data' must be a 'fitted_models' object.")
  }
  if (!is.null(mask) && !inherits(mask, c("SpatRaster", "SpatVector", "SpatExtent"))) {
    stop("Argument 'mask' must be a 'SpatRaster', 'SpatVector', or 'SpatExtent'.")
  }

  if (!inherits(type, "character") || length(type) != 1L) {
    stop("Argument 'type' must be a single 'character'.")
  }
  if (!type %in% c("basic", "simple", "detailed")) {
    stop("Invalid 'type' provided. Available options are: 'basic', 'simple', or 'detailed'.")
  }

  if (!inherits(where_distance, "character") || length(where_distance) != 1L) {
    stop("Argument 'where_distance' must be a single 'character'.")
  }
  if (!where_distance %in% c("in_range", "out_range", "all")) {
    stop("Invalid 'where_distance' provided. Available options are: 'in_range', 'out_range', or 'all'.")
  }

  if (!inherits(distance, "character") || length(distance) != 1L) {
    stop("Argument 'distance' must be a single 'character'.")
  }
  if (!distance %in% c("euclidean", "mahalanobis")) {
    stop("Invalid 'distance' provided. Available options are: 'euclidean' or 'mahalanobis'.")
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
  if (!inherits(calculate_distance, "logical")) {
    stop("Argument 'calculate_distance' must be 'logical'.")
  }
  if (!inherits(fix_NA, "logical")) {
    stop("Argument 'fix_NA' must be 'logical'.")
  }

  # Get calibration data
  m <- data$calibration_data[, -1, drop = FALSE]

  # Exclude categorical variables from M
  if (!is.null(data$categorical_variables)) {
    m <- m[, !colnames(m) %in% data$categorical_variables, drop = FALSE]
  }

  # Subset variables once before looping (supports both maxnet and glm)
  selected_vars <- NULL
  if (subset_variables) {
    coefs <- unique(unlist(lapply(data$Models, function(x) {
      fm <- x$Full_model
      if (inherits(fm, "glmnet")) {
        names(fm$betas)
      } else if (inherits(fm, "glm")) {
        names(stats::coef(fm)[-1L])
      } else {
        character(0)
      }
    }), use.names = FALSE))

    coefs <- ifelse(
      grepl("^(thresholds|hinge|categorical)\\(", coefs),
      sub(":.*$", "", coefs),
      coefs
    )
    coefs <- coefs[!grepl("^categorical\\(", coefs)]
    selected_vars <- unique(unlist(extract_variable_name(unique(coefs)), use.names = FALSE))

    if (!is.null(data$categorical_variables)) {
      selected_vars <- setdiff(selected_vars, data$categorical_variables)
    }
    m <- m[, selected_vars, drop = FALSE]
  }

  # Configure cores safely for CRAN (used inside mop::mop)
  if (parallel && is.null(ncores)) {
    ncores <- max(1L, parallel::detectCores() - 1L, na.rm = TRUE)
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      ncores <- min(ncores, 2L)
    }
  }

  # Get raster pattern and scenario paths
  raster_pattern <- projection_data$raster_pattern
  res_path <- check_pred_scenarios(projection_data = projection_data,
                                   out_dir = out_dir)
  n_scen <- nrow(res_path)

  # Setup progress bar safely for CRAN
  old_pbopt <- pbapply::pboptions(
    type = if (progress_bar) progress_type else "none"
  )
  on.exit(pbapply::pboptions(old_pbopt), add = TRUE)

  # Loop over projection scenarios using pblapply
  pbapply::pblapply(
    X = seq_len(n_scen),
    FUN = function(i) {
      # Read rasters
      r <- terra::rast(
        list.files(path = res_path$input_path[i], pattern = raster_pattern,
                   full.names = TRUE)
      )

      if (!is.null(mask)) {
        r <- terra::crop(x = r, y = mask, mask = TRUE)
      }

      # Exclude categorical variables from G
      if (!is.null(data$categorical_variables)) {
        nc <- setdiff(names(r), data$categorical_variables)
        r <- r[[nc]]
      }

      # Do PCA if present
      if (!is.null(data$pca)) {
        if (!("vars_out" %in% names(data$pca))) {
          r <- terra::predict(r[[data$pca$vars_in]], data$pca)
        } else {
          r <- c(terra::predict(r[[data$pca$vars_in]], data$pca),
                 r[[data$pca$vars_out]])
        }
      }

      # Subset variables if requested
      if (!is.null(selected_vars)) {
        r <- r[[selected_vars]]
      }

      # Check if there are variables in G absent in M
      g_out_m <- setdiff(names(r), names(m))
      if (length(g_out_m) > 0L) {
        r <- r[[names(m)]]
        warning(length(g_out_m),
                " variable(s) in 'G' were removed because they are not present in 'M'.")
      }

      # Run MOP
      mop_i <- mop::mop(m = m, g = r, type = type,
                        calculate_distance = calculate_distance,
                        where_distance = where_distance, distance = distance,
                        scale = scale, center = center, fix_NA = fix_NA,
                        percentage = percentage, comp_each = comp_each,
                        tol = tol, rescale_distance = rescale_distance,
                        parallel = parallel, n_cores = ncores,
                        progress_bar = FALSE)

      # Save results
      out_base <- res_path$output_path[i]
      if (!dir.exists(dirname(out_base))) {
        dir.create(dirname(out_base), recursive = TRUE, showWarnings = FALSE)
      }

      s <- data.frame(
        variables = paste(mop_i$summary$variables, collapse = ", "),
        type = type, scale = scale, center = center,
        calculate_distance = calculate_distance, distance = distance,
        percentage = percentage, rescale_distance = rescale_distance,
        fix_NA = fix_NA, N_m = mop_i$summary$N_m, N_g = mop_i$summary$N_g,
        m_min = paste(round(mop_i$summary$m_ranges[1, ], 1), collapse = ", "),
        m_max = paste(round(mop_i$summary$m_ranges[2, ], 1), collapse = ", "),
        stringsAsFactors = FALSE
      )
      utils::write.csv(s, paste0(out_base, "_summary.csv"), row.names = FALSE)

      if (!is.null(mop_i$mop_distances) &&
          all(!is.na(terra::minmax(mop_i$mop_distances)))) {
        if (!na_in_range) {
          mop_i$mop_distances[is.na(mop_i$mop_distances)] <- 0
          mop_i$mop_distances <- terra::mask(mop_i$mop_distances, r[[1]])
        }
        terra::writeRaster(mop_i$mop_distances,
                           paste0(out_base, "_mop_distances.tif"),
                           overwrite = overwrite)
      }

      if (!is.null(mop_i$mop_basic) &&
          all(!is.na(terra::minmax(mop_i$mop_basic)))) {
        if (!na_in_range) {
          mop_i$mop_basic[is.na(mop_i$mop_basic)] <- 0
          mop_i$mop_basic <- terra::mask(mop_i$mop_basic, r[[1]])
        }
        terra::writeRaster(mop_i$mop_basic,
                           paste0(out_base, "_mopbasic.tif"),
                           overwrite = overwrite)
      }

      if (!is.null(mop_i$mop_simple) &&
          all(!is.na(terra::minmax(mop_i$mop_simple)))) {
        if (!na_in_range) {
          mop_i$mop_simple[is.na(mop_i$mop_simple)] <- 0
          mop_i$mop_simple <- terra::mask(mop_i$mop_simple, r[[1]])
          levels(mop_i$mop_simple) <- NULL
        }
        terra::writeRaster(mop_i$mop_simple,
                           paste0(out_base, "_mopsimple.tif"),
                           overwrite = overwrite)
      }

      if (!is.null(mop_i$mop_detailed[[1]])) {
        lapply(names(mop_i$mop_detailed)[-1L], function(x) {
          res_x <- mop_i$mop_detailed[[x]]
          all_in_range <- sapply(res_x, function(z) {
            !any(is.na(terra::minmax(z)))
          })

          if (sum(all_in_range) > 0L) {
            res_x <- res_x[[all_in_range]]

            if (!na_in_range) {
              res_x[is.na(res_x)] <- 0
              res_x <- terra::mask(res_x, r[[1]])

              if (all(terra::is.factor(res_x))) {
                levels(res_x) <- rbind(
                  data.frame(id = 0, category = "Within ranges"),
                  levels(res_x)[[1]]
                )
              }
            }

            terra::writeRaster(res_x,
                               paste0(out_base, "_mop_", x, ".tif"),
                               overwrite = overwrite)
          }
        })

        # Save interpretation_combined
        if (!na_in_range) {
          interpretation_combined <- rbind(
            data.frame(values = 0, extrapolation_variables = "Within ranges"),
            mop_i$mop_detailed$interpretation_combined
          )
        } else {
          interpretation_combined <- mop_i$mop_detailed$interpretation_combined
        }
        utils::write.csv(interpretation_combined,
                         paste0(out_base, "_interpretation_combined.csv"),
                         row.names = FALSE)
      }

      return(invisible(NULL))
    },
    cl = NULL
  )

  # Return root directory and output path
  res_path$output_path <- dirname(res_path$output_path)
  res_final <- list("paths" = unique(res_path),
                    "root_directory" = out_dir)
  class(res_final) <- "mop_projections"

  # Save mop_projections object
  saveRDS(res_final, file.path(out_dir, "mop_projections.rds"))

  return(res_final)
}
