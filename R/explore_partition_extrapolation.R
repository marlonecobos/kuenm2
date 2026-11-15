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
#'                                        progress_bar = FALSE,
#'                                        progress_type = "timer", ...)
#'
#' @param data an object of class `prepared_data` returned by the
#' [prepare_data()] function.
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
#' available are `"basic"`, `"simple"`, and `"detailed"`. Default is
#' `"detailed"`. See \code{\link{projection_mop}}() for more details.
#' @param calculate_distance (logical) whether to calculate distances
#' (dissimilarities) between train and test data. Default is TRUE.
#' @param where_distance (character) specifies which values in train data should
#' be used to calculate distances. Options are: `"in_range"` (only conditions
#' within the train range), `"out_range"` (only conditions outside the
#' train range), and `"all"` (all conditions). Default is `"all"`.
#' @param progress_bar (logical) whether to display a progress bar during
#' processing. Default is FALSE.
#' @param progress_type (character) type of the progress bar. See
#' `pbapply::pbtypes()` for available progress bar types depending on operating
#' system. Default is `"timer"`.
#' @param ... additional arguments passed to \code{\link[mop]{mop}()}.
#'
#' @return
#' An object of class `explore_partition` containing:
#' - `Mop_results`: a `data.frame` with MOP distances (if
#'   `calculate_distance = TRUE`), indicators of non-analogous conditions,
#'   number and names of variables outside the training range, and novel
#'   categories (if categorical variables are present) for each partition;
#' - `calibration_data`: a `data.frame` with the calibration data ordered by
#'   partition;
#' - `categorical_variables`: character vector with categorical variable names
#'   (if present);
#' - `pca`: PCA object (if `do_pca = TRUE` was used in `prepare_data()`).
#'
#' @export
#' @importFrom mop mop
#' @importFrom pbapply pboptions pblapply
#' @examples
#' # Prepare data
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
#' res <- explore_partition_extrapolation(data = sp_swd)
#'
explore_partition_extrapolation <- function(data,
                                            include_train_background = TRUE,
                                            include_test_background = FALSE,
                                            variables = NULL,
                                            mop_type = "detailed",
                                            calculate_distance = TRUE,
                                            where_distance = "all",
                                            progress_bar = FALSE,
                                            progress_type = "timer", ...) {
  # Check data
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (!inherits(data, "prepared_data")) {
    stop("'data' must be a 'prepared_data' object.")
  }
  if (!inherits(include_train_background, "logical") ||
      length(include_train_background) != 1L) {
    stop("'include_train_background' must be a single 'logical'.")
  }
  if (!inherits(include_test_background, "logical") ||
      length(include_test_background) != 1L) {
    stop("'include_test_background' must be a single 'logical'.")
  }
  if (!is.null(variables)) {
    if (!inherits(variables, "character")) {
      stop("'variables' must be 'character' or NULL.")
    }
    var_out <- setdiff(variables, colnames(data$calibration_data))
    if (length(var_out) > 0L) {
      stop("Some 'variables' provided are absent from 'data': ",
           paste(var_out, collapse = ", "))
    }
  }

  if (!inherits(mop_type, "character") || length(mop_type) != 1L ||
      !(mop_type %in% c("basic", "simple", "detailed"))) {
    stop("Invalid 'mop_type' provided. Available options are: 'basic', 'simple', or 'detailed'.")
  }
  if (!inherits(calculate_distance, "logical") ||
      length(calculate_distance) != 1L) {
    stop("Argument 'calculate_distance' must be a single 'logical'.")
  }
  if (!inherits(where_distance, "character") || length(where_distance) != 1L ||
      !(where_distance %in% c("in_range", "out_range", "all"))) {
    stop("Invalid 'where_distance' provided. Available options are: 'in_range', 'out_range', and 'all'.")
  }
  if (!inherits(progress_bar, "logical") || length(progress_bar) != 1L) {
    stop("Argument 'progress_bar' must be a single 'logical'.")
  }

  # Always exclude 'pr_bg' and categorical variables from continuous MOP variables
  v <- setdiff(colnames(data$calibration_data),
               c("pr_bg", data$categorical_variables))

  # Subset variables if requested
  if (!is.null(variables)) {
    v <- intersect(variables, v)
  }

  if (length(v) == 0L) {
    stop("No continuous variables available to perform MOP analysis.")
  }

  part_names <- names(data$part_data)
  if (is.null(part_names)) {
    part_names <- as.character(seq_along(data$part_data))
  }

  has_xy <- !is.null(data$data_xy)
  xy_cols <- if (has_xy) colnames(data$data_xy) else NULL
  cat_vars <- data$categorical_variables
  calib_df <- data$calibration_data

  # Setup progress bar options safely for CRAN
  old_pbopt <- pbapply::pboptions(
    type = if (progress_bar) progress_type else "none"
  )
  on.exit(pbapply::pboptions(old_pbopt), add = TRUE)

  # Evaluate partitions with pblapply
  all_res_list <- pbapply::pblapply(
    X = seq_along(data$part_data),
    FUN = function(i) {
      i_name <- part_names[i]
      partition_i <- data$part_data[[i]]

      test_i <- calib_df[partition_i, , drop = FALSE]
      train_i <- calib_df[-partition_i, , drop = FALSE]

      if (!include_train_background) {
        train_i <- train_i[train_i$pr_bg == 1, , drop = FALSE]
      }
      if (!include_test_background) {
        test_keep <- which(test_i$pr_bg == 1)
        test_i <- test_i[test_keep, , drop = FALSE]
      } else {
        test_keep <- seq_len(nrow(test_i))
      }

      # Run MOP analysis
      res_i <- mop_with_records(train_data = train_i,
                                test_data = test_i,
                                variables = v,
                                categorical_variables = cat_vars,
                                mop_type = mop_type,
                                calculate_distance = calculate_distance,
                                where_distance = where_distance,
                                progress_bar = FALSE, ...)

      df_out <- res_i$mop_records
      df_out$pr_bg <- test_i$pr_bg

      # Append coordinates if available
      if (has_xy) {
        xy <- data$data_xy[partition_i[test_keep], xy_cols, drop = FALSE]
        df_out <- cbind(xy, df_out)
      }

      df_out$Partition <- i_name
      return(df_out)
    },
    cl = NULL
  )

  all_res <- do.call("rbind", all_res_list)
  rownames(all_res) <- NULL

  # Organize columns cleanly
  first_cols <- intersect(c("Partition", "pr_bg", xy_cols), colnames(all_res))
  other_cols <- setdiff(colnames(all_res), first_cols)
  all_res <- all_res[, c(first_cols, other_cols), drop = FALSE]

  # Reorganize calibration data by partition order
  partition_order <- unlist(data$part_data, use.names = FALSE)
  calib_data <- calib_df[partition_order, , drop = FALSE]

  # Return final results
  res_final <- list("Mop_results" = all_res,
                    "calibration_data" = calib_data,
                    "categorical_variables" = cat_vars,
                    "pca" = data$pca)
  class(res_final) <- "explore_partition"
  return(res_final)
}
