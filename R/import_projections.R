#' Import rasters resulting from projection functions
#'
#' @description
#' This function facilitates the import of rasters that have been generated and
#' written to disk by the `project_selected()`, `projection_changes()`,
#' `variability_projections()`, and `projection_mop()` functions. Users can
#' select specific periods (past/future), emission scenarios, General Circulation
#' Models (GCMs), and result types for import.
#'
#' @param projection an object of class `model_projections`,
#' `changes_projections`, `variability_projections`, or `mop_projections`. This
#' object is the direct output from one of the projection functions listed in
#' the description.
#' @param consensus (character) consensus measures to import. Available options
#' are: 'median', 'range', 'mean' and 'stdev' (standard deviation). Default is
#' c("median", "range", "mean", "stdev"), which imports all options. Only
#' applicable if `projection` is a `model_projections` object.
#' @param present (logical) wheter to import present-day projections. Default is
#' TRUE. Not applicable if projection is a  `changes_projections` object.
#' @param past_period (character) names of specific past periods (e.g., 'LGM' or
#' 'MID') to import. Default is NULL, meaning all available past periods will be
#' imported.
#' @param past_gcm (character) names of specific General Circulation Models
#' (GCMs) from the past to import. Default is NULL, meaning all available past
#' GCMs will be imported.
#' @param future_period (character) names of specific future periods (e.g.,
#' '2041-2060' or '2081-2100') to import. Default is NULL, meaning all available
#' future periods will be imported.
#' @param future_pscen (character) names of specific future emission scenarios
#' (e.g., 'ssp126' or 'ssp585') to import. Default is NULL, meaning all
#' available future scenarios will be imported.
#' @param future_gcm (character) names of specific General Circulation Models
#' (GCMs) from the future to import. Default is NULL, meaning all available
#' future GCMs will be imported.
#' @param change_types (character) names of the type of computed changes to
#' import. Available options are: 'summary', 'by_gcm', and 'by_change'. Default
#' is c("summary", "by_gcm", "by_change"), importing all types. Only applicable
#' if projection is a `changes_projections` object.
#' @param mop_types (character) type(s) of MOP to import. Available options are:
#' basic', 'simple', 'towards_high_combined', 'towards_low_combined',
#' towards_high_end', and 'towards_low_end'. Default is NULL, meaning all
#' available MOPs will be imported. Only applicable if projection is a
#' `mop_projections` object.
#'
#' @return A SpatRaster or a list of SpatRasters, structured according to the
#' input `projection` class:
#'   \itemize{
#'     \item If `projection` is `model_projections`: A stacked `SpatRaster`
#'           containing all selected projections.
#'     \item If `projection` is `changes_projections`: A list of `SpatRaster`s,
#'           organized by the selected `change_types` (e.g., 'summary', 'by_gcm', and/or 'by_change').
#'     \item If `projection` is `mop_projections`: A list of `SpatRaster`s,
#'           organized by the selected `mop_types` (e.g., 'simple' and 'basic').
#'     \item If `projection` is `variability_projections`: A list of `SpatRaster`s,
#'           containing the computed variability.
#'   }
#'
#' @export
#'
#' @importFrom terra rast
#'
#' @seealso
#' [prepare_projection()], [projection_changes()], [projection_variability()],
#' [projection_mop()]
#'
#' @examples
#' # Load packages
#' library(terra)
#' # Step 1: Organize variables for current projection
#' ## Import current variables (used to fit models)
#' var <- terra::rast(system.file("extdata", "Current_variables.tif",
#'                                package = "kuenm2"))
#'
#' ## Create a folder in a temporary directory to copy the variables
#' out_dir_current <- file.path(tempdir(), "Current_raw2")
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
#' out_dir_future <- file.path(tempdir(), "Future_raw2")
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
#' # Step 4: Run multiple model projections
#' ## A folder to save projection results
#' out_dir <- file.path(tempdir(), "Projection_results/maxnet")
#' dir.create(out_dir, recursive = TRUE)
#'
#' ## Project selected models to multiple scenarios
#' p <- project_selected(models = fitted_model_maxnet, projection_data = pr,
#'                       out_dir = out_dir)
#'
#' # Use import_projections to import results:
#' raster_p <- import_projections(projection = p, consensus = "mean")
#' plot(raster_p)
#'
#' # Step 5: Identify areas of change in projections
#' ## Contraction, expansion and stability
#' changes <- projection_changes(model_projections = p, output_dir = out_dir,
#'                               overwrite = TRUE)
#' # Use import_projections to import results:
#' raster_changes <- import_projections(projection = changes,
#'                                      change_type = c("summary", "by_gcm"))
#'
#' plot(raster_changes$by_gcm)
#' plot(raster_changes$Summary)
#'
#' # Step 6: Perform MOP for all projection scenarios
#' ## Create a folder to save MOP results
#' out_dir <- file.path(tempdir(), "MOP_results")
#' dir.create(out_dir, recursive = TRUE)
#'
#' #Import prepared data to serve as a base for MOP comparisons
#' data(sp_swd_cat, package = "kuenm2")
#'
#' ## Run MOP
#' kmop <- projection_mop(data = sp_swd_cat, projection_data = pr,
#'                        out_dir = out_dir, fitted_models = fitted_model_maxnet,
#'                        type = "detailed")
#' # Use import_projections to import results:
#' raster_mop <- import_projections(projection = kmop,
#'                                  mop_types = c("simple", "basic",
#'                                                "towards_high_combined",
#'                                                "towards_low_combined"))
#' plot(raster_mop$simple)
#' plot(raster_mop$basic)
#' plot(raster_mop$towards_high_combined)
#' plot(raster_mop$towards_low_combined)
#'
#' # Step 7: Compute variance from distinct sources
#' ## Set folder to save results
#' out_dir <- file.path(tempdir())
#' v <- projection_variability(model_projections = p, by_replicate = FALSE,
#'                             write_files = TRUE, output_dir = out_dir,
#'                             overwrite=TRUE)
#' v
#' raster_variability <- import_projections(projection = v,
#'                                          future_period = "2041-2060")
#' plot(raster_variability)
#'
import_projections <- function(projection,
                               consensus = c("median", "range", "mean", "stdev"),
                               present = TRUE,
                               past_period = NULL,
                               past_gcm = NULL,
                               future_period = NULL,
                               future_pscen = NULL,
                               future_gcm = NULL,
                               change_types = c("summary", "by_gcm", "by_change"),
                               mop_types = c("simple", "basic",
                                             "towards_high_combined",
                                             "towards_low_combined",
                                             "towards_high_end",
                                             "towards_low_end")){

  #Check general data
  if(!class(projection) %in% c("model_projections", "mop_projections",
                              "variability_projections", "changes_projections"))
    stop("'projection' must be an object of class 'model_projections',
    'mop_projections', 'variability_projections', or 'changes_projections'")

  if (!inherits(present, "logical")) {
    stop("Argument 'present' must be 'logical'.")
  }
  if (!is.null(past_period) & !inherits(past_period, "character")) {
    stop("Argument 'past_period' must be NULL or a 'character'.")
  }
  if (!is.null(past_gcm) & !inherits(past_gcm, "character")) {
    stop("Argument 'past_gcm' must be NULL or a 'character'.")
  }
  if (!is.null(future_period) & !inherits(future_period, "character")) {
    stop("Argument 'future_period' must be NULL or a 'character'.")
  }
  if (!is.null(future_pscen) & !inherits(future_pscen, "character")) {
    stop("Argument 'future_pscen' must be NULL or a 'character'.")
  }
  if (!is.null(future_gcm) & !inherits(future_gcm, "character")) {
    stop("Argument 'future_gcm' must be NULL or a 'character'.")
  }

  #Get paths and check data linked to model_projections
  if(inherits(projection, "model_projections")){
    #Check consensus
    if (!inherits(consensus, "character")) {
      stop("Argument 'consensus' must be a 'character'.")
    }
    consensus_out <- setdiff(consensus, c("median", "range", "mean", "stdev"))
    if (length(consensus_out) > 0) {
      stop("Invalid 'consensus' provided.",
           "\nAvailable options are: 'median', 'range', 'mean' and 'stdev'.")
    }

    #Check past_period
    if(!is.null(past_period)){
      available_past_period <- na.omit(unique(projection$paths$Period[projection$paths$Time == "Past"]))
      past_period_out <- setdiff(past_period, available_past_period)
      if (length(past_period_out) > 0) {
        stop("Invalid 'past_period' provided.",
             "\nAvailable options are: ", paste(available_past_period,
                                                collapse = "; "))
      }
    }

    #Check past_gcm
    if(!is.null(past_gcm)){
      available_past_gcm <- na.omit(unique(projection$paths$GCM[projection$paths$Time == "Past"]))
      past_gcm_out <- setdiff(past_gcm, available_past_gcm)
      if (length(past_gcm_out) > 0) {
        stop("Invalid 'past_gcm' provided.",
             "\nAvailable options are: ", paste(available_past_gcm,
                                                collapse = "; "))
      }
    }

    #Check future_period
    if(!is.null(future_period)){
      available_future_period <- na.omit(unique(projection$paths$Period[projection$paths$Time == "Future"]))
      future_period_out <- setdiff(future_period, available_future_period)
      if (length(future_period_out) > 0) {
        stop("Invalid 'future_period' provided.",
             "\nAvailable options are: ", paste(available_future_period,
                                                collapse = "; "))
      }
    }

    #Check future_gcm
    if(!is.null(future_gcm)){
      available_future_gcm <- na.omit(unique(projection$paths$Period[projection$paths$GCM == "Future"]))
      future_gcm_out <- setdiff(future_gcm, available_future_gcm)
      if (length(future_gcm_out) > 0) {
        stop("Invalid 'future_gcm' provided.",
             "\nAvailable options are: ", paste(available_future_gcm,
                                                collapse = "; "))
      }
    }

    #Check ssps
    if(!is.null(future_pscen)){
      available_future_pscen <- na.omit(unique(projection$paths$ssp))
      future_pscen_out <- setdiff(future_pscen, available_future_pscen)
      if (length(future_pscen_out) > 0) {
        stop("Invalid 'future_pscen' provided.",
             "\nAvailable options are: ", paste(available_future_pscen,
                                                collapse = "; "))
      }
    }

    #Get paths
    p_paths <- projection$paths}


  if(inherits(projection, "mop_projections")){
    #Check data
    mop_out <- setdiff(mop_types, c("simple", "basic",
                                    "towards_high_combined",
                                    "towards_low_combined",
                                    "towards_high_end",
                                    "towards_low_end"))
    if(length(mop_out) > 0)
      stop("One of the mop_types provided is not valid. The options are 'simple',
           'basic', 'towards_high_combined', 'towards_low_combined',
           'towards_high_end', and 'towards_low_end'.")

    #Check if root_directory exists
    if(is.null(projection$root_directory)){
      stop("'projection' does not have a root directory.
      Make sure you ran 'projection_changes()' with write_results set to 'TRUE'")
    }

    #Check if root_directory exists:
    if(!file.exists(projection$root_directory)){
      stop("The following root directory does not exist on this system:\n",
           projection$root_directory)
    }

    #Get paths
    p_paths <- projection$paths
    root_directory = projection$root_directory}


  #Subset scenarios from mop and model projections
  if(inherits(projection, "mop_projections") |
     inherits(projection, "model_projections")){
  #Extract periods and gcms
  Times <- unique(projection$paths$Time)
  # periods <- unique(p$paths$Period)
  # scenarios <- na.omit(unique(p$paths$ssp))
  # gcms <- na.omit(unique(p$paths$GCM))


  #Subset scenarios
  if(present){
    present_paths <- p_paths[p_paths$Time == "Present",]} else {
      present_paths <- NULL
    }

  if("Past" %in% Times){
    past_paths <- p_paths[p_paths$Time == "Past",]
    if(!is.null(past_period))
      past_paths <- p_paths[p_paths$Period %in% past_period,]

    if(!is.null(past_gcm))
      past_paths <- past_paths[past_paths$GCM %in% past_gcm,]
    } else {
      past_paths <- NULL }#End of Past in times


  if("Future" %in% Times){
    future_paths <- p_paths[p_paths$Time == "Future",]
    if(!is.null(future_period))
      future_paths <- future_paths[future_paths$Period %in% future_period,]

    if(!is.null(future_pscen))
      future_paths <- future_paths[future_paths$ssp %in% future_pscen,]

    if(!is.null(future_gcm))
      future_paths <- future_paths[future_paths$GCM %in% future_gcm,]
  } else {
    future_paths <- NULL }#End of Future in Times

  #Join data
  p_paths <- do.call("rbind", list(present_paths, past_paths, future_paths))
  }

  #Import model_projections ####
  if(inherits(projection, "model_projections")){
    output_rasters <- terra::rast(lapply(1:nrow(p_paths), function(i){
      r_i <- terra::rast(paste0(p_paths$output_path[i], "/General_consensus.tif"))
      #Extract consensus
      r_i <- r_i[[consensus]]
      #Rename
      new_name <- paste(na.omit(c(p_paths$Time[i], p_paths$Period[i],
                                  p_paths$Scenario[i], p_paths$ssp[i],
                                  p_paths$GCM[i])),
                        collapse = "_")
      if(new_name == "Present_Present_Present")
        new_name <- "Present"
      names(r_i) <- paste(new_name, names(r_i), sep = "_")
      return(r_i)
    }))
    }

  #### Import mop ####
  if(inherits(projection, "mop_projections")){
  # 1. List all raster files in the root directory
  all_files_in_root <- list.files(path = root_directory, recursive = TRUE,
                                  full.names = TRUE, pattern = "\\.tif$")

  # Initialize the list to store the results
  output_rasters <- list()

  # Iterate over each desired MOP type
  for (mop_type in mop_types) {
    # Initialize a temporary list to store the rasters for the current MOP type
    current_mop_rasters <- list()

    # Iterate over each row of the p_paths dataframe
    for (i in 1:nrow(p_paths)) {
      row_info <- p_paths[i, ]

      # Construct the filename pattern based on the dataframe information
      if (mop_type == "towards_high_combined" || mop_type == "towards_low_combined" ||
          mop_type == "towards_high_end" || mop_type == "towards_low_end") {
        mop_suffix <- paste0("_mop_", mop_type, ".tif$")
      } else {
        mop_suffix <- paste0("_mop", mop_type, ".tif$")
      }

      # Construct the expected file path
      if (row_info$Time == "Present") {
        # For "Present", the path is simpler (Present/Present_mop_type.tif)
        file_partial_path <- paste0("Present/Present", mop_suffix)
        # GCM might be NA for "Present", so the raster name will just be "Present"
        raster_name_id <- "Present"
      } else {
        # For other cases, use Time/Period/ssp/GCM_mop_type.tif
        file_partial_path <- paste0(
          row_info$Time, "/",
          row_info$Period, "/",
          row_info$ssp, "/",
          row_info$GCM, mop_suffix
        )
        # Construct the name for the raster (e.g., "Future_2081-2100_ssp126_ACCESS-CM2")
        raster_name_id <- paste(row_info$Time, row_info$Period, row_info$ssp,
                                row_info$GCM, sep = "_")
      }

      # Find the corresponding file in the list of all files
      escaped_file_partial_path <- gsub("\\.", "\\\\.", file_partial_path)

      # Filter files using grepl
      matching_files <- all_files_in_root[grepl(escaped_file_partial_path,
                                                all_files_in_root, fixed = FALSE)]


      # Check if exactly one file was found
      if (length(matching_files) == 1) {
        # Read the raster
        raster_obj <- terra::rast(matching_files)

        # Add the raster to the temporary list
        current_mop_rasters[[raster_name_id]] <- raster_obj
      } else if (length(matching_files) == 0) {
        warning(paste("No file found for pattern:", file_partial_path))
      } else {
        warning(paste("Multiple files found for pattern:", file_partial_path, ". Using the first one."))
        raster_obj <- terra::rast(matching_files[1])
        current_mop_rasters[[raster_name_id]] <- raster_obj
      }
    }
    # Add the list of rasters for the current MOP type to the main output list
    output_rasters[[mop_type]] <- current_mop_rasters

    # Combine if mop_type is simple, basic, or a combined type
    if (mop_type %in% c("basic", "simple", "towards_high_combined",
                        "towards_low_combined")) {
      # If the temporary list is not empty, combine the rasters
      if (length(output_rasters[[mop_type]]) > 0) {
        # Combine all SpatRasters into a single multi-layer SpatRaster
        output_rasters[[mop_type]] <- terra::rast(output_rasters[[mop_type]])
      } else {
        output_rasters[[mop_type]] <- NULL # Or an empty SpatRaster
        warning(paste("No rasters found to combine for mop_type:", mop_type))
      }
    }
  }
  }

  #### projection_changes ####
  if(inherits(projection, "changes_projections")){
    #Check change types
    change_out <- setdiff(change_types, c("summary", "by_gcm", "by_change"))
    if(length(change_out) > 0){
      stop("One or more of the 'change_types' provided are not valid.
      Available options are: 'summary', 'by_gcm', and 'by_change'.")
    }

    #Check if root_directory exists
    if(is.null(projection$root_directory)){
      stop("'projection' does not have a root directory.
      Make sure you ran 'projection_changes()' with write_results set to 'TRUE'")
    }

    #Check if root_directory exists:
    if(!file.exists(projection$root_directory)){
      stop("The following root directory does not exist on this system:\n",
       projection$root_directory)
    }

    #List files
    all_files <- list.files(path = projection$root_directory, pattern = ".tif$",
                            recursive = TRUE)
    #Create list to store results
    output_rasters <- list()

    #If by_gcm
    if("by_gcm" %in% change_types){
      #Get file
      f_gcm <- file.path(projection$root_directory, "Changes_by_GCM.tif")
      #Check if file exists
      if(file.exists(f_gcm)){
        output_rasters[["by_gcm"]] <- terra::rast(f_gcm)
      } else {
        warning("Changes_by_GCM.tif not found in root directory")
      }
    }

    #If by_change
    if("by_change" %in% change_types){
      #Get path
      path_change <- file.path(projection$root_directory, "Results_by_change")
      #Check if file exists
      if(file.exists(path_change)){
        #List files
        f_change <- list.files(path_change,pattern = ".tif$")
        #Split
        f_split <- strsplit(f_change, "_")

        #Convert do dataframe to subset periods
        df <- as.data.frame(t(sapply(f_split, `[`, 1:3))) #Convert to dataframe
        colnames(df) <- c("Time", "Period", "ssp") #Rename columns
        df$ssp <- gsub("\\.tif", "", df$ssp) #Remove tif

        #Subset if necessary
        #By period
        all_periods <- c(past_period, future_period)
        if(!is.null(all_periods))
          df <- df[df$Period %in% all_periods,]

        #By scenario
        if(!is.null(future_pscen)){
          df <- df[df$ssp %in% future_pscen | is.na(df$ssp), ]
        }

        #Get files to read
        files_to_read <- paste0(paste(df$Time, df$Period, df$ssp, sep = "_"),
                                ".tif")
        files_to_read <- file.path(projection$root_directory,
                                   "Results_by_change", files_to_read)
        #Read
        r <- lapply(files_to_read, terra::rast)
        names(r) <- paste(df$Time, df$Period, df$ssp, sep = "_")

        #Store results
        output_rasters[["by_change"]] <- r
      } else {
        warning("Changes_by_GCM.tif not found in root directory")
      }
    }

    #If by_gcm
    if("summary" %in% change_types){
      #Get file
      f_summary <- file.path(projection$root_directory, "Changes_summary.tif")
      #Check if file exists
      if(file.exists(f_summary)){
        output_rasters[["Summary"]] <- terra::rast(f_summary)
      } else {
        warning("Changes_summary.tif not found in root directory")
      }
    }

    #Append root directory
    output_rasters[["root_directory"]] <- projection$root_directory

  } #End of changes


  #### Variability ####
  if(inherits(projection, "variability_projections")){

    #Check if root_directory exists
    if(is.null(projection$root_directory)){
      stop("'projection' does not have a root directory.
      Make sure you run 'projection_variability()' with write_files set to 'TRUE'")
    }

    #Check if root_directory exists:
    if(!file.exists(projection$root_directory)){
      stop("The following root directory does not exist on this system:\n",
           projection$root_directory)
    }

    #Get root directory
    rd <- projection$root_directory
    #If present
    if(present){
    present_file <- file.path(rd, "Present.tif")
    if(!file.exists(present_file)){
      warning("Variance from present time not found in root directory")
      v_present <- NULL } else {
      v_present <- terra::rast(present_file)
      names(v_present) <- paste("Present", names(v_present), sep = "_")
    }
    } else {
      v_present <- NULL
    }

    #Subset other times
    other_times <- list.files(rd, pattern = ".tif$")
    #Remove Present
    other_times <- other_times[other_times != "Present.tif"]
    #Split
    f_split <- strsplit(other_times, "_")

    #Convert do dataframe to subset periods
    df <- as.data.frame(t(sapply(f_split, `[`, 1:3))) #Convert to dataframe
    colnames(df) <- c("Time", "Period", "ssp") #Rename columns
    df$ssp <- gsub("\\.tif", "", df$ssp) #Remove tif

    #Subset if necessary
    #By period
    all_periods <- c(past_period, future_period)
    if(!is.null(all_periods))
      df <- df[df$Period %in% all_periods,]

    #By scenario
    if(!is.null(future_pscen)){
      df <- df[df$ssp %in% future_pscen | is.na(df$ssp), ]
    }

    #Get files to read
    files_to_read <- paste0(paste(df$Time, df$Period, df$ssp, sep = "_"),
                            ".tif")
    files_to_read <- file.path(rd, files_to_read)
    #Read
    v_other_times <- terra::rast(lapply(seq_along(files_to_read), function(i){
      r_i <- terra::rast(files_to_read[i])
      names(r_i) <- paste(df$Time[i], df$Period[i], df$ssp[i], names(r_i),
                          sep = "_")
      return(r_i)
      }))

    #Append results
    if(!is.null(v_present)){
      output_rasters <- c(v_present, v_other_times)} else {
        output_rasters <- v_other_times
      }
    } #End of variability_projections

  #Final results
  return(output_rasters)
}
