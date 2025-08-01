#' Explore the Distribution of Partitions in Environmental Space
#'
#' @description
#' Plots training and testing data (presences and backgrounds) in a
#' two-dimensional environmental space. This space can be defined either by
#' performing a PCA on all environmental variables or by specifying two
#' environmental variables manually.
#'
#' @usage explore_partition_env(data, show_unused_data = FALSE,
#'                              raster_variables = NULL, mask = NULL,
#'                              variables = NULL, type_of_plot = "combined",
#'                              use_pca = TRUE, pcs = c("PC1", "PC2"),
#'                              partition_palette = "cols25",
#'                              custom_partition_palette = NULL,
#'                              include_test_background = TRUE,
#'                              pr_train_col = "#009E73",
#'                              pr_test_col = "#D55E00",
#'                              bg_train_col = "grey",
#'                              bg_test_col = "#56B4E9", pr_transparency = 0.75,
#'                              bg_transparency = 0.4, pch = 19, cex_plot = 1.2,
#'                              size_text_legend = 1, ...)
#'
#' @param data an object of class `prepared_data` returned by the prepare_data() function
#' @param show_unused_data (logical) whether to plot the distribution of
#' environmental conditions that are not represented by the background points.
#' If set to TRUE, the `raster_variables` must be provided. Only applicable when
#' `type_of_plot = "combined"`.
#' @param raster_variables a `SpatRaster` object representing the predictor
#' variables used to calibrate the models. Preferably the same object used in
#' `prepare_data`. Required only when `show_unused_data = TRUE`. Default is NULL.
#' @param mask (SpatRaster, SpatVector, or SpatExtent) spatial object used to
#' mask `raster_variables` to the area where the model will be calibrated.
#' Preferably the same object used in `prepare_data` (if applicable). Only used
#' when `show_unused_data = TRUE`. Default is NULL.
#' @param variables (character) names of the variables in `data` to define the
#' two-dimensional environmental space. If `use_pca = TRUE`, these variables
#' will be used to perform the PCA. If `use_pca = FALSE`, this must be a
#' character vector with exactly two variable names (e.g.,
#' `c("bio_1", "bio_12")`). Default is NULL, meaning all variables in `data`
#' will be used.
#' @param type_of_plot (character) the type of plot. Options are "combined" and
#' "individual". See details. Default is "combined".
#' @param use_pca (logical) whether to use PCA variables to define the
#' environmental space. If TRUE, a PCA will be performed on the variables,
#' unless `data` already includes a PCA object from using
#' `prepare_data(do_pca = TRUE)`. Default is TRUE.
#' @param pcs (character) the two PCA axes to use to define the two-dimensional
#' environmental space. Default is `c("PC1", "PC2")`, meaning the first two axes
#' will be used. Only applicable if `use_pca = TRUE`.
#' @param partition_palette (character) the color palette used to color the
#' different partitions. See `?kuenm2_discrete_palettes` to check available
#' options. Default is `"cols25"`.
#' @param custom_partition_palette (character) a character vector defining
#' custom colors for the different partitions. The number of values must match
#' the number of partitions in `data`. Default is NULL, meaning the palette
#' defined in `partition_palette` will be used.
#' @param include_test_background (logical) whether to include background points
#' that were not used for training when plotting individual partition plots.
#' Default is TRUE.
#' @param pr_train_col (character) the color used for train records in the
#' individual plots. Default is "009E73".
#' @param pr_test_col (character) the color used for test records in the
#' individual plots. Default is "D55E00".
#' @param bg_train_col (character) the color used for train backgrounds in the
#' individual plots. Default is "56B4E9".
#' @param bg_test_col (character) the color used for test backgrounds in the
#' individual plots. Default is "gray". Only applicable if
#' `include_test_background = TRUE`.
#' @param pr_transparency (numeric) a value between 0 (fully transparent) and 1
#' (fully opaque) defining the transparency of the points representing presences.
#' Default is 0.75.
#' @param bg_transparency (numeric) a value between 0 (fully transparent) and 1
#' (fully opaque) defining the transparency of the points representing
#' background points. Default is 0.4.
#' @param pch (numeric) a value between 1 and 25 to specify the point shape. See
#' `?pch` for details. Default is `19` (solid circle).
#' @param cex_plot (numeric) specify the size of the points in the plot. Default
#' is `1.2`.
#' @param size_text_legend (numeric) specify the size of the text of the legend.
#' Default is `1`.
#' @param ... additional arguments passed to `plot`.
#'
#' @details
#' The function provides two types of plots:
#'
#' - **combined**: two plots side by side, one showing the presences and another
#' showing the background points. The colors of the points represent the
#' partitions. This is the default option.
#'
#' - **individual**: one plot per partition. In each plot, the colors of the
#' points represent those used as train records, test records, train background,
#' or test background (i.e., not used during training in the specified
#' partition).
#'
#' To obtain both types of plots, set:
#' `type_of_plot = c("combined", "individual")`.
#'
#' @importFrom terra extract as.data.frame crop
#' @importFrom stats setNames
#' @importFrom grDevices rgb
#' @importFrom stats predict prcomp
#' @importFrom graphics points par plot.new legend
#' @return
#' Plots showing the training and testing data in a two-dimensional
#' environmental space.
#' @export
#'
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
#' # Explore the Distribution of Partitions in Environmental Space
#' explore_partition_env(data = sp_swd, show_unused_data = TRUE,
#'                       raster_variables = var,
#'                       type_of_plot = c("combined", "individual"))

explore_partition_env <- function(data,
                                  show_unused_data = FALSE,
                                  raster_variables = NULL,
                                  mask = NULL,
                                  variables = NULL,
                                  type_of_plot = "combined",
                                  use_pca = TRUE,
                                  pcs = c("PC1", "PC2"),
                                  partition_palette = "cols25",
                                  custom_partition_palette = NULL,
                                  include_test_background = TRUE,
                                  pr_train_col = "#009E73",
                                  pr_test_col = "#D55E00",
                                  bg_train_col = "grey",
                                  bg_test_col = "#56B4E9",
                                  pr_transparency = 0.75,
                                  bg_transparency = 0.4,
                                  pch = 19,
                                  cex_plot = 1.2,
                                  size_text_legend = 1,
                                  ...){
  if (!is.null(mask) &
      !inherits(mask, c("SpatRaster", "SpatVector", "SpatExtent"))) {
    stop("Argument 'mask' must be a 'SpatVector', 'SpatExtent' or 'SpatRaster'.")
  }

  #Check data
  if (missing(data)) {
    stop("Argument 'data' must be defined.")
  }
  if (!is.null(raster_variables) && !inherits(raster_variables, "SpatRaster")) {
    stop("Argument 'raster_variables' must be a 'SpatRaster' or 'NULL'")
  }

  if(!inherits(data, "prepared_data")){
    stop("'data' must be a 'prepared_data' object.")
  }

  if(show_unused_data){
    if(is.null(raster_variables)){
      stop("If 'show_unused_data = TRUE', you must provide the 'raster_variables'")
    }
    if(is.null(data$data_xy)){
      stop("If 'show_unused_data = TRUE', the 'data' must have the 'xy' element.
Make sure you run 'prepare_data()' with 'include_xy' set to TRUE")
    }
  }

  if(!any(type_of_plot %in% "combined")){
    show_unused_data <- FALSE
    warning("'show_unused_data = TRUE' is only applycable when type_of_plot = 'combined'")
  }

  if(!inherits(type_of_plot, "character"))
    stop("'type_of_plot' must be a 'character'")
  type_out <- setdiff(type_of_plot, c("combined", "individual"))
  if(length(type_out) > 0)
    stop("The 'type_of_plot' is not valid.
Available options are 'combined' or 'individual'")

  if(length(pr_transparency) > 1 || !inherits(pr_transparency, "numeric") || pr_transparency < 0 ||
     pr_transparency > 1)
    stop("'pr_transparency' must be a single numeric value between 0 and 1")

  if(length(bg_transparency) > 1 || !inherits(bg_transparency, "numeric") || bg_transparency < 0 ||
     bg_transparency > 1)
    stop("'bg_transparency' must be a single numeric value between 0 and 1")

  if(length(pch) > 1 || !inherits(pch, "numeric") || pch < 0 ||
     pch > 25)
    stop("'pch' must be a single numeric value between 0 and 25. Check '?pch' for details")

  if(!inherits(cex_plot, "numeric")){
    stop("'cex_plot' must be 'numeric'")
  }

  if(!inherits(size_text_legend, "numeric")){
    stop("'size_text_legend' must be 'numeric'")
  }

  #Check variables
  if(!is.null(variables)){
    if(!inherits(variables, "character"))
      stop("'variables' must be 'character' or 'NULL'")
    var_out <- setdiff(variables, colnames(data$calibration_data))
    if(length(var_out) > 0)
      stop("Some 'variables' provided are not present in the 'data': ",
           paste(var_out, collapse = " ,"))
  }

  if(!use_pca && length(variables) != 2)
    stop("If 'use_pca = FALSE', you must provide 2 'variables'")

  if(is.null(custom_partition_palette) &&
     !(partition_palette %in% names(kuenm2::kuenm2_discrete_palletes)))
    stop("Invalid 'partition_palette'. Check the available options in '?kuenm2_palletes'")

  if(!is.null(custom_partition_palette) &&
     length(custom_partition_palette) < length(data$part_data))
    stop("Insufficient number of colors in 'custom_partition_palette'.
Provide at least ", length(data$part_data), " colors.")

  if(use_pca && length(pcs) > 2){
    stop("If 'use_pca = TRUE', you must provide two strings to select the PCA-axis (e.g., c('PC1', 'Pcd')")
  }

  #Remove categorical
  if(!is.null(data$categorical_variables)){
    cd <- data$calibration_data[, setdiff(colnames(data$calibration_data),
                                        c(data$categorical_variables, "pr_bg"))]
  } else {
    cd <- data$calibration_data[, setdiff(colnames(data$calibration_data),
                                          "pr_bg")]
  }

  #If show_unused_data = TRUE...
  if(show_unused_data){
    id_rasters <- terra::extract(raster_variables[[1]], data$data_xy,
                                 cells = TRUE)[["cell"]]
    #Convert to daframe only the part of spatraster not covevered by xy
    df_r <- na.omit(terra::as.data.frame(raster_variables[-id_rasters],
                                         na.rm = TRUE))
  }


  # Give warning if data has already been PCA-transformed
  if (!is.null(data$pca) & use_pca) {
    warning("Variables have already been transformed using PCA; it is not necessary to run PCA again.")
  }


  #Subset variables, if necessary
  if(is.null(data$pca) && !is.null(variables)){
    cd <- cd[, variables]
  } else {
    variables <- colnames(cd)
  }


  #Run PCA
   if(use_pca && is.null(data$pca)){
    cd <- stats::predict(stats::prcomp(cd, center = TRUE, scale = TRUE), cd)
    #Get variables to plot
    v1 <- pcs[1]
    v2 <- pcs[2]
    #If show unused data
    if(show_unused_data){
    df_r <- stats::predict(stats::prcomp(data$calibration_data[, variables],
                                         center = TRUE, scale = TRUE), df_r)}
    } else if (!is.null(data$pca)) {
     v1 <- pcs[1]
     v2 <- pcs[2]
     if(show_unused_data){
       df_r <- stats::predict(data$pca, df_r)}
     } else if (is.null(data$pca) && !use_pca){
      v1 <- variables[1]
      v2 <- variables[2]
   }


  # Get partitions
  np <- names(data$part_data)
  n_plots <- length(np)


  #If combined plot#

  if(any(type_of_plot == "combined")){
    #Get partitions data
    res_test <- data.frame()
    for(i in np){
      test_i <- data.frame("Partition" = i, "test_train" = "test",
                           "pr_bg" = data$calibration_data[data$part_data[[i]],
                                                           "pr_bg"],
                           cd[data$part_data[[i]], ])
      res_test <- rbind(res_test, test_i) }
    #Split in background and presences
    pr_data <- res_test[res_test$pr_bg == 1, ]
    bg_data <- res_test[res_test$pr_bg == 0, ]
    if(show_unused_data){
      df_r <- data.frame(cbind("Partition" = "Unused data",
                               "test_train" = "Unused data",
                               "pr_bg" = "Unused data", df_r))
    }

    #Adjust layout
    h <- 0.13 * ceiling(n_plots / 5)

    layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, h))
    par(mar = c(4, 4, 2, 1)) # Standard plot margins for the individual plots

    #Get x and y limits
    xlimits <- range(res_test[[v1]])
    ylimits <- range(res_test[[v2]])

    #Adjust colors
    if(is.null(custom_partition_palette)){
      cores <- kuenm2::kuenm2_discrete_palletes[[partition_palette]]
    } else {
      cores <- custom_partition_palette
    }

    #Adjust colors if pch if filled
    if(pch %in% 21:25){
      col_plot <- "black"
      bg_plot <- alpha_color(cores[factor(pr_data$Partition)], pr_transparency)
      bg_legend <- unique(stats::setNames(cores[factor(pr_data$Partition)],
                                          factor(pr_data$Partition)))
    } else {
      col_plot <- alpha_color(cores[factor(pr_data$Partition)], pr_transparency)
      bg_plot <- NULL
      bg_legend <- NA
    }

    #Plot
    plot(pr_data[[v1]], pr_data[[v2]], cex = cex_plot,
         col = col_plot, bg = bg_plot,
         pch = pch, main = "Presences", xlab = v1, ylab = v2,
         xlim = xlimits, ylim = ylimits, ...)
    if(show_unused_data){
      plot(df_r[[v1]], df_r[[v2]], cex = cex_plot,
           col = alpha_color("gray50", bg_transparency),
           pch = 19, main = "Background", xlab = v1, ylab = v2, ...)
      graphics::points(bg_data[[v1]], bg_data[[v2]], cex = cex_plot,
             col = col_plot, bg = bg_plot,
             pch = pch, main = "Background", xlab = v1, ylab = v2)
      # Plot legend
      graphics::par(mar = c(0, 0, 0, 0))  # remove margins
      graphics::plot.new()
      graphics::legend("center",
             legend = unique(c("Unused data", unique(bg_data$Partition))),
             col = c("gray50", col_plot),
             #pt.bg = unique(bg_legend),
             pch = 19,
             bty = "n",
             cex = size_text_legend, ncol = ifelse(n_plots + 1 < 5, n_plots, 5))
    } else {
      #If do not show unused data
    plot(bg_data[[v1]], bg_data[[v2]], cex = cex_plot,
         col = col_plot, bg = bg_plot,
         pch = pch, main = "Background", xlab = v1, ylab = v2, ...)
    # Plot legend
    graphics::par(mar = c(0, 0, 0, 0))  # remove margins
    graphics::plot.new()
    graphics::legend("center", legend = unique(factor(bg_data$Partition)),
           col = col_plot, pt.bg = bg_legend,
           pch = pch,
           bty = "n",
           cex = size_text_legend, ncol = ifelse(n_plots < 5, n_plots, 5))
    }
    #Reset grid
    graphics::par(mfrow = c(1, 1),
        oma = c(0, 0, 0, 0),
        mar = c(5.1, 4.1, 4.1, 2.1)
    )
  }

  if(any(type_of_plot == "individual")){
    #### Adjust plot layout ####
    # Define a square-ish layout for the plots
    n_cols <- ceiling(sqrt(n_plots))
    n_rows <- ceiling(n_plots / n_cols)

    # Calculate the total number of cells in the main plot grid
    total_grid_cells <- n_rows * n_cols

    # Initialize a vector for all cells
    # Assign 0 initially to all cells
    cell_ids <- rep(0, total_grid_cells)

    # Assign plot IDs (1 to n_plots) to the first n_plots cells
    cell_ids[1:n_plots] <- 1:n_plots

    # Define the ID for the legend
    legend_cell_id <- n_plots + 1

    # Assign the legend ID to any remaining cells in the main grid that are 0
    # These are the "empty" cells AFTER all plots have been assigned their spots
    cell_ids[cell_ids == 0] <- legend_cell_id

    # Reshape this sequence into the main plot grid matrix
    # Use byrow = TRUE to fill row by row, matching the original intention
    plot_layout_matrix <- matrix(cell_ids, nrow = n_rows, ncol = n_cols, byrow = TRUE)
    plot_layout_matrix[plot_layout_matrix > n_plots] <- 0

    # Add an extra row specifically for the legend (occupying all columns)
    # This guarantees the legend gets its own dedicated line below the plots.
    layout_matrix_final <- rbind(plot_layout_matrix, rep(legend_cell_id, n_cols))

    # Adjust heights: the last row will be smaller (20% the height of the plots)
    heights <- c(rep(1, n_rows), 0.2)

    # Apply the layout
    layout(mat = layout_matrix_final, heights = heights)
    par(mar = c(4, 4, 2, 1)) # Standard plot margins for the individual plots

    levels_ordenados <- c("train_background", "test_background",
                          "train_presence", "test_presence")
    col_to_use <- c(bg_train_col, bg_test_col, pr_train_col, pr_test_col)
    col_to_use <- stats::setNames(col_to_use, levels_ordenados)

    for(i in np){

      test_i <- data.frame("Partition" = i, "test_train" = "test",
                           "pr_bg" = data$calibration_data[data$part_data[[i]],
                                                           "pr_bg"],
                           cd[data$part_data[[i]], ])

      train_i <- data.frame("Partition" = i, "test_train" = "train",
                            "pr_bg" = data$calibration_data[-data$part_data[[i]],
                                                            "pr_bg"],
                            cd[-data$part_data[[i]], ])
      # Merge data
      test_train_i <- rbind(train_i, test_i)
      test_train_i$pr_bg <- ifelse(test_train_i$pr_bg == 1, "presence",
                                   "background")
      test_train_i$category <- factor(paste(test_train_i$test_train,
                                            test_train_i$pr_bg, sep = "_"),
                                      levels = levels_ordenados)
      # Define levels of categories
      test_train_i$category <- factor(test_train_i$category,
                                      levels = levels_ordenados)
      # Sort levels
      test_train_i <- test_train_i[order(test_train_i$category), ]

      # Remove test background, if necessary
      if(!include_test_background){
        test_train_i <- test_train_i[test_train_i$category != "test_background",]
        test_train_i$category <- droplevels(test_train_i$category)
        col_to_use <- col_to_use[names(col_to_use) != "test_background"]
      }

      #Set transparency
      transp <- ifelse(test_train_i$pr_bg == "presence", pr_transparency,
                       bg_transparency)


      #Adjust colors if pch is filled
      if(pch %in% 21:25){
        col_plot <- "black"
        bg_plot <- alpha_color(col_to_use[as.numeric(test_train_i$category)],
                               transp)
        bg_legend <- stats::setNames(cores[factor(pr_data$Partition)],
                                     factor(pr_data$Partition))
      } else {
        col_plot <- alpha_color(col_to_use[as.numeric(test_train_i$category)],
                                transp)
        bg_plot <- NULL
        bg_legend <- NA
      }


      plot(test_train_i[[v1]], test_train_i[[v2]],
           col = col_plot, bg = bg_plot,
           pch = pch,
           cex = cex_plot,
           xlab = v1,
           ylab = v2,
           main = i, ...)
    }

    # Plot legend
    graphics::par(mar = c(0, 0, 0, 0))  # remove margins
    graphics::plot.new()
    graphics::legend("center",
           legend = levels(test_train_i$category),
           col = unique(col_plot), pt.bg = bg_legend,
           pch = pch,
           horiz = TRUE,
           bty = "n",
           cex = size_text_legend, pt.cex = cex_plot)
    #Reset grid
    graphics::par(mfrow = c(1, 1),
        oma = c(0, 0, 0, 0),
        mar = c(5.1, 4.1, 4.1, 2.1)
    )
  }

  return(invisible(NULL))

}
