plot_explore_partition <- function(explore_partition,
                                   data = NULL,
                                   space = "G",
                                   type_of_plot = "mop_distance",
                                   variables = NULL,
                                   use_pca = TRUE,
                                   pcs = c("PC1", "PC2"),
                                   calibration_area = NULL,
                                   include_test_background = FALSE,
                                   distance_pallete = pals::brewer.rdylgn(5),
                                   custom_pallete = NULL,
                                   break_type = "pretty",
                                   in_range_color = "#009E73",
                                   out_range_color = "#D55E00",
                                   calibration_area_col = "gray90",
                                   pr_transparency = 1,
                                   bg_transparency = 0.4,
                                   pch_in_range = 19,
                                   pch_out_range = 18,
                                   cex_plot = 1.4,
                                   size_text_legend = 1,
                                   ncols = NULL){

  #Get spatial data
  d <- explore_partition$Mop_results #Mop results
  cd <- explore_partition$calibration_data
  #Append results
  d <- cbind(d, cd)

  #Get partitions
  partitions <- unique(d$Partition)

  #Make PCA, if necessary
  if(space == "E"){
    # Give warning if data has already been PCA-transformed
    if (!is.null(explore_partition$pca) & use_pca) {
      warning("Variables have already been transformed using PCA; it is not necessary to run PCA again.")
    }

    #Subset variables, if necessary
    if(is.null(explore_partition$pca) && !is.null(variables)){
      cd <- cd[, variables]
    } else {
      variables <- setdiff(colnames(explore_partition$calibration_data),
                           c("pr_bg", explore_partition$categorical_variables))
      cd <- cd[, variables]
    }

    #Run pca
    if(use_pca && is.null(data$pca)){
      cd <- stats::predict(stats::prcomp(cd, center = TRUE, scale = TRUE), cd)
      #Get variables to plot
      v1 <- pcs[1]
      v2 <- pcs[2] } else if (!is.null(explore_partition$pca)) {
        v1 <- pcs[1]
        v2 <- pcs[2]
      } else if (is.null(explore_partition$pca) && !use_pca){
        v1 <- variables[1]
        v2 <- variables[2]
      }

    #Get ranges
    r_v1 <- range(cd[,v1])
    r_v2 <- range(cd[,v2])

    #Merge PCAs
    d <- cbind(d, cd[,c(v1, v2)])
  }



  #Order dataframe
  if(type_of_plot == "mop_distance"){
    d <- d[order(d$mop_distance),]
  } else {
    d <- d[order(d$n_var_out),]
  }

  if(!include_test_background){
    d <- d[d$pr_bg == 1, ]
  }

  #Get n of plots
  n_plots <- length(partitions)

  #### Adjust plot layout ####
  # Define a square-ish layout for the plots
  if(is.null(ncols)){
  n_cols <- ceiling(sqrt(n_plots)) } else {
    n_cols = ncols}
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

  #Define colors for mop distance
  if(is.null(custom_pallete) & type_of_plot == "mop_distance"){
    cores <- distance_pallete
    pal_cores <- grDevices::colorRampPalette(cores)
  } else if (!is.null(custom_pallete) && type_of_plot == "mop_distance"){
    pal_cores <- grDevices::colorRampPalette(custom_pallete)
  }

  #Adjust colors  and legend
  if(type_of_plot == "mop_distance"){
    normalized_values <- (d$mop_distance - min(d$mop_distance)) /
      (max(d$mop_distance) - min(d$mop_distance))
    point_colors <- pal_func(100)[round(normalized_values * 99) + 1]
    point_colors <- stats::setNames(point_colors, d$Partition)
    #Get min and max values
    min_val <- min(d$mop_distance, na.rm = TRUE)
    max_val <- max(d$mop_distance, na.rm = TRUE)
    #Set legend labels
    if(break_type == "pretty"){
      legend_labels <- pretty(min_val:max_val)} else {
        legend_labels <- quantile(d$mop_distance,
                                  probs = c(0, 0.25, 0.5, 0.75, 1))
      }
  } #End of adjust colors if mop_distance

  if(space == "G"){

    #Define calibration area
    if(inherits(calibration_area, "SpatRaster")){
      calibration_area <- calibration_area[[1]]
      calibration_area[!is.na(calibration_area)] <- 0
    }

    #### Make plots for mop_distance ####
    #Plot mop_distance
    if(type_of_plot == "mop_distance"){
      for(i in partitions){
        partition_i <- d[d$Partition == i,]
        #Adjust colors
        cores_i <- point_colors[names(point_colors) == i]
        #Adjust pch
        pch_i <- ifelse(d$inside_range, pch_in_range, pch_out_range)
        #Adjust transparency
        transp <- ifelse(d$pr_bg == 1, pr_transparency, bg_transparency)
        cores_i <- alpha_color(cores_i, transp)

        #Adjust colors if pch is filled
        if(any(pch_i %in% 21:25)){
          bg_plot <- cores_i
          cores_i <- "black"
          bg_legend <- c(in_range_color, out_range_color)
        } else {
          bg_plot <- NULL
          bg_legend <- NA
        }

        #Start plotting partitions
        terra::plot(calibration_area, col = calibration_area_col, legend = FALSE,
                    main = i, ...)
        terra::points(partition_i[, c("x", "y")],
                      col = cores_i, bg = bg_plot, pch = pch_i, cex = cex_plot)

    }
    # Plot legend
    graphics::par(mar = c(0, 0, 0, 0))  # remove margins
    graphics::plot.new()
    SpectrumLegend("top",                             # Legend position
      horiz = TRUE,
      palette = pal_cores,                     # Display our chosen palette
      legend = legend_labels,  # Annotate positions on legend
      title = "Distance",
      bty = "n", cex = size_text_legend)

    }

    #### Make plots for mop_simple ####
    #Adjust colors  and legend
    if(type_of_plot == "mop_simples"){
      for(i in partitions){
        partition_i <- d[d$Partition == i,]

        #Set colors
        cores_i <- ifelse(partition_i$inside_range, in_range_color,
                          out_range_color)
        #Adjust transparency
        transp <- ifelse(d$pr_bg == 1, pr_transparency, bg_transparency)
        cores_i <- alpha_color(cores_i, transp)

        #Adjust pch
        pch_i <- ifelse(partition_i$inside_range, pch_in_range, pch_out_range)

        #Adjust colors if pch is filled
        if(any(pch_i %in% 21:25)){
          bg_plot <- cores_i
          cores_i <- "black"
          bg_legend <- c(in_range_color, out_range_color)
          col_legend <- "black"
        } else {
          bg_plot <- NULL
          bg_legend <- NA
          col_legend <- c(in_range_color, out_range_color)
        }

        #Start plotting
        terra::plot(calibration_area, col = calibration_area_col, legend = FALSE,
                    main = i, ...)
        terra::points(partition_i[, c("x", "y")],
                      col = cores_i, bg = bg_plot, pch = pch_i, cex = cex_plot)

      }
      # Plot legend
      graphics::par(mar = c(0, 0, 0, 0))  # remove margins
      graphics::plot.new()
      legend("top",
             horiz = TRUE,
             col = col_legend,
             pt.bg = bg_legend,
             pch = c(pch_in_range, pch_out_range),
             legend = c("Within range", "Out of range"),
             title = "Within range of training data?",
             bty = "n", cex = size_text_legend)

    }
  } #End of if plot in G

  if(space == "E"){

    #####E - Mop distance####

    #Start plot partitions
    for(i in partitions){
      #Subset calibration data
      d_i <- d[d$Partition == i,]
      #Adjust colors
      cores_i <- point_colors[names(point_colors) == i]
      #Adjust pch
      pch_i <- ifelse(d_i$inside_range, pch_in_range, pch_out_range)
      #Adjust transparency
      transp <- ifelse(d_i$pr_bg == 1, pr_transparency, bg_transparency)
      cores_i <- alpha_color(cores_i, transp)

      #Adjust colors if pch is filled
      if(any(pch_i %in% 21:25)){
        bg_plot <- cores_i
        cores_i <- "black"
        bg_legend <- c(in_range_color, out_range_color)
      } else {
        bg_plot <- NULL
        bg_legend <- NA
      }

      #Plot
      plot(d_i[[v1]], d_i[[v2]], cex = cex_plot,
           col = cores_i, bg = bg_plot,
           pch = pch_i, main = i, xlab = v1, ylab = v2,
           xlim = r_v1, ylim = r_v2, ...)
    }
    # Plot legend
    graphics::par(mar = c(0, 0, 0, 0))  # remove margins
    graphics::plot.new()
    SpectrumLegend("top",                             # Legend position
                   horiz = TRUE,
                   palette = pal_cores,                     # Display our chosen palette
                   legend = legend_labels,  # Annotate positions on legend
                   title = "Distance",
                   bty = "n", cex = size_text_legend)

    } #End of space == E

  #Reset grid
  graphics::par(mfrow = c(1, 1),
                oma = c(0, 0, 0, 0),
                mar = c(5.1, 4.1, 4.1, 2.1)
  )

}

# ?explore_partition_extrapolation
#
# #Prepare data
# # Import occurrences
# data(occ_data, package = "kuenm2")
#
# # Import raster layers
# var <- terra::rast(system.file("extdata", "Current_variables.tif",
#                                package = "kuenm2"))
#
# # Prepare data for maxnet model
# sp_swd <- prepare_data(algorithm = "maxnet", occ = occ_data,
#                        x = "x", y = "y",
#                        raster_variables = var,
#                        species = occ_data[1, 1],
#                        n_background = 100,
#                        categorical_variables = "SoilType",
#                        features = c("l", "lq"),
#                        r_multiplier = 1,
#                        partition_method = "kfolds")
#
# # Analysis of extrapolation risks in partitions
# res <- explore_partition_extrapolation(data = sp_swd,
#                                        include_test_background = TRUE)
#
# explore_partition = res
# space = "G"
# type_of_plot = "mop_distance"
# calibration_area = NULL
# include_test_background = FALSE
# range_pallete = "cols25"
# distance_pallete = pals::brewer.rdylgn(8)
# n_distances = 9
# pr_transparency = 0.75
# bg_transparency = 0.4
# pch = 19
# cex_plot = 1.2
# size_text_legend = 1
