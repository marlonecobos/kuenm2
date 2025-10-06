#' Plot extrapolation risks for partitions
#'
#' @description
#' Visualize data from an `explore_partition` object generated with the
#' `explore_partition_extrapolation` function.
#'
#' @param explore_partition an object of class `explore_partition` returned by
#' the `explore_partition_extrapolation()` function.
#' @param space (character) vector specifying the space to plot. Available
#' options are  'G' for geographical space and E' for environmental space.
#' Default is c("G","E"), meaning both spaces are plotted.
#' @param type_of_plot (character) vector specifying the type(s) of plot.
#' Options are "simple", which shows whether the record in a partition is within
#' the range of the other partitions, and "distance", which shows the Euclidean
#' distance of the record to the set of conditions in the other partitions.
#' Default is c("simple", "distance"), meaning both plots are produced.
#' @param variables (character) A pair of variables used to define the axes of
#' the environmental space. Default is NULL, meaning the first two continuous
#' variables available in `explore_partition` are used to define the E space.
#' @param calibration_area (SpatRaster, SpatVector, or SpatExtent) A spatial
#' object representing the calibration area. Preferably, the same object
#' provided to `prepare_data`. Required only when "G" is included in
#' `type_of_plot`. Default is NULL, meaning the coordinates of the records are
#' used to define the calibration area.
#' @param show_limits (logical) whether to draw a box representing the lower and
#' upper limits of the variables, considering the other partitions (i.e., in
#' Partition 1, the box represents the limits considering Partitions 2, 3, and
#' 4. Only applicable when "E" is included in `type_of_plot`. Default is TRUE.
#' @param include_background (logical) whether to plot background points
#' together with presence records. Only applicable if `explore_partition` was
#' obtained using presence and background points (i.e., with
#' `include_test_background = TRUE` in `explore_partition_extrapolation`).
#' Default is FALSE.
#' @param distance_palette (character) a vector of valid colors used to
#' interpolate a palette for representing distance. Default is NULL, meaning a
#' built-in palette is used (green for lower distances and red for higher
#' distances). Only applicable if "distance" is included in `type_of_plot`.
#' @param break_type (character) specifies the method used to define distance
#' breaks. Options are "pretty" or "quantile". Default is "pretty", which uses
#' the `pretty()` function to set the breaks. Only applicable if "distance" is
#' included in `type_of_plot`.
#' @param in_range_color (character) a color used to represent records that
#' fall within the range of the other partitions. Default is "#009E73" (Seafoam
#' Green).
#' @param out_range_color (character) A color used to represent records that
#' fall outside the range of the other partitions. Default is "#D55E00"
#' (reddish-orange).
#' @param calibration_area_col (character) A color used to represent the
#' calibration area. Default is "gray90".
#' @param pr_alpha (numeric) specifies the transparency of presence records.
#' Default is 1, meaning fully opaque.
#' @param bg_alpha (numeric) specifies the transparency of background points.
#' Default is 0.4. Only applicable if `include_background` is set to TRUE.
#' @param pch_in_range (numeric) specifies the symbol used for points that fall
#' within the range of the other partitions. Default is 21 (filled circle).
#' See `?pch` for other available options.
#' @param pch_out_range (numeric) specifies the symbol used for points that fall
#' outside the range of the other partitions. Default is 24 (filled triangle).
#' See `?pch` for other available options.
#' @param cex_plot (numeric) specifies the size of points in the plot. Default
#' is 1.4
#' @param size_text_legend (numeric) specifies the size of the text in the
#' legend. Default  is 1.
#' @param legend.margin (numeric) specifies the height of the row in the layout
#' that contains the legend. Default is 0.4, meaning the row will be 40% the
#' height of the other rows in the layout.
#' @param lwd_legend (numeric) specifies the width of the legend bar
#' representing distance. Default is 12. Applicable only if "distance" is
#' included in `type_of_plot`.
#' @param ncols (numeric) specifies the number of columns in the plot layout.
#' Default is NULL, meaning the number of columns is determined automatically
#' based on the number of partitions.
#' @param ... additional arguments passed to `plot()`.
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom stats setNames quantile
#' @importFrom terra vect crop plot points
#' @importFrom graphics layout par plot.new plot rect
#' @export
#'
#' @examples
#' # Load prepared_data with spatial blocks as the partitioning method (from ENMeval)
#' data(swd_spatial_block, package = "kuenm2")
#' # Analyze extrapolation risks across partitions
#' res <- explore_partition_extrapolation(data = swd_spatial_block,
#'                                        include_test_background = TRUE)
#' # Plot partition distribution in Geographic Space (Distance and Simple MOP)
#' plot_explore_partition(explore_partition = res, space = "G",
#'                        variables = c("bio_7", "bio_15"))
#'
#' # Plot partition distribution in Environmental Space (Distance and Simple MOP)
#' plot_explore_partition(explore_partition = res, space = "E",
#'                        variables = c("bio_7", "bio_15"))
#'
plot_explore_partition <- function(explore_partition,
                                   space = c("G", "E"),
                                   type_of_plot = c("distance", "simple"),
                                   variables = NULL,
                                   calibration_area = NULL,
                                   show_limits = TRUE,
                                   include_background = FALSE,
                                   distance_palette = NULL,
                                   break_type = "pretty",
                                   in_range_color = "#009E73",
                                   out_range_color = "#D55E00",
                                   calibration_area_col = "gray90",
                                   pr_alpha = 1,
                                   bg_alpha = 0.4,
                                   pch_in_range = 21,
                                   pch_out_range = 24,
                                   cex_plot = 1.4,
                                   size_text_legend = 1,
                                   legend.margin = 0.4,
                                   lwd_legend = 12,
                                   ncols = NULL, ...){
  # Check data #
  if (missing(explore_partition)) {
    stop("Argument 'explore_partition' must be defined.")
  }
  if (!inherits(explore_partition, "explore_partition")) {
    stop("'explore_partition' must be a 'explore_partition' object.")
  }

  space_out <- setdiff(space, c("G", "E"))
  if (length(space_out) > 0) {
    stop("Invalid 'space' provided.",
         "\nAvailable options are: 'E' or 'G'")
  }

  type_out <- setdiff(type_of_plot, c("simple", "distance"))
  if (length(type_out) > 0) {
    stop("Invalid 'type_of_plot' provided.",
         "\nAvailable options are: 'simple' or 'distance'")
  }

  if(!is.null(variables)){
    variables_out <- setdiff(variables,
                             colnames(explore_partition$calibration_data[, -1]))
    if (length(variables_out) > 0) {
      stop("Invalid 'variables' provided.",
           "\nAvailable options are: ",
           paste(colnames(explore_partition$calibration_data[, -1]), collapse = ", "))
    }
    if(any(variables %in% explore_partition$categorical_variables)){
      stop("Categorical variable not supported.",
           "\nPlease provide a continuous variable")
    }
    if(length(variables) != 2){
      stop("Incorrect number of variables. Provide two variables.")
    }

    }

  if(!is.null(calibration_area)){
    if(!inherits(calibration_area, c("SpatVector", "SpatRaster", "SpatExtent"))){
      stop("Argument 'calibration_area' must be a 'SpatVector', 'SpatRaster', 'SpatExtent' or 'NULL'.")
    }
  }

  if (!inherits(include_background, "logical")) {
    stop("'include_background' must be 'logical'.")
  }

  # Check colors
  if(!is.null(distance_palette)){
    valid_palette <- sapply(distance_palette, valid_color)
    if(any(!valid_palette)){
      stop("Invalid colors provided in 'distance_palette'")
    }
  }

  if(!valid_color(in_range_color) | !inherits(in_range_color, "character")){
    stop("Invalid 'in_range_color' provided.",
         "Provide a valid color")
  }
  if(!valid_color(out_range_color) | !inherits(out_range_color, "character")){
    stop("Invalid 'out_range_color' provided.",
         "Provide a valid color")
  }
  if(!valid_color(calibration_area_col) | !inherits(calibration_area_col, "character")){
    stop("Invalid 'calibration_area_col' provided.",
         "Provide a valid color")
  }

  break_out <- setdiff(break_type, c("pretty", "quantile"))
  if (length(break_out) > 0) {
    stop("Invalid 'break_type' provided.",
         "\nAvailable options are: 'pretty' or 'quantile'")
  }

  if(pr_alpha <= 0 | pr_alpha > 1 |
     !inherits(pr_alpha, "numeric")){
    stop("'pr_alpha' must be a positive number between 0 and 1")
  }
  if(bg_alpha <= 0 | bg_alpha > 1 |
     !inherits(bg_alpha, "numeric")){
    stop("'bg_alpha' must be a positive number between 0 and 1")
  }

  if(pch_in_range < 0 | pch_in_range > 25 |
     !inherits(pch_in_range, "numeric") | pch_in_range != trunc(pch_in_range)){
    stop("'pch_in_range' must be a positive integer number between 0 and 25.")
  }
  if(pch_out_range < 0 | pch_out_range > 25 |
     !inherits(pch_out_range, "numeric") | pch_out_range != trunc(pch_out_range)){
    stop("'pch_out_range' must be a positive integer number between 0 and 25.")
  }
  if(pch_in_range %in% 21:25 & !(pch_out_range %in% 21:25))
    stop("You provided a filled symbol for 'pch_in_range'.
Please also provide a filled symbol for 'pch_out_range' (choose from 21, 22, 23, 24, or 25) or set both to non-filled symbols")
  if(pch_out_range %in% 21:25 & !(pch_in_range %in% 21:25))
    stop("You provided a filled symbol for 'pch_out_range'.
Please also provide a filled symbol for 'pch_in_range' (choose from 21, 22, 23, 24, or 25) or set both to non-filled symbols.")

  if(cex_plot <= 0 | !inherits(cex_plot, "numeric")){
    stop("'cex_plot' must be a positive number")
  }
  if(size_text_legend <= 0 | !inherits(size_text_legend, "numeric")){
    stop("'size_text_legend' must be a positive number")
  }
  if(legend.margin <= 0 | !inherits(legend.margin, "numeric")){
    stop("'legend.margin' must be a positive number")
  }
  if(lwd_legend <= 0 | !inherits(lwd_legend, "numeric")){
    stop("'lwd_legend' must be a positive number")
  }
  if(!is.null(ncols)){
    if(ncols <= 0 | !inherits(ncols, "numeric") | ncols != trunc(ncols)){
      stop("'ncols' must be a positive integer number")
  }}

  #Get spatial data
  d <- explore_partition$Mop_results #Mop results
  cd <- explore_partition$calibration_data[
    explore_partition$calibration_data$pr_bg %in% unique(d$pr_bg),]
  cd$pr_bg <- NULL

  #Append results
  d_cd <- cbind(d, cd)

  # Get variables if necessary
  if(is.null(variables)){
    variables <- colnames(cd[-1,])
    # Remove categorical variables if necessary
    if(!is.null(explore_partition$categorical_variables))
      variables <- setdiff(variables, explore_partition$categorical_variables)
    # Get first 2 variables
    variables <- variables[1:2]
  }

  # Get variables
  v1 <- variables[1]
  v2 <- variables[2]

  #Get ranges
  r_v1 <- range(cd[,v1])
  r_v2 <- range(cd[,v2])

  #Get partitions
  partitions <- unique(d_cd$Partition)

  # Include background in the plot?
  if(!include_background){
    d <- d_cd[d_cd$pr_bg == 1, ]
  } else {
    d <- d_cd
  }

  #Get n of plots
  n_plots <- length(partitions)

  #### Adjust plot layout ####
  # Define a square-ish layout for the plots
  if(is.null(ncols) & n_plots > 5){
  n_cols <- ceiling(sqrt(n_plots)) } else if (is.null(ncols) & n_plots <= 5){
    n_cols <- n_plots
  } else {
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

  # Adjust heights: the last row will be smaller (i.e, 20% the height of the plots)
  heights <- c(rep(1, n_rows), legend.margin)

  #Define colors for mop distance
  if(is.null(distance_palette) & "distance" %in% type_of_plot){
    cores <- c("#1A9641", "#A6D96A", "#FFFFBF", "#FDAE61", "#D7191C")
    pal_cores <- grDevices::colorRampPalette(cores)
  } else if (!is.null(distance_palette) && "distance" %in% type_of_plot){
    cores <- distance_palette
    pal_cores <- grDevices::colorRampPalette(distance_palette)
  }

  #Adjust colors  and legend
  if("distance" %in% type_of_plot){
    normalized_values <- (d$mop_distance - min(d$mop_distance)) /
      (max(d$mop_distance) - min(d$mop_distance))
    pal_func <- grDevices::colorRampPalette(cores)
    point_colors <- pal_func(100)[round(normalized_values * 99) + 1]
    point_colors <- stats::setNames(point_colors, d$Partition)
    #Get min and max values
    min_val <- min(d$mop_distance, na.rm = TRUE)
    max_val <- max(d$mop_distance, na.rm = TRUE)
    #Set legend labels
    if(break_type == "pretty"){
      legend_labels <- pretty(min_val:max_val)} else {
        legend_labels <- stats::quantile(d$mop_distance,
                                  probs = c(0, 0.25, 0.5, 0.75, 1))
      }
  } #End of adjust colors if mop_distance

  if("G" %in% space){

    #Define calibration area
    if(inherits(calibration_area, "SpatRaster")){
      calibration_area <- calibration_area[[1]]
      calibration_area[!is.na(calibration_area)] <- 0
    }

    if(is.null(calibration_area)){
      calibration_area <- terra::vect(system.file("extdata", "world.gpkg",
                                                  package = "kuenm2"))
      xy <- explore_partition$Mop_results[, c("x", "y")]
      ext_xy <- c(min(xy$x) - 1, max(xy$x) + 1, min(xy$y) - 1, max(xy$y) + 1)
      calibration_area <- terra::crop(calibration_area, ext_xy)
    }

    #### Make plots for mop_distance ####
    #Plot mop_distance
    if("distance" %in% type_of_plot){
      # Apply the layout
      graphics::layout(mat = layout_matrix_final, heights = heights)
      graphics::par(mar = c(4, 4, 2, 1)) # Standard plot margins for the individual plots

      for(i in partitions){
        partition_i <- d[d$Partition == i,]
        partition_i <- partition_i[order(partition_i$mop_distance),]
        #Adjust colors
        cores_i <- point_colors[names(point_colors) == i]
        #Adjust pch
        pch_i <- ifelse(partition_i$inside_range, pch_in_range, pch_out_range)
        #Adjust transparency
        transp <- ifelse(partition_i$pr_bg == 1, pr_alpha, bg_alpha)
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
    SpectrumLegend("bottom",  # Legend position
      horiz = T,
      palette = pal_cores, # Display our chosen palette
      legend = legend_labels,  # Annotate positions on legend
      title = "MOP Distance", lwd = lwd_legend,
      bty = "n", cex = size_text_legend, adj = c(0.5, 1))

    }

    #### Make plots for mop_simple ####
    #Adjust colors  and legend
    if("simple" %in% type_of_plot){
      # Apply the layout
      graphics::layout(mat = layout_matrix_final, heights = heights)
      graphics::par(mar = c(4, 4, 2, 1)) # Standard plot margins for the individual plots
      for(i in partitions){
        partition_i <- d[d$Partition == i,]
        partition_i <- partition_i[order(partition_i$n_var_out),]
        #Set colors
        cores_i <- ifelse(partition_i$inside_range, in_range_color,
                          out_range_color)
        #Adjust transparency
        transp <- ifelse(partition_i$pr_bg == 1, pr_alpha, bg_alpha)
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

  if("E" %in% space){

    if("distance" %in% type_of_plot){
      # Apply the layout
      graphics::layout(mat = layout_matrix_final, heights = heights)
      graphics::par(mar = c(4, 4, 2, 1)) # Standard plot margins for the individual plots
      #Start plot partitions
      for(i in partitions){
        #Subset calibration data
        d_i <- d[d$Partition == i,]
        d_i <- d_i[order(d_i$mop_distance),]
        # If show limits...
        if(show_limits){
          l1 <- c(min(d_cd[d_cd$Partition != i, v1]),
                  max(d_cd[d_cd$Partition != i, v1]))
          l2 <- c(min(d_cd[d_cd$Partition != i, v2]),
                  max(d_cd[d_cd$Partition != i, v2]))
        }
        #Adjust colors
        cores_i <- point_colors[names(point_colors) == i]
        #Adjust pch
        pch_i <- ifelse(d_i$inside_range, pch_in_range, pch_out_range)
        #Adjust transparency
        transp <- ifelse(d_i$pr_bg == 1, pr_alpha, bg_alpha)
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
        if(show_limits){
          graphics::rect(xleft = l1[1], ybottom = l2[1],
               xright = l1[2], ytop = l2[2],
               border = "red", lty = "dashed",
               lwd = 1.5)
        }
      }
      # Plot legend
      graphics::par(mar = c(0, 0, 0, 0))  # remove margins
      graphics::plot.new()
      SpectrumLegend("top",  # Legend position
                              horiz = T,
                              palette = pal_cores, # Display our chosen palette
                              legend = legend_labels,  # Annotate positions on legend
                              title = "MOP Distance", lwd = lwd_legend,
                              bty = "n", cex = size_text_legend, adj = c(0.5, 1))
    }

    if("simple" %in% type_of_plot){
      # Apply the layout
      graphics::layout(mat = layout_matrix_final, heights = heights)
      graphics::par(mar = c(4, 4, 2, 1)) # Standard plot margins for the individual plots
      # Order dataframe
      for(i in partitions){
        partition_i <- d[d$Partition == i,]
        partition_i <- partition_i[order(partition_i$n_var_out),]
        # If show limits...
        if(show_limits){
          l1 <- c(min(d_cd[d_cd$Partition != i, v1]),
                  max(d_cd[d_cd$Partition != i, v1]))
          l2 <- c(min(d_cd[d_cd$Partition != i, v2]),
                  max(d_cd[d_cd$Partition != i, v2]))
        }
        #Set colors
        cores_i <- ifelse(partition_i$inside_range, in_range_color,
                          out_range_color)
        #Adjust transparency
        transp <- ifelse(partition_i$pr_bg == 1, pr_alpha, bg_alpha)
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

        #Plot
        graphics::plot(partition_i[[v1]], partition_i[[v2]], cex = cex_plot,
             col = cores_i, bg = bg_plot,
             pch = pch_i, main = i, xlab = v1, ylab = v2,
             xlim = r_v1, ylim = r_v2, ...)
        if(show_limits){
          graphics::rect(xleft = l1[1], ybottom = l2[1],
               xright = l1[2], ytop = l2[2],
               border = "red", lty = "dashed",
               lwd = 1.5)
        }

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

    } #End of space == E

  #Reset grid
  graphics::par(mfrow = c(1, 1),
                oma = c(0, 0, 0, 0),
                mar = c(5.1, 4.1, 4.1, 2.1)
  )

}
