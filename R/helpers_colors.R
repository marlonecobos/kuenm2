#Set rgb colors
set_rgb_adjusted <- function(l){
  # Apply the function to each row of the dataframe
  sapply(1:nrow(l), function(x){
    col_x <- l$color[x]    # Original color (hexcode or color name)
    alpha_x <- l$alpha[x]  # Original color's alpha (0-1)

    # Convert the base color to RGB (0-255).
    col_rgb_base <- grDevices::col2rgb(col = col_x)

    # Normalize the base RGB components to 0-1
    r_base <- col_rgb_base[1] / 255
    g_base <- col_rgb_base[2] / 255
    b_base <- col_rgb_base[3] / 255

    # Define the white background (normalized to 0-1)
    # Calculate the new RGB components adjusted for perceived opacity.
    # This uses the same logic as the hex2RGB_custom function we created earlier.
    # C' = C * alpha + F * (1 - alpha)
    r_adjusted <- r_base * alpha_x + 1 * (1 - alpha_x)
    g_adjusted <- g_base * alpha_x + 1 * (1 - alpha_x)
    b_adjusted <- b_base * alpha_x + 1 * (1 - alpha_x)

    # Convert the adjusted RGB components (0-1) back to 0-255,
    # round them, and format them to two-digit hexadecimal.
    # Then, create the new color with full alpha (255).
    new_col_hex <- grDevices::rgb(
      red = round(r_adjusted * 255),
      green = round(g_adjusted * 255),
      blue = round(b_adjusted * 255),
      alpha = 255, # Fixed alpha at 255 (fully opaque)
      maxColorValue = 255,
      names = NULL
    )
    return(new_col_hex)
  })
}


#Helper to ser colors by resul change####
set_colors_by_change <- function(x_i, change, color, min_alpha, max_alpha){
  #Get levels
  l_i <- terra::levels(x_i)[[1]]
  colnames(l_i) <- "value"
  #Get number of gcms
  l_i$n_gcms <- gsub("\\D+", "", l_i[,2])
  #Set alphas
  a <- data.frame(n_gcms = 0:max(l_i$n_gcms),
                  alpha = seq(from = min_alpha, to = max_alpha,
                              length.out = max(as.numeric(l_i$n_gcms)) + 1))
  l_i <- merge(l_i, a)
  #Set colors
  l_i$color <- color
  #Adjust alpha: 0.1 if 0 gcms
  l_i$alpha[l_i$n_gcms == 0] <- 0.09
  #Set colors based on alpha
  l_i$rgb <- set_rgb_adjusted(l_i)
  terra::coltab(x_i) <- data.frame(value = l_i$value,
                                   col = l_i$rgb)
  return(x_i)
}

alpha_color <- function(colors, alpha) {
  # 1. Validate the alpha value
  # Ensures that alpha is between 0 and 1
  if (any(alpha < 0 | alpha > 1)) {
    stop("The 'alpha' value must be between 0 and 1.")
  }

  # 2. Convert input colors to 6-digit hexadecimal RGB format
  # col2rgb() returns a matrix with R, G, B components.
  # rgb() with maxColorValue=255 and no alpha channel.
  # We use t(col2rgb(colors)) to transpose it for easier iteration.
  rgb_matrix <- t(col2rgb(colors))
  hex_rgb <- apply(rgb_matrix, 1, function(color_row) {
    grDevices::rgb(color_row[1], color_row[2], color_row[3], maxColorValue = 255)
  })

  # 3. Convert the alpha value to hexadecimal format (00 to FF)
  # The alpha value (0-1) is scaled to 0-255 and rounded to the nearest integer.
  # We use sprintf to format the hexadecimal value, ensuring two digits and zero-padding.
  alpha_hex <- sapply(alpha * 255, function(alpha_val) {
    sprintf("%02X", round(alpha_val))
  })

  # 4. Combine the RGB hexadecimal color with the alpha hexadecimal
  # The rgb() function already adds the '#' prefix, so we just need to append the alpha.
  # We'll need to remove the alpha part if it already exists from any input colors,
  # but for standard color names/hex_rgb from rgb(), this won't be an issue.
  colors_with_alpha <- paste0(substr(hex_rgb, 1, 7), alpha_hex)

  return(colors_with_alpha)
}

#Get hex colors from pals package
# kuenm2_discrete_palletes <- list(
#   alphabet = pals::alphabet(),
#   alphabet2 = pals::alphabet2(),
#   cols25 = pals::cols25(n = 25),
#   glasbey = pals::glasbey(n = 32),
#   kelly = pals::kelly(n = 22),
#   polychrome = pals::polychrome(n = 36),
#   stepped = pals::stepped(n = 24),
#   stepped2 = pals::stepped2(n = 20),
#   stepped3 = pals::stepped3(n = 20),
#   okabe = pals::okabe(n = 8),
#   tableau20 = pals::tableau20(n = 20),
#   tol = pals::tol(n = 12),
#   tol.groundcover = pals::tol.groundcover(n = 14),
#   trubetskoy = pals::trubetskoy(n = 22),
#   watlington = pals::watlington(n = 16)
# )
# usethis::use_data(kuenm2_discrete_palletes)

#Add spectrum legend from PlotTools r Package
SpectrumLegend <- function(
    x = "topright", ...,
    palette,
    legend,
    lty = 1, lwd = 4,
    bty = "o",
    adj = if (horiz) c(0.5, 0.5) else c(0, 0.5),
    horiz = FALSE,
    lend = "butt",
    cex = 1,
    seg.len = 1
) {
  if (is.function(palette)) {
    palette <- palette(256)
  }
  nCol <- length(palette)
  if (nCol < 1) {
    stop("palette has length zero")
  }

  lgd <- legend(x = x,
                legend = legend,
                horiz = horiz,
                adj = adj,
                cex = cex,
                bty = ifelse(horiz, "n", bty),
                lty = 0, ncol = 1,
                seg.len = seg.len,
                ...)
  textXY <- lgd[["text"]]

  Cex <- cex * par("cex")
  xyc <- xyinch(par("cin"), warn.log = FALSE)

  if (horiz) {
    xEnds <- range(textXY[["x"]])
    yc <- Cex * xyc[[2]]
    barSpace <- yc
    yEnds <- textXY[["y"]][c(1, 1)] - barSpace

    # as not plotting lines:
    lgd[[c("rect", "left")]] <-  lgd[[c("rect", "left")]] + (barSpace / 2)
    lgd[[c("rect", "h")]] <-  lgd[[c("rect", "h")]] + barSpace

    if (bty == "o") {
      .DrawBox(lgd[["rect"]], ...)
    }
  } else {
    xc <- Cex * xyc[[1]]
    xEnds <- textXY[["x"]][c(1, 1)] - xc - (xc * seg.len / 2)
    yEnds <- range(textXY[["y"]])
  }

  .DrawLegend(xEnds, yEnds, nCol, palette, lwd, lty, lend)

  # Return:
  invisible(lgd)
}
