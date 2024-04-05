#### Prepare SWD ####
#library(terra)
prepare_data <- function(occ,
                        species = NULL,
                        x,
                        y,
                        spat_variables,
                        categorical_variables = NULL,
                        nbg = 10000,
                        kfolds = 4,
                        weights = NULL,
                        include_xy = TRUE,
                        write_files = F,
                        file_name = NULL,
                        seed = 1,
                        verbose = TRUE){
  #Check if weights has the same lenght of occ
  if(!is.null(weights)){
    if(nrow(occ) != length(weights)){
  stop("length of weights does not match the number of occurrences in occ")
    }}


  #Extract name of the specie
  if(!is.null(species)) {
    sp_name <- as.character(occ[1, "species"])
  } else {
    sp_name <- NULL
  }


  #Extract variables to occurrences
  xy <- as.matrix(occ[,c(x,y)])
  colnames(xy) <- c("x", "y")
  occ_var <- cbind(xy, terra::extract(x = spat_variables, y = xy))
  occ_var$pr_bg <- 1

  #Get background points
  cell_samp <- terra::as.data.frame(spat_variables[[1]], na.rm = TRUE,
                                    cells = TRUE)[, "cell"]
  cell_samp <- sample(cell_samp, size = nbg, replace = FALSE,
                      prob = NULL)
  bg_var <- terra::extract(x = spat_variables, y = cell_samp, xy = TRUE)
  bg_var$pr_bg <- 0
  #Join occ and background
  occ_bg <- rbind(occ_var, bg_var)
  #Reorder columns
  occ_bg <- occ_bg[,c(1, 2, which(names(occ_bg) == "pr_bg"),
                      (3:(ncol(occ_bg)-1))
                      [-which(names(occ_bg) == "pr_bg")])]

  #Remove NA from data
  #Identify rows with NA
  na_rows <- which(!complete.cases(occ_bg))
  if(length(na_rows) > 0){
    occ_bg <- occ_bg[-na_rows,] #Remove NAs from calibration data
    #Remove NAs from weights
    if(!is.null(weights)){
      weights <- weights[-na_rows]
    }
    if(verbose){
      message(length(na_rows), " rows were excluded from database because NAs were found")
    }
  }

  #Check if weights and calibration data have the same lenght
  if(!is.null(weights)) {
    if(nrow(occ_bg) != length(weights)) {
    stop("length of weights does not match number of rows in calibration_data")
  }}


  if(include_xy) {
  occ_bg_xy <- occ_bg[, 1:2]
  } else {
    occ_bg_xy <- NULL
  }

  occ_bg <- occ_bg[, -(1:2)]

  #Convert categorical variables to factor
  if(!is.null(categorical_variables)) {
    occ_bg[categorical_variables] <- lapply(occ_bg[categorical_variables],
                                            factor)
      }

  #Split in kfolds?
  k_f <- enmpa::kfold_partition(data = occ_bg,
                                  dependent = "pr_bg",
                                  k = kfolds, seed = seed)

  data <- list(species = sp_name,
               calibration_data = occ_bg,
               kfolds = k_f,
               data_xy = occ_bg_xy,
               categorical_variables = categorical_variables,
               weights = weights)

  #Save results?
  if(write_files) {
    saveRDS(data, paste0(file_name, ".RDS"))
  }
  return(data)
}
