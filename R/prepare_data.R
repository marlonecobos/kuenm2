#' Prepare data for model calibration
#'
#' @importFrom terra as.data.frame extract crop prcomp predict
#' @export

prepare_data <- function(occ,
                        species = NULL,
                        x,
                        y,
                        spat_variables,
                        mask = NULL,
                        categorical_variables = NULL,
                        do_pca = FALSE,
                        deviance_explained = 95,
                        min_explainded = 5,
                        exclude_from_pca = NULL,
                        center = TRUE,
                        scale = TRUE,
                        write_pca = FALSE,
                        output_pca = NULL,
                        nbg = 10000,
                        kfolds = 4,
                        weights = NULL,
                        include_xy = TRUE,
                        write_files = FALSE,
                        file_name = NULL,
                        seed = 1){

  if(write_files & is.null(file_name)){
    stop("If write_files = TRUE, you must specify a file_name")
  }

  #Check if weights has the same lenght of occ
  if(!is.null(weights)){
    if(nrow(occ) != length(weights)){
  stop("length of weights does not match the number of occurrences in occ")
    }}


  #Extract name of the specie
  if(is.null(species)) {
    sp_name <- as.character(occ[1, "species"])
  } else {
    sp_name <- species
  }

  #Mask variables?
  if(!is.null(mask)){
    spat_variables <- terra::crop(spat_variables, mask, mask = TRUE)
  }

  #### DO PCA? ####
  if(do_pca){
    #Select variables
    if(!is.null(exclude_from_pca)){
      var_to_pca <- setdiff(names(spat_variables), exclude_from_pca)
    } else {var_to_pca <- names(spat_variables)}

    #Do pca
    pca <- terra::prcomp(spat_variables[[var_to_pca]], center = center,
                         scale = scale)
    pca$x <- NULL #Remove matrix: unnecessary
    pca$vars_in <- var_to_pca #Vars included in PCA
    pca$vars_out <- exclude_from_pca #Vars to do not include in PCA
    #Get deviance explained
    d_exp <- cumsum(pca$sdev/sum(pca$sdev)) * 100
    if(!is.null(min_explainded)){
      d_exp <- d_exp[(pca$sdev/sum(pca$sdev) * 100) > min_explainded]
    }
    #N axis
    if(max(d_exp) > deviance_explained){
    ind_exp <- min(which(d_exp >= deviance_explained))} else {
      ind_exp <- length(d_exp)
    }
    env <- predict(spat_variables[[var_to_pca]], pca, index = 1:ind_exp)
    if(!is.null(exclude_from_pca)){
    env <- c(env, spat_variables[[exclude_from_pca]])}

    #write PCA?
    if(write_pca){
      terra::writeRaster(env, paste0(output_pca, "/", "pca_var.tif"))
    }

  } else {
    env <- spat_variables
    pca <- NULL}


  #Extract variables to occurrences
  xy <- as.matrix(occ[,c(x,y)])
  colnames(xy) <- c("x", "y")
  occ_var <- cbind(xy, terra::extract(x = env, y = xy))
  occ_var$pr_bg <- 1

  #Get background points
  cell_samp <- terra::as.data.frame(env[[1]], na.rm = TRUE,
                                    cells = TRUE)[, "cell"]
  if(length(cell_samp) < nbg){
    nbg <- length(cell_samp)
  }
  cell_samp <- sample(cell_samp, size = nbg, replace = FALSE,
                      prob = NULL)
  bg_var <- terra::extract(x = env, y = cell_samp, xy = TRUE)
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
    warning(length(na_rows), " rows were excluded from database because NAs were found")
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
               weights = weights,
               pca = pca)

  #Save results?
  if(write_files) {
    saveRDS(data, paste0(file_name, ".RDS"))
  }
  return(data)
}
