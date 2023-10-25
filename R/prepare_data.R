#### Prepare SWD ####
library(terra)
occ <- read.csv("Models/Myrcia_hatschbachii/Occurrences.csv")
species <- "Myrcia_hatschbachii"
x = "x"
y = "y"
spat_variables <- rast("Models/Myrcia_hatschbachii/PCA_var.tiff")
nbg = 10000
kfolds = 4
include_xy = TRUE
writeFiles = F
out_dir = NULL

#library(terra)
prepare_data <- function(occ,
                        species = NULL,
                        x,
                        y,
                        spat_variables,
                        categorical_variables = NULL,
                        nbg = 10000,
                        kfolds = 4,
                        include_xy = TRUE,
                        write_files = F,
                        file_name = NULL,
                        seed = 1,
                        verbose = TRUE){
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
  n_before <- nrow(occ_bg)
  occ_bg <- na.omit(occ_bg )
  n_after <- n_before - nrow(occ_bg)
  if(verbose & n_after > 0){
    message(n_after, " rows were excluded from database because NAs were found")
  }


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
                                  k = 4, seed = seed)

  data <- list(species = sp_name,
               calibration_data = occ_bg,
               kfolds = k_f,
               data_xy = occ_bg_xy,
               categorical_variables = categorical_variables)

  #Save results?
  if(write_files) {
    saveRDS(data, paste0(file_name, ".RDS"))
  }
  return(data)
}

# #Test
tt <- prepare_data(occ = occ,
                   species = "species",
                   x = "x",
                   y = "y",
                   spat_variables = rast("Models/Myrcia_hatschbachii/PCA_var.tiff"),
                   categorical_variables = "SoilType",
                   nbg = 10000,
                   kfolds = 4,
                   include_xy = TRUE,
                   write_files = T,
                   file_name =  "Models/Myrcia_hatschbachii/Data")

tt2 <- readRDS("Models/Myrcia_hatschbachii/Data.RDS")
tt2$calibration_data %>% str()

any(is.na(tt2$calibration_data))
