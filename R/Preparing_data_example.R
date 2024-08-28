# #### Arrange environmental data ####
# library(terra)
# library(geodata)
# library(dplyr)
#
# #Worldclim data
# worldclim_global(var = "bio", res = 10,
#                  path = "../test_kuenm2/WC_Current_raw_10")
# var <- rast(list.files("../test_kuenm2/WC_Current_raw_10/wc2.1_10m/", full.names = TRUE))
# names(var)
# #Rename
# names(var) <- gsub("wc2.1_10m_", "", names(var))
#
# #Select variables
# var <- var[[c("bio_1", "bio_7", "bio_12", "bio_15")]]
# #Cut variables
# #Occurrences
# load("data/occ_data.rda")
# #Add convex ploygon and buffer around records
# pts <- vect(occ_data, geom = c(x = "x", y = "y"), crs = crs(var))
# m <- terra::convHull(pts) %>% buffer(., width = 200*1000)
# var_m <- crop(var, m, mask =TRUE)
# var_m %>% values() %>% na.omit() %>% length()
# plot(var_m$bio_1)
# plot(pts, add = TRUE)
# names(var_m)
#
# #Append soil type
# soil <- rast("../Models_to_run/Current_Neotropic/Variables.tiff")
# soil <- soil$SoilType
# #Mask
# soil <- crop(soil, m, mask = TRUE)
# plot(soil)
# #Aggregate
# f_to_aggreg <- res(var_m)/res(soil)
# #Resample
# soil10 <- terra::resample(soil, var_m$bio_1, method = "near")
# plot(soil10)
# names(soil10)
# var_m <- c(var_m, soil10)
#
# #Save as new data
# writeRaster(var_m, "inst/extdata/Current_variables.tif")
#
# #Projections
# #Get data
# dir.create("../test_kuenm2/WC-Future_Raw_10")
# #GCMS
# gcms <- c("ACCESS-CM2", "MIROC6")
# ssps <- c("ssp126", "ssp585")
# periods <- c("2041-2060", "2081-2100")
# #Create grid
# g <- expand.grid(gcms = gcms, ssps = ssps, periods = periods)
# #Create url
# g$url <- paste0("https://geodata.ucdavis.edu/cmip6/10m/",
#                 g$gcms, "/", g$ssps, "/wc2.1_10m_bioc_", g$gcms, "_", g$ssps, "_", g$periods, ".tif")
# #Create dist file
# g$path <- paste0("../test_kuenm2/WC-Future_Raw_10/",
#                  "wc2.1_10m_bioc_", g$gcms, "_", g$ssps, "_", g$periods, ".tif")
#
# #Download
# lapply(1:nrow(g), function(i){
#   utils::download.file(url = g$url[i],
#                        destfile = g$path[i],
#                        mode = "wb")
# })
#
#
# #Import, fix and save variables
# f <- list.files("../test_kuenm2/WC-Future_Raw_10/", full.names = TRUE)
# f
#
#
# #Read rasters, cut and save variables, using same name
# lapply(f, function(i){
#   var_x <- rast(i)
#   #Rename variables
#   names(var_x) <- gsub(".*_(.*?)$", "bio_\\1", names(var_x))
#   names(var_x) <- gsub("bio0", "bio", names(var_x))
#   names(var_x) <- gsub("bio", "bio_", names(var_x))
#   #Select and cut variables
#   var_x <- crop(var_x[[names(var_m)[-5]]], m, mask = TRUE)
#   #Get name to save
#   s <- fs::path_ext_remove(basename(i))
#   #Save
#   writeRaster(var_x,
#               paste0("inst/extdata/", s, ".tif"), overwrite = TRUE)
# })
#
# # #Past
# # p <- list.dirs("../test_kuenm2/WC_Past_raw/", recursive = FALSE)
# # #get lgm and mid of two gcms
# # p <- p[c(1,2,5,6)]
# #
# # #Read, rename and save variables
# # lapply(p, function(x){
# #   var_x <- list.files(x, full.names = TRUE)
# #   var_x <- rast(lapply(var_x, rast))
# #   #Rename variables
# #   names(var_x) <- as.numeric(gsub("\\D", "", names(var_x)))
# #   names(var_x) <- paste0("bio_", names(var_x))
# #   #Select and mask variables
# #   var_x <- crop(var_x[[names(v)[-5]]], m, mask = TRUE)
# #   #Fix variables
# #   var_x$bio_1 <- var_x$bio_1/10
# #   var_x$bio_7 <- var_x$bio_7/10
# #   #Get name to save
# #   s <- fs::path_ext_remove(basename(x))
# #   #Save
# #   writeRaster(var_x,
# #               paste0("inst/extdata/", s, ".tif"), overwrite = TRUE)
# # })
# #
