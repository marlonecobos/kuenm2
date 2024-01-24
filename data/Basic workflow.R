####Basic workflow to run pre-KUENM2####
#Load packages
library(terra)
library(dplyr)
library(foreach)


#### 1. Prepare SWD ####
#Load function
source("R/prepare_data.R")

#Import data
  #PCA variables
var <- rast("data/PCA_Variables.tiff")
  #Occurrences
occ <- readRDS("data/Occurrences.RDS")


#Crate folder outside package to save results
mydir <- "../test_kuenm2"
dir.create(mydir)

#Prepare swd
sp_swd <- prepare_data(occ = occ,
                       species = occ[1,1],
                       x = "x",
                       y = "y",
                       spat_variables = var,
                       categorical_variables = NULL,
                       nbg = 10000,
                       kfolds = 4,
                       include_xy = TRUE,
                       write_files = T,
                       file_name = "../test_kuenm2/Myrcia",
                       seed = 1,
                       verbose = TRUE)
#Clean environment
rm(list = ls())

#### 2. Create grid of formulas ####
#Load functions
source("R/calibration_grid_glmnetmx.R")
#Load swd data
data <- readRDS("../test_kuenm2/Myrcia.RDS")
#Creta grid
g <- calibration_grid_glmnetmx(swd = data$calibration_data, x = NULL, y = NULL,
      fold_column = "folds",
      min.number = 3,
      categorical_var = NULL,
      features = c("l", "q", "lq", "lqp", "p"),
      min_continuous = NULL,
      regm = c(0.1, 1, 3))
#Save grid of formulas
saveRDS(g, "../test_kuenm2/formulas.RDS")
#Clean environment
rm(list = ls())


#### 3. Fit candidate models ####
#Load functions
source("R/calibration_glmnetmx.R")
source("R/glmnet_mx.R")
source("R/helpers_calibration_glmnetmx.R")
source("R/helpers_glmnetmx.R")
source("R/evaluation_functions.R")
source("R/sel_best_models.R")

#Load data
data <- readRDS("../test_kuenm2/Myrcia.RDS") #SWD data
g <- readRDS("../test_kuenm2/formulas.RDS") #SWD data

#Use 70% of the available cores
ncores <- round(parallel::detectCores()* 0.7, 0)

#Fit models
  #In parallel
m <- calibration_glmnetmx(data = data, #Data in **CLASS??** format
                          #pr_bg, #Column name with presence (1) or background (0)
                          formula_grid = g, #Grid with formulas
                          test_concave = TRUE, #Test concave curves in quadratic models?
                          #folds = 4, #Columns name with k_folds or vector indicating k_folds
                          parallel = TRUE,
                          ncores = 8,
                          progress_bar = TRUE, #Show progress bar? Only works if parallel_type = "doSNOW"
                          write_summary = FALSE, #Write summary of each candidate evaluations?
                          out_dir = NULL, #Name of the folder to write candidate evaluations
                          parallel_type = "doSNOW",
                          return_replicate = TRUE,
                          omrat_thr = 5,
                          omrat_threshold = 5,
                          allow_tolerance = TRUE, #If omission rate is higher than set, select the model with minimum omission rate
                          tolerance = 0.01,
                          AIC = "ws",
                          delta_aic = 2,
                          skip_existing_models = FALSE, #Only works if write_summary = TRUE
                          verbose = TRUE)

#Not in parallel
m2 <- calibration_glmnetmx(data = data, #Data in **CLASS??** format
                          #pr_bg, #Column name with presence (1) or background (0)
                          formula_grid = g, #Grid with formulas
                          test_concave = TRUE, #Test concave curves in quadratic models?
                          #folds = 4, #Columns name with k_folds or vector indicating k_folds
                          parallel = FALSE,
                          ncores = 1,
                          progress_bar = TRUE, #Show progress bar? If parallel = TRUE, only works if parallel_type = "doSNOW"
                          write_summary = FALSE, #Write summary of each candidate evaluations?
                          out_dir = NULL, #Name of the folder to write candidate evaluations
                          parallel_type = "doSNOW",
                          return_replicate = TRUE,
                          omrat_thr = c(5, 10),
                          omrat_threshold = 10,
                          allow_tolerance = TRUE, #If omission rate is higher than set, select the model with minimum omission rate
                          tolerance = 0.01,
                          AIC = "ws",
                          delta_aic = 2,
                          skip_existing_models = FALSE, #Only works if write_summary = TRUE
                          verbose = TRUE)


#Save candidate models
saveRDS(m2, "../test_kuenm2/Candidate_and_Best_models.RDS")


#### Fit best models ####
source("R/fit_selected_glmnetmx.R")
source("R/helpers_calibration_glmnetmx.R")
source("R/helpers_glmnetmx.R")
source("R/calibration_glmnetmx.R")
source("R/evaluation_functions.R")
source("R/glmnet_mx.R")
source("R/part_data.R")

#Load selectes models
m <- readRDS("../test_kuenm2/Candidate_and_Best_models.RDS")

#In parallel
fm <- fit_selected_glmnetmx(calibration_results = m,
                            n_replicates = 10,
                            rep_type = "kfold",
                            train_portion = 0.7,
                            write_models = TRUE, #Write files?
                            file_name = "../test_kuenm2/Best_Models", #Name of the folder to write final models
                            parallel = TRUE,
                            ncores = 4,
                            parallelType = "doSNOW",
                            progress_bar = TRUE,
                            verbose = TRUE)

#Not parallel
fm2 <- fit_selected_glmnetmx(calibration_results = m,
                            n_replicates = 10,
                            rep_type = "kfold",
                            train_portion = 0.7,
                            write_models = FALSE, #Write files?
                            file_name = NULL, #Name of the folder to write final models
                            parallel = F,
                            ncores = 4,
                            parallelType = "doSNOW",
                            progress_bar = TRUE,
                            verbose = TRUE)




####Predict best models####
#Load functions
source("R/predict_selected_glmnetmx.R")
source("R/helpers_glmnetmx.R")
source("R/helpers_calibration_glmnetmx.R")

#Load fitted best models
bm <- readRDS("../test_kuenm2/Best_Models.RDS")
#Load spatrasters
r <- rast("data/PCA_Variables.tiff")

#Predict without clamp
p <- predict_selected_glmnetmx(models = bm,
                               spat_var = r,
                               write_files = FALSE,
                               write_replicates = FALSE,
                               out_dir = NULL,
                               consensus_per_model = TRUE,
                               consensus_general = TRUE,
                               consensus = c("median", "range", "mean", "stdev"), #weighted mean
                               clamping = FALSE,
                               var_to_clamp = NULL,
                               type = "cloglog",
                               overwrite = FALSE)
#Predict with clamp (all variables)
p_clamp <- predict_selected_glmnetmx(models = bm,
                               spat_var = r,
                               write_files = FALSE,
                               write_replicates = FALSE,
                               out_dir = NULL,
                               consensus_per_model = TRUE,
                               consensus_general = TRUE,
                               consensus = c("median", "range", "mean", "stdev"), #weighted mean
                               clamping = TRUE,
                               var_to_clamp = NULL,
                               type = "cloglog",
                               overwrite = FALSE)
#Predict with clamp with only one variable
p_clamp_pc1 <- predict_selected_glmnetmx(models = bm,
                                     spat_var = r,
                                     write_files = FALSE,
                                     write_replicates = FALSE,
                                     out_dir = NULL,
                                     consensus_per_model = TRUE,
                                     consensus_general = TRUE,
                                     consensus = c("median", "range", "mean", "stdev"), #weighted mean
                                     clamping = TRUE,
                                     var_to_clamp = "PC1",
                                     type = "cloglog",
                                     overwrite = FALSE)
#Compare
plot(p$General_consensus$median)
plot(p_clamp$General_consensus$median)
plot(p$General_consensus$median - p_clamp$General_consensus$median)
plot(p$General_consensus$median - p_clamp_pc1$General_consensus$median)
plot(p_clamp$General_consensus$median - p_clamp_pc1$General_consensus$median)


# ####Compare time between lapply, pblapply and for with progress bar####
# #install.packages("rbenchmark")
# library(rbenchmark)
# library(pbapply)
#
# #Load functions
# source("R/calibration_glmnetmx.R")
# source("R/glmnet_mx.R")
# source("R/helpers_calibration_glmnetmx.R")
# source("R/helpers_glmnetmx.R")
# source("R/evaluation_functions.R")
# source("R/sel_best_models.R")
# source("R/helpers_lapply_with_progress_bar.R")
#
# #Load data
# data <- readRDS("../test_kuenm2/Myrcia.RDS") #SWD data
# g <- readRDS("../test_kuenm2/formulas.RDS") #SWD data
#
# #Set progress bar
# pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
#                      max = nrow(g), # Maximum value of the progress bar
#                      style = 3,    # Progress bar style (also available style = 1 and style = 2)
#                      width = 80,   # Progress bar width. Defaults to getOption("width")
#                      char = "=")
#
# #Benchmarck
# teste <- benchmark(
#   "pblapply" = {pblapply(X = 1:nrow(g), function(x)
#   fit_eval_models(x, formula_grid2 = g, data = data, formula_grid = g,
#                   omrat_thr = 10,
#                   write_summary = FALSE, return_replicate = TRUE))},
#                   "classic_for" = {
#                     results2 <- vector("list", length = nrow(g))
#
#                     # Loop for com barra de progresso manual
#                     for (x in 1:nrow(g)) {
#                       # Execute a função fit_eval_models
#                       results2[[x]] <- fit_eval_models(x, formula_grid2 = g, data = data, formula_grid = g,
#                                                        omrat_thr = 10,
#                                                        write_summary = FALSE, return_replicate = TRUE)
#
#                       # Sets the progress bar to the current state
#                       setTxtProgressBar(pb, x)
#                     }
#                   },
#   replications = 10,
#   columns = c("test", "replications", "elapsed",
#               "relative", "user.self", "sys.self"))
# flextable(teste)
# write.csv(teste, "../test_kuenm2/Time_test2.csv")
# #Read time test
# teste <- read.csv("../test_kuenm2/Time_test2.csv")
#


