#' Calibration Results (glm)
#'
#' @description
#' A `calibration_results` object resulted from \code{calibration()} using maxnet algorithm
#'
#' @usage data("calib_results_glm")
#'
#' @format A \code{calibration_results} with the following elements:
#' \describe{
#'   \item{species}{Species names}
#'   \item{calibration_data}{A \code{data.frame} with the variables extracted to presence and background points}
#'   \item{formula_grid}{A \code{data.frame} with the ID, formulas, and regularization multipliers of each candidate model}
#'   \item{part_data}{A \code{list} with the partition data, where each element corresponds to a replicate and contains the **indices of the test points** for that replicate}
#'   \item{partition_method}{A \code{character} indicating the partition method}
#'   \item{n_replicates}{A \code{numeric} value indicating the number of replicates or k-folds}
#'   \item{train_proportion}{A \code{numeric} value indicating the proportion of occurrences used as train points when the partition method is 'subsample' or 'boostrap'}
#'   \item{data_xy}{A \code{data.frame} with the coordinates of the occurrence and bakground points}
#'   \item{continuous_variables}{A \code{character} indicating the names of the continuous variables}
#'   \item{categorical_variables}{A \code{character} indicating the names of the categorical variables}
#'   \item{weights}{A \code{numeric} value specifying weights for the occurrence records. It's NULL, meaning it was not set weights.}
#'   \item{pca}{A \code{prcomp} object storing PCA information. Is NULL, meaning PCA was not performed}
#'   \item{algorithm}{A \code{character} indicanting the algorithm (glm)}
#'   \item{calibration_results}{A \code{list} containing the evaluation metrics for each candidate model}
#'   \item{omission_rate}{A \code{numeric} value indicating the omission rate used to evaluate the models (10%)}
#'   \item{addsamplestobackground}{A \code{logical} value indicating whether to add to the background any presence sample that is not already there.}
#'   \item{selected_models}{A \code{data.frame} with the formulas and evaluation metrics for each selected model}
#'   \item{summary}{A \code{list} with the number and the ID of the models removed and selected during selection procedure}
#' }
"calib_results_glm"


#' Calibration Results (Maxnet)
#'
#' @description
#' A `calibration_results` object resulted from \code{calibration()} using maxnet algorithm
#'
#' @usage data("calib_results_maxnet")
#'
#' @format A \code{calibration_results} with the following elements:
#' \describe{
#'   \item{species}{Species names}
#'   \item{calibration_data}{A \code{data.frame} with the variables extracted to presence and background points}
#'   \item{formula_grid}{A \code{data.frame} with the ID, formulas, and regularization multipliers of each candidate model}
#'   \item{part_data}{A \code{list} with the partition data, where each element corresponds to a replicate and contains the **indices of the test points** for that replicate}
#'   \item{partition_method}{A \code{character} indicating the partition method}
#'   \item{n_replicates}{A \code{numeric} value indicating the number of replicates or k-folds}
#'   \item{train_proportion}{A \code{numeric} value indicating the proportion of occurrences used as train points when the partition method is 'subsample' or 'boostrap'}
#'   \item{data_xy}{A \code{data.frame} with the coordinates of the occurrence and bakground points}
#'   \item{continuous_variables}{A \code{character} indicating the names of the continuous variables}
#'   \item{categorical_variables}{A \code{character} indicating the names of the categorical variables}
#'   \item{weights}{A \code{numeric} value specifying weights for the occurrence records. It's NULL, meaning it was not set weights.}
#'   \item{pca}{A \code{prcomp} object storing PCA information. Is NULL, meaning PCA was not performed}
#'   \item{algorithm}{A \code{character} indicanting the algorithm (maxnet)}
#'   \item{calibration_results}{A \code{list} containing the evaluation metrics for each candidate model}
#'   \item{omission_rate}{A \code{numeric} value indicating the omission rate used to evaluate the models (10%)}
#'   \item{addsamplestobackground}{A \code{logical} value indicating whether to add to the background any presence sample that is not already there.}
#'   \item{selected_models}{A \code{data.frame} with the formulas and evaluation metrics for each selected model}
#'   \item{summary}{A \code{list} with the number and the ID of the models removed and selected during selection procedure}
#' }
"calib_results_maxnet"


#' Spatial Blocks from ENMeval
#'
#' @description
#' A list resulting from \code{ENMeval::get.block()} to partition occurrence and background localities into bins for training and validation (or, evaluation and calibration). This object is used in the "Prepare Data for Model Calibration" vignette to demonstrate how to implement custom data partitions generated by \code{ENMeval} in \code{kuenm2}.
#'
#' @usage data("enmeval_block")
#'
#' @format A \code{list} with the following elements:
#' \describe{
#'    \item{occs.grp}{A \code{numeric} vector indicating the spatial group to which each occurrence belongs}
#'    \item{bg.grp}{A \code{numeric} vector indicating the spatial group to which each background point belongs}
#' }
"enmeval_block"


#' Fitted model with CHELSA variables
#'
#' @description
#' A `fitted_models` object resulting from \code{fit_selected()} using calibration data based on CHELSA variables.
#'
#' @usage data("fitted_model_chelsa")
#'
#' @format A \code{fitted_models} with the following elements:
#' \describe{
#'   \item{species}{Species names}
#'   \item{Models}{A \code{list} with the fitted maxnet models (replicates and full models)}
#'   \item{calibration_data}{A \code{data.frame} containing the variables extracted for presence and background points}
#'   \item{continuous_variables}{A \code{character} indicating the names of the continuous variables}
#'   \item{categorical_variables}{A \code{character} indicating the names of the categorical variables}
#'   \item{selected_models}{A \code{data.frame} with formulas and evaluation metrics for each selected model}
#'   \item{weights}{A \code{numeric} vector specifying weights for the occurrence records. \code{NULL} if no weights were set.}
#'   \item{pca}{A \code{prcomp} object containing PCA results. \code{NULL} if PCA was not performed.}
#'   \item{addsamplestobackground}{A \code{logical} value indicating whether to add any presence point not already included to the background.}
#'   \item{omission_rate}{A \code{numeric} value indicating the omission rate used to evaluate models.}
#'   \item{thresholds}{A \code{numeric} vector with thresholds used to binarize each replicate and the consensus (mean and median), calculated based on the omission rate defined in \code{calibration()}.}
#'   \item{algorithm}{A \code{character} string indicating the algorithm used (maxnet).}
#'   \item{partition_method}{A \code{character} string indicating the partitioning method used.}
#'   \item{n_replicates}{A \code{numeric} value indicating the number of replicates or folds.}
#'   \item{train_proportion}{A \code{numeric} value indicating the proportion of occurrences used for training when the partition method is 'subsample' or 'bootstrap'.}
#'   }
"fitted_model_chelsa"


#' Fitted model with glm algorithm
#'
#' @description
#' A glm `fitted_models` object resulting from \code{fit_selected()} using calibration data with based on WorldCLim variables.
#'
#' @usage data("fitted_model_glm")
#'
#' @format A \code{fitted_models} with the following elements:
#' \describe{
#'   \item{species}{Species names}
#'   \item{Models}{A \code{list} with the fitted maxnet models (replicates and full models)}
#'   \item{calibration_data}{A \code{data.frame} containing the variables extracted for presence and background points}
#'   \item{continuous_variables}{A \code{character} indicating the names of the continuous variables}
#'   \item{categorical_variables}{A \code{character} indicating the names of the categorical variables}
#'   \item{selected_models}{A \code{data.frame} with formulas and evaluation metrics for each selected model}
#'   \item{weights}{A \code{numeric} vector specifying weights for the occurrence records. \code{NULL} if no weights were set.}
#'   \item{pca}{A \code{prcomp} object containing PCA results. \code{NULL} if PCA was not performed.}
#'   \item{addsamplestobackground}{A \code{logical} value indicating whether to add any presence point not already included to the background.}
#'   \item{omission_rate}{A \code{numeric} value indicating the omission rate used to evaluate models.}
#'   \item{thresholds}{A \code{numeric} vector with thresholds used to binarize each replicate and the consensus (mean and median), calculated based on the omission rate defined in \code{calibration()}.}
#'   \item{algorithm}{A \code{character} string indicating the algorithm used (glm).}
#'   \item{partition_method}{A \code{character} string indicating the partitioning method used.}
#'   \item{n_replicates}{A \code{numeric} value indicating the number of replicates or folds.}
#'   \item{train_proportion}{A \code{numeric} value indicating the proportion of occurrences used for training when the partition method is 'subsample' or 'bootstrap'.}
#'   }
"fitted_model_glm"


#' Fitted model with maxnet algorithm
#'
#' @description
#' A maxnet `fitted_models` object resulting from \code{fit_selected()} using calibration data with based on WorldCLim variables.
#'
#' @usage data("fitted_model_maxnet")
#'
#' @format A \code{fitted_models} with the following elements:
#' \describe{
#'   \item{species}{Species names}
#'   \item{Models}{A \code{list} with the fitted maxnet models (replicates and full models)}
#'   \item{calibration_data}{A \code{data.frame} containing the variables extracted for presence and background points}
#'   \item{continuous_variables}{A \code{character} indicating the names of the continuous variables}
#'   \item{categorical_variables}{A \code{character} indicating the names of the categorical variables}
#'   \item{selected_models}{A \code{data.frame} with formulas and evaluation metrics for each selected model}
#'   \item{weights}{A \code{numeric} vector specifying weights for the occurrence records. \code{NULL} if no weights were set.}
#'   \item{pca}{A \code{prcomp} object containing PCA results. \code{NULL} if PCA was not performed.}
#'   \item{addsamplestobackground}{A \code{logical} value indicating whether to add any presence point not already included to the background.}
#'   \item{omission_rate}{A \code{numeric} value indicating the omission rate used to evaluate models.}
#'   \item{thresholds}{A \code{numeric} vector with thresholds used to binarize each replicate and the consensus (mean and median), calculated based on the omission rate defined in \code{calibration()}.}
#'   \item{algorithm}{A \code{character} string indicating the algorithm used (maxnet).}
#'   \item{partition_method}{A \code{character} string indicating the partitioning method used.}
#'   \item{n_replicates}{A \code{numeric} value indicating the number of replicates or folds.}
#'   \item{train_proportion}{A \code{numeric} value indicating the proportion of occurrences used for training when the partition method is 'subsample' or 'bootstrap'.}
#'   }
"fitted_model_maxnet"


#' Spatial Blocks from flexsdm
#'
#' @description
#' A list resulting from \code{flexsdm::part_sblock()}, used to partition occurrence and background localities into bins for training and evaluation. This object is used in the "Prepare Data for Model Calibration" vignette to demonstrate how to implement custom data partitions generated by \code{flexsdm} in \code{kuenm2}
#'
#' @usage data("enmeval_block")
#'
#' @format A \code{list} with the following elements:
#' \describe{
#'    \item{part}{A \code{tibble} object with information used in 'data' arguments and a additional column .part with partition group.}
#'    \item{best_part_info}{A \code{tibble} with information about the best partition.}
#' }
"flexsdm_block"

#' Discrete palettes based on pals R package
#'
#' @description
#' Color palettes designed for discrete, categorical data. Palettes retrived from pals R package
#'
#' @usage data("kuenm2_discrete_palletes")
#' @format A \code{list} with the following color palettes: "alphabet",
#' "alphabet2", "cols25", "glasbey", "kelly", "polychrome", "stepped",
#' "stepped2", "stepped3", "okabe", "tableau20", "tol", "tol.groundcover",
#' "trubetskoy", and "watlington"
#'
#' @references
#' Wright K (2023). pals: Color Palettes, Colormaps, and Tools to Evaluate
#' Them_. R package version 1.8, <https://CRAN.R-project.org/package=pals>.
"kuenm2_discrete_palletes"


#' Species Occurrence with Erroneous Records
#'
#' @description
#' A \code{data.frame} containing the coordinates of 51 valid occurrences of *Myrcia hatschbachii* (a tree endemic to southern Brazil), along with a set of erroneous records used to demonstrate data cleaning procedures. The valid occurrences were sourced from Trindade & Marques (2024).
#'
#' @usage data("occ_data_noclean")
#'
#' @format A \code{data.frame} with the following columns:
#' \describe{
#'   \item{species}{The species name.}
#'   \item{x}{Longitude.}
#'   \item{y}{Latitude.}
#' }
#' @references
#' Trindade, W.C.F., Marques, M.C.M., 2023. The Invisible Species: Big Data Unveil Coverage Gaps in the Atlantic Forest Hotspot. Diversity and Distributions 30, e13931. https://doi.org/10.1111/ddi.13931
"occ_data_noclean"


#' Species Occurrence
#'
#' @description
#' A \code{data.frame} containing the coordinates of 51 valid occurrences of *Myrcia hatschbachii* (a tree endemic to southern Brazil). The valid occurrences were sourced from Trindade & Marques (2024) and contains only the records retrieved
#' from GBIF and SpeciesLink.
#'
#' @usage data("occ_data")
#'
#' @format A \code{data.frame} with the following columns:
#' \describe{
#'   \item{species}{The species name.}
#'   \item{x}{Longitude.}
#'   \item{y}{Latitude.}
#' }
#' @references
#' Trindade, W.C.F., Marques, M.C.M., 2023. The Invisible Species: Big Data Unveil Coverage Gaps in the Atlantic Forest Hotspot. Diversity and Distributions 30, e13931. https://doi.org/10.1111/ddi.13931
"occ_data"

#' Prepared Data for glm models
#'
#' @description
#' A `prepared_data` object resulted from \code{prepare_data()} to calibrate models using 'glm' algorithm.
#'
#' @usage data("sp_swd_glm")
#'
#' @format A \code{prepared_data} object with the following elements:
#' \describe{
#'   \item{species}{Species names}
#'   \item{calibration_data}{A \code{data.frame} containing the variables extracted for presence and background points}
#'   \item{formula_grid}{A \code{data.frame} with the ID, formulas, and regularization multipliers of each candidate model}
#'   \item{part_data}{A \code{list} with the partition data, where each element corresponds to a replicate and contains the **indices of the test points** for that replicate}
#'   \item{partition_method}{A \code{character} indicating the partition method}
#'   \item{n_replicates}{A \code{numeric} value indicating the number of replicates or k-folds}
#'   \item{train_proportion}{A \code{numeric} value indicating the proportion of occurrences used as train points when the partition method is 'subsample' or 'boostrap'}
#'   \item{data_xy}{A \code{data.frame} with the coordinates of the occurrence and bakground points}
#'   \item{continuous_variables}{A \code{character} indicating the names of the continuous variables}
#'   \item{categorical_variables}{A \code{character} indicating the names of the categorical variables}
#'   \item{weights}{A \code{numeric} value specifying weights for the occurrence records. It's NULL, meaning it was not set weights.}
#'   \item{pca}{A \code{prcomp} object storing PCA information. Is NULL, meaning PCA was not performed}
#'   \item{algorithm}{A \code{character} indicanting the algorithm (glm)}
#'   }
"sp_swd_glm"


#' Prepared Data for maxnet models
#'
#' @description
#' A `prepared_data` object resulted from \code{prepare_data()} to calibrate models using 'glm' algorithm.
#'
#' @usage data("sp_swd")
#'
#' @format A \code{prepared_data} object with the following elements:
#' \describe{
#'   \item{species}{Species names}
#'   \item{calibration_data}{A \code{data.frame} containing the variables extracted for presence and background points}
#'   \item{formula_grid}{A \code{data.frame} with the ID, formulas, and regularization multipliers of each candidate model}
#'   \item{part_data}{A \code{list} with the partition data, where each element corresponds to a replicate and contains the **indices of the test points** for that replicate}
#'   \item{partition_method}{A \code{character} indicating the partition method}
#'   \item{n_replicates}{A \code{numeric} value indicating the number of replicates or k-folds}
#'   \item{train_proportion}{A \code{numeric} value indicating the proportion of occurrences used as train points when the partition method is 'subsample' or 'boostrap'}
#'   \item{data_xy}{A \code{data.frame} with the coordinates of the occurrence and bakground points}
#'   \item{continuous_variables}{A \code{character} indicating the names of the continuous variables}
#'   \item{categorical_variables}{A \code{character} indicating the names of the categorical variables}
#'   \item{weights}{A \code{numeric} value specifying weights for the occurrence records. It's NULL, meaning it was not set weights.}
#'   \item{pca}{A \code{prcomp} object storing PCA information. Is NULL, meaning PCA was not performed}
#'   \item{algorithm}{A \code{character} indicanting the algorithm (glm)}
#'   }
"sp_swd"

#' Prepared data with spatial blocks created with ENMeval
#'
#' @description
#' A `prepared_data` object resulted from \code{prepare_data()} to calibrate
#' models using 'glmnet' algorithm. In this object, the original partitioning
#' was replaced with spatial blocks generated using the \code{get.block()}
#' method from the ENMeval R package.
#'
#' @usage data("sp_swd")
"swd_spatial_block"

#' User Custom Calibration Data
#'
#' @description
#' A \code{data.frame} containing presence and background records along with environmental variables used to demonstrate data preparation with user-supplied data.
#'
#' @usage data("user_data")
#'
#' @format A \code{data.frame} with the following columns:
#' \describe{
#'   \item{pr_bg}{Column indicating presences (1) and background (0).}
#'   \item{bio_1}{The extracted values for the variable bio_1 at presence and background points.}
#'   \item{bio_7}{The extracted values for the variable bio_12 at presence and background points.}
#'   \item{bio_12}{The extracted values for the variable bio_12 at presence and background points.}
#'   \item{bio_15}{The extracted values for the variable bio_15 at presence and background points.}
#'   \item{bio_15}{The extracted values for the variable soilType at presence and background points.}
#' }
"user_data"

#' Fitted model with concave curves
#'
#' @description
#' A maxnet `fitted_models` object resulting from \code{fit_selected()} with a model presenting concave curves.
#'
#' @usage data("fitted_model_concave")
#'
#' @format A \code{fitted_models} with the following elements:
#' \describe{
#'   \item{species}{Species names}
#'   \item{Models}{A \code{list} with the fitted maxnet models (replicates and full models)}
#'   \item{calibration_data}{A \code{data.frame} containing the variables extracted for presence and background points}
#'   \item{continuous_variables}{A \code{character} indicating the names of the continuous variables}
#'   \item{categorical_variables}{A \code{character} indicating the names of the categorical variables}
#'   \item{selected_models}{A \code{data.frame} with formulas and evaluation metrics for each selected model}
#'   \item{weights}{A \code{numeric} vector specifying weights for the occurrence records. \code{NULL} if no weights were set.}
#'   \item{pca}{A \code{prcomp} object containing PCA results. \code{NULL} if PCA was not performed.}
#'   \item{addsamplestobackground}{A \code{logical} value indicating whether to add any presence point not already included to the background.}
#'   \item{omission_rate}{A \code{numeric} value indicating the omission rate used to evaluate models.}
#'   \item{thresholds}{A \code{numeric} vector with thresholds used to binarize each replicate and the consensus (mean and median), calculated based on the omission rate defined in \code{calibration()}.}
#'   \item{algorithm}{A \code{character} string indicating the algorithm used (maxnet).}
#'   \item{partition_method}{A \code{character} string indicating the partitioning method used.}
#'   \item{n_replicates}{A \code{numeric} value indicating the number of replicates or folds.}
#'   }
"fitted_model_concave"

#' Independent Species Occurrence
#'
#' @description
#' A \code{data.frame} containing the coordinates of 82 occurrences of *Myrcia hatschbachii* (a tree endemic to southern Brazil). The valid occurrences were
#' sourced from NeotropicTree (Oliveira-Filho, 2017) and were used as
#' independent data to test the models fitted with the `occ_data`.
#'
#' @usage data("new_occ")
#'
#' @format A \code{data.frame} with the following columns:
#' \describe{
#'   \item{species}{The species name.}
#'   \item{x}{Longitude.}
#'   \item{y}{Latitude.}
#' }
#' @references
#' Oliveira_Filho, A.T. 2017. NeoTropTree, Flora arbórea da Região Neotropical: Um banco de dados envolvendo biogeografia, diversidade e conservação. Universidade Federal de Minas Gerais. (http://www.neotroptree,info).
"new_occ"

#### Document ext data ####

#' Example Bias File
#'
#' A `SpatRaster` object representing a bias layer used for extracting
#' background points with the `prepare_data()` function.
#'
#' @format A `SpatRaster` object.
#'
#' @name bias
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' bias <- terra::rast(system.file("extdata", "bias_file.tif",
#'                                 package = "kuenm2"))
#'
#' terra::plot(bias)
NULL

#' SpatRaster Representing LGM Conditions (GCM: CCSM4)
#'
#' Raster layer containing bioclimatic variables representing Last Glacial
#' Maximum (LGM) climatic conditions based on the CCSM4 General Circulation
#' Model (GCM). The variables were resampled to 10arc-minutes and masked using
#' the `m` provided in the package. Data sourced from CHELSA:
#' \url{https://chelsa-climate.org/last-glacial-maximum-climate/}
#'
#' @format A `SpatRaster` object.
#'
#' @name chelsa_lgm_ccsm4
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' chelsa_lgm_ccsm4 <- terra::rast(system.file("extdata",
#'                                            "CHELSA_LGM_CCSM4.tif",
#'                                             package = "kuenm2"))
#'
#' terra::plot(chelsa_lgm_ccsm4)
NULL

#' SpatRaster Representing LGM Conditions (GCM: CNRM-CM5)
#'
#' Raster layer containing bioclimatic variables representing Last Glacial
#' Maximum (LGM) climatic conditions based on the CNRM-CM5 General Circulation
#' Model. The variables were resampled to 10arc-minutes and masked using the `m`
#' provided in the package. Data sourced from CHELSA:
#' \url{https://chelsa-climate.org/last-glacial-maximum-climate/}
#'
#' @format A `SpatRaster` object.
#'
#' @name chelsa_lgm_cnrm_cm5
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' chelsa_lgm_cnrm_cm5 <- terra::rast(system.file("extdata",
#'                                                "CHELSA_LGM_CNRM-CM5.tif",
#'                                                package = "kuenm2"))
#' terra::plot(chelsa_lgm_cnrm_cm5)
NULL

#' SpatRaster Representing LGM Conditions (GCM: FGOALS-g2)
#'
#' Raster layer containing bioclimatic variables representing Last Glacial
#' Maximum (LGM) climatic conditions based on the FGOALS-g2 General Circulation
#' Model. The variables were resampled to 10arc-minutes and masked using the `m`
#' provided in the package. Data sourced from CHELSA:
#' \url{https://chelsa-climate.org/last-glacial-maximum-climate/}
#'
#' @format A `SpatRaster` object.
#'
#' @name chelsa_lgm_fgoals_g2
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' chelsa_lgm_fgoals_g2 <- terra::rast(system.file("extdata",
#'                                                "CHELSA_LGM_FGOALS-g2.tif",
#'                                                package = "kuenm2"))
#' terra::plot(chelsa_lgm_fgoals_g2)
NULL

#' SpatRaster Representing LGM Conditions (GCM: IPSL-CM5A-LR)
#'
#' Raster layer containing bioclimatic variables representing Last Glacial
#' Maximum (LGM) climatic conditions based on the IPSL-CM5A-LR General
#' Circulation Model. The variables were resampled to 10arc-minutes and masked
#' using the `m` provided in the package. Data sourced from CHELSA:
#' \url{https://chelsa-climate.org/last-glacial-maximum-climate/}
#'
#' @format A `SpatRaster` object.
#'
#' @name chelsa_lgm_ipsl
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' chelsa_lgm_ipsl <- terra::rast(system.file("extdata",
#'                                            "CHELSA_LGM_IPSL-CM5A-LR.tif",
#'                                            package = "kuenm2"))
#' terra::plot(chelsa_lgm_ipsl)
NULL

#' SpatRaster Representing LGM Conditions (GCM: MIROC-ESM)
#'
#' Raster layer containing bioclimatic variables representing Last Glacial
#' Maximum (LGM) climatic conditions based on the MIROC-ESM General Circulation
#' Model. The variables were resampled to 10arc-minutes and masked using the `m`
#' provided in the package. Data sourced from CHELSA:
#' \url{https://chelsa-climate.org/last-glacial-maximum-climate/}
#'
#' @format A `SpatRaster` object.
#'
#' @name chelsa_lgm_miroc
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' chelsa_lgm_miroc <- terra::rast(system.file("extdata",
#'                                            "CHELSA_LGM_MIROC-ESM.tif",
#'                                            package = "kuenm2"))
#' terra::plot(chelsa_lgm_miroc)
NULL

#' SpatRaster Representing LGM Conditions (GCM: MPI-ESM-P)
#'
#' Raster layer containing bioclimatic variables representing Last Glacial
#' Maximum (LGM) climatic conditions based on the MPI-ESM-P General Circulation
#' Model. The variables were resampled to 10arc-minutes and masked using the `m`
#' provided in the package. Data sourced from CHELSA:
#' \url{https://chelsa-climate.org/last-glacial-maximum-climate/}
#'
#' @format A `SpatRaster` object.
#'
#' @name chelsa_lgm_mpi
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' chelsa_lgm_mpi <- terra::rast(system.file("extdata",
#'                                            "CHELSA_LGM_MPI-ESM-P.tif",
#'                                            package = "kuenm2"))
#' terra::plot(chelsa_lgm_mpi)
NULL

#' SpatRaster Representing LGM Conditions (GCM: MRI-CGCM3)
#'
#' Raster layer containing bioclimatic variables representing Last Glacial
#' Maximum (LGM) climatic conditions based on the MRI-CGCM3 General Circulation
#' Model. The variables were resampled to 10arc-minutes and masked using the `m`
#' provided in the package. Data sourced from CHELSA:
#' \url{https://chelsa-climate.org/last-glacial-maximum-climate/}
#'
#' @format A `SpatRaster` object.
#'
#' @name chelsa_lgm_mri
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' chelsa_lgm_mri <- terra::rast(system.file("extdata",
#'                                            "CHELSA_LGM_MRI-CGCM3.tif",
#'                                            package = "kuenm2"))
#' terra::plot(chelsa_lgm_mri)
NULL

#' SpatRaster Representing present-day Conditions (CHELSA)
#'
#' Raster layer containing bioclimatic variables representing present-day
#' climatic conditions. The variables were resampled to a 10 arc-minute
#' resolution and masked using the `m` region provided in the package. Data
#' sourced from CHELSA:
#' \url{https://chelsa-climate.org/}
#'
#' @format A `SpatRaster` object.
#'
#' @name chelsa_current
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' chelsa_current <- terra::rast(system.file("extdata",
#'                                            "Current_CHELSA.tif",
#'                                            package = "kuenm2"))
#' terra::plot(chelsa_current)
NULL

#' SpatRaster Representing present-day Conditions (WorldClim)
#'
#' Raster layer containing bioclimatic variables representing present-day
#' climatic conditions. The variables were obtained at a 10 arc-minute
#' resolution and masked using the `m` region provided in the package. Data
#' sourced from WorldClim:
#' \url{https://worldclim.org/data/worldclim21.html}
#'
#' @format A `SpatRaster` object.
#'
#' @name var
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' var <- terra::rast(system.file("extdata",
#'                                "Current_variables.tif",
#'                                 package = "kuenm2"))
#' terra::plot(var)
NULL

#' SpatRaster Representing Future Conditions (2041-2060, SSP126, GCM: ACCESS-CM2)
#'
#' A raster layer containing bioclimatic variables representing future climatic
#' conditions (2041-2060) based on the ACCESS-CM2 General Circulation Model
#' under the SSP126 scenario. The variables were obtained at a 10 arc-minute
#' resolution and masked using the `m` region provided in the package. Data
#' sourced from WorldClim: \url{https://worldclim.org/data/cmip6/cmip6climate.html}
#'
#' @format A `SpatRaster` object.
#'
#' @name future_2050_ssp126_access
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' future_2050_ssp126_access <- terra::rast(system.file("extdata",
#'                                     "wc2.1_10m_bioc_ACCESS-CM2_ssp126_2041-2060.tif",
#'                                      package = "kuenm2"))
#' terra::plot(future_2050_ssp126_access)
NULL

#' SpatRaster Representing Future Conditions (2081-2100, SSP126, GCM: ACCESS-CM2)
#'
#' A raster layer containing bioclimatic variables representing future climatic
#' conditions (2081-2100) based on the ACCESS-CM2 General Circulation Model
#' under the SSP126 scenario. The variables were obtained at a 10 arc-minute
#' resolution and masked using the `m` region provided in the package. Data
#' sourced from WorldClim: \url{https://worldclim.org/data/cmip6/cmip6climate.html}
#'
#' @format A `SpatRaster` object.
#'
#' @name future_2100_ssp126_access
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' future_2100_ssp126_access <- terra::rast(system.file("extdata",
#'                                     "wc2.1_10m_bioc_ACCESS-CM2_ssp126_2081-2100.tif",
#'                                      package = "kuenm2"))
#' terra::plot(future_2100_ssp126_access)
NULL

#' SpatRaster Representing Future Conditions (2041-2060, SSP585, GCM: ACCESS-CM2)
#'
#' A raster layer containing bioclimatic variables representing future climatic
#' conditions (2041-2060) based on the ACCESS-CM2 General Circulation Model
#' under the SSP585 scenario. The variables were obtained at a 10 arc-minute
#' resolution and masked using the `m` region provided in the package. Data
#' sourced from WorldClim: \url{https://worldclim.org/data/cmip6/cmip6climate.html}
#'
#' @format A `SpatRaster` object.
#'
#' @name future_2050_ssp585_access
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' future_2050_ssp585_access <- terra::rast(system.file("extdata",
#'                                     "wc2.1_10m_bioc_ACCESS-CM2_ssp585_2041-2060.tif",
#'                                      package = "kuenm2"))
#' terra::plot(future_2050_ssp585_access)
NULL

#' SpatRaster Representing Future Conditions (2081-2100, SSP585, GCM: ACCESS-CM2)
#'
#' A raster layer containing bioclimatic variables representing future climatic
#' conditions (2081-2100) based on the ACCESS-CM2 General Circulation Model
#' under the SSP585 scenario. The variables were obtained at a 10 arc-minute
#' resolution and masked using the `m` region provided in the package. Data
#' sourced from WorldClim: \url{https://worldclim.org/data/cmip6/cmip6climate.html}
#'
#' @format A `SpatRaster` object.
#'
#' @name future_2100_ssp585_access
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' future_2100_ssp585_access <- terra::rast(system.file("extdata",
#'                                     "wc2.1_10m_bioc_ACCESS-CM2_ssp585_2081-2100.tif",
#'                                      package = "kuenm2"))
#' terra::plot(future_2100_ssp585_access)
NULL

#' SpatRaster Representing Future Conditions (2041-2060, SSP126, GCM: MIROC6)
#'
#' A raster layer containing bioclimatic variables representing future climatic
#' conditions (2041-2060) based on the MIROC6 General Circulation Model
#' under the SSP126 scenario. The variables were obtained at a 10 arc-minute
#' resolution and masked using the `m` region provided in the package. Data
#' sourced from WorldClim: \url{https://worldclim.org/data/cmip6/cmip6climate.html}
#'
#' @format A `SpatRaster` object.
#'
#' @name future_2050_ssp126_miroc
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' future_2050_ssp126_miroc <- terra::rast(system.file("extdata",
#'                                     "wc2.1_10m_bioc_MIROC6_ssp126_2041-2060.tif",
#'                                      package = "kuenm2"))
#' terra::plot(future_2050_ssp126_miroc)
NULL

#' SpatRaster Representing Future Conditions (2081-2100, SSP126, GCM: MIROC6)
#'
#' A raster layer containing bioclimatic variables representing future climatic
#' conditions (2081-2100) based on the MIROC6 General Circulation Model
#' under the SSP126 scenario. The variables were obtained at a 10 arc-minute
#' resolution and masked using the `m` region provided in the package. Data
#' sourced from WorldClim: \url{https://worldclim.org/data/cmip6/cmip6climate.html}
#'
#' @format A `SpatRaster` object.
#'
#' @name future_2100_ssp126_miroc
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' future_2100_ssp126_miroc <- terra::rast(system.file("extdata",
#'                                     "wc2.1_10m_bioc_MIROC6_ssp126_2081-2100.tif",
#'                                      package = "kuenm2"))
#' terra::plot(future_2100_ssp126_miroc)
NULL

#' SpatRaster Representing Future Conditions (2041-2060, SSP585, GCM: MIROC6)
#'
#' A raster layer containing bioclimatic variables representing future climatic
#' conditions (2041-2060) based on the MIROC6 General Circulation Model
#' under the SSP585 scenario. The variables were obtained at a 10 arc-minute
#' resolution and masked using the `m` region provided in the package. Data
#' sourced from WorldClim: \url{https://worldclim.org/data/cmip6/cmip6climate.html}
#'
#' @format A `SpatRaster` object.
#'
#' @name future_2050_ssp585_miroc
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' future_2050_ssp585_miroc <- terra::rast(system.file("extdata",
#'                                     "wc2.1_10m_bioc_MIROC6_ssp585_2041-2060.tif",
#'                                      package = "kuenm2"))
#' terra::plot(future_2050_ssp585_miroc)
NULL

#' SpatRaster Representing Future Conditions (2081-2100, SSP585, GCM: MIROC6)
#'
#' A raster layer containing bioclimatic variables representing future climatic
#' conditions (2081-2100) based on the MIROC6 General Circulation Model
#' under the SSP585 scenario. The variables were obtained at a 10 arc-minute
#' resolution and masked using the `m` region provided in the package. Data
#' sourced from WorldClim: \url{https://worldclim.org/data/cmip6/cmip6climate.html}
#'
#' @format A `SpatRaster` object.
#'
#' @name future_2100_ssp585_miroc
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' bring raster variables to analysis.
#'
#' @examples
#' future_2100_ssp585_miroc <- terra::rast(system.file("extdata",
#'                                     "wc2.1_10m_bioc_MIROC6_ssp585_2081-2100.tif",
#'                                      package = "kuenm2"))
#' terra::plot(future_2100_ssp585_miroc)
NULL

#' SpatVector Representing Calibration Area for *Myrcia hatschbachii*
#'
#' A spatial vector defining the calibration area used to extract background
#' points for fitting models of *Myrcia hatschbachii*. The area was generated by
#' creating a minimum convex polygon around presence records (`occ_data`), then
#' applying a 300 km buffer.
#'
#' @format A `Spatvector` object.
#'
#' @name m
#'
#' @return No return value. Used with function \code{\link[terra]{vect}} to
#' bring raster variables to analysis.
#'
#' @examples
#' m <- terra::vect(system.file("extdata",
#'                              "m.gpkg",
#'                               package = "kuenm2"))
#' terra::plot(m)
NULL

#' World country polygons from Natural Earth
#'
#' A spatial vector of the world countries. This is a simplified version of the `countries110` from rnaturalearth R package.
#'
#' @format A `Spatvector` object.
#'
#' @name world
#'
#' @return No return value. Used with function \code{\link[terra]{vect}} to
#' bring vector variables to analysis.
#'
#' @examples
#' m <- terra::vect(system.file("extdata",
#'                              "world.gpkg",
#'                               package = "kuenm2"))
#' terra::plot(m)
NULL
