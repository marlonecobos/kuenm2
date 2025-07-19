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
#' A \code{data.frame} containing the coordinates of 51 valid occurrences of *Myrcia hatschbachii* (a tree endemic to southern Brazil). The valid occurrences were sourced from Trindade & Marques (2024).
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
