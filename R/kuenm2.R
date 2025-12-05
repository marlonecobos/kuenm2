#' kuenm2: Detailed Development of Ecological Niche Models
#'
#' kuenm2 A new set of tools to help with the development of detailed ecological
#' niche models using multiple algorithms, at the moment Maxnet and GLM.
#' Pre-modeling analyses and explorations can be done to prepare data. Model
#' calibration (model selection) can be done by creating and testing several
#' candidate models, that are later selected based on a multicriteria approach.
#' Handy options for producing final models with transfers are included. Other
#' tools to assess extrapolation risks and variability in model transfers are
#' also available.
#'
#' @section Main functions by stage in the ENM process:
#'
#' ## Pre-modeling steps
#'
#' - Data preparation: [initial_cleaning()], [advanced_cleaning()],
#' [prepare_data()], [prepare_user_data()]
#' - Data exploration: [explore_calibration_hist()], [explore_partition_env()],
#' [explore_partition_geo()], [explore_partition_extrapolation()],
#' [plot_calibration_hist()], [plot_explore_partition()]
#'
#' ## Modeling process
#' - Model calibration: [calibration()], [select_models()]
#' - Model exploration: [fit_selected()], [variable_importance()],
#' [plot_importance()], [response_curve()], [all_response_curves()],
#' [bivariate_response()], [partition_response_curves()]
#' - Model projection: [predict_selected()], [organize_for_projection()],
#' [organize_future_worldclim()], [prepare_projection()], [project_selected()]
#'
#' ## Post-modeling analysis
#' - Variability: [projection_changes()], [projection_variability()]
#' - Uncertainty: [projection_mop()]
#'
"_PACKAGE"
