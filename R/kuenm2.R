#' kuenm2: Detailed Development of Ecological Niche Models
#'
#' kuenm2 A new set of tools to help with the development of detailed ecological
#' niche models using multiple algorithms. Pre-modeling analyses and
#' explorations can be done to prepare data. Model calibration (model selection)
#' can be done by creating and testing several candidate models. Handy options
#' for producing final models with transfers are included. Other tools to
#' assess extrapolation risks and variability in model transfers are also
#' available.
#'
#' @section Functions by stage in the ENM process:
#'
#' - Data preparation: [initial_cleaning()], [advanced_cleaning()],
#' [prepare_data()], [prepare_user_data()]
#' - Model calibration: [calibration()], [select_models()]
#' - Model exploration: [fit_selected()], [variable_importance()],
#' [plot_importance()], [response_curve()], [bivariate_response()]
#' - Model projection: [predict_selected()], [project_selected()]
#' - Post-modeling: [projection_variability()], [projection_changes()],
#' [projection_mop()]
#'
#'
"_PACKAGE"
