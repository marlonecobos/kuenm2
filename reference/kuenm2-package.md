# kuenm2: Detailed Development of Ecological Niche Models

kuenm2 A new set of tools to help with the development of detailed
ecological niche models using multiple algorithms, at the moment Maxnet
and GLM. Pre-modeling analyses and explorations can be done to prepare
data. Model calibration (model selection) can be done by creating and
testing several candidate models, that are later selected based on a
multicriteria approach. Handy options for producing final models with
transfers are included. Other tools to assess extrapolation risks and
variability in model transfers are also available.

## Main functions by stage in the ENM process

### Pre-modeling steps

- Data preparation:
  [`initial_cleaning()`](https://marlonecobos.github.io/kuenm2/reference/initial_cleaning.md),
  [`advanced_cleaning()`](https://marlonecobos.github.io/kuenm2/reference/advanced_cleaning.md),
  [`prepare_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md),
  [`prepare_user_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_user_data.md)

- Data exploration:
  [`explore_calibration_hist()`](https://marlonecobos.github.io/kuenm2/reference/explore_calibration_hist.md),
  [`explore_partition_env()`](https://marlonecobos.github.io/kuenm2/reference/explore_partition_env.md),
  [`explore_partition_geo()`](https://marlonecobos.github.io/kuenm2/reference/explore_partition_geo.md),
  [`explore_partition_extrapolation()`](https://marlonecobos.github.io/kuenm2/reference/explore_partition_extrapolation.md),
  [`plot_calibration_hist()`](https://marlonecobos.github.io/kuenm2/reference/plot_calibration_hist.md),
  [`plot_explore_partition()`](https://marlonecobos.github.io/kuenm2/reference/plot_explore_partition.md)

### Modeling process

- Model calibration:
  [`calibration()`](https://marlonecobos.github.io/kuenm2/reference/calibration.md),
  [`select_models()`](https://marlonecobos.github.io/kuenm2/reference/select_models.md)

- Model exploration:
  [`fit_selected()`](https://marlonecobos.github.io/kuenm2/reference/fit_selected.md),
  [`variable_importance()`](https://marlonecobos.github.io/kuenm2/reference/variable_importance.md),
  [`plot_importance()`](https://marlonecobos.github.io/kuenm2/reference/plot_importance.md),
  [`response_curve()`](https://marlonecobos.github.io/kuenm2/reference/response_curve.md),
  [`all_response_curves()`](https://marlonecobos.github.io/kuenm2/reference/response_curve.md),
  [`bivariate_response()`](https://marlonecobos.github.io/kuenm2/reference/bivariate_response.md),
  [`partition_response_curves()`](https://marlonecobos.github.io/kuenm2/reference/partition_response_curves.md)

- Model projection:
  [`predict_selected()`](https://marlonecobos.github.io/kuenm2/reference/predict_selected.md),
  [`organize_for_projection()`](https://marlonecobos.github.io/kuenm2/reference/organize_for_projection.md),
  [`organize_future_worldclim()`](https://marlonecobos.github.io/kuenm2/reference/organize_future_worldclim.md),
  [`prepare_projection()`](https://marlonecobos.github.io/kuenm2/reference/prepare_projection.md),
  [`project_selected()`](https://marlonecobos.github.io/kuenm2/reference/project_selected.md)

### Post-modeling analysis

- Variability:
  [`projection_changes()`](https://marlonecobos.github.io/kuenm2/reference/projection_changes.md),
  [`projection_variability()`](https://marlonecobos.github.io/kuenm2/reference/projection_variability.md)

- Uncertainty:
  [`projection_mop()`](https://marlonecobos.github.io/kuenm2/reference/projection_mop.md)

## See also

Useful links:

- <https://marlonecobos.github.io/kuenm2/>

- Report bugs at <https://github.com/marlonecobos/kuenm2/issues>

## Author

**Maintainer**: Weverton C. F. Trindade <wevertonf1993@gmail.com>
([ORCID](https://orcid.org/0000-0003-2045-4555))

Authors:

- Luis F. Arias-Giraldo <lfarias.giraldo@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-4861-8064))

- Luis Osorio-Olvera <luismurao@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-0701-5398))

- A. Townsend Peterson <town@ku.edu>
  ([ORCID](https://orcid.org/0000-0003-0243-2379))

- Marlon E. Cobos <manubio13@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-2611-1767))
