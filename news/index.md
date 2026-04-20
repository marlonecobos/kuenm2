# Changelog

## kuenm2 0.1.3

- In
  [`bivariate_response()`](https://marlonecobos.github.io/kuenm2/reference/bivariate_response.md),
  store and restore graphical parameters only when plotting area is
  modified.
- In `response_curves()`, omit graphical parameter storing and resetting
  to allow plotting in multipanel figures.
- In
  [`projection_changes()`](https://marlonecobos.github.io/kuenm2/reference/projection_changes.md),
  create `data.frame` to set levels representing changes outside of the
  loop. This ensures it works correctly even when `by_gcm` is set to
  `FALSE`.

## kuenm2 0.1.2

CRAN release: 2026-03-29

- Initial CRAN submission.
