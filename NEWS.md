# kuenm2 0.1.3

* In `bivariate_response()`, store and restore graphical parameters only when plotting area is modified.
* In `response_curves()`, omit graphical paramater storing and resetting to allow plotting in multipanel figures.
* In `projection_changes()`, create `data.frame` to set levels representing changes outside of the loop. This ensures it works correctly even when `by_gcm` is set to `FALSE`.

# kuenm2 0.1.2

* Initial CRAN submission.
