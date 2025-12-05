# Extract predictor names from formulas

Extract predictor names from formulas

## Usage

``` r
extract_var_from_formulas(formulas, ...)
```

## Arguments

- formulas:

  (character or formula) model formulas.

- ...:

  Arguments to pass to
  [`all.vars()`](https://rdrr.io/r/base/allnames.html)

## Value

A character vector or a list of the same length as `formulas`,
containing the names of the predictors each formula.

## Examples

``` r
# Import an example of calibration results
data(calib_results_maxnet, package = "kuenm2")

# Extract predictor names
vars <- extract_var_from_formulas(calib_results_maxnet$formula_grid$Formulas)
```
