# Maxent-like Generalized Linear Models (GLM)

This function fits a Generalized Linear Model (GLM) to binary
presence-background data. It allows for the specification of custom
weights, with a default in which presences have a weight of 1 and
background 100.

## Usage

``` r
glm_mx(formula, family = binomial(link = "cloglog"), data,
       weights = NULL, ...)
```

## Arguments

- formula:

  A formula specifying the model to be fitted, in the format used by
  [`glm`](https://rdrr.io/r/stats/glm.html).

- family:

  A description of the error distribution and link function to be used
  in the model. Defaults to `binomial(link = "cloglog")`, which is
  commonly used for presence-background data.

- data:

  A `data.frame` containing the variables in the model. Must include a
  column named `pr_bg` that indicates whether a record is a presence (1)
  or background (0), and at least another column with an independent
  variable (predictor).

- weights:

  Optional. A numeric vector of weights for each observation. If not
  provided, default weights of 1 for presences and 100 for background
  are used.

- ...:

  Additional arguments to be passed to
  [`glm`](https://rdrr.io/r/stats/glm.html).

## Value

A fitted [`glm`](https://rdrr.io/r/stats/glm.html) object. The model
object includes the minimum and maximum values of the non-factor
variables in the dataset, stored as `model$varmin` and `model$varmax`.

## Details

For more details about glms using presence and background emulating what
Maxent does, see Fithian and Hastie (2013) <doi:10.1214/13-AOAS667>.
