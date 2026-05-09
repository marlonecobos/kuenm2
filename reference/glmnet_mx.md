# Maxent-like glmnet models

This function fits Maxent-like models using the `glmnet` package,
designed for presence-background data.

## Usage

``` r
glmnet_mx(p, data, f, regmult = 1.0, regfun = maxnet.default.regularization,
          addsamplestobackground = TRUE, weights = NULL, ...)
```

## Arguments

- p:

  A vector of binary presence-background labels, where 1 indicates
  presence and 0 indicates background.

- data:

  A `data.frame` containing the predictor variables for the model. This
  must include the same number of rows as the length of `p`.

- f:

  A formula specifying the model to be fitted, in the format used by
  [`model.matrix`](https://rdrr.io/r/stats/model.matrix.html).

- regmult:

  (numeric) Regularization multiplier, default is 1.0.

- regfun:

  A function that calculates regularization penalties. Default is
  `maxnet.default.regularization`.

- addsamplestobackground:

  (logical) Whether to add presence points not in the background to the
  background data. Default is `TRUE`.

- weights:

  (numeric) A numeric vector of weights for each observation. Default is
  `NULL`, which sets weights to 1 for presence points and 100 for
  background points.

- ...:

  Additional arguments to pass to
  [`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html).

## Value

A fitted Maxent-like model object of class `glmnet_mx`, which includes
model coefficients, AIC (if requested), and other elements such as
feature mins and maxes, sample means, and entropy.

## Details

This function is modified from the package maxnet and fits a Maxent-like
model using regularization to avoid over-fitting. Regularization weights
are computed using a provided function (which can be changed) and can be
multiplied by a regularization multiplier (`regmult`). The function also
includes an option to calculate AIC.
