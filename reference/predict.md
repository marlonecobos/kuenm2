# Predict method for glmnet_mx (maxnet) models

Predict method for glmnet_mx (maxnet) models

## Usage

``` r
predict.glmnet_mx(object, newdata, clamp = FALSE,
                  type = c("link", "exponential", "cloglog", "logistic",
                  "cumulative"))
```

## Arguments

- object:

  a glmnet_mx object.

- newdata:

  data to predict on.

- clamp:

  (logical) whether to clamp predictions. Default = FALSE.

- type:

  (character) type of prediction to be performed. Options are: "link",
  "exponential", "cloglog", "logistic", and cumulative. Defaults to
  "link" if not defined.

## Value

A glmnet_mx (maxnet) prediction.
