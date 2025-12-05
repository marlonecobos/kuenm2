# Summary plot for variable importance in models

See details in
[`plot_importance`](https://rdrr.io/pkg/enmpa/man/plot_importance.html)

## Usage

``` r
plot_importance(x, xlab = NULL, ylab = "Relative contribution",
                main = "Variable importance", extra_info = TRUE, ...)
```

## Arguments

- x:

  data.frame output from
  [`variable_importance`](https://marlonecobos.github.io/kuenm2/reference/variable_importance.md)().

- xlab:

  (character) a label for the x axis.

- ylab:

  (character) a label for the y axis.

- main:

  (character) main title for the plot.

- extra_info:

  (logical) when results are from more than one model, it adds
  information about the number of models using each predictor and the
  mean contribution found.

- ...:

  additional arguments passed to barplot or boxplot.

  Value A barplot or boxplot depending on the number of models
  considered.
