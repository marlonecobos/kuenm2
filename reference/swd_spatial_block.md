# Prepared data with spatial blocks created with ENMeval

A `prepared_data` object resulted from
[`prepare_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)
to calibrate models using 'glmnet' algorithm. In this object, the
original partitioning was replaced with spatial blocks generated using
the `get.block()` method from the ENMeval R package.

## Usage

``` r
data("sp_swd")
```

## Format

An object of class `prepared_data` of length 13.
