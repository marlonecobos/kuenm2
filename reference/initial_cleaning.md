# Initial occurrence data cleaning steps

Simple occurrence data cleaning procedures.

## Usage

``` r
initial_cleaning(data, species, x, y,
                 other_columns = NULL, keep_all_columns = TRUE,
                 sort_columns = TRUE, remove_na = TRUE, remove_empty = TRUE,
                 remove_duplicates = TRUE, by_decimal_precision = FALSE,
                 decimal_precision = 0, longitude_precision = NULL,
                 latitude_precision = NULL)

sort_columns(data, species, x, y, keep_all_columns = FALSE)

remove_missing(data, columns = NULL, remove_na = TRUE,
               remove_empty = TRUE, keep_all_columns = TRUE)

remove_duplicates(data, columns = NULL, keep_all_columns = TRUE)

remove_corrdinates_00(data, x, y)

filter_decimal_precision(data, x,
                         y, decimal_precision = 0,
                         longitude_precision = NULL,
                         latitude_precision = NULL)
```

## Arguments

- data:

  data.frame with occurrence records.

- species:

  (character) name of the column in `data` containing species name.

- x:

  (character) name of the column in `data` containing longitude values.

- y:

  (character) name of the column in `data` containing latitude values.

- other_columns:

  (character) vector of other column name(s) in `data` to be considered
  while performing cleaning steps, default = NULL.

- keep_all_columns:

  (logical) whether to keep all columns in `data`. Default = TRUE.

- sort_columns:

  (logical) whether to sort species, longitude, and latitude columns in
  `data`. Default = TRUE.

- remove_na:

  (logical) whether to remove NA values in the columns considered.
  Default = TRUE.

- remove_empty:

  (logical) whether to remove empty (missing) values in the columns
  considered. Default = TRUE.

- remove_duplicates:

  (logical) whether to remove duplicates in the columns considered.
  Default = TRUE.

- by_decimal_precision:

  (logical) whether to remove certain records with coordinate precision
  lower than that of the following three parameters. Default = FALSE

- decimal_precision:

  (numeric) decimal precision threshold for coordinates. Default = 0.
  Ignored if the following two parameters are defined.

- longitude_precision:

  (numeric) decimal precision threshold for longitude. Default = NULL.

- latitude_precision:

  (numeric) decimal precision threshold for latitude. Default = NULL.

- columns:

  (character) vector of additional column name(s) in `data` to be
  considered while removing missing or duplicate records, default =
  NULL.

## Value

A data.frame with resulting occurrence records.

## Details

Function `initial_cleaning` helps to perform all simple steps of data
cleaning.

## See also

[`advanced_cleaning`](https://marlonecobos.github.io/kuenm2/reference/advanced_cleaning.md)

## Examples

``` r
# Import occurrences
data(occ_data_noclean, package = "kuenm2")

# remove missing data
mis <- remove_missing(data = occ_data_noclean, columns = NULL, remove_na = TRUE,
                      remove_empty = TRUE)

# remove exact duplicates
mis_dup <- remove_duplicates(data = mis, columns = NULL, keep_all_columns = TRUE)

# remove records with 0 for x and y coordinates
mis_dup_00 <- remove_corrdinates_00(data = mis_dup, x = "x", y = "y")

# remove coordinates with low decimal precision.
mis_dup_00_dec <- filter_decimal_precision(data = mis_dup_00, x = "x", y = "y",
                                           decimal_precision = 2)

# all basic cleaning steps
clean_init <- initial_cleaning(data = occ_data_noclean, species = "species",
                               x = "x", y = "y", remove_na = TRUE,
                               remove_empty = TRUE, remove_duplicates = TRUE,
                               by_decimal_precision = TRUE,
                               decimal_precision = 2)
```
