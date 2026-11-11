# Print a prepared detection matrix

Display matrix entries and row/column names without preparation
metadata.

## Usage

``` r
# S3 method for class 'modgee_matrix'
print(x, ...)
```

## Arguments

- x:

  A matrix returned by
  [`read_matrix()`](https://zhizhouhe2006.github.io/modgee/reference/read_matrix.md).

- ...:

  Additional arguments passed to
  [`print.default()`](https://rdrr.io/r/base/print.default.html).

## Value

The original object, invisibly. Preparation metadata are retained.
