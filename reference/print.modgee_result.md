# Print a MODGEE result

Display the overall summary, group summaries and k comparison when
present. The k comparison includes an AIC interpretation and, where
supported, its direction. All other result components remain available.

## Usage

``` r
# S3 method for class 'modgee_result'
print(x, ...)
```

## Arguments

- x:

  A result returned by
  [`fit_modgee()`](https://zhizhouhe2006.github.io/modgee/reference/fit_modgee.md).

- ...:

  Additional arguments passed to
  [`print()`](https://rdrr.io/r/base/print.html).

## Value

The complete object, invisibly.
