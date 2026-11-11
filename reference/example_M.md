# Prepared bird detection matrix

A numeric 0/1 matrix prepared with read_matrix(example_file, timestep =
10, baseline_min_species = 10). Species names are row names and interval
labels are column names. The preparation attribute records aggregation
and coverage; the final interval contains records only through 2020.

## Usage

``` r
example_M
```

## Format

A numeric matrix.

## Source

Maintainer-provided birds_matrix.csv and birds_groups.csv. Verify
original data attribution and redistribution terms before public
release.

## Usage Examples


    dim(example_M)

## See also

**Functions:**
[`read_matrix`](https://zhizhouhe2006.github.io/modgee/reference/read_matrix.md),
[`fit_modgee`](https://zhizhouhe2006.github.io/modgee/reference/fit_modgee.md),
[`plot_modgee`](https://zhizhouhe2006.github.io/modgee/reference/plot_modgee.md)

**Examples:**
[`example_file`](https://zhizhouhe2006.github.io/modgee/reference/example_file.md),
[`example_groups`](https://zhizhouhe2006.github.io/modgee/reference/example_groups.md),
[`example_result`](https://zhizhouhe2006.github.io/modgee/reference/example_result.md),
[`example_grouped_result`](https://zhizhouhe2006.github.io/modgee/reference/example_grouped_result.md)

## Display

Printing example_M displays only the matrix. Preparation metadata remain
available through attr(example_M, "preparation") and are retained for
fitting. The object inherits from modgee_matrix.
