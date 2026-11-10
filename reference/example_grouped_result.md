# Fitted grouped example model

A complete result list computed with fit_modgee(M = example_M, groups =
example_groups, reference_group = "primary + old secondary", do_k =
TRUE) after set.seed(123). The stored object includes the fixed-k and
free-k comparison. Loading it does not refit the model.

## Usage

``` r
example_grouped_result
```

## Format

A fitted MODGEE result list.

## Source

Maintainer-provided birds_matrix.csv and birds_groups.csv. Verify
original data attribution and redistribution terms before public
release.

## Usage Examples


    example_grouped_result$group_results$summary
    example_grouped_result$do_k

## See also

**Functions:**
[`read_matrix`](https://zhizhouhe2006.github.io/modgee/reference/read_matrix.md),
[`fit_modgee`](https://zhizhouhe2006.github.io/modgee/reference/fit_modgee.md),
[`plot_modgee`](https://zhizhouhe2006.github.io/modgee/reference/plot_modgee.md)

**Examples:**
[`example_file`](https://zhizhouhe2006.github.io/modgee/reference/example_file.md),
[`example_M`](https://zhizhouhe2006.github.io/modgee/reference/example_M.md),
[`example_groups`](https://zhizhouhe2006.github.io/modgee/reference/example_groups.md),
[`example_result`](https://zhizhouhe2006.github.io/modgee/reference/example_result.md)
