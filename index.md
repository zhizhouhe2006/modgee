# modgee

Estimate extinction from species detection histories with MODGEE
(Matrix-Of-Detections-Gives-Extinction-Estimates).

## Installation

Install the development version from GitHub:

``` r

install.packages("remotes")
remotes::install_github("zhizhouhe2006/modgee")
```

The package contains C++ code. Installing from source on Windows
requires Rtools compatible with your R version.

## Quick start

``` r

library(modgee)

# Prepare detections without overwriting the installed examples.
M <- read_matrix(example_file, timestep = 10, baseline_min_species = 10)
set.seed(123)
result <- fit_modgee(M)
result$summary
plot_modgee(result)
save_result(result, path = "results", prefix = "analysis")
```

Input rows represent species and columns represent equally spaced
numeric time labels. Entries must be 0 or 1. Row and column names are
required.

## Grouped models and k comparison

``` r

set.seed(123)
result_grouped <- fit_modgee(
  M = example_M, groups = example_groups,
  reference_group = "primary + old secondary", do_k = TRUE
)
result_grouped$group_results$summary
result_grouped$do_k
plot_modgee(result_grouped)
```

## Bootstrap confidence intervals

``` r

set.seed(123)
result_bootstrap <- fit_modgee(
  M = example_M, groups = example_groups,
  reference_group = "primary + old secondary", bootstrap = 200, CI = 0.95
)
plot_modgee(result_bootstrap)
save_result(result_bootstrap, path = "results", prefix = "bootstrap")
```

Bootstrap refits models and may take substantial time. Checkpoints are
stored under `modgee_bootstrap` in the working directory. Rerun the same
call there to resume. Inspect convergence and valid replicate counts
before interpreting the intervals.

## Precomputed examples

The package provides `example_file`, `example_M`, `example_groups`,
`example_result`, `example_grouped_result` and
`example_bootstrap_result`. The result objects are precomputed; viewing
or plotting them does not refit models. The bootstrap example uses only
10 replicates for demonstration.

``` r

plot_modgee(example_bootstrap_result)
save_result(example_bootstrap_result, path = "results", prefix = "example")
```

## CSV exports

[`save_result()`](https://zhizhouhe2006.github.io/modgee/reference/save_result.md)
writes `summary.csv`, `time_result.csv` and `species_result.csv`, with
the requested prefix. It additionally writes `k_comparison.csv` when a k
comparison is present and `CI.csv` when bootstrap curve intervals are
present. Overall and group results share the same files. Existing files
are protected unless `overwrite = TRUE`.

## Documentation and development

See
[`?read_matrix`](https://zhizhouhe2006.github.io/modgee/reference/read_matrix.md),
[`?fit_modgee`](https://zhizhouhe2006.github.io/modgee/reference/fit_modgee.md),
[`?plot_modgee`](https://zhizhouhe2006.github.io/modgee/reference/plot_modgee.md)
and
[`?save_result`](https://zhizhouhe2006.github.io/modgee/reference/save_result.md)
for argument and output definitions.

``` r

devtools::test()
devtools::check()
```

Scripts in `data-raw/` regenerate example datasets; ordinary
installation does not need to run them.

## License and attribution

See `LICENSE`, `LICENSE.md` and `inst/COPYRIGHTS` for copyright and
licensing information. Report bugs through the repository Issues page.
