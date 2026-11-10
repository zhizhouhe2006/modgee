# MODGEE estimation from species detection histories

Fit the original MODGEE model to estimate detection parameters,
cumulative extinction and dark extinctions.

## Details

Use
[`fit_modgee`](https://zhizhouhe2006.github.io/modgee/reference/fit_modgee.md)
to fit a prepared species-by-interval binary matrix. Use
[`example_M`](https://zhizhouhe2006.github.io/modgee/reference/example_M.md)
for a documented development example.

Time aggregation and baseline selection are performed before fitting.
Grouped models share a temporal extinction baseline and may compare
fixed-k and free-k models. Bootstrap confidence intervals are not
computed by these functions.

## See also

**Functions:**
[`read_matrix`](https://zhizhouhe2006.github.io/modgee/reference/read_matrix.md),
[`fit_modgee`](https://zhizhouhe2006.github.io/modgee/reference/fit_modgee.md),
[`plot_modgee`](https://zhizhouhe2006.github.io/modgee/reference/plot_modgee.md)

**Examples:**
[`example_file`](https://zhizhouhe2006.github.io/modgee/reference/example_file.md),
[`example_M`](https://zhizhouhe2006.github.io/modgee/reference/example_M.md),
[`example_groups`](https://zhizhouhe2006.github.io/modgee/reference/example_groups.md),
[`example_result`](https://zhizhouhe2006.github.io/modgee/reference/example_result.md),
[`example_grouped_result`](https://zhizhouhe2006.github.io/modgee/reference/example_grouped_result.md)

## Author

**Maintainer**: Zhizhou He <hezhizhou2006@gmail.com>

Authors:

- Zhizhou He <hezhizhou2006@gmail.com>
