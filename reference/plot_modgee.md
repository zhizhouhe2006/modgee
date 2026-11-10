# Plot MODGEE result

Visualize cumulative extinction, species extinction probabilities, or
conditional expected extinction years from an existing MODGEE fit.

## Usage

``` r
plot_modgee(result, type = "curve")
```

## Arguments

- result:

  A result list returned by
  [`fit_modgee()`](https://zhizhouhe2006.github.io/modgee/reference/fit_modgee.md)
  using the current output format. Older objects must be refitted or
  updated to that format.

- type:

  A single string: `"curve"` (default), `"probability"`, or `"time"`.
  Full names are required.

## Value

A `ggplot` object. The function neither saves a file nor explicitly
prints the plot. Use [`print()`](https://rdrr.io/r/base/print.html)
inside loops and
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
to save it.

## Details

All plots use a white background, a closed black panel border, no grid,
larger axis text, and a two-column legend below the panel. Groups use
green, orange, blue and purple first; additional groups extend the
palette.

**Cumulative extinction curves**

`type = "curve"` plots `result$curve$time` against
`result$curve$cumulative_extinction_rate`. The axes are labelled *Year*
and *Cumulative extinction rate*. Rates are shown as percentages. The
overall curve is a thick black line. For grouped fits, coloured group
curves from `result$group_results$curve` are also shown, and the overall
curve is labelled *Overall* in the legend. The fitted overall curve is
used directly, not an unweighted average of group curves. No smoothing
is added. When bootstrap intervals are available, a grey band shows the
overall pointwise confidence interval from `result$bootstrap$curve`, at
the confidence level used in fitting. Group curves retain their original
colours without additional bands. Missing interval bounds leave gaps in
the band.

**Extinction probability histograms**

`type = "probability"` counts recorded species using
`result$species$p_extinct`. Bins have width 0.1, with boundaries
`[0, 0.1]`, `(0.1, 0.2]`, ..., `(0.9, 1]`. The axes are *Extinction
probability* and *Number of species*. Grouped fits use stacked colours,
so each bar's total is the total number of species in that interval.
Counts are not normalized by group size.

**Expected extinction years**

`type = "time"` plots `estimated_extinction_year` against `p_extinct`,
with one point per eligible species. The axes are *Expected extinction
year* and *Extinction probability*. Estimated years are conditional
means given extinction, not confirmed extinction dates. They are read
directly from the species table and are already mapped to the numeric
time labels supplied during fitting. If interval indices were supplied
instead of years, the plot retains those indices despite the default
year labels; labels can be changed with
[`ggplot2::labs()`](https://ggplot2.tidyverse.org/reference/labs.html).
No jitter or regression line is added.

**Groups and missing estimates**

Group membership comes from `result$species$belonged_group`. Group order
and the colour mapping use `result$group_results$summary`. The same
group has the same colour across the three plot types. All-zero species
are excluded from species-level plots. Groups containing only all-zero
species have no fitted curve and no recorded species to plot. Dark
extinction estimates are not converted into artificial species records.

For probability histograms, missing or non-finite estimates are omitted
with a warning giving their count. Time plots silently omit species with
missing or non-finite estimates or zero extinction probability, because
their conditional extinction time is undefined. Missing curve estimates,
unknown groups, duplicate species identifiers, invalid probabilities, or
no eligible observations result in an error. Only the fields needed by
the requested plot are required.

## Usage Examples


    plot_modgee(example_result)
    plot_modgee(example_result, type = "probability")
    plot_modgee(example_result, type = "time")

    plot_modgee(example_grouped_result)
    plot_modgee(example_grouped_result, type = "probability")
    plot_modgee(example_grouped_result, type = "time")

    # Precomputed bootstrap fit with an overall confidence band.
    plot_modgee(example_bootstrap_result)

    # Save a plot to a PNG file.
    p <- plot_modgee(example_bootstrap_result)
    ggplot2::ggsave("modgee_curve.png", plot = p,
                    width = 8, height = 6, units = "in", dpi = 300)

## See also

**Functions:**
[`read_matrix`](https://zhizhouhe2006.github.io/modgee/reference/read_matrix.md),
[`fit_modgee`](https://zhizhouhe2006.github.io/modgee/reference/fit_modgee.md),
[`save_result`](https://zhizhouhe2006.github.io/modgee/reference/save_result.md)

**Examples:**
[`example_file`](https://zhizhouhe2006.github.io/modgee/reference/example_file.md),
[`example_M`](https://zhizhouhe2006.github.io/modgee/reference/example_M.md),
[`example_groups`](https://zhizhouhe2006.github.io/modgee/reference/example_groups.md),
[`example_result`](https://zhizhouhe2006.github.io/modgee/reference/example_result.md),
[`example_grouped_result`](https://zhizhouhe2006.github.io/modgee/reference/example_grouped_result.md),
[`example_bootstrap_result`](https://zhizhouhe2006.github.io/modgee/reference/example_bootstrap_result.md)
