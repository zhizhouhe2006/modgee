# Save MODGEE results as CSV files

Export existing estimates without refitting or recomputing intervals.

## Usage

``` r
save_result(result, path = ".", prefix = "modgee", overwrite = FALSE)
```

## Arguments

- result:

  A complete result returned by
  [`fit_modgee()`](https://zhizhouhe2006.github.io/modgee/reference/fit_modgee.md).

- path:

  Output directory, created recursively if needed. Defaults to the
  current working directory.

- prefix:

  Filename prefix. An underscore is appended unless the prefix is empty
  or already ends with one. Defaults to `"modgee"`.

- overwrite:

  Allow replacement of existing target files. Defaults to FALSE.

## Value

Writes three to five CSV files to `path`. Each filename starts with
`prefix` followed by an underscore; an empty prefix adds nothing. For
example, `prefix = "hk"` produces `hk_summary.csv` and, when bootstrap
results are available, `hk_CI.csv`.

**CSV files**

- `summary.csv`:

  One overall row, followed by one row per supplied group when groups
  are present. Columns `scope` ("overall" or "group") and `group`
  identify the rows; the overall group label is `NA`. The remaining
  overall columns are `baseline_timestep`, `final_timestep`,
  `excluded_species`, `known_species`, `singleton_species`,
  `known_extant`, `known_extinct`, `dark_extinct`, `total_extinct`,
  `original_richness`, `cumulative_extinction_rate` and `AIC`, copied
  from `result$summary`. Group rows contain the corresponding group
  estimates. With groups, `l` is also included: it is the group's hazard
  multiplier, with 1 for the reference group and NA for the overall row
  and unfitted groups. AIC is an overall model statistic and is NA on
  group rows. Groups containing only all-zero species remain in this
  table, with undefined estimates as NA.

- `time_result.csv`:

  A long-format table: one row per time step for the overall model,
  followed by one row per time step for each fitted group. Columns are
  `scope`, `group`, `time`, `h`, `m` and `cumulative_extinction_rate`.
  The overall rows contain the shared survey parameter h, the reference
  hazard m and the overall cumulative extinction rate. Group rows
  contain their own hazard m and cumulative extinction rate; their h is
  NA because the shared value is already recorded in the overall rows.
  Time labels are copied from the fitted result.

- `species_result.csv`:

  One row per input species, in the original input order. Copies
  `result$species`, including `species`, `belonged_group` (when
  grouped), `detected_timesteps`, `is_singleton`,
  `first_detection_timestep`, `last_detection_timestep`, `d`,
  `p_extinct` and `estimated_extinction_year`. All-zero species are
  retained with NA estimates. The extinction year is conditional on
  extinction, as defined in
  [`fit_modgee`](https://zhizhouhe2006.github.io/modgee/reference/fit_modgee.md).

- `k_comparison.csv`:

  Written only when `result$do_k` is not NULL. Copies that table:
  fixed-k and free-k model rows, with `model`, `k`, group multipliers
  `l_<group>` when present, `AIC` and `converged`. No new models are
  fitted.

- `CI.csv`:

  Written only when `result$bootstrap` is not NULL. Combines
  `result$bootstrap$curve` and, when present,
  `result$bootstrap$group_results$curve`. Overall rows come first,
  followed by group rows. Columns `scope`, `group` and `time` identify
  each curve point. `estimate` is the original fitted cumulative
  extinction rate; `lower` and `upper` are its pointwise bootstrap
  confidence limits at the confidence level used in
  [`fit_modgee()`](https://zhizhouhe2006.github.io/modgee/reference/fit_modgee.md).
  `n_valid` is the number of valid bootstrap estimates used for that
  point. Undefined limits remain NA. Only curve intervals are exported,
  not intervals for summary metrics, parameters or individual species.

## Details

Always writes summary.csv, time_result.csv and species_result.csv, with
the requested prefix. Summary and time tables contain overall rows
first, followed by group rows. The scope column is "overall" or "group";
group is NA for overall rows, so a real group named "overall" remains
unambiguous. Missing columns are filled with NA. Group summaries retain
l; their AIC is NA. Group time rows contain the group's own m and
cumulative extinction rate; h is NA because h is shared and recorded on
overall rows only. All-zero-only groups remain in summary but have no
fitted curve rows. Species results retain all input rows, columns and
their original order.

If result\$do_k is non-NULL, k_comparison.csv contains that table
unchanged apart from omission of row names. If result\$bootstrap is
non-NULL, CI.csv combines bootstrap\$curve and
bootstrap\$group_results\$curve, with scope and group columns. These are
pointwise cumulative-extinction intervals: time, estimate, lower, upper
and n_valid. No summary-metric intervals or replicate tables are
exported.

CSV files use UTF-8, comma separators, NA for missing values and no row
numbers. Estimates are not rounded, clipped or recomputed. All tables
and filename conflicts are checked before writing. CSVs are staged
before installation; a failed installation attempts to restore previous
files. This is not a concurrent-writer transaction; use a separate
directory for simultaneous exports. Optional files left by earlier calls
are not deleted; use distinct prefixes for different analyses. Use
saveRDS to preserve the complete R object, including bootstrap settings
and recovery information.

## Usage Examples


    # Save ungrouped results.
    save_result(example_result, path = "results")

    # Save grouped results.
    save_result(example_grouped_result, path = "results", prefix = "grouped")

    # Save bootstrap results.
    save_result(example_bootstrap_result, path = "results", prefix = "bootstrap")

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
[`example_grouped_result`](https://zhizhouhe2006.github.io/modgee/reference/example_grouped_result.md),
[`example_bootstrap_result`](https://zhizhouhe2006.github.io/modgee/reference/example_bootstrap_result.md)
