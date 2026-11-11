# Read and prepare a species detection matrix

Read and prepare a species detection matrix

## Usage

``` r
read_matrix(
  data,
  timestep = 1,
  time_range = NULL,
  baseline = NULL,
  baseline_min_species = NULL
)
```

## Arguments

- data:

  A CSV file path or a numeric or logical detection matrix. In a CSV,
  the first column contains species names and remaining column names
  contain starting years. In a matrix, species names must be row names
  and starting years must be column names; do not include a species
  column inside the matrix. Detection values must be 0/1 or, for a
  logical matrix, FALSE/TRUE. Missing values are not allowed.

- timestep:

  Positive output interval width, in years. Defaults to 1.

- time_range:

  Optional numeric vector of length two giving inclusive limits on
  output interval labels.

- baseline:

  Optional starting interval label within the selected range.

- baseline_min_species:

  Optional positive integer. Select the first interval in the selected
  range with at least this many detected species.

## Value

A numeric 0/1 matrix with species row names and time column names. A
`preparation` attribute records the input type, input interval width,
output interval width, retained labels, baseline selection and coverage.
Coverage endpoints use an exclusive upper bound. The returned matrix has
class `modgee_matrix`; printing shows only matrix entries and row/column
names. Preparation metadata remain accessible through
`attr(M, "preparation")`. Existing informational messages are retained.

## Details

At least two time columns are required. Column names must represent
finite positive integers (zero, negative and fractional years are
rejected). Integer-valued numeric forms such as 2000.0 are accepted. The
names represent the starting years of consecutive, equally spaced
observation intervals. The input interval width is inferred from the
spacing between column names. Missing observation intervals must not be
omitted from the input. Equally spaced labels are interpreted as
consecutive intervals, not isolated years.

The output width must be an integer multiple of the input width. An
invalid width error reports both the supplied timestep and input width.
Input intervals must fit within calendar-aligned output intervals
without crossing their boundaries. Output labels follow
`floor(time / timestep) * timestep`. For input labels 1900, 1910, 1920,
use `timestep = 10` to retain the original interval width. Such data
cannot be split into annual records.

Processing proceeds in this order: aggregate intervals, apply
`time_range`, then select the baseline. Aggregation records whether a
species was detected at least once within each output interval. The
upper bound of `time_range` refers to an interval label, not its ending
year. If the requested range extends beyond the available output labels
but overlaps them, the overlap is retained with an informational message
reporting the retained label range before baseline selection. A range
with no overlapping labels is an error.

Invalid time labels are reported with their original names and
positions. Duplicate years are compared numerically, so 1900 and 1900.0
are duplicates. Time-column positions exclude the species-name column in
a CSV.

`baseline` and `baseline_min_species` are mutually exclusive. If neither
is supplied, the first selected interval is used. All-zero species rows
and the original species order are retained. Incomplete boundary
intervals are retained with an informational message. For annual input
ending in 2026 and `timestep = 10`, the column labelled 2020 summarizes
only 2020–2026. Unobserved years are not added.

Existing custom attributes on matrix input are not carried forward.
Preparation information is recalculated from the supplied time columns.
For already aggregated data, the original coverage of a partial interval
cannot be recovered from its starting-year label alone; each input
column is interpreted as covering one full inferred input interval.

## Usage Examples


    # Read annual records from the installed CSV.
    M_annual <- read_matrix(data = example_file)

    # Aggregate into ten-year intervals and select the baseline.
    M <- read_matrix(
      data = example_file,
      timestep = 10,
      baseline_min_species = 10
    )
    dim(M)
    attr(M, "preparation")

    # Prepare an existing matrix instead of a CSV file.
    M_decadal <- read_matrix(data = M_annual, timestep = 10)

    # Keep a specified range of interval labels.
    M_range <- read_matrix(
      data = example_file, timestep = 10,
      time_range = c(1900, 2000)
    )

    # Specify the baseline directly.
    M_baseline <- read_matrix(
      data = example_file, timestep = 10,
      baseline = 1900
    )

    # Apply a detection threshold within a selected time range.
    M_selected <- read_matrix(
      data = example_file, timestep = 10,
      time_range = c(1900, 2020),
      baseline_min_species = 10
    )

## See also

**Functions:**
[`fit_modgee`](https://zhizhouhe2006.github.io/modgee/reference/fit_modgee.md),
[`plot_modgee`](https://zhizhouhe2006.github.io/modgee/reference/plot_modgee.md),
[`save_result`](https://zhizhouhe2006.github.io/modgee/reference/save_result.md)

**Examples:**
[`example_file`](https://zhizhouhe2006.github.io/modgee/reference/example_file.md),
[`example_M`](https://zhizhouhe2006.github.io/modgee/reference/example_M.md),
[`example_groups`](https://zhizhouhe2006.github.io/modgee/reference/example_groups.md),
[`example_result`](https://zhizhouhe2006.github.io/modgee/reference/example_result.md),
[`example_grouped_result`](https://zhizhouhe2006.github.io/modgee/reference/example_grouped_result.md),
[`example_bootstrap_result`](https://zhizhouhe2006.github.io/modgee/reference/example_bootstrap_result.md)

## Console display

Every successful call prints the prepared matrix, including calls
assigned to a variable. The complete matrix is returned invisibly to
avoid duplicate automatic printing; preparation metadata are retained
but not displayed.
