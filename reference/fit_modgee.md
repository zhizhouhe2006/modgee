# Fit the original MODGEE model

Estimate detection parameters, cumulative extinction and dark
extinctions from species detection histories. This implementation
supports mutually exclusive species groups sharing a temporal extinction
baseline.

## Usage

``` r
fit_modgee(
  M,
  groups = NULL,
  reference_group = NULL,
  do_k = FALSE,
  bootstrap = NULL,
  CI = 0.95
)
```

## Arguments

- M:

  A numeric or logical matrix with species in rows and equally spaced
  time intervals in columns. Entries must be 0 or 1, with no missing
  values. A value of 1 indicates at least one detection during an
  interval; 0 indicates no recorded detection, not confirmed absence or
  extinction. Row names are required and must identify species. Leading
  and trailing whitespace is removed before checking that names are
  nonempty and unique. Column names are required and must convert to
  finite, unique, strictly increasing, equally spaced numeric time
  labels. They are the only source of time labels; there is no separate
  `times` argument. Character and complex matrices are not accepted.
  Logical entries are converted to numeric 0/1 values.

- groups:

  Optional CSV path, two-column data frame (species, group), character
  vector, factor, or numeric/logical 0/1 membership matrix. Names are
  trimmed and matched to all input species, including all-zero rows.
  Unnamed vectors and matrices follow input row order. Extra named
  species are ignored with a warning. Each species must have one group.
  CSV and character labels use radix sorting; factors use level order;
  membership matrices use column order. The first active group is the
  reference unless reference_group is specified. Empty groups are
  removed and reported. Every active group needs at least two
  repeated-detection species (a minimum safeguard, not a guarantee of
  statistical identifiability).

- reference_group:

  Optional single group name. Leading and trailing whitespace is
  removed. When NULL, the first active group under the default ordering
  is the reference. An explicit name overrides factor levels and
  membership-matrix column order; other groups retain their relative
  order. The named group must contain recorded species and meet the
  repeated-detection requirement. Supplying a reference without groups
  is an error. The reference log multiplier is fixed at zero, so its
  hazard multiplier is one. The selected name is returned in
  `groups$reference_group`.

- do_k:

  Logical scalar. If TRUE, compare the fixed-k model with a model that
  estimates a single k shared by all species. The two models share the
  first-stage detection estimates and use the same conditional
  likelihood in the extinction stage. For the alternative, the
  extinction hazard is `m[t] * d[i]^k * l[group[i]]`, with the group
  multiplier equal to 1 when groups are absent. The alternative starts
  at k = 0. Main estimates use k = 0.

- bootstrap:

  NULL (default), or an integer of at least 3 giving the number of
  scheduled full-pipeline species bootstrap replicates. Failed and
  nonconverged replicates count toward this number and are not
  resampled.

- CI:

  Confidence level strictly between 0 and 1; default 0.95. Intervals are
  pointwise percentile intervals (quantile type 7).

## Value

A list containing the following components:

- `summary`:

  A named list with the following entries, in order:

  `baseline_timestep`, `final_timestep`

  :   First and last numeric column labels of the supplied matrix
      (calendar years when year labels are supplied).

  `excluded_species`

  :   Number of all-zero rows in the supplied matrix.

  `known_species`

  :   Number of nonzero rows.

  `singleton_species`

  :   Number of species detected in exactly one timestep.

  `known_extant`

  :   Expected surviving recorded species:
      `known_species - known_extinct`. This is not the number detected
      in the final timestep.

  `known_extinct`

  :   Sum of extinction probabilities over recorded species.

  `dark_extinct`

  :   Estimated unrecorded extinctions; negative estimates are retained.
      This is not the number of all-zero rows.

  `total_extinct`

  :   `known_extinct + dark_extinct`.

  `original_richness`

  :   `known_species + dark_extinct`.

  `cumulative_extinction_rate`

  :   Final cumulative extinction proportion, on the 0–1 scale; equals
      `total_extinct / original_richness`.

  `AIC`

  :   Fixed-k main-model AIC, using the second-stage conditional log
      likelihood and the original two-stage parameter-count convention:
      `2 * p - 2 * logLik`, where `p = (R - 1) + T + (T - 3) + (G - 1)`;
      logLik is the internal second-stage conditional log likelihood; R
      is the number of repeated-detection species, T the number of
      timesteps, and G the number of active groups (one without groups).
      Returns `NA` if either main fitting stage did not converge. When
      `do_k = TRUE`, matches the fixed-k comparison row.

- `parameters`:

  Named vectors `d`, `h` and `m`. Grouped models also have `l`, the
  actual hazard multipliers, with reference value one. The starting
  values `starts$ls_start` remain log multipliers.

- `curve`:

  Data frame: `time`, `h`, `m`, `cumulative_extinction_rate`. With
  groups, m is the reference hazard, and the overall rate is weighted by
  inferred baseline richness.

- `species`:

  Data frame in original input order, with columns: `species`,
  `belonged_group` (only with groups), `detected_timesteps`,
  `is_singleton`, `first_detection_timestep`, `last_detection_timestep`,
  `d`, `p_extinct`, `estimated_extinction_year`. Detection timesteps are
  counted, not listed. First and last detections use actual numeric
  column labels. Estimated extinction year is the conditional mean
  extinction timestep mapped onto those labels, not a confirmed
  extinction date. All-zero rows have zero detected timesteps, FALSE for
  is_singleton, and NA for detection years and estimates. Undefined
  extinction times are NA.

- `groups`:

  NULL without groups; otherwise matching information, `species_group`,
  membership matrix `H`, active `levels`, `reference_group`,
  `removed_groups`, and `matching`.

- `group_results`:

  NULL without groups; otherwise a list of two data frames:

  `curve`

  :   `group`, `time`, `m`, `cumulative_extinction_rate`, for active
      groups only. The hazard m is the group's own hazard.

  `summary`

  :   `group`, then `l` (the actual hazard multiplier, not its
      logarithm), followed by the same columns as the overall summary,
      except AIC. Counts and estimates refer to each group. Supplied
      groups containing only all-zero species are included for counting;
      their dark_extinct, total_extinct, original_richness and
      cumulative_extinction_rate are NA because these groups are not
      fitted. The reference group has l = 1; groups excluded from
      fitting have l = NA. Group multipliers match those in
      `parameters$l`.

- `do_k`:

  NULL unless requested; otherwise the existing fixed-k versus free-k
  comparison table. Main outputs always use the fixed-k model.

- `status`, `converged`, `convergence`:

  Completion flags, optimizer diagnostics and convergence information.

- `bootstrap`:

  NULL unless requested; bootstrap settings, status, interval tables,
  diagnostics and compact replicate estimates.

- `M`, `times`:

  Validated input matrix and numeric time labels.

- `starts`:

  Starting values, with species entries in input order.

## Details

Active groups have at least one recorded species in the supplied
analysis matrix. Groups with only all-zero rows, or unused factor
levels, are removed from fitting. Each active group must additionally
have at least two species detected in two or more timesteps; otherwise
fitting stops.

Prepare the analysis matrix before calling this function. In particular,
aggregate annual records if necessary and select the analysis baseline
using
[`read_matrix()`](https://zhizhouhe2006.github.io/modgee/reference/read_matrix.md)
or other preparation code. This function does not aggregate, subset,
relabel or reorder time columns. It does not generate missing species
names or time labels.

Use numeric column names such as `"1870"` and `"1880"`, not `"1870s"` or
`"X1870"`. To use interval indices, explicitly assign column names such
as `"1"`, `"2"`, and so on. Invalid or duplicate time labels are
reported with their original column positions and names. Adjacent
spacings must differ from the first spacing by no more than an absolute
tolerance of 1e-8, as in
[`read_matrix()`](https://zhizhouhe2006.github.io/modgee/reference/read_matrix.md).
The optional preparation attribute is not required for fitting and is
not used instead of validating the matrix itself.

The current implementation supports 4 to 215 time intervals and requires
at least two species detected in at least two intervals each. These are
minimum computational requirements, not guarantees that the data contain
enough information for reliable estimation.

All-zero species are retained in the species table but excluded from
fitting. Their parameter and extinction estimates are `NA`. Species
detected in exactly one interval are treated as singletons using the
original procedure.

The reference species is moved to the first row internally. The public
species table and `parameters$d` retain the original input order.
Internal reordering does not change public row order.

Use [`set.seed()`](https://rdrr.io/r/base/Random.html) before fitting to
reproduce the random starting values. Optimization can still depend on
starting values and the numerical environment. A normal convergence code
does not establish a global optimum or validate the ecological
assumptions.

## Interpreting the estimates

The parameter `m` is not an extinction proportion. The corresponding
interval extinction probability is `1 - exp(-m)`; cumulative extinction
is `1 - exp(-cumsum(m))`. The first, second and final interval
extinction parameters are fixed at zero by the original implementation,
not independently estimated as zero. With groups, replace m by m \* l
for each group. Overall cumulative extinction is not calculated from the
reference m alone and is not an unweighted average of group curves.

`p_extinct` is a species-level probability of extinction by the final
interval, conditional on its detection history and the fitted
parameters.

Expected extinction times and their variances are conditional on
extinction. The variance is in squared interval units.
`estimated_extinction_year` maps interval indices onto the supplied time
labels; it is not an exact date of extinction. Conditional times may be
undefined (`NaN`) when the fitted probability of extinction is zero.

Dark extinctions represent losses without a recorded detection in the
analysis period. The estimate is not a count of supplied all-zero rows.
A negative estimate is not a literal negative number of species; it
indicates that the estimate needs interpretation in light of the data
and model assumptions.

## Bootstrap and recovery

Nonzero species rows are sampled globally with replacement. Every
replicate refits detection, extinction, group and singleton parameters
at fixed k = 0. The reference group and active group set remain fixed;
invalid samples are saved as failures. All-zero rows are not sampled.
Negative dark extinction estimates are retained. No species-level or
parameter intervals are computed.

Checkpoints are written below `getwd()/modgee_bootstrap`, in numbered
run directories. An unfinished matching run resumes in the same working
directory. After all requested replicates have been attempted, intervals
generated and final output printed, a completion marker is committed. A
subsequent call starts a new run, refitting the original model and all
requested replicates. Old runs are retained. Failed or nonconverged
replicates count as attempted; they do not prevent completion. Changing
bootstrap or CI after completion also starts a new run.

Sampling plans and random states are saved before fitting. An
interrupted replicate restarts from its saved sample and random state,
not its optimizer iteration. Within an unfinished run, saved failures
and nonconverged fits are not retried. A changed bootstrap count uses
the requested prefix of saved replicates or extends that unfinished run;
changed CI recomputes intervals. An unfinished run owns its saved seed.
For a fresh run, the caller's current random state is used; repeating
the same seed can reproduce the same numerical results despite fresh
computation.

Changing the input, model implementation or numerical dependency
versions starts a new task. Task locks cover run selection, fitting and
completion. Invalid plans or completion markers stop recovery; corrupt
result files are quarantined and refitted. Write failures leave a run
unfinished. Legacy flat checkpoints from older versions are preserved
but are not resumed by this run-directory format.

Bootstrap restores the caller's random state. Reproducibility assumes
the same numerical environment. Only fully converged replicates with
finite values for the requested metric enter its interval. Fewer than
two values yield NA endpoints; fewer than 100 successful replicates
produce a warning. Ten replicates are suitable only for demonstrations.

The bootstrap component contains settings, status, summary, curve,
group_results, diagnostics and replicates. Summary tables use metric,
estimate, lower, upper and n_valid; curve tables replace metric with
time. Group tables also contain group. estimate is always the original
point estimate. Six summary metrics receive intervals: known_extant,
known_extinct, dark_extinct, total_extinct, original_richness and
cumulative_extinction_rate. The compact replicate tables have a
replicate index and a value column. Diagnostics retain warnings, errors
and convergence codes for each attempt. The original do_k comparison is
cached separately and is never repeated for bootstrap samples.

## Usage Examples


    # Run or resume a bootstrap task in the current working directory.
    result_bootstrap <- fit_modgee(
      M = example_M, groups = example_groups,
      reference_group = "primary + old secondary", bootstrap = 10
    )
    result_bootstrap$bootstrap$status
    result_bootstrap$bootstrap$curve

    # Fit without groups using automatic starting values.
    set.seed(123)
    result <- fit_modgee(M = example_M)
    result$summary

    # Fit groups using the default reference group.
    set.seed(123)
    result_grouped <- fit_modgee(
      M = example_M,
      groups = example_groups
    )
    result_grouped$groups$reference_group
    result_grouped$group_results$summary

    # Select the reference group explicitly.
    set.seed(123)
    result_reference <- fit_modgee(
      M = example_M,
      groups = example_groups,
      reference_group = "primary + old secondary"
    )

    # Compare fixed-k and free-k models for grouped data.
    set.seed(123)
    result_grouped_k <- fit_modgee(
      M = example_M,
      groups = example_groups,
      reference_group = "primary + old secondary",
      do_k = TRUE
    )
    result_grouped_k$do_k

    # The k comparison is also available without groups.
    set.seed(123)
    result_k <- fit_modgee(M = example_M, do_k = TRUE)
    result_k$do_k

    # Supply named group labels instead of a two-column data frame.
    group_labels <- setNames(example_groups[[2]], example_groups[[1]])
    set.seed(123)
    result_labels <- fit_modgee(M = example_M, groups = group_labels)

## See also

**Functions:**
[`read_matrix`](https://zhizhouhe2006.github.io/modgee/reference/read_matrix.md),
[`plot_modgee`](https://zhizhouhe2006.github.io/modgee/reference/plot_modgee.md),
[`save_result`](https://zhizhouhe2006.github.io/modgee/reference/save_result.md)

**Examples:**
[`example_file`](https://zhizhouhe2006.github.io/modgee/reference/example_file.md),
[`example_M`](https://zhizhouhe2006.github.io/modgee/reference/example_M.md),
[`example_groups`](https://zhizhouhe2006.github.io/modgee/reference/example_groups.md),
[`example_result`](https://zhizhouhe2006.github.io/modgee/reference/example_result.md),
[`example_grouped_result`](https://zhizhouhe2006.github.io/modgee/reference/example_grouped_result.md),
[`example_bootstrap_result`](https://zhizhouhe2006.github.io/modgee/reference/example_bootstrap_result.md)

## Console messages

When bootstrap is NULL, each completed call reports convergence of the
main model's detection, extinction and singleton fits. Unsuccessful
stages are listed separately; singleton entries include the species
name. Missing convergence codes are reported as unverified. The
convergence table includes `species`, which is NA for the two main
stages. Code 999 records an optimizer error; a singleton may still have
a fallback estimate.

If do_k is TRUE, the free-k comparison status is reported separately. An
unsuccessful comparison may have been skipped or failed its validity
checks. It does not change the main model's convergence flag.

Negative overall or group dark-extinction estimates produce an
explanatory message and are retained unchanged. Undetected species
surviving at the final timestep are one possible explanation, not a
confirmed cause.

These summaries use [`message()`](https://rdrr.io/r/base/message.html)
and can be hidden with
[`suppressMessages()`](https://rdrr.io/r/base/message.html). Calls with
bootstrap enabled, including the original fit and all replicates, do not
emit per-fit convergence summaries. Existing optimizer warnings and
bootstrap diagnostics remain available. Calls that stop with an error do
not print a completion summary.

## Default display

The returned object has class `modgee_result` and retains all result
components. Printing shows `summary`, followed by
`group_results$summary` and `do_k` when present. Every completed call
explicitly prints these summaries, including calls assigned to a
variable. The complete object is returned invisibly to avoid duplicate
automatic printing. This display also applies once to completed
bootstrap results, not each replicate. Detection, extinction and
singleton-stage progress is printed as fitting proceeds, including the
original bootstrap fit and each newly fitted replicate. The
singleton-stage line is also printed when there are no singleton
species; in that case no singleton optimization is needed. Failed fits
only report stages reached. Reused checkpoint fits are not fitted again
and do not emit these stage lines. Internal do_k and do_l flags are not
printed. For non-bootstrap calls, the generic duplicate convergence
warning is suppressed. Other optimizer warnings are retained.
Convergence messages omit numerical status codes; the original codes
remain available in `result$convergence`.

## Singleton numerical errors

If a singleton log-likelihood is undefined during the initial-value
search, fitting stops with an error identifying the species and search
stage. Sparse detection histories or numerical instability are possible
causes, not confirmed diagnoses. The failed candidate is not silently
skipped.

## Interpretation messages

Negative dark-extinction messages clarify that a negative estimate does
not, by itself, invalidate the cumulative extinction rate or fitted
parameters. For do_k, the printed comparison is followed by an AIC
interpretation. A higher free-k AIC indicates no clear evidence of a
relationship; a lower AIC indicates evidence of a relationship. With
lower AIC, negative k indicates higher risk among less detectable
species and possible underestimation by the fixed-k model; positive k
indicates the reverse and possible overestimation. These are model-based
interpretations, not proof of independence or a quantified bias. Equal
AIC and unavailable comparisons receive separate messages. Main
estimates remain fixed-k. Bootstrap prints these interpretations once
for the final result, not each replicate.
