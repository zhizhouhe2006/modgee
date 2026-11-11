#' Fit the original MODGEE model
#'
#' Estimate detection parameters, cumulative extinction and dark
#' extinctions from species detection histories. This implementation
#' supports mutually exclusive species groups sharing a temporal extinction
#' baseline.
#'
#' @param M A numeric or logical matrix with species in rows
#'   and equally spaced time intervals in columns. Entries must be
#'   0 or 1, with no missing values. A value of 1 indicates at least
#'   one detection during an interval; 0 indicates no recorded
#'   detection, not confirmed absence or extinction.
#'   Row names are required and must identify species. Leading and trailing
#'   whitespace is removed before checking that names are nonempty and unique.
#'   Column names are required and must convert to finite, unique, strictly
#'   increasing, equally spaced numeric time labels. They are the only source
#'   of time labels; there is no separate \code{times} argument.
#'   Character and complex matrices are not accepted. Logical entries are
#'   converted to numeric 0/1 values.
#'
#' @param groups Optional CSV path, two-column data frame (species, group),
#'   character vector, factor, or numeric/logical 0/1 membership matrix.
#'   Names are trimmed and matched to all input species, including all-zero
#'   rows. Unnamed vectors and matrices follow input row order. Extra named
#'   species are ignored with a warning. Each species must have one group.
#'   CSV and character labels use radix sorting; factors use level order;
#'   membership matrices use column order. The first active group is the
#'   reference unless reference_group is specified. Empty groups are removed and reported. Every active group
#'   needs at least two repeated-detection species (a minimum safeguard,
#'   not a guarantee of statistical identifiability).
#'
#' @param reference_group Optional single group name. Leading and trailing
#'   whitespace is removed. When NULL, the first active group under the
#'   default ordering is the reference. An explicit name overrides factor
#'   levels and membership-matrix column order; other groups retain their
#'   relative order. The named group must contain recorded species and meet
#'   the repeated-detection requirement. Supplying a reference without groups
#'   is an error. The reference log multiplier is fixed at zero, so its
#'   hazard multiplier is one. The selected name is returned in
#'   \code{groups$reference_group}.
#'
#' @param do_k Logical scalar. If TRUE, compare the fixed-k model
#'   with a model that estimates a single k shared by all species.
#'   The two models share the first-stage detection estimates and use the
#'   same conditional likelihood in the extinction stage. For the alternative,
#'   the extinction hazard is \code{m[t] * d[i]^k * l[group[i]]}, with the
#'   group multiplier equal to 1 when groups are absent. The alternative
#'   starts at k = 0. Main estimates use k = 0.
#'
#' @param bootstrap NULL (default), or an integer of at least 3 giving the number
#'   of scheduled full-pipeline species bootstrap replicates. Failed and
#'   nonconverged replicates count toward this number and are not resampled.
#' @param CI Confidence level strictly between 0 and 1; default 0.95.
#'   Intervals are pointwise percentile intervals (quantile type 7).
#'
#' @details
#' Active groups have at least one recorded species in the supplied analysis
#' matrix. Groups with only all-zero rows, or unused factor levels, are
#' removed from fitting. Each active group must additionally have at least
#' two species detected in two or more timesteps; otherwise fitting stops.
#'
#' Prepare the analysis matrix before calling this function.
#' In particular, aggregate annual records if necessary and select
#' the analysis baseline using \code{read_matrix()} or other preparation code.
#' This function does not aggregate, subset, relabel or reorder time columns.
#' It does not generate missing species names or time labels.
#'
#' Use numeric column names such as \code{"1870"} and \code{"1880"}, not
#' \code{"1870s"} or \code{"X1870"}. To use interval indices, explicitly
#' assign column names such as \code{"1"}, \code{"2"}, and so on.
#' Invalid or duplicate time labels are reported with their original column
#' positions and names. Adjacent spacings must differ from the first spacing
#' by no more than an absolute tolerance of 1e-8, as in \code{read_matrix()}.
#' The optional preparation attribute is not required for fitting and is
#' not used instead of validating the matrix itself.
#'
#' The current implementation supports 4 to 215 time intervals
#' and requires at least two species detected in at least two
#' intervals each. These are minimum computational requirements,
#' not guarantees that the data contain enough information for
#' reliable estimation.
#'
#' All-zero species are retained in the species table but excluded
#' from fitting. Their parameter and extinction estimates are
#' \code{NA}. Species detected in exactly one interval are treated
#' as singletons using the original procedure.
#'
#' The reference species is moved to the first row internally.
#' The public species table and \code{parameters$d} retain the
#' original input order. Internal reordering does not change public row order.
#'
#' Use \code{set.seed()} before fitting
#' to reproduce the random starting values. Optimization can still
#' depend on starting values and the numerical environment.
#' A normal convergence code does not establish a global optimum
#' or validate the ecological assumptions.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{\code{summary}}{A named list with the following entries, in order:
#'     \describe{
#'       \item{\code{baseline_timestep}, \code{final_timestep}}{First and last numeric
#'         column labels of the supplied matrix (calendar years when year labels are supplied).}
#'       \item{\code{excluded_species}}{Number of all-zero rows in the supplied matrix.}
#'       \item{\code{known_species}}{Number of nonzero rows.}
#'       \item{\code{singleton_species}}{Number of species detected in exactly one timestep.}
#'       \item{\code{known_extant}}{Expected surviving recorded species:
#'         \code{known_species - known_extinct}. This is not the number detected in the final timestep.}
#'       \item{\code{known_extinct}}{Sum of extinction probabilities over recorded species.}
#'       \item{\code{dark_extinct}}{Estimated unrecorded extinctions; negative estimates are retained.
#'         This is not the number of all-zero rows.}
#'       \item{\code{total_extinct}}{\code{known_extinct + dark_extinct}.}
#'       \item{\code{original_richness}}{\code{known_species + dark_extinct}.}
#'       \item{\code{cumulative_extinction_rate}}{Final cumulative extinction proportion,
#'         on the 0--1 scale; equals \code{total_extinct / original_richness}.}
#'       \item{\code{AIC}}{Fixed-k main-model AIC, using the second-stage conditional
#'         log likelihood and the original two-stage parameter-count convention:
#'         \code{2 * p - 2 * logLik}, where
#'         \code{p = (R - 1) + T + (T - 3) + (G - 1)};
#'         logLik is the internal second-stage conditional log likelihood;
#'         R is the number of repeated-detection species, T the number of timesteps,
#'         and G the number of active groups (one without groups).
#'         Returns \code{NA} if either main fitting stage did not converge.
#'         When \code{do_k = TRUE}, matches the fixed-k comparison row.}
#'     }
#'   }
#'   \item{\code{parameters}}{Named vectors \code{d}, \code{h} and \code{m}.
#'     Grouped models also have \code{l}, the actual hazard multipliers, with
#'     reference value one. The starting values \code{starts$ls_start}
#'     remain log multipliers.}
#'   \item{\code{curve}}{Data frame: \code{time}, \code{h}, \code{m},
#'     \code{cumulative_extinction_rate}. With groups, m is the reference hazard,
#'     and the overall rate is weighted by inferred baseline richness.}
#'   \item{\code{species}}{Data frame in original input order, with columns:
#'     \code{species}, \code{belonged_group} (only with groups),
#'     \code{detected_timesteps}, \code{is_singleton},
#'     \code{first_detection_timestep}, \code{last_detection_timestep},
#'     \code{d}, \code{p_extinct}, \code{estimated_extinction_year}.
#'     Detection timesteps are counted, not listed. First and last detections use
#'     actual numeric column labels. Estimated extinction year is the conditional
#'     mean extinction timestep mapped onto those labels, not a confirmed extinction
#'     date. All-zero rows have zero detected timesteps, FALSE for is_singleton,
#'     and NA for detection years and estimates. Undefined extinction times are NA.}
#'   \item{\code{groups}}{NULL without groups; otherwise matching information,
#'     \code{species_group}, membership matrix \code{H}, active \code{levels},
#'     \code{reference_group}, \code{removed_groups}, and \code{matching}.}
#'   \item{\code{group_results}}{NULL without groups; otherwise a list of two data frames:
#'     \describe{
#'       \item{\code{curve}}{\code{group}, \code{time}, \code{m},
#'         \code{cumulative_extinction_rate}, for active groups only.
#'         The hazard m is the group's own hazard.}
#'       \item{\code{summary}}{\code{group}, then \code{l} (the actual hazard multiplier, not its logarithm), followed by the same columns as the
#'         overall summary, except AIC. Counts and estimates refer to each group.
#'         Supplied groups containing only all-zero species are included for counting;
#'         their dark_extinct, total_extinct, original_richness and
#'         cumulative_extinction_rate are NA because these groups are not fitted.
#'         The reference group has l = 1; groups excluded from fitting have l = NA.
#'         Group multipliers match those in \code{parameters$l}.}
#'     }
#'   }
#'   \item{\code{do_k}}{NULL unless requested; otherwise the existing fixed-k versus
#'     free-k comparison table. Main outputs always use the fixed-k model.}
#'   \item{\code{status}, \code{converged}, \code{convergence}}{Completion flags,
#'     optimizer diagnostics and convergence information.}
#'   \item{\code{bootstrap}}{NULL unless requested; bootstrap settings, status,
#'     interval tables, diagnostics and compact replicate estimates.}
#'   \item{\code{M}, \code{times}}{Validated input matrix and numeric time labels.}
#'   \item{\code{starts}}{Starting values, with species entries in input order.}
#' }
#'
#' @section Interpreting the estimates:
#' The parameter \code{m} is not an extinction proportion.
#' The corresponding interval extinction probability is
#' \code{1 - exp(-m)}; cumulative extinction is
#' \code{1 - exp(-cumsum(m))}. The first, second and final
#' interval extinction parameters are fixed at zero by the
#' original implementation, not independently estimated as zero.
#' With groups, replace m by m * l for each group. Overall cumulative
#' extinction is not calculated from the reference m alone and is not an
#' unweighted average of group curves.
#'
#' \code{p_extinct} is a species-level probability of extinction
#' by the final interval, conditional on its detection history
#' and the fitted parameters.
#'
#' Expected extinction times and their variances are conditional
#' on extinction. The variance is in squared interval units.
#' \code{estimated_extinction_year} maps interval indices onto
#' the supplied time labels; it is not an exact date of extinction.
#' Conditional times may be undefined (\code{NaN}) when the
#' fitted probability of extinction is zero.
#'
#' Dark extinctions represent losses without a recorded detection
#' in the analysis period. The estimate is not a count of supplied
#' all-zero rows. A negative estimate is not a literal negative
#' number of species; it indicates that the estimate needs
#' interpretation in light of the data and model assumptions.
#'
#' @section Console messages:
#' When bootstrap is NULL, each completed call reports convergence of the main
#' model's detection, extinction and singleton fits. Unsuccessful stages are
#' listed separately; singleton entries include the species name. Missing
#' convergence codes are reported as unverified. The convergence table includes
#' \code{species}, which is NA for the two main stages. Code 999 records an
#' optimizer error; a singleton may still have a fallback estimate.
#'
#' If do_k is TRUE, the free-k comparison status is reported separately. An
#' unsuccessful comparison may have been skipped or failed its validity checks.
#' It does not change the main model's convergence flag.
#'
#' Negative overall or group dark-extinction estimates produce an explanatory
#' message and are retained unchanged. Undetected species surviving at the final
#' timestep are one possible explanation, not a confirmed cause.
#'
#' These summaries use \code{message()} and can be hidden with
#' \code{suppressMessages()}. Calls with bootstrap enabled, including the original
#' fit and all replicates, do not emit per-fit convergence summaries. Existing
#' optimizer warnings and bootstrap diagnostics remain available. Calls that
#' stop with an error do not print a completion summary.
#'
#' @section Bootstrap and recovery:
#' Nonzero species rows are sampled globally with replacement. Every replicate
#' refits detection, extinction, group and singleton parameters at fixed k = 0.
#' The reference group and active group set remain fixed; invalid samples are
#' saved as failures. All-zero rows are not sampled. Negative dark extinction
#' estimates are retained. No species-level or parameter intervals are computed.
#'
#' Checkpoints are written below \code{getwd()/modgee_bootstrap}, in numbered
#' run directories. An unfinished matching run resumes in the same working directory.
#' After all requested replicates have been attempted, intervals generated and final
#' output printed, a completion marker is committed. A subsequent call starts a new
#' run, refitting the original model and all requested replicates. Old runs are retained.
#' Failed or nonconverged replicates count as attempted; they do not prevent completion.
#' Changing bootstrap or CI after completion also starts a new run.
#'
#' Sampling plans and random states are saved before fitting. An interrupted replicate
#' restarts from its saved sample and random state, not its optimizer iteration.
#' Within an unfinished run, saved failures and nonconverged fits are not retried.
#' A changed bootstrap count uses the requested prefix of saved replicates or extends
#' that unfinished run; changed CI recomputes intervals. An unfinished run owns its
#' saved seed. For a fresh run, the caller's current random state is used; repeating
#' the same seed can reproduce the same numerical results despite fresh computation.
#'
#' Changing the input, model implementation or numerical dependency versions starts
#' a new task. Task locks cover run selection, fitting and completion. Invalid plans
#' or completion markers stop recovery; corrupt result files are quarantined and
#' refitted. Write failures leave a run unfinished. Legacy flat checkpoints from
#' older versions are preserved but are not resumed by this run-directory format.
#'
#' Bootstrap restores the caller's random state. Reproducibility assumes the
#' same numerical environment. Only fully converged replicates with finite
#' values for the requested metric enter its interval. Fewer than two values
#' yield NA endpoints; fewer than 100 successful replicates produce a warning.
#' Ten replicates are suitable only for demonstrations.
#'
#' The bootstrap component contains settings, status, summary, curve,
#' group_results, diagnostics and replicates. Summary tables use metric,
#' estimate, lower, upper and n_valid; curve tables replace metric with time.
#' Group tables also contain group. estimate is always the original point
#' estimate. Six summary metrics receive intervals: known_extant, known_extinct,
#' dark_extinct, total_extinct, original_richness and cumulative_extinction_rate.
#' The compact replicate tables have a replicate index and a value column.
#' Diagnostics retain warnings, errors and convergence codes for each attempt.
#' The original do_k comparison is cached separately and is never repeated
#' for bootstrap samples.
#'
#' @section Default display:
#' The returned object has class \code{modgee_result} and retains all result
#' components. Printing shows \code{summary}, followed by
#' \code{group_results$summary} and \code{do_k} when present. Every completed
#' call explicitly prints these summaries, including calls assigned to a variable.
#' The complete object is returned invisibly to avoid duplicate automatic printing.
#' This display also applies once to completed bootstrap results, not each replicate.
#' Detection, extinction and singleton-stage progress is printed as fitting
#' proceeds, including the original bootstrap fit and each newly fitted replicate.
#' The singleton-stage line is also printed when there are no singleton species;
#' in that case no singleton optimization is needed. Failed fits only report
#' stages reached. Reused checkpoint fits are not fitted again and do not emit
#' these stage lines. Internal do_k and do_l flags are not printed.
#' For non-bootstrap calls, the generic duplicate convergence warning is
#' suppressed. Other optimizer warnings are retained.
#' Convergence messages omit numerical status codes; the original codes remain
#' available in \code{result$convergence}.
#'
#' @section Singleton numerical errors:
#' If a singleton log-likelihood is undefined during the initial-value search,
#' fitting stops with an error identifying the species and search stage.
#' Sparse detection histories or numerical instability are possible causes,
#' not confirmed diagnoses. The failed candidate is not silently skipped.
#'
#' @section Interpretation messages:
#' Negative dark-extinction messages clarify that a negative estimate does not,
#' by itself, invalidate the cumulative extinction rate or fitted parameters.
#' For do_k, the printed comparison is followed by an AIC interpretation. A higher
#' free-k AIC indicates no clear evidence of a relationship; a lower AIC indicates
#' evidence of a relationship. With lower AIC, negative k indicates higher risk
#' among less detectable species and possible underestimation by the fixed-k model;
#' positive k indicates the reverse and possible overestimation. These are model-based
#' interpretations, not proof of independence or a quantified bias. Equal AIC and
#' unavailable comparisons receive separate messages. Main estimates remain fixed-k.
#' Bootstrap prints these interpretations once for the final result, not each replicate.
#' @export
#' @seealso
#' \strong{Functions:}
#' \code{\link{read_matrix}}, \code{\link{plot_modgee}}, \code{\link{save_result}}
#'
#' \strong{Examples:}
#' \code{\link{example_file}}, \code{\link{example_M}}, \code{\link{example_groups}}, \code{\link{example_result}}, \code{\link{example_grouped_result}}, \code{\link{example_bootstrap_result}}
#' @section Usage Examples:
#' \preformatted{
#' # Run or resume a bootstrap task in the current working directory.
#' result_bootstrap <- fit_modgee(
#'   M = example_M, groups = example_groups,
#'   reference_group = "primary + old secondary", bootstrap = 10
#' )
#' result_bootstrap$bootstrap$status
#' result_bootstrap$bootstrap$curve
#'
#' # Fit without groups using automatic starting values.
#' set.seed(123)
#' result <- fit_modgee(M = example_M)
#' result$summary
#'
#' # Fit groups using the default reference group.
#' set.seed(123)
#' result_grouped <- fit_modgee(
#'   M = example_M,
#'   groups = example_groups
#' )
#' result_grouped$groups$reference_group
#' result_grouped$group_results$summary
#'
#' # Select the reference group explicitly.
#' set.seed(123)
#' result_reference <- fit_modgee(
#'   M = example_M,
#'   groups = example_groups,
#'   reference_group = "primary + old secondary"
#' )
#'
#' # Compare fixed-k and free-k models for grouped data.
#' set.seed(123)
#' result_grouped_k <- fit_modgee(
#'   M = example_M,
#'   groups = example_groups,
#'   reference_group = "primary + old secondary",
#'   do_k = TRUE
#' )
#' result_grouped_k$do_k
#'
#' # The k comparison is also available without groups.
#' set.seed(123)
#' result_k <- fit_modgee(M = example_M, do_k = TRUE)
#' result_k$do_k
#'
#' # Supply named group labels instead of a two-column data frame.
#' group_labels <- setNames(example_groups[[2]], example_groups[[1]])
#' set.seed(123)
#' result_labels <- fit_modgee(M = example_M, groups = group_labels)
#' }
fit_modgee <- function(M, groups = NULL, reference_group = NULL, do_k = FALSE,
                       bootstrap = NULL, CI = 0.95) {
  if (missing(M)) {
    stop("Provide a detection matrix with species row names and numeric time column names.",
         call. = FALSE)
  }
  .bootstrap_validate(bootstrap, CI)
  if (!is.null(bootstrap)) {
    result <- .bootstrap_modgee(M, groups, reference_group, do_k,
      as.integer(bootstrap), CI, .finalize = function(result) {
        class(result) <- c("modgee_result", "list")
        .modgee_report_dark(result)
        print.modgee_result(result)
        result
      })
    return(invisible(result))
  }
  result <- .modgee_quiet_fit(.fit_modgee_impl(
    M = M, groups = groups, reference_group = reference_group, do_k = do_k))
  result["bootstrap"] <- list(NULL)
  .modgee_report_fit(result)
  class(result) <- c("modgee_result", "list")
  print.modgee_result(result)
  invisible(result)
}

.fit_modgee_impl <- function(M, starts = NULL, groups = NULL,
                       reference_group = NULL, do_k = FALSE, .capture = NULL) {
  if (!is.logical(do_k) || length(do_k) != 1L ||
      !is.null(dim(do_k)) || is.na(do_k)) {
    stop("do_k must be TRUE or FALSE.", call. = FALSE)
  }
  k_start <- if (do_k && is.list(starts)) starts[["k_start"]] else NULL
  if (!is.null(k_start) &&
      (!(is.double(k_start) || is.integer(k_start)) ||
       is.object(k_start) || !is.null(dim(k_start)) ||
       length(k_start) != 1L || !is.finite(k_start))) {
    stop("starts$k_start must be a finite numeric scalar.", call. = FALSE)
  }
  if (missing(M)) {
    stop("Provide a detection matrix with species row names and numeric time column names.",
         call. = FALSE)
  }
  input <- .validate_modgee_matrix(M)
  M <- input$M
  times <- input$times
  species_names <- rownames(M)
  n <- nrow(M)
  TT <- ncol(M)
  counts <- rowSums(M)
  spacing <- diff(times)
  group_info <- .normalize_modgee_groups(groups, species_names, counts, reference_group)

  recorded <- which(counts > 0)
  repeated <- which(counts >= 2)

  if (length(repeated) < 2L) {
    stop("At least two species must be detected in two or more intervals.")
  }

  # Put a repeated-detection species first; preserve all other order.
  anchor <- repeated[1]
  fit_rows <- c(anchor, recorded[recorded != anchor])
  fit_M <- M[fit_rows, , drop = FALSE]

  if (is.null(starts)) {
    initial <- est_detect_parms(fit_M)
    fit_starts <- list(
      ds_start = initial$d_est,
      hs_start = initial$h_est,
      ms_start = c(0, 0, stats::runif(TT - 2, 0, 1))
    )

    # Store starting values in the original input species order.
    input_d <- rep(NA_real_, n)
    input_d[fit_rows] <- fit_starts$ds_start
    starts <- list(
      ds_start = input_d,
      hs_start = fit_starts$hs_start,
      ms_start = fit_starts$ms_start
    )
  } else {
    required <- c("ds_start", "hs_start", "ms_start")

    if (!is.list(starts) || !all(required %in% names(starts))) {
      stop("starts must contain ds_start, hs_start and ms_start.")
    }
    group_start <- starts$ls_start
    starts <- starts[required]

    expected_lengths <- c(n, TT, TT)
    for (i in seq_along(required)) {
      if (
        !(is.double(starts[[i]]) || is.integer(starts[[i]])) ||
        !is.null(dim(starts[[i]])) ||
        length(starts[[i]]) != expected_lengths[i]
      ) {
        stop("Invalid starting-value vector: ", required[i])
      }
    }

    if (
      !is.finite(starts$ds_start[anchor]) ||
      abs(starts$ds_start[anchor] - 1) > 1e-8
    ) {
      stop(
        "ds_start must be 1 for the reference species: ",
        species_names[anchor]
      )
    }

    fit_starts <- starts
    fit_starts$ds_start <- starts$ds_start[fit_rows]
  }

  used <- c(
    fit_starts$ds_start[rowSums(fit_M) >= 2],
    fit_starts$hs_start,
    fit_starts$ms_start[3:(TT - 1)]
  )
  if (any(!is.finite(used)) || any(used <= 0)) {
    stop(
      "Starting values used by the optimizer must be positive and finite. ",
      "The data may be too sparse for the original initialization."
    )
  }

  fit_H <- NULL
  ls_start <- NULL
  if (!is.null(group_info)) {
    fit_H <- group_info$H[fit_rows, , drop = FALSE]
    ls_start <- .modgee_group_starts(
      if (exists("group_start", inherits = FALSE)) group_start else NULL,
      group_info$levels
    )
    starts$ls_start <- ls_start
  }
  # A single group is the standard model with an explicit label.
  grouped <- !is.null(group_info) && length(group_info$levels) > 1L
  fitter <- if (grouped) .est_modgee_grouped else est_detect_extinct_parms_ML
  fit <- fitter(
    M = fit_M,
    H = if (grouped) fit_H else NULL,
    ds_start = fit_starts$ds_start,
    hs_start = fit_starts$hs_start,
    ms_start = fit_starts$ms_start,
    k_start = NULL,
    ls_start = if (grouped) ls_start else NULL
  )

  if (is.null(fit)) {
    stop("MODGEE fitting failed; no valid fit was returned.")
  }
  if (any(!is.finite(c(
    fit$ds_est, fit$hs_est, fit$ms_est, fit$ls_est, fit$LL
  )))) {
    stop("The fit contains non-finite parameter estimates or likelihood.")
  }

  if (!is.null(group_info)) {
    fit$H <- fit_H
    if (!grouped) fit$ls_est <- 0
    names(fit$ls_est) <- group_info$levels
  }
  statistics <- if (is.null(group_info)) compute_summary_stats(fit) else
    .summary_modgee_grouped(fit, times, group_info)

  if (any(!is.finite(c(
    statistics$Mus_overall,
    statistics$p_exts,
    statistics$X_T_hat
  )))) {
    stop("The extinction summary contains non-finite results.")
  }

  singleton_names <- rownames(fit_M)[rowSums(fit_M) == 1L]
  convergence <- do.call(rbind, lapply(
    seq_along(fit$optimizer_runs), function(i) {
      x <- fit$optimizer_runs[[i]]
      data.frame(
        call = i,
        stage = if (i == 1L) "detection" else
          if (i == 2L) "extinction" else "singleton",
        species = if (i > 2L) singleton_names[i - 2L] else NA_character_,
        method = rownames(x),
        convcode = x$convcode,
        objective = x$value,
        row.names = NULL
      )
    }
  ))

  converged <- all(
    !is.na(convergence$convcode) & convergence$convcode == 0
  )

  curve <- data.frame(
    time = times,
    h = as.numeric(fit$hs_est),
    m = as.numeric(fit$ms_est),
    cumulative_extinction = as.numeric(statistics$Mus_overall)
  )

  # Restore the original input order and retain all-zero species.
  species <- data.frame(
    species = species_names,
    detected_intervals = counts,
    included_in_fit = counts > 0,
    exclusion_reason = ifelse(
      counts == 0, "No M in supplied intervals", NA_character_
    ),
    d = rep(NA_real_, n),
    p_extinct = rep(NA_real_, n),
    expected_extinction_interval = rep(NA_real_, n),
    extinction_interval_variance = rep(NA_real_, n),
    expected_extinction_time = rep(NA_real_, n),
    row.names = NULL
  )

  species$d[fit_rows] <- as.numeric(fit$ds_est)
  species$p_extinct[fit_rows] <- as.numeric(statistics$p_exts)
  species$expected_extinction_interval[fit_rows] <-
    as.numeric(statistics$E_t_exts)
  species$extinction_interval_variance[fit_rows] <-
    as.numeric(statistics$var_t_exts)
  species$expected_extinction_time[fit_rows] <-
    times[1] + (as.numeric(statistics$E_t_exts) - 1) * spacing[1]

  parameters <- list(
    d = stats::setNames(species$d, species_names),
    h = stats::setNames(as.numeric(fit$hs_est), times),
    m = stats::setNames(as.numeric(fit$ms_est), times)
  )

  if (!is.null(group_info)) {
    parameters$l <- exp(fit$ls_est)
    species$group <- unname(group_info$species_group)
  }

  comparison <- if (do_k) {
    .modgee_compare_k(fit, grouped, k_start)
  } else {
    NULL
  }
  if (!is.null(k_start)) starts$k_start <- k_start
  public <- .modgee_public_output(M, times, fit, statistics, species, group_info)

  if (!is.null(.capture)) .capture(fit, grouped)

  list(
    summary = public$summary,
    parameters = parameters,
    groups = group_info,
    group_results = public$group_results,
    curve = public$curve,
    species = public$species,
    do_k = comparison,
    status = list(
      completed = TRUE,
      converged = converged,
      message = if (converged) {
        "Optimization completed with normal convergence codes."
      } else {
        "Results returned, but optimizer convergence requires review."
      },
      input_species = n,
      fitted_species = length(fit_rows),
      excluded_zero_species = sum(counts == 0)
    ),
    M = M,
    times = times,
    starts = starts,
    convergence = convergence,
    converged = converged
  )
}


.validate_modgee_matrix <- function(M) {
  if (!is.matrix(M) ||
      !(is.double(M) || is.integer(M) || is.logical(M))) {
    stop("M must be a real numeric or logical matrix.", call. = FALSE)
  }
  if (nrow(M) == 0L || ncol(M) < 4L || ncol(M) > 215L) {
    stop("Provide at least one species and 4 to 215 time intervals.",
         call. = FALSE)
  }

  species_names <- rownames(M)
  if (is.null(species_names)) {
    stop("Species row names are required; automatic names are not generated.",
         call. = FALSE)
  }
  invalid_species <- which(is.na(species_names) |
                             !nzchar(trimws(species_names)))
  if (length(invalid_species)) {
    stop("Species names must be nonempty and not missing. Invalid rows: ",
         paste(invalid_species, collapse = ", "), ".", call. = FALSE)
  }
  species_names <- trimws(species_names)
  duplicate_species <- unique(species_names[duplicated(species_names)])
  if (length(duplicate_species)) {
    stop("Species names must be unique after trimming whitespace. Duplicates: ",
         paste(encodeString(duplicate_species, quote = '"'), collapse = ", "),
         ".", call. = FALSE)
  }

  time_names <- colnames(M)
  if (is.null(time_names)) {
    stop("Time column names are required; interval indices are not generated.",
         call. = FALSE)
  }
  label_details <- function(index) {
    labels <- encodeString(time_names[index], quote = '"')
    labels[is.na(time_names[index])] <- "<NA>"
    paste0("time column ", index, " (", labels, ")", collapse = "; ")
  }
  missing_labels <- which(is.na(time_names) | !nzchar(trimws(time_names)))
  if (length(missing_labels)) {
    stop("Time column names must be present and non-empty. Invalid: ",
         label_details(missing_labels), ".", call. = FALSE)
  }
  times <- suppressWarnings(as.numeric(time_names))
  invalid_labels <- which(!is.finite(times))
  if (length(invalid_labels)) {
    stop("Time column names must be finite numbers. Invalid: ",
         label_details(invalid_labels), ".", call. = FALSE)
  }
  duplicate_years <- unique(times[duplicated(times)])
  if (length(duplicate_years)) {
    details <- vapply(duplicate_years, function(year) {
      paste0("year ", as.character(year), ": ",
             label_details(which(times == year)))
    }, character(1))
    stop("Time column names must be unique. Duplicate years: ",
         paste(details, collapse = "; "), ".", call. = FALSE)
  }
  spacing <- diff(times)
  reversed <- which(spacing <= 0)
  if (length(reversed)) {
    stop("Time column names must be strictly increasing. Out-of-order pairs: ",
         paste0(reversed, " -> ", reversed + 1L, collapse = ", "),
         ". Prepare the time order before fitting.", call. = FALSE)
  }
  # Use the same finite-value check and absolute tolerance as read_matrix().
  if (!all(.matrix_near(spacing, spacing[1L]))) {
    stop("Time intervals must be finite and equally spaced (absolute tolerance 1e-8).",
         call. = FALSE)
  }

  invalid_values <- is.na(M) | !(M == 0 | M == 1)
  if (any(invalid_values)) {
    position <- which(invalid_values, arr.ind = TRUE)[1, ]
    stop("Detection must be 0 or 1: species '", species_names[position[1]],
         "', time '", time_names[position[2]], "'. Missing values are not allowed.",
         call. = FALSE)
  }
  storage.mode(M) <- "double"
  rownames(M) <- species_names
  colnames(M) <- as.character(times)
  list(M = M, times = times)
}
