#' Read and prepare a species detection matrix
#'
#' @param data A CSV file path or a numeric or logical detection matrix.
#'   In a CSV, the first column contains species names and remaining column
#'   names contain starting years. In a matrix, species names must be row
#'   names and starting years must be column names; do not include a species
#'   column inside the matrix. Detection values must be 0/1 or, for a logical
#'   matrix, FALSE/TRUE. Missing values are not allowed.
#' @param timestep Positive output interval width, in years. Defaults to 1.
#' @param time_range Optional numeric vector of length two giving inclusive
#'   limits on output interval labels.
#' @param baseline Optional starting interval label within the selected range.
#' @param baseline_min_species Optional positive integer. Select the first
#'   interval in the selected range with at least this many detected species.
#'
#' @details
#' At least two time columns are required. Numeric column names represent the
#' starting years of consecutive, equally spaced observation intervals.
#' The input interval width is inferred from the spacing between column names.
#' Missing observation intervals must not be omitted from the input. Equally
#' spaced labels are interpreted as consecutive intervals, not isolated years.
#'
#' The output width must be an integer multiple of the input width.
#' Input intervals must fit within calendar-aligned output intervals without
#' crossing their boundaries. Output labels follow
#' \code{floor(time / timestep) * timestep}.
#' For input labels 1900, 1910, 1920, use \code{timestep = 10} to retain
#' the original interval width. Such data cannot be split into annual records.
#'
#' Processing proceeds in this order: aggregate intervals, apply
#' \code{time_range}, then select the baseline. Aggregation records whether
#' a species was detected at least once within each output interval.
#' The upper bound of \code{time_range} refers to an interval label, not its
#' ending year. If the requested range extends beyond the available output
#' labels but overlaps them, the overlap is retained with an informational
#' message reporting the retained label range before baseline selection.
#' A range with no overlapping labels is an error.
#'
#' Invalid time labels are reported with their original names and positions.
#' Duplicate years are compared numerically, so 1900 and 1900.0 are duplicates.
#' Time-column positions exclude the species-name column in a CSV.
#'
#' \code{baseline} and \code{baseline_min_species} are mutually exclusive.
#' If neither is supplied, the first selected interval is used.
#' All-zero species rows and the original species order are retained.
#' Incomplete boundary intervals are retained with an informational message.
#' For annual input ending in 2026 and \code{timestep = 10}, the column
#' labelled 2020 summarizes only 2020--2026. Unobserved years are not added.
#'
#' Existing custom attributes on matrix input are not carried forward.
#' Preparation information is recalculated from the supplied time columns.
#' For already aggregated data, the original coverage of a partial interval
#' cannot be recovered from its starting-year label alone; each input column
#' is interpreted as covering one full inferred input interval.
#'
#' @return A numeric 0/1 matrix with species row names and time column names.
#'   A \code{preparation} attribute records the input type, input interval width,
#'   output interval width, retained labels, baseline selection and coverage.
#'   Coverage endpoints use an exclusive upper bound.
#'
#' @export
#' @seealso
#' \strong{Functions:}
#' \code{\link{fit_modgee}}, \code{\link{plot_modgee}}, \code{\link{save_result}}
#'
#' \strong{Examples:}
#' \code{\link{example_file}}, \code{\link{example_M}}, \code{\link{example_groups}}, \code{\link{example_result}}, \code{\link{example_grouped_result}}, \code{\link{example_bootstrap_result}}
#' @section Usage Examples:
#' \preformatted{
#' # Read annual records from the installed CSV.
#' M_annual <- read_matrix(data = example_file)
#'
#' # Aggregate into ten-year intervals and select the baseline.
#' M <- read_matrix(
#'   data = example_file,
#'   timestep = 10,
#'   baseline_min_species = 10
#' )
#' dim(M)
#' attr(M, "preparation")
#'
#' # Prepare an existing matrix instead of a CSV file.
#' M_decadal <- read_matrix(data = M_annual, timestep = 10)
#'
#' # Keep a specified range of interval labels.
#' M_range <- read_matrix(
#'   data = example_file, timestep = 10,
#'   time_range = c(1900, 2000)
#' )
#'
#' # Specify the baseline directly.
#' M_baseline <- read_matrix(
#'   data = example_file, timestep = 10,
#'   baseline = 1900
#' )
#'
#' # Apply a detection threshold within a selected time range.
#' M_selected <- read_matrix(
#'   data = example_file, timestep = 10,
#'   time_range = c(1900, 2020),
#'   baseline_min_species = 10
#' )
#' }
read_matrix <- function(data, timestep = 1, time_range = NULL,
                        baseline = NULL, baseline_min_species = NULL) {
  scalar_number <- function(x) {
    (is.double(x) || is.integer(x)) &&
      length(x) == 1L && is.finite(x)
  }

  if (missing(data)) {
    stop("'data' must be a CSV file path or a detection matrix.", call. = FALSE)
  }
  if (!scalar_number(timestep) || timestep <= 0) {
    stop("'timestep' must be a positive number.", call. = FALSE)
  }
  if (!is.null(time_range)) {
    if (!(is.double(time_range) || is.integer(time_range)) ||
        length(time_range) != 2L || any(!is.finite(time_range)) ||
        time_range[1] > time_range[2]) {
      stop("'time_range' must contain two finite, ordered numbers.",
           call. = FALSE)
    }
  }
  if (!is.null(baseline) && !scalar_number(baseline)) {
    stop("'baseline' must be a finite number.", call. = FALSE)
  }
  if (!is.null(baseline_min_species)) {
    if (!scalar_number(baseline_min_species) ||
        baseline_min_species < 1 ||
        baseline_min_species != floor(baseline_min_species)) {
      stop("'baseline_min_species' must be a positive integer.",
           call. = FALSE)
    }
  }
  if (!is.null(baseline) && !is.null(baseline_min_species)) {
    stop("Use either 'baseline' or 'baseline_min_species', not both.",
         call. = FALSE)
  }

  input <- .read_detection_input(data)
  times <- input$times
  M <- input$M
  input_width <- times[2] - times[1]

  if (!all(.matrix_near(diff(times), input_width))) {
    stop("Time columns must be equally spaced; missing intervals are not allowed.",
         call. = FALSE)
  }
  ratio <- timestep / input_width
  if (!is.finite(ratio) || ratio < 1 ||
      !.matrix_near(ratio, round(ratio))) {
    stop("'timestep' must be an integer multiple of the input interval width.",
         call. = FALSE)
  }
  if (!all(.matrix_near(times / input_width,
                        round(times / input_width)))) {
    stop("Input intervals are not aligned with calendar-based output intervals.",
         call. = FALSE)
  }

  scaled <- times / timestep
  near_integer <- .matrix_near(scaled, round(scaled))
  scaled[near_integer] <- round(scaled[near_integer])
  bins <- floor(scaled)
  unique_bins <- unique(bins)
  labels <- unique_bins * timestep
  if (length(unique_bins) > 1L && any(diff(unique_bins) != 1)) {
    stop("An internal output interval is missing.", call. = FALSE)
  }

  aggregated <- matrix(
    0, nrow = nrow(M), ncol = length(labels),
    dimnames = list(rownames(M), as.character(labels))
  )
  coverage <- data.frame(
    label = labels,
    observed_start = numeric(length(labels)),
    observed_end_exclusive = numeric(length(labels)),
    complete = logical(length(labels))
  )
  for (j in seq_along(labels)) {
    index <- which(bins == unique_bins[j])
    aggregated[, j] <- as.numeric(rowSums(M[, index, drop = FALSE]) > 0)
    coverage$observed_start[j] <- min(times[index])
    coverage$observed_end_exclusive[j] <- max(times[index]) + input_width
    coverage$complete[j] <-
      .matrix_near(coverage$observed_start[j], labels[j]) &&
      .matrix_near(coverage$observed_end_exclusive[j], labels[j] + timestep)
  }

  keep <- rep(TRUE, length(labels))
  if (!is.null(time_range)) {
    keep <- labels >= time_range[1] & labels <= time_range[2]
  }
  if (!any(keep)) {
    stop("No intervals fall within 'time_range'.", call. = FALSE)
  }
  if (!is.null(time_range) &&
      (time_range[1] < min(labels) || time_range[2] > max(labels))) {
    message(
      "Requested time_range [", time_range[1], ", ", time_range[2],
      "] extends beyond available interval labels [",
      min(labels), ", ", max(labels),
      "]. Retained interval labels before baseline selection: [",
      min(labels[keep]), ", ", max(labels[keep]), "]."
    )
  }
  M <- aggregated[, keep, drop = FALSE]
  labels <- labels[keep]
  coverage <- coverage[keep, , drop = FALSE]

  start <- .select_matrix_baseline(
    M, labels, baseline, baseline_min_species
  )
  keep <- seq.int(start, ncol(M))
  M <- M[, keep, drop = FALSE]
  labels <- labels[keep]
  coverage <- coverage[keep, , drop = FALSE]
  rownames(coverage) <- NULL

  if (any(!coverage$complete)) {
    message(
      "Incomplete boundary intervals retained: ",
      paste(coverage$label[!coverage$complete], collapse = ", "),
      ". Detections use only the observed time coverage."
    )
  }
  if (ncol(M) < 4L) {
    message("Fewer than four intervals remain; check fit_modgee() requirements.")
  }
  if (sum(rowSums(M) >= 2) < 2L) {
    message("Fewer than two species have detections in multiple intervals.")
  }

  baseline_method <- if (!is.null(baseline)) {
    "specified"
  } else if (!is.null(baseline_min_species)) {
    "minimum_species"
  } else {
    "first_interval"
  }
  attr(M, "preparation") <- list(
    input_type = input$type,
    input_timestep = input_width,
    timestep = timestep,
    time_range = time_range,
    retained_labels = labels,
    baseline_method = baseline_method,
    requested_baseline = baseline,
    baseline_min_species = baseline_min_species,
    baseline = labels[1],
    coverage = coverage
  )
  M
}


.read_detection_input <- function(data) {
  if (is.matrix(data)) {
    if (!(is.double(data) || is.integer(data) || is.logical(data))) {
      stop("Matrix input must be numeric or logical, with 0/1 values.",
           call. = FALSE)
    }
    species <- rownames(data)
    time_names <- colnames(data)
    M <- data
    input_type <- "matrix"
  } else if (is.character(data) && length(data) == 1L && !is.na(data)) {
    if (!file.exists(data) || dir.exists(data)) {
      stop("'data' must refer to an existing CSV file.", call. = FALSE)
    }
    x <- utils::read.csv(
      data, check.names = FALSE, colClasses = "character",
      row.names = NULL, na.strings = c("", "NA"),
      strip.white = TRUE, fill = FALSE
    )
    if (nrow(x) == 0L || ncol(x) < 3L) {
      stop("The CSV must contain species rows and at least two time columns.",
           call. = FALSE)
    }
    species <- x[[1]]
    time_names <- names(x)[-1]
    M <- matrix(
      suppressWarnings(as.numeric(as.matrix(x[-1]))),
      nrow = nrow(x), ncol = ncol(x) - 1L
    )
    input_type <- "csv"
  } else {
    stop("'data' must be a CSV file path or a detection matrix.", call. = FALSE)
  }

  if (nrow(M) == 0L || ncol(M) < 2L) {
    stop("Input must contain species rows and at least two time columns.",
         call. = FALSE)
  }
  if (is.null(species) || anyNA(species) || any(!nzchar(trimws(species)))) {
    stop("Species names must not be missing or empty; matrices need row names.",
         call. = FALSE)
  }
  species <- trimws(species)
  if (anyDuplicated(species)) {
    stop("Duplicate species name: ", species[duplicated(species)][1],
         call. = FALSE)
  }
  if (is.null(time_names)) {
    stop("Time column names must be present and non-empty: all column names are missing.",
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

  # Rebuild a plain numeric matrix without inherited custom attributes.
  M <- matrix(
    as.numeric(M), nrow = nrow(M), ncol = ncol(M),
    dimnames = list(species, as.character(times))
  )
  invalid <- is.na(M) | !(M == 0 | M == 1)
  if (any(invalid)) {
    position <- which(invalid, arr.ind = TRUE)[1, ]
    stop("Detection must be 0 or 1: species '", species[position[1]],
         "', time '", times[position[2]], "'.", call. = FALSE)
  }
  order_time <- order(times)
  if (!identical(order_time, seq_along(times))) {
    message("Time columns were reordered into increasing order.")
    times <- times[order_time]
    M <- M[, order_time, drop = FALSE]
  }
  list(M = M, times = times, type = input_type)
}


.matrix_near <- function(x, y) {
  is.finite(x) & is.finite(y) & abs(x - y) <= 1e-8
}


.select_matrix_baseline <- function(M, times, baseline,
                                    baseline_min_species) {
  if (!is.null(baseline)) {
    index <- which(.matrix_near(times, baseline))
    if (length(index) != 1L) {
      stop("'baseline' must match an interval label within the selected range.",
           call. = FALSE)
    }
    return(index)
  }
  if (!is.null(baseline_min_species)) {
    index <- which(colSums(M) >= baseline_min_species)
    if (!length(index)) {
      stop("No interval meets 'baseline_min_species' within the selected range.",
           call. = FALSE)
    }
    return(index[1])
  }
  1L
}
