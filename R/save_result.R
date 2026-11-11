#' Save MODGEE results as CSV and RDS files
#'
#' Export existing estimates without refitting or recomputing intervals.
#' @param result A complete result returned by \code{fit_modgee()}.
#' @param path Output directory, created recursively if needed. Defaults to
#'   the current working directory.
#' @param prefix Filename prefix. An underscore is appended unless the prefix
#'   is empty or already ends with one. Defaults to \code{"modgee"}.
#' @param overwrite Allow replacement of existing target files. Defaults to FALSE.
#' @return
#' Writes three to five CSV files and one RDS file to \code{path}. Each filename starts with
#' \code{prefix} followed by an underscore; an empty prefix adds nothing.
#' For example, \code{prefix = "hk"} produces \file{hk_summary.csv} and,
#' when bootstrap results are available, \file{hk_CI.csv}.
#'
#' \strong{CSV files}
#' \describe{
#'   \item{\file{summary.csv}}{One overall row, followed by one row per supplied
#' group when groups are present. Columns \code{scope} ("overall" or "group")
#' and \code{group} identify the rows; the overall group label is \code{NA}.
#' The remaining overall columns are \code{baseline_timestep},
#' \code{final_timestep}, \code{excluded_species}, \code{known_species},
#' \code{singleton_species}, \code{known_extant}, \code{known_extinct},
#' \code{dark_extinct}, \code{total_extinct}, \code{original_richness},
#' \code{cumulative_extinction_rate} and \code{AIC}, copied from
#' \code{result$summary}. Group rows contain the corresponding group estimates.
#' With groups, \code{l} is also included: it is the group's hazard multiplier,
#' with 1 for the reference group and NA for the overall row and unfitted groups.
#' AIC is an overall model statistic and is NA on group rows. Groups containing
#' only all-zero species remain in this table, with undefined estimates as NA.}
#'   \item{\file{time_result.csv}}{A long-format table: one row per time step
#' for the overall model, followed by one row per time step for each fitted group.
#' Columns are \code{scope}, \code{group}, \code{time}, \code{h}, \code{m}
#' and \code{cumulative_extinction_rate}. The overall rows contain the shared
#' survey parameter h, the reference hazard m and the overall cumulative
#' extinction rate. Group rows contain their own hazard m and cumulative
#' extinction rate; their h is NA because the shared value is already recorded
#' in the overall rows. Time labels are copied from the fitted result.}
#'   \item{\file{species_result.csv}}{One row per input species, in the original
#' input order. Copies \code{result$species}, including \code{species},
#' \code{belonged_group} (when grouped), \code{detected_timesteps},
#' \code{is_singleton}, \code{first_detection_timestep},
#' \code{last_detection_timestep}, \code{d}, \code{p_extinct} and
#' \code{estimated_extinction_year}. All-zero species are retained with NA
#' estimates. The extinction year is conditional on extinction, as defined
#' in \code{\link{fit_modgee}}.}
#'   \item{\file{k_comparison.csv}}{Written only when \code{result$do_k} is
#' not NULL. Copies that table: fixed-k and free-k model rows, with \code{model},
#' \code{k}, group multipliers \code{l_<group>} when present, \code{AIC}
#' and \code{converged}. No new models are fitted.}
#'   \item{\file{CI.csv}}{Written only when \code{result$bootstrap} is not NULL.
#' Combines \code{result$bootstrap$curve} and, when present,
#' \code{result$bootstrap$group_results$curve}. Overall rows come first,
#' followed by group rows. Columns \code{scope}, \code{group} and \code{time}
#' identify each curve point. \code{estimate} is the original fitted cumulative
#' extinction rate; \code{lower} and \code{upper} are its pointwise bootstrap
#' confidence limits at the confidence level used in \code{fit_modgee()}.
#' \code{n_valid} is the number of valid bootstrap estimates used for that point.
#' Undefined limits remain NA. Only curve intervals are exported, not intervals
#' for summary metrics, parameters or individual species.}
#' }
#'
#' @details
#' Always writes summary.csv, time_result.csv and species_result.csv, with
#' the requested prefix. Summary and time tables contain overall rows first,
#' followed by group rows. The scope column is "overall" or "group"; group is
#' NA for overall rows, so a real group named "overall" remains unambiguous.
#' Missing columns are filled with NA. Group summaries retain l; their AIC
#' is NA. Group time rows contain the group's own m and cumulative extinction
#' rate; h is NA because h is shared and recorded on overall rows only.
#' All-zero-only groups remain in summary but have no fitted curve rows.
#' Species results retain all input rows, columns and their original order.
#'
#' If result$do_k is non-NULL, k_comparison.csv contains that table unchanged
#' apart from omission of row names. If result$bootstrap is non-NULL,
#' CI.csv combines bootstrap$curve and
#' bootstrap$group_results$curve, with scope and group columns. These are
#' pointwise cumulative-extinction intervals: time, estimate, lower, upper
#' and n_valid. No summary-metric intervals or replicate tables are exported to CSV.
#'
#' CSV files use UTF-8, comma separators, NA for missing values and no row
#' numbers. Estimates are not rounded, clipped or recomputed. All tables and
#' filename conflicts are checked before writing. All files are staged before
#' installation; a failed installation attempts to restore previous files.
#' This is not a concurrent-writer transaction; use a separate directory for
#' simultaneous exports. Optional files left by earlier calls are not deleted;
#' use distinct prefixes for different analyses. The prefixed result.rds file
#' preserves the complete R object, including diagnostics and bootstrap information.
#' @section Complete R object:
#' The RDS filename is \file{modgee_result.rds} with the default prefix,
#' \file{hk_result.rds} with \code{prefix = "hk"}, or \file{result.rds} with
#' an empty prefix. Restore the complete object with \code{readRDS()}.
#' RDS exports use the same overwrite checks, staging and rollback as CSV files.
#' The function invisibly returns a data frame with columns \code{table},
#' \code{file}, \code{nrow} and \code{ncol}. The RDS row has table = "result"
#' and NA for both dimensions.
#' @export
#' @seealso
#' \strong{Functions:}
#' \code{\link{read_matrix}}, \code{\link{fit_modgee}}, \code{\link{plot_modgee}}
#'
#' \strong{Examples:}
#' \code{\link{example_file}}, \code{\link{example_M}}, \code{\link{example_groups}}, \code{\link{example_result}}, \code{\link{example_grouped_result}}, \code{\link{example_bootstrap_result}}
#' @section Usage Examples:
#' \preformatted{
#' # Save ungrouped results.
#' save_result(example_result, path = "results")
#'
#' # Save grouped results.
#' save_result(example_grouped_result, path = "results", prefix = "grouped")
#'
#' # Save bootstrap results.
#' save_result(example_bootstrap_result, path = "results", prefix = "bootstrap")
#' }
save_result <- function(result, path = ".", prefix = "modgee", overwrite = FALSE) {
  scalar_text <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && is.null(dim(x))
  }
  if (!scalar_text(path) || !nzchar(trimws(path))) {
    stop("path must be one nonempty directory path.", call. = FALSE)
  }
  if (!scalar_text(prefix) || grepl('[<>:"/\\\\|?*[:cntrl:]]', prefix)) {
    stop("prefix must be one filename prefix without path separators or reserved characters.",
         call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L ||
      is.na(overwrite) || !is.null(dim(overwrite))) {
    stop("overwrite must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.list(result) || !all(c("summary", "curve", "species") %in% names(result))) {
    stop("result must be a complete fit_modgee result.", call. = FALSE)
  }
  check_table <- function(x, label, required = character()) {
    if (!is.data.frame(x) || !nrow(x) || is.null(names(x)) ||
        anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x)) ||
        !all(required %in% names(x)) ||
        any(!vapply(x, function(v) is.atomic(v) && is.null(dim(v)), logical(1)))) {
      stop("Invalid result table: ", label, ".", call. = FALSE)
    }
    x
  }
  tag <- function(x, scope) {
    if ("scope" %in% names(x) || (scope == "overall" && "group" %in% names(x))) {
      stop("Reserved scope/group columns in overall or untagged table.", call. = FALSE)
    }
    if (scope == "group") {
      if (!"group" %in% names(x) || anyNA(x$group) ||
          any(!nzchar(trimws(as.character(x$group))))) {
        stop("Group tables require nonempty group labels.", call. = FALSE)
      }
      x$group <- as.character(x$group)
    } else {
      x$group <- NA_character_
    }
    x$scope <- scope
    x[c("scope", "group", setdiff(names(x), c("scope", "group")))]
  }
  bind_tables <- function(x, y) {
    columns <- union(names(x), names(y))
    for (nm in setdiff(columns, names(x))) x[[nm]] <- NA
    for (nm in setdiff(columns, names(y))) y[[nm]] <- NA
    out <- rbind(x[columns], y[columns])
    rownames(out) <- NULL
    out
  }
  check_time <- function(x, times, label, grouped = FALSE) {
    if (!is.numeric(x$time) || any(!is.finite(x$time))) {
      stop("Invalid time labels in ", label, ".", call. = FALSE)
    }
    blocks <- if (grouped) split(x$time, as.character(x$group)) else list(x$time)
    if (any(vapply(blocks, function(t) anyDuplicated(t) > 0L ||
                   !setequal(t, times), logical(1)))) {
      stop("Duplicated or mismatched time labels in ", label, ".", call. = FALSE)
    }
  }

  summary <- result$summary
  if (!is.data.frame(summary)) {
    if (!is.list(summary) || !length(summary) || is.null(names(summary)) ||
        anyDuplicated(names(summary)) || anyNA(names(summary)) ||
        any(!nzchar(names(summary))) ||
        any(!vapply(summary, function(v) is.atomic(v) && length(v) == 1L &&
                      is.null(dim(v)), logical(1)))) {
      stop("result$summary must be a named scalar list or one-row data frame.",
           call. = FALSE)
    }
    summary <- as.data.frame(summary, check.names = FALSE, stringsAsFactors = FALSE)
  }
  summary <- check_table(summary, "summary", "cumulative_extinction_rate")
  if (nrow(summary) != 1L) stop("Overall summary must have one row.", call. = FALSE)
  curve <- check_table(result$curve, "curve",
                       c("time", "h", "m", "cumulative_extinction_rate"))
  check_time(curve, curve$time, "curve")
  species <- check_table(result$species, "species", "species")
  if (anyNA(species$species) || anyDuplicated(species$species) ||
      any(!nzchar(trimws(as.character(species$species))))) {
    stop("Species identifiers must be nonempty and unique.", call. = FALSE)
  }
  tables <- list(summary = tag(summary, "overall"),
                 time_result = tag(curve, "overall"), species_result = species)
  groups <- result$group_results
  if (!is.null(groups)) {
    gs <- check_table(groups$summary, "group_results$summary", c("group", "l"))
    gc <- check_table(groups$curve, "group_results$curve",
                      c("group", "time", "m", "cumulative_extinction_rate"))
    gs <- tag(gs, "group")
    gc <- tag(gc, "group")
    if (anyDuplicated(gs$group) || !all(gc$group %in% gs$group)) {
      stop("Inconsistent group summary and curve labels.", call. = FALSE)
    }
    check_time(gc, curve$time, "group_results$curve", TRUE)
    tables$summary <- bind_tables(tables$summary, gs)
    tables$time_result <- bind_tables(tables$time_result, gc)
  }
  if (!is.null(result$do_k)) {
    tables$k_comparison <- check_table(result$do_k, "do_k")
  }
  boot <- result$bootstrap
  if (!is.null(boot)) {
    required <- c("time", "estimate", "lower", "upper", "n_valid")
    bc <- check_table(boot$curve, "bootstrap$curve", required)
    check_time(bc, curve$time, "bootstrap$curve")
    bs <- tag(bc, "overall")
    bg <- boot$group_results
    if (!is.null(bg) && !is.null(bg$curve) && nrow(bg$curve) > 0L) {
      bgc <- check_table(bg$curve, "bootstrap$group_results$curve", c("group", required))
      bgc <- tag(bgc, "group")
      if (is.null(groups) || !setequal(bgc$group, groups$curve$group)) {
        stop("Bootstrap curve groups do not match fitted curve groups.", call. = FALSE)
      }
      check_time(bgc, curve$time, "bootstrap$group_results$curve", TRUE)
      bs <- bind_tables(bs, bgc)
    } else if (!is.null(groups)) {
      stop("Grouped bootstrap result is missing its group curve intervals.", call. = FALSE)
    }
    tables$CI <- bs
  }

  stem <- if (!nzchar(prefix) || endsWith(prefix, "_")) prefix else paste0(prefix, "_")
  path <- path.expand(path)
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create output directory: ", path, call. = FALSE)
  }
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  targets <- file.path(path, c(paste0(stem, names(tables), ".csv"),
                              paste0(stem, "result.rds")))
  existing <- file.exists(targets)
  if (any(dir.exists(targets))) stop("A target file path is a directory.", call. = FALSE)
  if (!overwrite && any(existing)) {
    stop("Files already exist; use overwrite = TRUE: ",
         paste(basename(targets[existing]), collapse = ", "), call. = FALSE)
  }
  stage <- tempfile(".modgee-export-", tmpdir = path)
  if (!dir.create(stage)) stop("Cannot stage result files in: ", path, call. = FALSE)
  keep_stage <- FALSE
  on.exit(if (!keep_stage) unlink(stage, recursive = TRUE), add = TRUE)
  staged <- file.path(stage, basename(targets))
  backups <- file.path(stage, paste0("backup-", seq_along(targets)))
  for (i in seq_along(tables)) {
    utils::write.csv(tables[[i]], staged[i], row.names = FALSE,
                     na = "NA", fileEncoding = "UTF-8")
  }
  saveRDS(result, staged[length(staged)])
  if (any(existing) && !all(file.copy(targets[existing], backups[existing]))) {
    stop("Cannot back up existing result files; no files were replaced.", call. = FALSE)
  }
  installed <- integer()
  tryCatch({
    for (i in seq_along(targets)) {
      installed <- c(installed, i)
      if (!file.copy(staged[i], targets[i], overwrite = overwrite)) {
        stop("Cannot write file: ", targets[i], call. = FALSE)
      }
    }
  }, error = function(e) {
    restored <- TRUE
    for (i in installed) {
      if (existing[i]) {
        restored <- isTRUE(file.copy(backups[i], targets[i], overwrite = TRUE)) && restored
      } else {
        restored <- (unlink(targets[i]) == 0L) && restored
      }
    }
    if (!restored) {
      # Keep backups available if an operating-system failure prevents rollback.
      keep_stage <<- TRUE
      warning("Some files could not be restored; backups remain at: ", stage,
              call. = FALSE)
    }
    stop(e)
  })
  files <- data.frame(table = c(names(tables), "result"), file = targets,
                      nrow = c(vapply(tables, nrow, integer(1)), NA_integer_),
                      ncol = c(vapply(tables, ncol, integer(1)), NA_integer_),
                      row.names = NULL)
  message("Saved ", length(tables), " CSV files and 1 RDS file to ", path)
  invisible(files)
}
