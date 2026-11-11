.bootstrap_validate <- function(bootstrap, CI) {
  scalar <- function(x) (is.double(x) || is.integer(x)) && !is.object(x) &&
    is.null(dim(x)) && length(x) == 1L && is.finite(x)
  if (!is.null(bootstrap) &&
      (!scalar(bootstrap) || bootstrap < 3 || bootstrap != floor(bootstrap) ||
       bootstrap > .Machine$integer.max))
    stop("bootstrap must be NULL or a single integer greater than or equal to 3.", call. = FALSE)
  if (!scalar(CI) || CI <= 0 || CI >= 1)
    stop("CI must be a single finite number strictly between 0 and 1.", call. = FALSE)
}

.bootstrap_metrics <- function() c("known_extant", "known_extinct", "dark_extinct",
  "total_extinct", "original_richness", "cumulative_extinction_rate")

.bootstrap_values <- function(result) {
  metrics <- .bootstrap_metrics()
  summary <- data.frame(metric = metrics,
    value = vapply(metrics, function(m) as.numeric(result$summary[[m]]), numeric(1)))
  curve <- data.frame(time = result$curve$time,
                       value = result$curve$cumulative_extinction_rate)
  gs <- gc <- NULL
  if (!is.null(result$group_results)) {
    x <- result$group_results$summary
    gs <- do.call(rbind, lapply(seq_len(nrow(x)), function(i)
      data.frame(group = x$group[i], metric = metrics,
                 value = as.numeric(x[i, metrics]))))
    x <- result$group_results$curve
    gc <- data.frame(group = x$group, time = x$time,
                      value = x$cumulative_extinction_rate)
  }
  list(summary = summary, curve = curve, group_summary = gs, group_curve = gc)
}

.bootstrap_fit_replicate <- function(plan, M, info) {
  idx <- plan$indices
  Mb <- M[idx, , drop = FALSE]
  rownames(Mb) <- sprintf("boot_%06d", seq_along(idx))
  gb <- ref <- NULL
  if (!is.null(info)) {
    labels <- unname(info$species_group[idx])
    counts <- rowSums(Mb)
    repeats <- vapply(info$levels, function(g)
      sum(labels == g & counts >= 2L), integer(1))
    if (any(repeats < 2L))
      stop("Bootstrap sample lacks two repeated-detection species in: ",
           paste(info$levels[repeats < 2L], collapse = ", "), call. = FALSE)
    gb <- stats::setNames(factor(labels, levels = info$levels), rownames(Mb))
    ref <- info$reference_group
  }
  .bootstrap_set_rng(plan$fit_rng)
  .fit_modgee_impl(Mb, groups = gb, reference_group = ref, do_k = FALSE)
}

.bootstrap_run_replicate <- function(plan, M, info) {
  warnings <- character()
  started <- proc.time()[[3L]]
  # Interrupt conditions are deliberately not caught as ordinary errors.
  fit <- tryCatch(withCallingHandlers(
    .bootstrap_fit_replicate(plan, M, info),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  failed <- inherits(fit, "error")
  state <- if (failed) "failed" else if (isTRUE(fit$converged))
    "success" else "nonconverged"
  list(state = state, warnings = unique(warnings),
       error = if (failed) conditionMessage(fit) else NA_character_,
       elapsed = unname(proc.time()[[3L]] - started),
       convergence = if (failed) NULL else fit$convergence,
       values = if (failed) NULL else .bootstrap_values(fit),
       plan_hash = .bootstrap_hash(plan))
}

.bootstrap_modgee <- function(M, groups, reference_group, do_k, B, CI,
                              .finalize = base::identity) {
  if (!is.logical(do_k) || length(do_k) != 1L || !is.null(dim(do_k)) || is.na(do_k))
    stop("do_k must be TRUE or FALSE.", call. = FALSE)
  validated <- .validate_modgee_matrix(M)
  M <- validated$M
  counts <- rowSums(M)
  if (sum(counts >= 2L) < 2L)
    stop("At least two species must be detected in two or more intervals.", call. = FALSE)
  info <- .normalize_modgee_groups(groups, rownames(M), counts, reference_group)
  identity <- .bootstrap_identity(M, info)
  task <- .bootstrap_hash(identity)
  task_dir <- file.path(getwd(), "modgee_bootstrap", task)
  if (!dir.exists(task_dir) && !dir.create(task_dir, recursive = TRUE))
    stop("Cannot create bootstrap directory: ", task_dir, call. = FALSE)
  task_dir <- normalizePath(task_dir, winslash = "/", mustWork = TRUE)
  lock <- filelock::lock(file.path(task_dir, "task.lock"), timeout = 0)
  if (is.null(lock)) stop("Bootstrap task is locked by another process: ", task_dir,
                          call. = FALSE)
  on.exit(filelock::unlock(lock), add = TRUE)
  directory <- .bootstrap_select_run(task_dir, task)
  message("Bootstrap checkpoints: ", directory)
  old_rng <- .bootstrap_rng()
  on.exit(.bootstrap_set_rng(old_rng), add = TRUE)
  manifest_path <- file.path(directory, "manifest.rds")
  resumed <- file.exists(manifest_path)
  if (!resumed) {
    # Materialize an initial seed even when the caller has never set one.
    if (is.null(old_rng$seed)) stats::runif(1)
    .bootstrap_write(list(identity = identity, rng = .bootstrap_rng()),
                     manifest_path, task, "manifest")
  }
  manifest <- .bootstrap_read(manifest_path, task, "manifest")
  if (!identical(manifest$identity, identity)) stop("Bootstrap task identity mismatch.")
  original_path <- file.path(directory, "original_fit.rds")
  if (file.exists(original_path)) {
    original <- .bootstrap_read(original_path, task, "original")
    message("Reusing original fit.")
  } else {
    .bootstrap_set_rng(manifest$rng)
    private <- NULL
    result <- .fit_modgee_impl(M, groups = groups, reference_group = reference_group,
      do_k = FALSE, .capture = function(fit, grouped) {
        private <<- list(fit = fit, grouped = grouped)
      })
    # Independent per-replicate streams do not depend on optimizer consumption.
    seed <- sample.int(.Machine$integer.max, 1L)
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    set.seed(seed)
    original <- list(result = result, private = private, stream = .bootstrap_rng())
    .bootstrap_write(original, original_path, task, "original")
  }
  result <- original$result
  if (!isTRUE(result$converged))
    stop("Original fit did not converge; diagnostics saved in ", original_path,
         ". Review before starting bootstrap in a new task directory.", call. = FALSE)
  # Compute the comparison once from the saved fixed-k internal fit.
  if (do_k) {
    comparison_path <- file.path(directory, "k_comparison.rds")
    if (!file.exists(comparison_path)) {
      comparison <- .modgee_compare_k(original$private$fit, original$private$grouped)
      .bootstrap_write(comparison, comparison_path, task, "k_comparison")
    }
    result$do_k <- .bootstrap_read(comparison_path, task, "k_comparison")
  } else result["do_k"] <- list(NULL)
  for (subdir in c("plans", "results")) {
    p <- file.path(directory, subdir)
    if (!dir.exists(p) && !dir.create(p)) stop("Cannot create checkpoint directory: ", p)
  }
  recorded <- which(counts > 0L)
  stream <- original$stream
  records <- vector("list", B)
  for (b in seq_len(B)) {
    name <- sprintf("rep_%08d.rds", b)
    plan_path <- file.path(directory, "plans", name)
    result_path <- file.path(directory, "results", name)
    if (file.exists(plan_path)) {
      plan <- .bootstrap_read(plan_path, task, "plan", b)
      if (!identical(plan$stream, stream)) stop("Bootstrap random-stream mismatch: ", plan_path)
    } else {
      if (file.exists(result_path)) stop("Result exists without its sampling plan: ", result_path)
      .bootstrap_set_rng(stream)
      idx <- recorded[sample.int(length(recorded), length(recorded), replace = TRUE)]
      plan <- list(stream = stream, indices = idx,
                   source_species = rownames(M)[idx],
                   internal_species = sprintf("boot_%06d", seq_along(idx)),
                   fit_rng = .bootstrap_rng())
      .bootstrap_write(plan, plan_path, task, "plan", b)
    }
    stream$seed <- parallel::nextRNGStream(stream$seed)
    saved <- NULL
    if (file.exists(result_path)) {
      saved <- tryCatch(.bootstrap_read(result_path, task, "result", b),
                        error = function(e) NULL)
      if (!is.null(saved) && (!identical(saved$plan_hash, .bootstrap_hash(plan)) ||
          !saved$state %in% c("success", "failed", "nonconverged"))) saved <- NULL
      if (is.null(saved)) {
        quarantine <- tempfile(paste0(name, ".corrupt-"), tmpdir = dirname(result_path))
        if (!file.rename(result_path, quarantine)) stop("Cannot quarantine corrupt checkpoint.")
        warning("Corrupt result quarantined; rerunning the saved sample for replicate ", b,
                call. = FALSE)
      }
    }
    if (is.null(saved)) {
      message("Bootstrap ", b, "/", B, ": fitting")
      saved <- .bootstrap_run_replicate(plan, M, info)
      .bootstrap_write(saved, result_path, task, "result", b)
      message("Bootstrap ", b, "/", B, ": ", saved$state, " (saved)")
    } else message("Bootstrap ", b, "/", B, ": reused ", saved$state)
    records[[b]] <- saved
  }
  result$bootstrap <- .bootstrap_collect(result, records, CI, directory, resumed)
  # Finish public output before committing completion; errors remain resumable.
  result <- .finalize(result)
  .bootstrap_write(list(requested = B, CI = CI, do_k = do_k,
                        finished_at = as.character(Sys.time())),
                   file.path(directory, "completed.rds"), task, "completed")
  result
}

.bootstrap_collect <- function(result, records, CI, directory, resumed) {
  template <- .bootstrap_values(result)
  state <- vapply(records, `[[`, character(1), "state")
  diagnostics <- data.frame(replicate = seq_along(records), state = state,
    used_for_ci = FALSE,
    error = vapply(records, `[[`, character(1), "error"),
    elapsed = vapply(records, `[[`, numeric(1), "elapsed"))
  diagnostics$warnings <- I(lapply(records, `[[`, "warnings"))
  diagnostics$convergence <- I(lapply(records, `[[`, "convergence"))
  replicas <- intervals <- vector("list", length(template))
  names(replicas) <- names(intervals) <- names(template)
  for (name in names(template)) {
    target <- template[[name]]
    if (is.null(target)) next
    keys <- setdiff(names(target), "value")
    values <- matrix(NA_real_, length(records), nrow(target))
    for (b in seq_along(records)) {
      source <- records[[b]]$values[[name]]
      if (is.null(source)) next
      # Match structured keys without delimiter-dependent string concatenation.
      for (j in seq_len(nrow(target))) {
        match <- rep(TRUE, nrow(source))
        for (key in keys) match <- match & source[[key]] == target[[key]][j]
        if (sum(match) == 1L) values[b, j] <- source$value[match]
      }
    }
    replicas[[name]] <- do.call(rbind, lapply(seq_along(records), function(b) {
      x <- target; x$value <- values[b, ]; cbind(replicate = b, x)
    }))
    out <- target; names(out)[names(out) == "value"] <- "estimate"
    out$lower <- out$upper <- NA_real_
    out$n_valid <- integer(nrow(out))
    for (j in seq_len(nrow(out))) {
      eligible <- state == "success" & is.finite(values[, j])
      diagnostics$used_for_ci <- diagnostics$used_for_ci | eligible
      out$n_valid[j] <- sum(eligible)
      if (sum(eligible) >= 2L) {
        q <- stats::quantile(values[eligible, j],
          probs = c((1 - CI) / 2, (1 + CI) / 2), type = 7, names = FALSE)
        out$lower[j] <- q[1]; out$upper[j] <- q[2]
      }
    }
    intervals[[name]] <- out
  }
  B <- length(records)
  n_success <- sum(state == "success")
  if (n_success < 100L || any(state != "success"))
    warning("Bootstrap: ", n_success, "/", B,
      " successful fits. Inspect diagnostics and per-interval n_valid; small samples give unstable intervals.",
      call. = FALSE)
  if (any(vapply(intervals, function(x) !is.null(x) &&
      any(x$n_valid < 2L & is.finite(x$estimate)), logical(1))))
    warning("Some intervals have fewer than two valid estimates and are NA.", call. = FALSE)
  list(settings = list(requested = B, CI = CI,
    resampling = "global species rows with replacement; full refit; fixed k = 0",
    interval_method = "pointwise percentile (quantile type 7)",
    checkpoint_dir = directory, schema_version = 1L),
    status = list(requested = B, completed = B, successful = n_success,
      failed = sum(state == "failed"), nonconverged = sum(state == "nonconverged"),
      pending = 0L, resumed = resumed),
    summary = intervals$summary, curve = intervals$curve,
    group_results = if (is.null(result$group_results)) NULL else
      list(summary = intervals$group_summary, curve = intervals$group_curve),
    diagnostics = diagnostics, replicates = replicas)
}

# Called while holding the task-level lock. Legacy flat checkpoints stay untouched.
.bootstrap_select_run <- function(task_dir, task) {
  runs <- list.dirs(task_dir, full.names = TRUE, recursive = FALSE)
  runs <- sort(runs[grepl("^run_[0-9]{8}$", basename(runs))])
  pending <- character()
  for (run in runs) {
    marker <- file.path(run, "completed.rds")
    if (file.exists(marker)) {
      .bootstrap_read(marker, task, "completed")
    } else pending <- c(pending, run)
  }
  if (length(pending) > 1L)
    stop("Multiple unfinished bootstrap runs found in: ", task_dir, call. = FALSE)
  if (length(pending)) return(pending[[1L]])
  index <- if (length(runs)) max(as.integer(sub("run_", "", basename(runs)))) + 1L else 1L
  if (is.na(index) || index > 99999999L) stop("Bootstrap run limit reached.")
  directory <- file.path(task_dir, sprintf("run_%08d", index))
  if (!dir.create(directory)) stop("Cannot create bootstrap run: ", directory, call. = FALSE)
  directory
}
