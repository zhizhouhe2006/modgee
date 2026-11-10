# Checkpoints are immutable envelopes, committed on the same filesystem.
.bootstrap_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2)

.bootstrap_read <- function(path, task, kind, index = 0L) {
  x <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.list(x) || !identical(x$schema, 1L) ||
      !identical(x$task, task) || !identical(x$kind, kind) ||
      !identical(x$index, as.integer(index)) ||
      !identical(x$checksum, .bootstrap_hash(x$payload)))
    stop("Invalid bootstrap checkpoint: ", path, call. = FALSE)
  x$payload
}

.bootstrap_write <- function(payload, path, task, kind, index = 0L) {
  if (file.exists(path)) stop("Refusing to overwrite checkpoint: ", path, call. = FALSE)
  tmp <- tempfile(".pending-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  x <- list(schema = 1L, task = task, kind = kind, index = as.integer(index),
            checksum = .bootstrap_hash(payload), payload = payload)
  saveRDS(x, tmp, version = 2)
  .bootstrap_read(tmp, task, kind, index)
  if (!file.rename(tmp, path)) stop("Could not commit checkpoint: ", path, call. = FALSE)
  invisible(NULL)
}

.bootstrap_rng <- function() {
  list(kind = RNGkind(), seed = if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL)
}

.bootstrap_set_rng <- function(x) {
  do.call(RNGkind, as.list(x$kind))
  if (is.null(x$seed)) {
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", x$seed, .GlobalEnv)
  invisible(NULL)
}

.bootstrap_identity <- function(M, info) {
  # Hash numerical R bodies too, so development edits cannot reuse stale fits.
  env <- environment(.fit_modgee_impl)
  names <- c(".fit_modgee_impl", ".normalize_modgee_groups",
    ".modgee_group_starts", ".summary_modgee_grouped", ".modgee_public_output",
    ".est_modgee_grouped", ".compute_neg_LL_extinct_param_groups",
    ".compute_neg_LL_detect_param_groups", ".compute_neg_LL_detect_param_1_obs_groups",
    "est_detect_extinct_parms_ML", "est_detect_parms", "compute_summary_stats",
    "compute_P_extinct_i", "compute_neg_LL_detect_param",
    "compute_neg_LL_extinct_param", "compute_neg_LL_detect_param_1_obs")
  bodies <- lapply(names, function(n)
    paste(deparse(body(get(n, envir = env)), control = c("keepInteger", "keepNA")),
          collapse = "\n"))
  # Increment engine when numerical C++ source or fitting policy changes.
  list(schema = 1L, engine = "modgee-bootstrap-full-k0-v1", bodies = bodies,
       native_source = "5dff6e9cd227edd4d16a4c8c08c84b8a3523464018dafeda836516815dae2820",
       R = as.character(getRversion()),
       dependencies = vapply(c("Rcpp", "optimx"), function(p)
         as.character(utils::packageVersion(p)), character(1)),
       M = matrix(as.numeric(M), nrow(M), dimnames = dimnames(M)),
       groups = if (is.null(info)) NULL else list(
         labels = info$species_group, levels = info$levels,
         removed = info$removed_groups, reference = info$reference_group),
       method = "global_species_rows_full_refit_k0")
}
