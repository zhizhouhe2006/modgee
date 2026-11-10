# Normalize mutually exclusive groups before filtering or reordering species.
.normalize_modgee_groups <- function(groups, species_names, counts,
                                     reference_group = NULL) {
  if (!is.null(reference_group)) {
    if (!is.character(reference_group) || !is.null(dim(reference_group)) ||
        length(reference_group) != 1L || is.na(reference_group) ||
        !nzchar(trimws(reference_group)))
      stop("reference_group must be NULL or a single nonempty character string.",
           call. = FALSE)
    reference_group <- unname(trimws(reference_group))
    if (is.null(groups))
      stop("reference_group requires groups.", call. = FALSE)
  }
  if (is.null(groups)) return(NULL)
  clean <- function(x, what, unique_required = FALSE) {
    if (!is.character(x) && !is.factor(x))
      stop(what, " must be character or factor.", call. = FALSE)
    x <- trimws(as.character(x))
    if (anyNA(x) || any(!nzchar(x)))
      stop(what, " must be nonempty and not missing.", call. = FALSE)
    if (unique_required && anyDuplicated(x))
      stop(what, " must be unique after trimming whitespace.", call. = FALSE)
    x
  }
  source <- "vector"
  ids <- NULL
  level_order <- NULL
  if (is.character(groups) && length(groups) == 1L && is.null(dim(groups))) {
    source <- "csv"
    if (is.na(groups) || !file.exists(groups))
      stop("groups CSV file does not exist.", call. = FALSE)
    groups <- utils::read.csv(groups, check.names = FALSE,
                              stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  }
  if (is.data.frame(groups)) {
    if (source != "csv") source <- "data.frame"
    if (ncol(groups) != 2L) stop("groups must have exactly two columns.", call. = FALSE)
    ids <- clean(groups[[1]], "Group species names", TRUE)
    g <- groups[[2]]
  } else if (is.matrix(groups)) {
    source <- "matrix"
    if (!(is.double(groups) || is.integer(groups) || is.logical(groups)) ||
        anyNA(groups) || any(!(groups == 0 | groups == 1)) ||
        ncol(groups) < 1L || any(rowSums(groups) != 1))
      stop("Membership matrix must be 0/1 with exactly one 1 per row.", call. = FALSE)
    if (is.null(colnames(groups))) stop("Membership matrix needs group column names.", call. = FALSE)
    level_order <- clean(colnames(groups), "Group column names", TRUE)
    ids <- rownames(groups)
    g <- level_order[max.col(groups, ties.method = "first")]
  } else {
    if (!(is.character(groups) || is.factor(groups)) || !is.null(dim(groups)))
      stop("Unsupported groups input.", call. = FALSE)
    ids <- names(groups)
    g <- groups
  }
  if (is.factor(g)) level_order <- clean(levels(g), "Factor levels", TRUE)
  g <- clean(g, "Group labels")
  ignored <- character()
  if (!is.null(ids)) {
    ids <- clean(ids, "Group species names", TRUE)
    index <- match(species_names, ids)
    if (anyNA(index)) stop("Missing groups for: ",
      paste(species_names[is.na(index)], collapse = ", "), call. = FALSE)
    ignored <- setdiff(ids, species_names)
    if (length(ignored)) warning(length(ignored), " extra species in groups ignored.", call. = FALSE)
    g <- g[index]
  } else if (length(g) != length(species_names)) {
    stop("Unnamed groups must have one entry per input species.", call. = FALSE)
  }
  if (is.null(level_order)) level_order <- sort(unique(g), method = "radix")
  active <- level_order[level_order %in% g[counts > 0]]
  removed <- setdiff(level_order, active)
  if (!is.null(reference_group)) {
    if (!reference_group %in% active) {
      if (reference_group %in% removed)
        stop("reference_group '", reference_group,
             "' has no recorded species in the supplied detection matrix.",
             call. = FALSE)
      stop("Unknown reference_group '", reference_group,
           "'. Available active groups: ", paste(active, collapse = ", "),
           call. = FALSE)
    }
    active <- c(reference_group, active[active != reference_group])
  }
  if (length(removed)) warning("Groups without recorded species removed from fitting: ",
                              paste(removed, collapse = ", "), call. = FALSE)
  repeated <- vapply(active, function(x) sum(g == x & counts >= 2), integer(1))
  if (any(repeated < 2L)) stop("Each active group needs at least two species detected in two or more intervals: ",
                              paste(active[repeated < 2L], collapse = ", "), call. = FALSE)
  H <- outer(g, active, "==") * 1
  dimnames(H) <- list(species_names, active)
  list(species_group = stats::setNames(g, species_names), H = H,
       levels = active, reference_group = active[1], removed_groups = removed,
       matching = list(source = source, matched_by_name = !is.null(ids),
                       ignored_species = ignored, trimmed_whitespace = TRUE))
}

.modgee_group_starts <- function(x, levels) {
  if (is.null(x)) return(stats::setNames(rep(0, length(levels)), levels))
  if (!(is.double(x) || is.integer(x)) || !is.null(dim(x)) ||
      length(x) != length(levels) || any(!is.finite(x)))
    stop("ls_start must be a finite numeric vector with one value per active group.", call. = FALSE)
  if (!is.null(names(x))) {
    nm <- trimws(names(x))
    if (anyNA(nm) || anyDuplicated(nm) || !setequal(nm, levels))
      stop("ls_start names must match the active groups.", call. = FALSE)
    x <- x[match(levels, nm)]
  }
  if (abs(x[1]) > 1e-8) stop("Reference group ls_start must be zero.", call. = FALSE)
  x[1] <- 0
  stats::setNames(x, levels)
}

# Summarize only observed groups, using the original baseline-richness weighting.
.summary_modgee_grouped <- function(fit, times, group_info) {
  M <- fit$M
  n <- nrow(M)
  TT <- ncol(M)
  levels <- group_info$levels
  G <- length(levels)
  membership <- max.col(fit$H, ties.method = "first")
  hazards <- exp(fit$ls_est) %o% fit$ms_est
  if (any(!is.finite(hazards))) stop("Non-finite group hazards.")
  Mus <- t(vapply(seq_len(G), function(j) -expm1(-cumsum(hazards[j, ])), numeric(TT)))
  p <- matrix(NA_real_, n, TT)
  E_t <- var_t <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    z <- compute_P_extinct_i(M[i, ], fit$ds_est[i], fit$hs_est, hazards[membership[i], ])
    p[i, ] <- z[[1]]
    E_t[i] <- z[[2]]
    var_t[i] <- z[[3]]
  }
  if (any(!is.finite(p)) || any(p < -1e-8 | p > 1 + 1e-8))
    stop("Invalid species extinction probabilities in grouped summary.")
  S <- E <- matrix(0, G, TT, dimnames = list(levels, as.character(times)))
  for (j in seq_len(G)) {
    ii <- membership == j
    S[j, ] <- colSums(1 - p[ii, , drop = FALSE])
    E[j, ] <- colSums(p[ii, , drop = FALSE])
  }
  # exp(-sum(m)) avoids subtracting a rounded cumulative probability from 1.
  survival <- exp(-rowSums(hazards))
  if (any(!is.finite(survival) | survival <= 0))
    stop("Group survival underflow; baseline richness cannot be estimated.")
  richness <- S[, TT] / survival
  if (any(!is.finite(richness)) || sum(richness) <= 0)
    stop("Invalid inferred baseline richness.")
  total_extinct <- Mus * richness
  X <- total_extinct - E
  overall <- colSums(total_extinct) / sum(richness)
  curves <- do.call(rbind, lapply(seq_len(G), function(j)
    data.frame(group = levels[j], time = times, m = hazards[j, ],
               cumulative_extinction = Mus[j, ], row.names = NULL)))
  input_g <- group_info$species_group
  summary <- data.frame(
    group = levels,
    input_species = as.integer(table(factor(input_g, levels = levels))),
    recorded_species = tabulate(membership, nbins = G),
    repeated_species = tabulate(membership[rowSums(M) >= 2], nbins = G),
    l = exp(as.numeric(fit$ls_est)), hazard_multiplier = exp(as.numeric(fit$ls_est)),
    surviving_species = S[, TT], recorded_extinctions = E[, TT],
    cumulative_extinction = Mus[, TT], dark_extinction = X[, TT],
    original_richness = richness, row.names = NULL)
  list(Mus_overall = overall, p_exts = p[, TT], E_t_exts = E_t,
       var_t_exts = var_t, S_T_hat = sum(S[, TT]), E_T_hat = sum(E[, TT]),
       p_T_hat = overall[TT], X_T_hat = sum(X[, TT]),
       X_hat_by_subgroup = X,
       group_results = list(curve = curves, summary = summary))
}
