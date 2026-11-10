# Build the public output without changing the fitted model or diagnostics.
.modgee_public_output <- function(M, times, fit, statistics, species, group_info) {
  counts <- rowSums(M)
  TT <- ncol(M)
  known <- counts > 0
  n_known <- sum(known)
  extinct <- sum(species$p_extinct[known])
  dark <- as.numeric(statistics$X_T_hat)
  repeated <- sum(counts >= 2)
  n_groups <- if (is.null(group_info)) 1L else length(group_info$levels)
  n_parameters <- (repeated - 1L) + TT + (TT - 3L) + (n_groups - 1L)
  runs <- fit$optimizer_runs
  fit_ok <- length(runs) >= 2L && all(vapply(runs[1:2], function(x) {
    isTRUE(x$convcode[1L] == 0)
  }, logical(1L)))
  summary <- list(
    baseline_timestep = times[1L],
    final_timestep = times[TT],
    excluded_species = sum(!known),
    known_species = n_known,
    singleton_species = sum(counts == 1),
    known_extant = n_known - extinct,
    known_extinct = extinct,
    dark_extinct = dark,
    total_extinct = extinct + dark,
    original_richness = n_known + dark,
    cumulative_extinction_rate = as.numeric(statistics$Mus_overall[TT]),
    AIC = if (fit_ok) 2 * n_parameters - 2 * fit$LL else NA_real_
  )
  first <- last <- rep(NA_real_, nrow(M))
  for (i in which(known)) {
    detected <- which(M[i, ] == 1)
    first[i] <- times[detected[1L]]
    last[i] <- times[utils::tail(detected, 1L)]
  }
  output_species <- data.frame(species = rownames(M), row.names = NULL)
  if (!is.null(group_info))
    output_species$belonged_group <- unname(group_info$species_group)
  output_species$detected_timesteps <- as.integer(counts)
  output_species$is_singleton <- counts == 1
  output_species$first_detection_timestep <- first
  output_species$last_detection_timestep <- last
  output_species$d <- species$d
  output_species$p_extinct <- species$p_extinct
  output_species$estimated_extinction_year <- species$expected_extinction_time
  output_species$estimated_extinction_year[
    !is.finite(output_species$estimated_extinction_year)] <- NA_real_
  curve <- data.frame(time = times, h = as.numeric(fit$hs_est),
    m = as.numeric(fit$ms_est),
    cumulative_extinction_rate = as.numeric(statistics$Mus_overall))
  group_results <- NULL
  if (!is.null(group_info)) {
    old <- statistics$group_results$summary
    labels <- unname(group_info$species_group)
    # Include supplied all-zero-only groups, with unavailable model estimates.
    levels <- c(group_info$levels, setdiff(unique(labels), group_info$levels))
    group_summary <- do.call(rbind, lapply(levels, function(g) {
      rows <- labels == g
      recorded <- rows & known
      ng <- sum(recorded)
      eg <- sum(species$p_extinct[recorded])
      j <- match(g, old$group)
      dg <- if (is.na(j)) NA_real_ else old$dark_extinction[j]
      rate <- if (is.na(j)) NA_real_ else old$cumulative_extinction[j]
      data.frame(group = g, l = unname(exp(fit$ls_est[match(g, group_info$levels)])),
        baseline_timestep = times[1L],
        final_timestep = times[TT], excluded_species = sum(rows & !known),
        known_species = ng, singleton_species = sum(rows & counts == 1),
        known_extant = ng - eg, known_extinct = eg, dark_extinct = dg,
        total_extinct = eg + dg, original_richness = ng + dg,
        cumulative_extinction_rate = rate, row.names = NULL)
    }))
    group_curve <- statistics$group_results$curve
    names(group_curve)[names(group_curve) == "cumulative_extinction"] <-
      "cumulative_extinction_rate"
    group_results <- list(curve = group_curve, summary = group_summary)
  }
  list(summary = summary, curve = curve, species = output_species,
       group_results = group_results)
}
