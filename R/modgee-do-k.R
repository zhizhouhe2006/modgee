.modgee_compare_k <- function(fit, grouped, k_start = NULL) {
  M <- fit$M
  repeated <- rowSums(M) >= 2L
  TT <- ncol(M)
  log_l <- if (is.null(fit$ls_est)) numeric() else fit$ls_est
  reference_l <- exp(log_l)
  n_parameters <- (sum(repeated) - 1L) + TT + (TT - 3L) +
    length(log_l[-1L])

  # Only the shared detection stage and the extinction stage enter AIC.
  runs <- fit$optimizer_runs
  fixed_ok <- length(runs) >= 2L &&
    all(vapply(runs[1:2], function(run) {
      isTRUE(run$convcode[1L] == 0)
    }, logical(1L)))
  fixed_AIC <- if (fixed_ok) 2 * n_parameters - 2 * fit$LL else NA_real_
  free_k <- free_AIC <- NA_real_
  free_l <- rep(NA_real_, length(log_l))
  if (length(free_l)) free_l[1L] <- 1
  free_ok <- FALSE

  if (fixed_ok) {
    likelihood <- if (grouped) .compute_neg_LL_extinct_param_groups else
      compute_neg_LL_extinct_param
    fit_M <- M[repeated, , drop = FALSE]
    fit_H <- if (grouped) fit$H[repeated, , drop = FALSE] else NULL
    d <- fit$ds_est[repeated]
    h <- fit$hs_est[-TT]
    k_index <- TT - 1L
    initial <- c(log(fit$hs_est[TT]), log(fit$ms_est[3:(TT - 1L)]),
                 0, log_l[-1L])

    objective <- function(par) {
      value <- tryCatch(
        likelihood(par, l_M = fit_M, l_H = fit_H, l_ds = d, l_hs = h,
                   do_k = TRUE, do_l = grouped),
        error = function(e) Inf
      )
      if (is.finite(value)) value else Inf
    }
    tolerance <- 1e-7 * max(1, abs(fit$LL))
    nested_value <- objective(initial)
    if (any(d != 1) && is.finite(nested_value) &&
        abs(nested_value + fit$LL) <= tolerance) {
      if (!is.null(k_start)) initial[k_index] <- k_start
      alternative <- tryCatch(
        stats::optim(initial, objective, method = "BFGS",
                     control = list(maxit = 500L)),
        error = function(e) NULL
      )
      if (!is.null(alternative) &&
          isTRUE(alternative$convergence == 0L) &&
          all(is.finite(alternative$par)) &&
          is.finite(alternative$value) &&
          alternative$value <= nested_value + tolerance) {
        candidate_l <- if (length(log_l)) {
          c(1, exp(alternative$par[k_index + seq_along(log_l[-1L])]))
        } else {
          numeric()
        }
        if (all(is.finite(candidate_l) & candidate_l > 0)) {
          free_ok <- TRUE
          free_k <- alternative$par[k_index]
          free_l <- candidate_l
          free_AIC <- 2 * (n_parameters + 1L) + 2 * alternative$value
        }
      }
    }
  }

  comparison <- data.frame(model = c("fixed_k", "free_k"), k = c(0, free_k))
  for (i in seq_along(log_l)) {
    comparison[[paste0("l_", names(log_l)[i])]] <- c(reference_l[i], free_l[i])
  }
  comparison$AIC <- c(fixed_AIC, free_AIC)
  comparison$converged <- c(fixed_ok, free_ok)
  comparison
}

