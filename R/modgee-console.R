.modgee_report_fit <- function(result) {
  diagnostics <- result$convergence
  bad <- is.na(diagnostics$convcode) | diagnostics$convcode != 0
  if (!any(bad)) {
    message("All fitting stages converged.")
  } else {
    message("Some fitting stages did not converge or could not be verified:")
    for (i in which(bad)) {
      label <- switch(diagnostics$stage[i],
        detection = "Stage 1 (detection)",
        extinction = "Stage 2 (extinction)",
        singleton = paste0("Singleton fit for '", diagnostics$species[i], "'"))
      code <- diagnostics$convcode[i]
      if (is.na(code)) {
        message("  ", label, ": convergence could not be verified.")
      } else if (code == 999L) {
        message("  ", label, " failed because of an optimizer error.")
      } else {
        message("  ", label, " did not converge.")
      }
    }
  }
  if (!is.null(result$do_k)) {
    free <- result$do_k[result$do_k$model == "free_k", , drop = FALSE]
    if (nrow(free) == 1L && isTRUE(free$converged)) {
      message("Free-k comparison fit converged.")
    } else {
      message("Free-k comparison: no valid converged fit was obtained.")
    }
  }
  .modgee_report_dark(result)
  invisible(NULL)
}

.modgee_report_dark <- function(result) {
  explanation <- paste0(
    "The estimated number of dark extinctions is negative. ",
    "One possible explanation is that some species surviving at the final ",
    "timestep remain undetected. ",
    "This negative estimate does not, by itself, invalidate the estimated ",
    "cumulative extinction rate or fitted parameters.")
  if (isTRUE(result$summary$dark_extinct < 0)) message(explanation)
  groups <- result$group_results$summary
  if (!is.null(groups)) {
    for (i in which(!is.na(groups$dark_extinct) & groups$dark_extinct < 0)) {
      message("Group '", groups$group[i], "': ", explanation)
    }
  }
  invisible(NULL)
}

# Suppress only the warning duplicated by the fit summary; stream progress.
.modgee_quiet_fit <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (identical(conditionMessage(w),
                    "An optimization did not report normal convergence.")) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

.modgee_report_k <- function(comparison) {
  fixed <- comparison[comparison$model == "fixed_k", , drop = FALSE]
  free <- comparison[comparison$model == "free_k", , drop = FALSE]
  valid <- nrow(fixed) == 1L && nrow(free) == 1L &&
    isTRUE(fixed$converged) && isTRUE(free$converged) &&
    is.finite(fixed$AIC) && is.finite(free$AIC) && is.finite(free$k)
  if (!valid) {
    cat("The AIC comparison is unavailable because one or both models did not yield a valid converged fit.\n")
  } else if (free$AIC > fixed$AIC) {
    cat("The free-k model has a higher AIC than the fixed-k model, indicating no clear evidence of a relationship between extinction risk and detectability.\n")
  } else if (free$AIC < fixed$AIC) {
    cat("The free-k model has a lower AIC than the fixed-k model, indicating evidence of a relationship between extinction risk and detectability.\n")
    if (free$k < 0) {
      cat("The estimated k is negative, indicating that less detectable species have higher extinction risk. This may cause the fixed-k model to underestimate the cumulative extinction rate.\n")
    } else if (free$k > 0) {
      cat("The estimated k is positive, indicating that more detectable species have higher extinction risk. This may cause the fixed-k model to overestimate the cumulative extinction rate.\n")
    } else {
      cat("The estimated k is zero; no direction of association is indicated.\n")
    }
  } else {
    cat("The two models have equal AIC values; no difference in model support is indicated.\n")
  }
  invisible(NULL)
}
