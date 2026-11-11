#' Print a MODGEE result
#'
#' Display the overall summary, group summaries and k comparison when present.
#' The k comparison includes an AIC interpretation and, where supported, its direction.
#' The complete object is retained, including species results and diagnostics.
#' @param x A result returned by \code{fit_modgee()}.
#' @param ... Additional arguments passed to \code{print()}.
#' @return The complete object, invisibly.
#' @export
print.modgee_result <- function(x, ...) {
  cat("$summary\n")
  print(x$summary, ...)
  if (!is.null(x$group_results$summary)) {
    cat("\n$group_results$summary\n")
    print(x$group_results$summary, ...)
  }
  if (!is.null(x$do_k)) {
    cat("\n$do_k\n")
    print(x$do_k, ...)
    .modgee_report_k(x$do_k)
  }
  invisible(x)
}
