.modgee_singleton_loglik <- function(numerator, denominator) {
  LL <- log(numerator / denominator)
  if (is.nan(LL)) {
    error <- simpleError("The singleton log-likelihood is undefined (NaN).", call = NULL)
    class(error) <- c("modgee_singleton_undefined", class(error))
    stop(error)
  }
  LL
}

.modgee_singleton_initial_value <- function(expr, species) {
  tryCatch(expr, modgee_singleton_undefined = function(e) {
    stop(
      "Singleton fitting failed for species '", species,
      "' during the initial-value search.\n",
      "The log-likelihood could not be evaluated numerically (NaN). ",
      "This may reflect sparse detection histories or numerical instability ",
      "in the fitted parameters.\n",
      "No fit was returned. Review the detection histories and earlier ",
      "optimization warnings before retrying.", call. = FALSE
    )
  })
}
