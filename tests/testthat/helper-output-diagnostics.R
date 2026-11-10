# Expose internals only in historical numerical regression tests.
.fit_modgee_with_diagnostics <- function(...) {
  f <- .fit_modgee_impl
  b <- body(f)
  result <- b[[length(b)]]
  stopifnot(is.call(result), identical(result[[1L]], as.name("list")))
  fields <- alist(fit_species_order = species_names[fit_rows],
                 reference_species = species_names[anchor],
                 fit = fit, statistics = statistics)
  result <- as.call(c(as.list(result), fields))
  b[[length(b)]] <- result
  body(f) <- b
  f(...)
}
