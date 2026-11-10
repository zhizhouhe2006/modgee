test_that("invalid inferred time labels are rejected", {
  M <- matrix(c(1, 1, 0, 0, 1, 0, 1, 1), nrow = 2, byrow = TRUE,
              dimnames = list(c("species_a", "species_b"), NULL))
  colnames(M) <- c("X1990", "X2000", "X2010", "X2020")
  expect_error(.fit_modgee_with_diagnostics(M), "must be finite numbers", fixed = TRUE)
  colnames(M) <- c("1990", "2000", "Inf", "2020")
  expect_error(.fit_modgee_with_diagnostics(M), "must be finite numbers", fixed = TRUE)
  colnames(M) <- c("2020", "2010", "2000", "1990")
  expect_error(.fit_modgee_with_diagnostics(M), "strictly increasing", fixed = TRUE)
  colnames(M) <- c("1990", "2000", "2000", "2020")
  expect_error(.fit_modgee_with_diagnostics(M), "must be unique", fixed = TRUE)
  colnames(M) <- c("1990", "2000", "2015", "2020")
  expect_error(.fit_modgee_with_diagnostics(M), "equally spaced", fixed = TRUE)
})

test_that("time labels must be supplied through column names", {
  M <- matrix(1, 2, 4, dimnames = list(c("species_a", "species_b"),
                                      c("1990", "2000", "2010", "2020")))
  # The unused-argument error is localized, so do not match English wording.
  expect_error(.fit_modgee_with_diagnostics(M, times = c(1990, 2000, 2010, 2020)))
  colnames(M) <- NULL
  expect_error(.fit_modgee_with_diagnostics(M), "Time column names are required", fixed = TRUE)
})

test_that("explicit interval-index column labels preserve model estimates", {
  ref <- readRDS(test_path("fixtures", "birds-reference.rds"))
  M <- ref$M
  years <- as.numeric(ref$times)
  colnames(M) <- as.character(years)
  inferred <- .fit_modgee_with_diagnostics(M, starts = ref$starts)
  expect_equal(inferred$times, years)

  renamed <- M
  colnames(renamed) <- paste0("interval_", seq_len(ncol(M)))
  expect_error(.fit_modgee_with_diagnostics(renamed, starts = ref$starts),
               "must be finite numbers", fixed = TRUE)

  # Interval indices must now be assigned explicitly, rather than inferred.
  indexed_M <- M
  colnames(indexed_M) <- as.character(seq_len(ncol(M)))
  indexed <- .fit_modgee_with_diagnostics(indexed_M, starts = ref$starts)
  expect_equal(indexed$times, seq_len(ncol(M)))
  expect_equal(unname(indexed$parameters$d), unname(inferred$parameters$d), tolerance = 1e-8)
  expect_equal(unname(indexed$parameters$h), unname(inferred$parameters$h), tolerance = 1e-8)
  expect_equal(unname(indexed$parameters$m), unname(inferred$parameters$m), tolerance = 1e-8)
  expect_equal(indexed$fit$LL, inferred$fit$LL, tolerance = 1e-8)
  expect_equal(indexed$curve$cumulative_extinction_rate,
               inferred$curve$cumulative_extinction_rate, tolerance = 1e-8)
  expect_equal(indexed$species$p_extinct, inferred$species$p_extinct, tolerance = 1e-8)
  expect_equal(indexed$statistics$E_t_exts[match(indexed$species$species, indexed$fit_species_order)],
               inferred$statistics$E_t_exts[match(inferred$species$species, inferred$fit_species_order)], tolerance = 1e-8)
  expect_equal(indexed$statistics$var_t_exts[match(indexed$species$species, indexed$fit_species_order)],
               inferred$statistics$var_t_exts[match(inferred$species$species, inferred$fit_species_order)], tolerance = 1e-8)
  expect_equal(indexed$species$estimated_extinction_year,
               indexed$statistics$E_t_exts[match(indexed$species$species, indexed$fit_species_order)], tolerance = 1e-8)
  expect_equal(inferred$species$estimated_extinction_year,
               years[1] + (inferred$statistics$E_t_exts[match(inferred$species$species, inferred$fit_species_order)] - 1) * diff(years)[1],
               tolerance = 1e-8)
  expect_equal(indexed$summary$dark_extinct, inferred$summary$dark_extinct, tolerance = 1e-8)
})
